#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/SynchronizingEDProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

#include "DataFormats/FEDRawData/interface/RawDataBuffer.h"
#include "DataFormats/FEDRawData/interface/SLinkRocketHeaders.h"

#include "DataFormats/FTLDigiSoA/interface/alpaka/BTLDigiDeviceCollection.h"
#include "CondFormats/DataRecord/interface/BTLReadoutMapRcd.h"
#include "EventFilter/MTDRawToDigi/interface/BTLElectronicsSpecs.h"

#include "BTLRawToDigiInputData.h"
#include "BTLRawToDigiAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using btldigi::BTLDigiDeviceCollection;

  // SynchronizingEDProducer: the work up to the count of the digis is enqueued in acquire(),
  // then the framework waits for the queue in a non-blocking way (ExternalWork) before calling
  // produce(), where nDigisTot is available on the host to allocate the output collection.
  class BTLRawToDigi : public stream::SynchronizingEDProducer<> {
  public:
    explicit BTLRawToDigi(const edm::ParameterSet& iConfig);
    ~BTLRawToDigi() override = default;
    static void fillDescriptions(edm::ConfigurationDescriptions&);

    void acquire(device::Event const&, device::EventSetup const&) override;
    void produce(device::Event&, device::EventSetup const&) override;

  private:
    edm::EDGetTokenT<RawDataBuffer> rawDataBufferToken_;
    // the device collection is put directly: host consumers of BTLDigiHostCollection (e.g. BTLBaseRecHitSoAProducer)
    // get the copy to host done (asynchronously) by the framework
    device::EDPutToken<BTLDigiDeviceCollection> digiPutToken_;
    device::ESGetToken<BTLElectronicsToDetIdMappingDevice, BTLReadoutMapRcd> elecToDetIdToken_;

    BTLRawToDigiAlgo algo_;
    BTLElectronicsIndexer indexer_;

    // reused across events - member state is fine in stream::EDProducer
    // (one instance per stream, not shared across streams)
    int32_t channelCapacity_ = 0;
    std::unique_ptr<BTLRawToDigiInputData> inputDataHost_;

    // filled in acquire(), read in produce()
    int32_t nChannels_ = 0;
    std::optional<cms::alpakatools::host_buffer<int32_t>> nDigisTot_h_;

    static constexpr int kWordsPerChannel = 2;
  };

  BTLRawToDigi::BTLRawToDigi(const edm::ParameterSet& iConfig)
      : SynchronizingEDProducer(iConfig),
        rawDataBufferToken_(consumes<RawDataBuffer>(iConfig.getParameter<edm::InputTag>("rawDataBufferTag"))),
        digiPutToken_(produces()),
        elecToDetIdToken_(esConsumes()) {
    indexer_.firstFedId = BTLElectronicsSpecs::kFirstFEDId;
    indexer_.nFeds = BTLElectronicsSpecs::kNumberOfFEDs;
    indexer_.nFEDsPerSerenity = BTLElectronicsSpecs::nFEDsPerSerenity;
    indexer_.nHSLinksPerSerenity = BTLElectronicsSpecs::nHSLinksPerSerenity;
    indexer_.nELinksPerHSLink = BTLElectronicsSpecs::nELinksPerHSLink;
    indexer_.nChannelsPerELink = BTLElectronicsSpecs::nChannelsPerELink;
  }

  void BTLRawToDigi::acquire(device::Event const& iEvent, device::EventSetup const& iSetup) {
    auto const& rawDataBuffer = iEvent.get(rawDataBufferToken_);
    auto const& elecToDetId = iSetup.getData(elecToDetIdToken_);
    auto& queue = iEvent.queue();

    static constexpr std::size_t kHeaderSize = sizeof(SLinkRocketHeader_v3);
    static constexpr std::size_t kTrailerSize = sizeof(SLinkRocketTrailer_v3);
    static constexpr std::size_t kChannelBytes = kWordsPerChannel * sizeof(uint64_t);  // 128 bit / channel

    // --- Make a first iteration over the FEDs to compute the total buffer size
    int32_t nChannels = 0;
    for (uint32_t fed = BTLElectronicsSpecs::kFirstFEDId;
         fed < BTLElectronicsSpecs::kFirstFEDId + BTLElectronicsSpecs::kNumberOfFEDs;
         ++fed) {
      const auto& fedData = rawDataBuffer.fragmentData(fed);
      if (fedData.size() == 0)
        continue;
      const std::size_t payloadBytes = fedData.size() - kHeaderSize - kTrailerSize;
      nChannels += static_cast<int32_t>(payloadBytes / kChannelBytes);
    }

    // --- if no channels in the payload, nothing to do on the device: produce() puts an empty collection
    nChannels_ = nChannels;
    if (nChannels == 0)
      return;

    // --- reuse a pinned staging buffer across events; only grow, never shrink/realloc every event
    if (nChannels > channelCapacity_) {
      channelCapacity_ = nChannels + nChannels / 4;  // headroom, avoid re-growing every event (+25% arbitrary)
      inputDataHost_ = std::make_unique<BTLRawToDigiInputData>(queue, channelCapacity_);
    }

    // --- pass 2:
    // Copy the FED payload (header/trailer removed) into the contiguous staging buffer.
    // Each channel occupies two uint64_t words, and fill channelFedId_h in the same loop

    int32_t offset = 0;  // == number of channels already copied
    for (uint32_t fed = BTLElectronicsSpecs::kFirstFEDId;
         fed < BTLElectronicsSpecs::kFirstFEDId + BTLElectronicsSpecs::kNumberOfFEDs;
         ++fed) {
      const auto& fedData = rawDataBuffer.fragmentData(fed);
      if (fedData.size() == 0)
        continue;
      const auto payload = fedData.payload(kHeaderSize, kTrailerSize);
      const std::size_t payloadBytes = fedData.size() - kHeaderSize - kTrailerSize;

      if (payloadBytes % kChannelBytes != 0) {
        throw cms::Exception("BTLRawToDigi")
            << "FED " << fed << ": payload size is not aligned to channel size. "
            << "payloadBytes = " << payloadBytes << ", kChannelBytes = " << kChannelBytes
            << ", remainder = " << payloadBytes % kChannelBytes;
      }
      // copy payloadBytes from fedData payload (header and trailer removed) in rawWords_h_,
      // with an offset kWordsPerChannel * offset corresponding to 2 x number of channels laready copied
      std::memcpy(inputDataHost_->rawWords.data() + kWordsPerChannel * offset, payload.data(), payloadBytes);
      const int32_t nChInFed = static_cast<int32_t>(payloadBytes / kChannelBytes);
      std::fill(
          inputDataHost_->channelFedId.data() + offset, inputDataHost_->channelFedId.data() + offset + nChInFed, fed);
      offset += nChInFed;
    }

    // --- hand off to the device algo: decode, pairing and count of the digis.
    // The total number of digis is copied asynchronously to nDigisTot_h_
    nDigisTot_h_.emplace(cms::alpakatools::make_host_buffer<int32_t>(queue));
    algo_.prepareDigis(queue,
                       inputDataHost_->rawWords.data(),
                       inputDataHost_->channelFedId.data(),
                       nChannels,
                       elecToDetId,
                       indexer_,
                       *nDigisTot_h_);
  }

  void BTLRawToDigi::produce(device::Event& iEvent, device::EventSetup const& iSetup) {
    // --- no channels in the payload: return empty collection
    if (nChannels_ == 0) {
      iEvent.emplace(digiPutToken_, BTLDigiDeviceCollection(iEvent.queue(), 0));
      return;
    }

    auto const& elecToDetId = iSetup.getData(elecToDetIdToken_);
    auto& queue = iEvent.queue();

    // the queue has been synchronized by the framework after acquire(): nDigisTot is available
    const int32_t nDigisTot = *nDigisTot_h_->data();
    nDigisTot_h_.reset();

    // no wait: the fill kernel may still be running when produce() returns,
    // the framework synchronizes the consumers with the queue
    iEvent.emplace(digiPutToken_, algo_.produceDigis(queue, elecToDetId, indexer_, nDigisTot));
  }

  void BTLRawToDigi::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("rawDataBufferTag", edm::InputTag("btlRaw"));  // ???
    descriptions.addWithDefaultLabel(desc);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(BTLRawToDigi);
