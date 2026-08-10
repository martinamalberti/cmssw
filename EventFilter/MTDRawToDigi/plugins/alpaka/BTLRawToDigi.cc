#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

#include "DataFormats/FEDRawData/interface/RawDataBuffer.h"
#include "DataFormats/FEDRawData/interface/SLinkRocketHeaders.h"

#include "DataFormats/FTLDigiSoA/interface/BTLDigiSoA.h"
#include "DataFormats/FTLDigiSoA/interface/alpaka/BTLDigiDeviceCollection.h"
#include "CondFormats/DataRecord/interface/BTLReadoutMapRcd.h"
#include "EventFilter/MTDRawToDigi/interface/BTLElectronicsSpecs.h"

#include "BTLRawToDigiInputData.h"
#include "BTLRawToDigiAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using btldigi::BTLDigiDeviceCollection;

  class BTLRawToDigi : public stream::EDProducer<> {
  public:
    explicit BTLRawToDigi(const edm::ParameterSet& iConfig);
    ~BTLRawToDigi() override = default;
    static void fillDescriptions(edm::ConfigurationDescriptions&);
    
    void produce(device::Event&, device::EventSetup const&) override;

  private:
    edm::EDGetTokenT<RawDataBuffer> rawDataBufferToken_;
    device::ESGetToken<BTLElectronicsToDetIdMappingDevice, BTLReadoutMapRcd> elecToDetIdToken_;
    edm::EDPutTokenT<BTLDigiDeviceCollection> digiPutToken_;
    
    BTLRawToDigiAlgo algo_;
    BTLElectronicsIndexer indexer_;
    
    // reused across events - member state is fine in stream::EDProducer
    // (one instance per stream, not shared across streams)
    int32_t channelCapacity_ = 0; ///??????????????????????
    //cms::alpakatools::host_buffer<uint64_t[]> rawWords_h_;
    //cms::alpakatools::host_buffer<int32_t[]> channelFedId_h_;
    std::unique_ptr<BTLRawToDigiInputData> inputDataHost_;
    
    static constexpr int kWordsPerChannel = 2;
  };    


  

  BTLRawToDigi::BTLRawToDigi(const edm::ParameterSet& iConfig)
    : EDProducer(iConfig),
      rawDataBufferToken_(consumes<RawDataBuffer>(iConfig.getParameter<edm::InputTag>("rawDataBufferTag"))),
      elecToDetIdToken_(esConsumes()) {
    indexer_.firstFedId = BTLElectronicsSpecs::kFirstFEDId;
    indexer_.nFeds = BTLElectronicsSpecs::kNumberOfFEDs;
    indexer_.nHsLinks = BTLElectronicsSpecs::kNumberOfHsLinks;
    indexer_.nELinks = BTLElectronicsSpecs::kNumberOfELinks;
  }


  void BTLRawToDigi::produce(device::Event& iEvent, device::EventSetup const& iSetup) {
    auto const& rawDataBuffer = iEvent.get(rawDataBufferToken_);
    auto const& elecToDetId = iSetup.getData(elecToDetIdToken_);
    auto& queue = iEvent.queue();

    static constexpr std::size_t kHeaderSize  = sizeof(SLinkRocketHeader_v3);
    static constexpr std::size_t kTrailerSize = sizeof(SLinkRocketTrailer_v3);
    static constexpr std::size_t kChannelBytes = kWordsPerChannel * sizeof(uint64_t);  // 128 bit / channel
    
    // --- Make a first iteration over the FEDs to compute the total buffer size
    int32_t nChannels = 0;
    for (uint32_t fed = BTLElectronicsSpecs::kFirstFEDId; fed < BTLElectronicsSpecs::kFirstFEDId + BTLElectronicsSpecs::kNumberOfFEDs; ++fed) {
      const auto& fedData = rawDataBuffer.fragmentData(fed);
      if (fedData.size() == 0)
	continue;
      const std::size_t payloadBytes = fedData.size() - kHeaderSize - kTrailerSize;
      nChannels += static_cast<int32_t>(payloadBytes / kChannelBytes);
    }
    
    // --- reuse a pinned staging buffer across events; only grow, never shrink/realloc every event 
    if (nChannels > channelCapacity_) {
      channelCapacity_ = nChannels + nChannels / 4;  // headroom, avoid re-growing every event (+25% arbitrary)
      inputDataHost_ = std::make_unique<BTLRawToDigiInputData>(queue, channelCapacity_);
      //rawWords_h_ = cms::alpakatools::make_host_buffer<uint64_t[]>(queue, 2 * channelCapacity_);
      //channelFedId_h_ = cms::alpakatools::make_host_buffer<int32_t[]>(queue, channelCapacity_);
    }
    

    // --- pass 2:
    // Copy the FED payload (header/trailer removed) into the contiguous staging buffer.
    // Each channel occupies two uint64_t words, and fill channelFedId_h in the same loop
        
    int32_t offset = 0; // == number of channels already copied
    for (uint32_t fed = BTLElectronicsSpecs::kFirstFEDId; fed < BTLElectronicsSpecs::kFirstFEDId + BTLElectronicsSpecs::kNumberOfFEDs; ++fed) {
      const auto& fedData = rawDataBuffer.fragmentData(fed);
      if (fedData.size() == 0)
	continue;
      const auto payload = fedData.payload(kHeaderSize, kTrailerSize);
      const std::size_t payloadBytes = fedData.size() - kHeaderSize - kTrailerSize;
      std::memcpy(inputDataHost_->rawWords.data() + kWordsPerChannel * offset, payload.data() + kHeaderSize, payloadBytes); // copia payloadBytes da fedData (skippando l'header), dentro rawWords_h_, con un offset 2 * channelOffset che corrisponde a 2 x numero di canali gia' copiati: rawWords_h_.data() + 2 * channelOffset e' un puntatore alla posizione di destinazione in rawWords_h_. Il fattore 2 serve perche' ogni canale sono 128 bit, quindi occupa 2 unit64. 
      const int32_t nChInFed = static_cast<int32_t>(payloadBytes / kChannelBytes);
      std::fill(inputDataHost_->channelFedId.data() + offset, inputDataHost_->channelFedId.data() + offset + nChInFed, fed);
      offset += nChInFed;
    }


    // --- hand off to the device algo
    //    auto digis = algo_.process(queue, rawWords_h_.data(), channelFedId_h_.data(), nChannels, indexer_, elecToDetId); // nCHannels serve per delimitare la parte valida del buffer.
    auto digis = algo_.process(queue, inputDataHost_->rawWords.data(), inputDataHost_->channelFedId.data(), nChannels, indexer_, elecToDetId); // nCHannels serve per delimitare la parte valida del buffer.
    
    iEvent.emplace(digiPutToken_, std::move(digis));
    
  }
  
  void BTLRawToDigi::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("rawDataBufferTag", edm::InputTag("rawDataBufferTag")); // ???
    descriptions.addWithDefaultLabel(desc);
  }
  
  
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(BTLRawToDigi);
