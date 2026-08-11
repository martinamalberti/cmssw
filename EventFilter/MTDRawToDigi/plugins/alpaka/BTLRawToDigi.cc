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

#include "DataFormats/FTLDigiSoA/interface/BTLDigiHostCollection.h"
#include "DataFormats/FTLDigiSoA/interface/alpaka/BTLDigiDeviceCollection.h"
#include "CondFormats/DataRecord/interface/BTLReadoutMapRcd.h"
#include "EventFilter/MTDRawToDigi/interface/BTLElectronicsSpecs.h"

#include "BTLRawToDigiInputData.h"
#include "BTLRawToDigiAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using ::btldigi::BTLDigiHostCollection;
  using btldigi::BTLDigiDeviceCollection;
  
  class BTLRawToDigi : public stream::EDProducer<> { // stream::EDProdicer vedi: https://twiki.cern.ch/twiki/bin/view/CMSPublic/FWMultithreadedFrameworkModuleTypes?utm_source=chatgpt.com#Comparing_Stream_and_Global_Modu
  public:
    explicit BTLRawToDigi(const edm::ParameterSet& iConfig);
    ~BTLRawToDigi() override = default;
    static void fillDescriptions(edm::ConfigurationDescriptions&);
    
    void produce(device::Event&, device::EventSetup const&) override;

  private:
    edm::EDGetTokenT<RawDataBuffer> rawDataBufferToken_;
    //edm::EDPutTokenT<btldigi::BTLDigiDeviceCollection> digiPutToken_;
    edm::EDPutTokenT<BTLDigiHostCollection> digiPutToken_; // BTLBaseRecHitSoAProducer consuma una BTLDigiHostCollection
    device::ESGetToken<BTLElectronicsToDetIdMappingDevice, BTLReadoutMapRcd> elecToDetIdToken_;
    
    BTLRawToDigiAlgo algo_;
    BTLElectronicsIndexer indexer_;
    
    // reused across events - member state is fine in stream::EDProducer
    // (one instance per stream, not shared across streams)
    int32_t channelCapacity_ = 0;
    std::unique_ptr<BTLRawToDigiInputData> inputDataHost_;
    
    static constexpr int kWordsPerChannel = 2;
  };    


  

  BTLRawToDigi::BTLRawToDigi(const edm::ParameterSet& iConfig)
    : EDProducer(iConfig),
      rawDataBufferToken_(consumes<RawDataBuffer>(iConfig.getParameter<edm::InputTag>("rawDataBufferTag"))),
      digiPutToken_(produces()),
      elecToDetIdToken_(esConsumes()) {
    indexer_.firstFedId = BTLElectronicsSpecs::kFirstFEDId;
    indexer_.hsLinkOffset = BTLElectronicsSpecs::kHSLinksOffset;
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

      if (payloadBytes % kChannelBytes != 0) {
	throw cms::Exception("BTLRawToDigi")
	  << "FED " << fed
	  << ": payload size is not aligned to channel size. "
	  << "payloadBytes = " << payloadBytes
	  << ", kChannelBytes = " << kChannelBytes
	  << ", remainder = " << payloadBytes % kChannelBytes;
      }
      std::memcpy(inputDataHost_->rawWords.data() + kWordsPerChannel * offset, payload.data(), payloadBytes); // copia payloadBytes da fedData (dopo aver tolto header e trailer) in rawWords_h_, con un offset kWordsPerChannel * offset che corrisponde a 2 x numero di canali gia' copiati: rawWords_h_.data() + 2 * channelOffset e' un puntatore alla posizione di destinazione in rawWords_h_. Il fattore 2 serve perche' ogni canale sono 128 bit, quindi occupa 2 unit64. 
      const int32_t nChInFed = static_cast<int32_t>(payloadBytes / kChannelBytes);
      std::fill(inputDataHost_->channelFedId.data() + offset, inputDataHost_->channelFedId.data() + offset + nChInFed, fed);
      offset += nChInFed;
    }


    // --- hand off to the device algo
    auto digis_d = algo_.process(queue, inputDataHost_->rawWords.data(), inputDataHost_->channelFedId.data(), nChannels, indexer_, elecToDetId); // nCHannels serve per delimitare la parte valida del buffer.

    // Device -> Host
    BTLDigiHostCollection digis_h(digis_d.size());

    alpaka::memcpy(queue, digis_h.buffer(), digis_d.buffer());

    alpaka::wait(queue);    
    
    iEvent.emplace(digiPutToken_, std::move(digis_h));
    
  }
  
  void BTLRawToDigi::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("rawDataBufferTag", edm::InputTag("btlRaw")); // ???
    descriptions.addWithDefaultLabel(desc);

  }
  
  
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(BTLRawToDigi);
