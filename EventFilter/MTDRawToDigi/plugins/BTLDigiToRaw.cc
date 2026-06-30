#include <vector>
#include <map>
#include <cstdint>
#include <cstring>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FEDRawData/interface/RawDataBuffer.h"
#include "DataFormats/FTLDigi/interface/BTLDigi.h"
#include "DataFormats/FTLDigi/interface/MTDDigiCollections.h"

#include "CondFormats/MTDObjects/interface/BTLReadoutMap.h"
#include "CondFormats/MTDObjects/interface/BTLElectronicsId.h"
#include "CondFormats/DataRecord/interface/BTLReadoutMapRcd.h"

#include "EventFilter/MTDRawToDigi/interface/BitStream.h"

namespace btldigitoraw {
  using ChannelStream = BitStream<2>;
}

class BTLDigiToRaw : public edm::one::EDProducer<> {
public:
  explicit BTLDigiToRaw(const edm::ParameterSet&);
  ~BTLDigiToRaw() override;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  
  // -------------------------------------------------------------------
  // Encode one TOFHIR channel into a raw word.
  // Returns a BitStream<2> (128 bit)
  // -------------------------------------------------------------------
  btldigitoraw::ChannelStream encodeChannelPayload(int fed, int hslink, int elink,
						   const btldigi::BTLDigi& digi,
						   bool isPlus) const;

  // -------------------------------------------------------------------
  // Write the accumulated raw words for one FED/source ID
  // into the RawDataBuffer, adding header and trailer (128 bit each)
  // -------------------------------------------------------------------
  void fillFEDBuffer(int fedId,
                     const std::vector<btldigitoraw::ChannelStream>& words,
                     RawDataBuffer& rawData) const;

  // -- Input tokens
  const edm::EDGetTokenT<BTLDigiContentCollection> digiToken_;
  const edm::ESGetToken<BTLReadoutMap, BTLReadoutMapRcd> readoutMapToken_;
};


BTLDigiToRaw::BTLDigiToRaw(const edm::ParameterSet& iConfig)
    : digiToken_(consumes<BTLDigiContentCollection>(iConfig.getParameter<edm::InputTag>("btlDigiLabel"))),
      readoutMapToken_(esConsumes<BTLReadoutMap, BTLReadoutMapRcd>()) {
  produces<RawDataBuffer>();
}

BTLDigiToRaw::~BTLDigiToRaw() {}


void BTLDigiToRaw::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // -- Retrieve the BTL digi collection (sorted by BTLDetId)
  const auto& digis = iEvent.get(digiToken_);

  // -- Retrieve readout map from EventSetup
  const BTLReadoutMap& readoutMap = iSetup.getData(readoutMapToken_);

  // -- Get EventID and RunID
  unsigned int eventId_ = iEvent.id().event();
  
  // -- Output
  auto rawDataBuffer = std::make_unique<RawDataBuffer>();


  // -- Linear scan: the digi collection is sorted by BTLDetId, so crystals
  // belonging to the same FED are contiguous. We accumulate channel
  // streams into a single buffer and flush it as soon as the FED id
  // changes
  std::vector<btldigitoraw::ChannelStream> currentChannelsStream;
  int currentFed = -1;
  
  for (const auto& digi : digis) {
    BTLDetId detId(digi.krawId());

    // -- Forward lookup: DetId -> electronics
    BTLElectronicsIdPair elecIds = readoutMap.getElectronicsId(detId);
    
    // -- Sanity check: both sides must belong to the same FED and must have same hs-link and e-link Ids
    if (elecIds.minus.fedId() != elecIds.plus.fedId()) {
      edm::LogError("BTLDigiToRaw") << "BTLDigiToRaw::produce(): "
				    << "minus and plus sides of crystal " << std::hex << detId.rawId() << std::dec
				    << " belong to different FEDs ("
				    << elecIds.minus.fedId() << " vs " << elecIds.plus.fedId()
				    << ") -- skipping crystal.";
      continue;
    }
    
    const int fedId = elecIds.minus.fedId();
    const int hslinkId = elecIds.minus.hsLinkId();
    const int elinkId = elecIds.minus.eLinkId();

    // Check if FED id changed and flush the buffer accumulated so far and start a new one.
    if (fedId != currentFed) {
      if (currentFed >= 0) {
        fillFEDBuffer(currentFed, currentChannelsStream, *rawDataBuffer); // TO-DO
      }
      currentChannelsStream.clear();
      currentFed = fedId;
    }
    
    // Encode minus and plus side
    currentChannelsStream.push_back(encodeChannelPayload(fedId, hslinkId, elinkId, digi, false));
    currentChannelsStream.push_back(encodeChannelPayload(fedId, hslinkId, elinkId, digi, true));
  }

  iEvent.put(std::move(rawDataBuffer));
}

// ------------------------------------------------------------
// encodeChannelPayload
// packs all fields (for the requested side) into a
// 128-bit ChannelStream.
// ------------------------------------------------------------
btldigitoraw::ChannelStream BTLDigiToRaw::encodeChannelPayload(int fed, int hslink, int elink,  
							       const btldigi::BTLDigi& digi,
							       bool isPlusSide) const {
  btldigitoraw::ChannelStream stream;
  
  // Get the fields for the requested side
  uint16_t BC0count = digi.kBC0count();
  bool status = digi.kstatus();
  uint32_t BCcount = digi.kBCcount();
  uint8_t chID = isPlusSide ? digi.kchIDPlus() : digi.kchIDMinus();
  uint16_t T1coarse = isPlusSide ? digi.kT1coarsePlus() : digi.kT1coarseMinus();
  uint16_t T2coarse = isPlusSide ? digi.kT2coarsePlus() : digi.kT2coarseMinus();
  uint16_t EOIcoarse = isPlusSide ? digi.kEOIcoarsePlus() : digi.kEOIcoarseMinus();
  uint16_t Charge = isPlusSide ? digi.kChargePlus() : digi.kChargeMinus();
  uint16_t T1fine = isPlusSide ? digi.kT1finePlus() : digi.kT1fineMinus();
  uint16_t T2fine = isPlusSide ? digi.kT2finePlus() : digi.kT2fineMinus();
  uint16_t IdleTime = isPlusSide ? digi.kIdleTimePlus() : digi.kIdleTimeMinus();
  uint8_t PrevTrigF = isPlusSide ? digi.kPrevTrigFPlus() : digi.kPrevTrigFMinus();
  uint8_t TACID = isPlusSide ? digi.kTACIDPlus() : digi.kTACIDMinus();

  int slink = fed; 
  
  stream.set_bits(127, 10, static_cast<uint64_t>(BC0count));
  stream.set_bits(117,  1, static_cast<uint64_t>(status));
  stream.set_bits(110,  7, static_cast<uint64_t>(slink)); // will be removed?  0-12, but FEDId will have a offset....
  stream.set_bits(104,  6, static_cast<uint64_t>(hslink)); 
  stream.set_bits( 99,  5, static_cast<uint64_t>(elink));
  stream.set_bits( 94,  5, static_cast<uint64_t>(chID));
  stream.set_bits( 82, 12, static_cast<uint64_t>(BCcount));
  stream.set_bits( 67, 15, static_cast<uint64_t>(T1coarse));
  stream.set_bits( 57, 10, static_cast<uint64_t>(T2coarse));
  stream.set_bits( 47, 10, static_cast<uint64_t>(EOIcoarse));
  stream.set_bits( 37, 10, static_cast<uint64_t>(Charge));
  stream.set_bits( 27, 10, static_cast<uint64_t>(T1fine));
  stream.set_bits( 17, 10, static_cast<uint64_t>(T2fine));
  stream.set_bits(  7, 10, static_cast<uint64_t>(IdleTime));
  stream.set_bits(  3,  4, static_cast<uint64_t>(PrevTrigF));
  stream.set_bits(  0,  3, static_cast<uint64_t>(TACID));
  
  return stream;
}

// ------------------------------------------------------------
// Fill FED buffer
// Writes header (128 bit) + channel streams (128 bit each) + trailer
// (128 bit) into the RawDataBuffer for the given source ID.
// ------------------------------------------------------------
void BTLDigiToRaw::fillFEDBuffer(int fedId,
				 const std::vector<btldigitoraw::ChannelStream>& channelsStream,
				 RawDataBuffer& rawDataBuffer) const {
  

  // TODO TODO
  
  // -- Fill header

  // -- payload

  // -- Fill trailer


  ///  rawDataBuffer->addSource(.....) ???
  
}


// ------------------------------------------------------------
// fillDescriptions
// ------------------------------------------------------------
void BTLDigiToRaw::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("btlDigiLabel", edm::InputTag("btlDigiLabel", "BTLDigi"));
  descriptions.addWithDefaultLabel(desc);
}

