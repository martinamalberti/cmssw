#include <vector>
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
#include "DataFormats/FEDRawData/interface/SLinkRocketHeaders.h"
#include "DataFormats/FTLDigi/interface/BTLDigi.h"
#include "DataFormats/FTLDigi/interface/MTDDigiCollections.h"

#include "CondFormats/MTDObjects/interface/BTLReadoutMap.h"
#include "CondFormats/MTDObjects/interface/BTLElectronicsId.h"
#include "CondFormats/DataRecord/interface/BTLReadoutMapRcd.h"

#include "EventFilter/MTDRawToDigi/interface/BitStream.h"
#include "EventFilter/MTDRawToDigi/interface/BTLElectronicsSpecs.h"

namespace btldigitoraw {

  using ChannelStream = BitStream<2>;
  using Word = ChannelStream::word_t;
  
  const int kWordsPerChannel = ChannelStream::BITS / ChannelStream::WORD_BITS;

  constexpr std::size_t kSlinkHeaderSize  = sizeof(SLinkRocketHeader_v3);
  constexpr std::size_t kSlinkTrailerSize = sizeof(SLinkRocketTrailer_v3);


  
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
  btldigitoraw::ChannelStream encodeChannelPayload(int fed, int hslink, int elink, const btldigi::BTLDigi& digi, bool isPlus) const;

  // -------------------------------------------------------------------
  // Write the accumulated raw words for one FED/source ID
  // into the RawDataBuffer, adding header and trailer (128 bit each)
  // -------------------------------------------------------------------
  void fillFEDBuffer(int fedId,  uint64_t eventId, const std::vector<btldigitoraw::ChannelStream>& channelsStream, RawDataBuffer& rawDataBuffer) const;
  
  std::vector<btldigitoraw::ChannelStream> currentChannelsStream_;
  
  // -- Input tokens
  const edm::EDGetTokenT<BTLDigiContentCollection> digiToken_;
  const edm::ESGetToken<BTLReadoutMap, BTLReadoutMapRcd> readoutMapToken_;
};


BTLDigiToRaw::BTLDigiToRaw(const edm::ParameterSet& iConfig)
    : digiToken_(consumes<BTLDigiContentCollection>(iConfig.getParameter<edm::InputTag>("btlDigiCollection"))),
      readoutMapToken_(esConsumes<BTLReadoutMap, BTLReadoutMapRcd>()) {
  produces<RawDataBuffer>();
}


BTLDigiToRaw::~BTLDigiToRaw() {}


void BTLDigiToRaw::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace btldigitoraw;
  
  // -- Retrieve the BTL digi collection (sorted by BTLDetId)
  const auto& digis = iEvent.get(digiToken_);
  
  // -- Retrieve readout map from EventSetup
  const BTLReadoutMap& readoutMap = iSetup.getData(readoutMapToken_);
  
  // -- Vars needed to fill SLinkRocket header and trailer
  uint64_t eventId = iEvent.id().event();
  
  // -- Output
  // pre-allocated RawDataBuffer size (in Bytes)
  const std::size_t maxTotalSize = (ChannelStream::BITS/8) * digis.size() * 2 + (kSlinkHeaderSize+kSlinkTrailerSize) * BTLElectronicsSpecs::kNumberOfFEDs;   /// 128 bit (16Bytes) x number of BTL digis x 2  + number of FEDs * size (header + trailer)
  auto rawDataBuffer = std::make_unique<RawDataBuffer>(maxTotalSize);
  LogDebug("BTLDigiToRaw") << "BTL number of digis = " << digis.size() << "\n";
  LogDebug("BTLDigiToRaw") << "Max RawDataBuffer pre-allocated size = " << int(maxTotalSize)<< " Bytes \n";
  
  // -- Linear scan: the digi collection is sorted by BTLDetId, so crystals
  // belonging to the same FED are contiguous. We accumulate channel
  // streams into a single buffer and flush it as soon as the FED id
  // changes
  int currentFed = -1;
  
  for (const auto& digi : digis) {
    BTLDetId detId(digi.krawId());

    // -- Forward lookup: DetId -> electronics
    BTLElectronicsIdPair elecIds = readoutMap.getElectronicsId(detId);
    
    // -- Sanity check: both sides must belong to the same FED and must have same hs-link and e-link Ids
    if (elecIds.minus.fedId() != elecIds.plus.fedId() ||
	elecIds.minus.hsLinkId() != elecIds.plus.hsLinkId() ||
	elecIds.minus.eLinkId() != elecIds.plus.eLinkId() )
      {
	edm::LogError("BTLDigiToRaw") << "BTLDigiToRaw::produce(): "
				      << "minus and plus sides of crystal " << std::hex << detId.rawId() << std::dec
				      << " belong to different FEDs/HS-link/E-link: " << "\n"
				      << " FED: " << elecIds.minus.fedId() << " vs " << elecIds.plus.fedId()
				      << " HL-link: " << elecIds.minus.hsLinkId() << " vs " << elecIds.plus.hsLinkId()
				      << " e-link: " << elecIds.minus.eLinkId() << " vs " << elecIds.plus.eLinkId()
				      << ") -- skipping crystal.";
	continue;
      }
    
    const int fedId = elecIds.minus.fedId();
    const int hslinkId = elecIds.minus.hsLinkId();
    const int elinkId = elecIds.minus.eLinkId();

    // Check if FED id changed and flush the buffer accumulated so far and start a new one.
    if (fedId != currentFed) {
      if (currentFed >= 0) {
	fillFEDBuffer(currentFed, eventId, currentChannelsStream_, *rawDataBuffer); 
      }
      currentChannelsStream_.clear();
      currentFed = fedId;
    }

    // Encode minus and plus side
    currentChannelsStream_.push_back(encodeChannelPayload(fedId, hslinkId, elinkId, digi, false));
    currentChannelsStream_.push_back(encodeChannelPayload(fedId, hslinkId, elinkId, digi, true));
  }

  // -- Flush the last accumulated FED buffer
  if (currentFed >= 0 && !currentChannelsStream_.empty()) {
    fillFEDBuffer(currentFed, eventId, currentChannelsStream_, *rawDataBuffer);
  }
  
  iEvent.put(std::move(rawDataBuffer));
}

// ------------------------------------------------------------
// encodeChannelPayload
// packs all fields (for the requested side) into a
// 128-bit ChannelStream.
// ------------------------------------------------------------
btldigitoraw::ChannelStream BTLDigiToRaw::encodeChannelPayload(int fed, int hslink, int elink, const btldigi::BTLDigi& digi, bool isPlusSide) const {

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

  int slink = fed - BTLElectronicsSpecs::kFirstFEDId;
  
  stream.set_bits(118, 10, static_cast<uint64_t>(BC0count));
  stream.set_bits(117,  1, static_cast<uint64_t>(status));
  stream.set_bits(110,  7, static_cast<uint64_t>(slink)); // will be removed from the channel payload and put only in the Slink header ?  
  stream.set_bits(104,  6, static_cast<uint64_t>(hslink)); // 6-bits !?!?!?!?!? should be 7 bits to cover 4-75 hs-link values range?
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

  LogDebug("BTLDigiToRaw")
    << "  BC0count="  << static_cast<int>(BC0count)
    << "  status="    << static_cast<int>(status)
    << "  slink="     << static_cast<int>(slink)
    << "  hslink="    << static_cast<int>(hslink)
    << "  elink="     << static_cast<int>(elink)
    << "  chID="      << static_cast<int>(chID)
    << "  BCcount="   << static_cast<int>(BCcount)
    << "  T1coarse="  << static_cast<int>(T1coarse)
    << "  T2coarse="  << static_cast<int>(T2coarse)
    << "  EOIcoarse=" << static_cast<int>(EOIcoarse)
    << "  Charge="    << static_cast<int>(Charge)
    << "  T1fine="    << static_cast<int>(T1fine)
    << "  T2fine="    << static_cast<int>(T2fine)
    << "  IdleTime="  << static_cast<int>(IdleTime)
    << "  PrevTrigF=" << static_cast<int>(PrevTrigF)
    << "  TACID="     << static_cast<int>(TACID)
    << "\n";
  
  return stream;
}

// ------------------------------------------------------------
// Fill FED buffer
// Writes:
//   - S-Link header  (128 bit)
//   - Channel payloads (128 bit each)
//   - S-Link trailer (128 bit)
// into the RawDataBuffer.
// ------------------------------------------------------------
void BTLDigiToRaw::fillFEDBuffer(int fedId, uint64_t eventId, const std::vector<btldigitoraw::ChannelStream>& channelsStream, RawDataBuffer& rawDataBuffer) const {

  using namespace btldigitoraw;
  
  const std::size_t payloadBytes  = channelsStream.size() * sizeof(Word) * kWordsPerChannel;
  const std::size_t fragmentBytes = kSlinkHeaderSize + payloadBytes + kSlinkTrailerSize;

  LogDebug("BTLDigiToRaw") << " FED id = " << fedId 
			   << " event id = " << eventId
			   << " kSlinkHeaderSize = " << kSlinkHeaderSize
			   << " kSlinkTrailerSize = " << kSlinkTrailerSize
			   << " payloadBytes = " << payloadBytes
			   << " fragmentBytes = " << fragmentBytes;
  
  std::vector<unsigned char> buffer(fragmentBytes, 0);

  // --------------------------------------------------------
  // S-Link header (128 bit)
  // --------------------------------------------------------
  // -- provisional values for emu_status, l1a_types, l1a_phys taken as in https://github.com/cms-sw/cmssw/blob/master/DataFormats/FEDRawData/test/TestWriteRawDataBuffer.cc
  uint8_t emu_status = 2;  //set 2 indicating fragment generated by DTH (emulator)
  uint16_t l1a_types = 1;  //set provisionally to 1, to be revised later
  uint8_t l1a_phys = 0;
  auto* slinkHeader = new ((void*)buffer.data()) SLinkRocketHeader_v3(static_cast<uint32_t>(fedId), l1a_types, l1a_phys, emu_status, eventId);
  assert(slinkHeader->verifyMarker());

  // --------------------------------------------------------
  // Channel payload
  // --------------------------------------------------------
  std::vector<Word> payloadWords;
  payloadWords.reserve(channelsStream.size() * kWordsPerChannel);

  for (const auto& channel : channelsStream) {
    for (std::size_t i = 0; i < kWordsPerChannel; ++i) {
      payloadWords.push_back(channel.raw_data()[i]);
    }
  }
  std::memcpy(buffer.data() + kSlinkHeaderSize, payloadWords.data(), payloadWords.size() * sizeof(Word));
    
  // --------------------------------------------------------
  // S-Link trailer (128 bit)
  // --------------------------------------------------------
  // provisional values
  uint16_t slt_status = 0;  
  uint16_t crc = 0;
  uint32_t orbitId = 0x3989;
  uint16_t bxId = 2200;
  uint16_t daq_crc = 0;
  auto* slinkTrailer =
    new ((void*)(buffer.data() + fragmentBytes - kSlinkTrailerSize)) SLinkRocketTrailer_v3(slt_status,
											   crc,
											   orbitId,
											   bxId,
											   static_cast<uint32_t>(fragmentBytes >> SLR_WORD_NUM_BYTES_SHIFT),
											   daq_crc);
  assert(slinkTrailer->verifyMarker());

  
  // --------------------------------------------------------
  // Store FED fragment
  // --------------------------------------------------------
  rawDataBuffer.addSource(static_cast<uint32_t>(fedId), buffer.data(), static_cast<uint32_t>(fragmentBytes));
}



// ------------------------------------------------------------
// fillDescriptions
// ------------------------------------------------------------
void BTLDigiToRaw::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("btlDigiCollection", edm::InputTag("mix", "MTDBarrel"));
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(BTLDigiToRaw);
