#include <cstring>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/FEDRawData/interface/RawDataBuffer.h"
#include "DataFormats/FEDRawData/interface/SLinkRocketHeaders.h"

#include "EventFilter/MTDRawToDigi/interface/BitStream.h"

namespace btldigitoraw {
  using ChannelStream = BitStream<2>;
  constexpr std::size_t kWordsPerChannel = ChannelStream::BITS / ChannelStream::WORD_BITS;
}  // namespace btldigitoraw

class DumpBTLRawDataBuffer : public edm::global::EDAnalyzer<> {
public:
  explicit DumpBTLRawDataBuffer(edm::ParameterSet const&);
  void analyze(edm::StreamID, edm::Event const&, edm::EventSetup const&) const override;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  const edm::EDGetTokenT<RawDataBuffer> rawDataBufferToken_;
};

// ------------------------------------------------------------
// Constructor
// ------------------------------------------------------------
DumpBTLRawDataBuffer::DumpBTLRawDataBuffer(edm::ParameterSet const& iPSet)
    : rawDataBufferToken_(consumes<RawDataBuffer>(iPSet.getParameter<edm::InputTag>("rawDataBufferTag"))) {}

// ------------------------------------------------------------
// analyze
// ------------------------------------------------------------
void DumpBTLRawDataBuffer::analyze(edm::StreamID,
                                   edm::Event const& iEvent,
                                   edm::EventSetup const&) const {
  const auto& rawDataBuffer = iEvent.get(rawDataBufferToken_);

  constexpr std::size_t kHeaderSize  = sizeof(SLinkRocketHeader_v3);
  constexpr std::size_t kTrailerSize = sizeof(SLinkRocketTrailer_v3);
  constexpr std::size_t kWordBytes   = 16;  // 128-bit words

  edm::LogInfo("DumpBTLRawDataBuffer")
      << "=== Event " << iEvent.id().event()
      << " -- number of FED fragments: " << rawDataBuffer.map().size();

  // -- Iterate over all source IDs actually present in the buffer
  for (const auto& entry : rawDataBuffer.map()) {
    const uint32_t sourceId = entry.first;
    const auto& fragData = rawDataBuffer.fragmentData(sourceId);

    if (!fragData.isValid() || fragData.size() == 0) {
      edm::LogWarning("DumpBTLRawDataBuffer") << "No data for source ID = " << sourceId;
      continue;
    }

    edm::LogInfo("DumpBTLRawDataBuffer") << "--- SOURCE ID = " << sourceId << "  fragment size = " << fragData.size() << " bytes";

    // --------------------------------------------------------
    // S-Link header
    // --------------------------------------------------------
    if (fragData.size() < kHeaderSize) {
      edm::LogWarning("DumpBTLRawDataBuffer") << "Fragment too small for S-Link header: " << fragData.size() << " bytes";
      continue;
    }

    const auto headerSpan = fragData.dataHeader(kHeaderSize);
    auto headerView = makeSLinkRocketHeaderView(headerSpan);

    edm::LogInfo("DumpBTLRawDataBuffer")
        << "  [Header]"
        << "  version="    << static_cast<int>(headerView->version())
        << "  sourceID="   << headerView->sourceID()
        << "  eventID="    << headerView->globalEventID()
        << "  l1aTypes="   << headerView->l1aTypes()
        << "  l1aPhys="    << static_cast<int>(headerView->l1aPhysType())
        << "  emuStatus="  << static_cast<int>(headerView->emuStatus())
        << "  markerOK="   << headerView->verifyMarker();

    // --------------------------------------------------------
    // Payload (channel words, 128 bit each)
    // --------------------------------------------------------
    if (fragData.size() < kHeaderSize + kTrailerSize) {
      edm::LogWarning("DumpBTLRawDataBuffer")
          << "Fragment too small for header + trailer: " << fragData.size() << " bytes";
      continue;
    }

    const auto payloadSpan = fragData.payload(kHeaderSize, kTrailerSize); // return payload between header and trailer
    const std::size_t nWords = payloadSpan.size() / kWordBytes;

    edm::LogInfo("DumpBTLRawDataBuffer")
        << "  [Payload]  " << nWords << " channel word(s) of 128 bit each";

    auto printStream = [&](const BitStream<2>& bs) {
			 std::stringstream ss;
			 for (int index = static_cast<int>(bs.size())-1; index >= 0; --index) {
			   ss << bs.get_bit(static_cast<size_t>(index));
			 }
			 return ss.str();
		       };
    
     
    for (std::size_t iWord = 0; iWord < nWords; ++iWord) {
      const unsigned char* wordPtr = payloadSpan.data() + iWord * kWordBytes;

      // Copy the 16 bytes of this channel word into a BitStream<2>
      // and extract each field by name, mirroring the layout defined
      // in BTLDigiToRaw::encodeChannelPayload().
      btldigitoraw::ChannelStream stream;
      std::memcpy(stream.raw_data(), wordPtr, kWordBytes);
   
      const uint64_t BC0count  = stream.get_bits(118, 10);
      const uint64_t status    = stream.get_bits(117,  1);
      const uint64_t slink     = stream.get_bits(110,  7);
      const uint64_t hslink    = stream.get_bits(104,  6);
      const uint64_t elink     = stream.get_bits( 99,  5);
      const uint64_t chID      = stream.get_bits( 94,  5);
      const uint64_t BCcount   = stream.get_bits( 82, 12);
      const uint64_t T1coarse  = stream.get_bits( 67, 15);
      const uint64_t T2coarse  = stream.get_bits( 57, 10);
      const uint64_t EOIcoarse = stream.get_bits( 47, 10);
      const uint64_t Charge    = stream.get_bits( 37, 10);
      const uint64_t T1fine    = stream.get_bits( 27, 10);
      const uint64_t T2fine    = stream.get_bits( 17, 10);
      const uint64_t IdleTime  = stream.get_bits(  7, 10);
      const uint64_t PrevTrigF = stream.get_bits(  3,  4);
      const uint64_t TACID     = stream.get_bits(  0,  3);

      edm::LogInfo("DumpBTLRawDataBuffer")
	<< "  word " << iWord << "\n" 
	<< "  bit content: " << printStream(stream) << "\n"
	<< "  BC0count="  << BC0count
	<< "  status="    << status
	<< "  slink="     << slink
	<< "  hslink="    << hslink
	<< "  elink="     << elink
	<< "  chID="      << chID
	<< "  BCcount="   << BCcount
	<< "  T1coarse="  << T1coarse
	<< "  T2coarse="  << T2coarse
	<< "  EOIcoarse=" << EOIcoarse
	<< "  Charge="    << Charge
	<< "  T1fine="    << T1fine
	<< "  T2fine="    << T2fine
	<< "  IdleTime="  << IdleTime
	<< "  PrevTrigF=" << PrevTrigF
	<< "  TACID="     << TACID;
    }

    // --------------------------------------------------------
    // S-Link trailer
    // --------------------------------------------------------
    const auto trailerSpan = fragData.dataTrailer(kTrailerSize);
    auto trailerView = makeSLinkRocketTrailerView(trailerSpan, headerView->version());

    edm::LogInfo("DumpBTLRawDataBuffer")
        << "  [Trailer]"
        << "  orbitID="       << trailerView->orbitID()
        << "  bxID="          << trailerView->bxID()
        << "  eventLenBytes=" << trailerView->eventLenBytes()
        << "  status="        << trailerView->status()
        << "  crc=0x"         << std::hex << trailerView->crc()    << std::dec
        << "  daqCRC=0x"      << std::hex << trailerView->daqCRC() << std::dec
        << "  markerOK="      << trailerView->verifyMarker();
  }
}

// ------------------------------------------------------------
// fillDescriptions
// ------------------------------------------------------------
void DumpBTLRawDataBuffer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("rawDataBufferTag");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(DumpBTLRawDataBuffer);
