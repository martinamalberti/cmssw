#ifndef EventFilter_MTDRawToDigi_plugins_alpaka_BTLRawToDigiAlgo_h
#define EventFilter_MTDRawToDigi_plugins_alpaka_BTLRawToDigiAlgo_h

#include "DataFormats/FTLDigi/interface/alpaka/BTLDigiSoACollection.h"
#include "CondFormats/MTDObjects/interface/alpaka/BTLElectronicsToDetIdSoACollection.h" // ANCORA DA IMPLEMENTARE
#include "CondFormats/MTDObjects/interface/BTLChannelMaps.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"


// ---
// Intermediate SoA produced by the decoding kernel.
// Each row corresponds to one 128-bit electronics channel payload.
// This is not yet a BTLDigi: two channel payloads belonging to the
// same ASIC will later be paired to build one BTLDigi entry.
namespace btlraw {
  GENERATE_SOA_LAYOUT(BTLChannelPayloadSoALayout,
		      SOA_COLUMN(uint16_t, bc0count),
		      SOA_COLUMN(bool, status),
		      SOA_COLUMN(uint32_t, bcCount),
		      SOA_COLUMN(uint8_t, chId),      // 0..31, as read from the payload
		      SOA_COLUMN(uint16_t, t1Coarse),
		      SOA_COLUMN(uint16_t, t2Coarse),
		      SOA_COLUMN(uint16_t, eoiCoarse),
		      SOA_COLUMN(uint16_t, charge),
		      SOA_COLUMN(uint16_t, t1Fine),
		      SOA_COLUMN(uint16_t, t2Fine),
		      SOA_COLUMN(uint16_t, idleTime),
		      SOA_COLUMN(uint8_t, prevTrigF),
		      SOA_COLUMN(uint8_t, tacId),
		      
		      SOA_COLUMN(int32_t, fedId),
		      SOA_COLUMN(int16_t, hsLinkId),
		      SOA_COLUMN(int16_t, eLinkId),

		      // Dense key encoding (fedId, hsLinkId, eLinkId).
		      // Used to group channels belonging to the same ASIC before pairing.
		      SOA_COLUMN(uint32_t, chipKey))

  using BTLChannelPayloadSoA = BTLChannelPayloadSoALayout<>;
  using BTLChannelPayloadHost = PortableHostCollection<BTLChannelPayloadSoA>;
} // namespace btlraw

namespace ALPAKA_ACCELERATOR_NAMESPACE::btlraw {
  using namespace ::btlraw;
  using BTLChannelPayloadDevice = PortableCollection<BTLChannelPayloadSoA>;
}  // namespace btlraw

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(btlraw::BTLChannelPayloadDevice, btlraw::BTLChannelPayloadHost);
// ---


namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // Device pipeline only: decode -> bucket by chip -> pair -> fill digis.
  class BTLRawToDigiAlgo {
  public:

    BTLDigiDeviceCollection process(Queue& queue,
				    const uint64_t* rawWords_h,
				    const int32_t* channelFedId_h,
				    int32_t nChannels,
				    BTLElectronicsIndexer const& indexer,
				    BTLElectronicsToDetIdDeviceCollection const& elecToDetId) const;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
