#ifndef CondFormats_MTDObjects_BTLElectronicsToDetIdMappingSoA_h
#define CondFormats_MTDObjects_BTLElectronicsToDetIdMappingSoA_h

#include <alpaka/alpaka.hpp>
#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

// -----------------------------------------------------------------------
// Dense electronics -> crystal rawId lookup
// Flat array used by BTL unpacker on GPU to go from
// (fedId, hsLinkId, eLinkId, channelId) to the crystal's BTLDetId::rawId().
// -----------------------------------------------------------------------

GENERATE_SOA_LAYOUT(BTLElectronicsToDetIdMappingSoALayout,
                    SOA_COLUMN(uint32_t, rawId),
                    SOA_COLUMN(bool, valid),
                    SOA_COLUMN(uint8_t, partnerChannelId),
                    SOA_COLUMN(int8_t, side))

using BTLElectronicsToDetIdMappingSoA = BTLElectronicsToDetIdMappingSoALayout<>;

// The table is indexed by the dense channel index
// BTLElectronicsIndexing::denseChannelIndex(fedId, hsLinkId, eLinkId, channelId).

#endif
