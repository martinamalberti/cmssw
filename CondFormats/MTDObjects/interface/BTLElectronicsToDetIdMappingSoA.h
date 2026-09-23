#ifndef CondFormats_MTDObjects_BTLElectronicsToDetIdMappingSoA_h
#define CondFormats_MTDObjects_BTLElectronicsToDetIdMappingSoA_h

//#include <cstdint>
#include <alpaka/alpaka.hpp>
#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

// -----------------------------------------------------------------------
// Dense electronics -> crystal rawId lookup
// Flat array used by BTL unpacker on GPU to go from (fedId, hsLinkId, eLinkId, channel id) to the crystal's BTLDetId::rawId().
// -----------------------------------------------------------------------
GENERATE_SOA_LAYOUT(BTLElectronicsToDetIdMappingSoALayout,
		    SOA_COLUMN(uint32_t, rawId),
		    SOA_COLUMN(bool, valid),
		    SOA_COLUMN(uint8_t, partnerChannelId),
		    SOA_COLUMN(int8_t, side),
		    SOA_COLUMN(int8_t, pairId))
		    //SOA_COLUMN(uint8_t, minusChannelId),
                    //OA_COLUMN(uint8_t, plusxxsChannelId))

using BTLElectronicsToDetIdMappingSoA = BTLElectronicsToDetIdMappingSoALayout<>;

// Flat-index helper shared between the ESProducer (host, fills the table)
// and the pairing kernel (device, reads it).
struct BTLElectronicsIndexer {
  int32_t firstFedId;
  int32_t hsLinkOffset;
  int32_t nFeds;
  int32_t nHsLinks;
  int32_t nELinks;
  //  static constexpr int32_t nPairs = 16;  // canali o coppie ?????
  static constexpr int32_t nChannelsPerChip = 32;  // number of channels in one chip 
  static constexpr int32_t nPairsPerChip = 16;  // max number of pairs within one chip
  
  //ALPAKA_FN_ACC inline int32_t flatIndex(int32_t fedId, int32_t hs, int32_t el, int32_t pairIdx) const {
  ALPAKA_FN_ACC inline int32_t flatIndex(int32_t fedId, int32_t hs, int32_t el, int32_t channel) const {
    int32_t f = fedId - firstFedId;
    int32_t h = hs - hsLinkOffset;
    return ((f * nHsLinks + h) * nELinks + el) * nChannelsPerChip + channel;
    //return ((f * nHsLinks + h) * nELinks + el) * nPairs + pairIdx;
  }
  //inline int32_t size() const { return nFeds * nHsLinks * nELinks * nPairs; }
  inline int32_t size() const { return nFeds * nHsLinks * nELinks * nChannelsPerChip; }
};

#endif
