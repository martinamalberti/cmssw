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

// Flat-index helper shared between the ESProducer (host, fills the table)
// and the pairing kernel (device, reads it).
//
// BTL electronics hierarchy:
//   12 Serenity boards in the detector
//   Each Serenity drives 12 S-links = 12 FED IDs per Serenity --> 144 FED IDs total
//   Each Serenity has 36 HS-links total (12 S-links with 3 HS-links/S-link)
//   Each HS-link (1 RU) carries 24 elinks
//   Each elink = 1 TOFHIR2C chip = 32 channels
struct BTLElectronicsIndexer {
  int32_t firstFedId;        // first FED ID in the BTL (inclusive)
  int32_t nFeds;             // total number of BTL FED IDs (= 144)
  int32_t nFEDsPerSerenity;  // 12 Slink/FED ids per each serenity
  int32_t nHSLinksPerSerenity;
  int32_t nELinksPerHSLink;
  int32_t nChannelsPerELink;

  // Aliases kept for use in the kernel code
  static constexpr int32_t nChannelsPerChip = 32;
  static constexpr int32_t nPairsPerChip = 16;  // max pairs within one chip

  // Returns the flat array index for the entry (fedId, hs, el, channel).
  //
  // Parameters:
  //   fedId   raw FED ID from the DAQ fragment
  //   hs      absolute HS-link id within the Serenity, as encoded in the payload (range 0-35);
  //   el      elink id within the HS-link (0-23)
  //   channel channel id within the elink/chip (0-31)
  ALPAKA_FN_ACC inline int32_t flatIndex(int32_t fedId, int32_t hs, int32_t el, int32_t channel) const {
    int32_t serenityIdx = (fedId - firstFedId) / nFEDsPerSerenity;  // 0-11
    return ((serenityIdx * nHSLinksPerSerenity + hs) * nELinksPerHSLink + el) * nChannelsPerELink + channel;
  }

  // Checks the (hs, el) ids read from the payload before using them in flatIndex():
  //   - hs < nHSLinksPerSerenity and el < nELinksPerHSLink: otherwise flatIndex() points to another chip
  //     (or outside the mapping table), since hs has 6 bits (0-63) and el has 5 bits (0-31) in the payload;
  //   - hs belongs to this FED: FED k of a Serenity carries the HS-links [k * nHSLinksPerFED, (k + 1) * nHSLinksPerFED).
  //     The ESProducer of the mapping checks that the readout map respects this rule.
  // fedId is assumed in [firstFedId, firstFedId + nFeds), since it comes from the FED loop on the host.
  ALPAKA_FN_HOST_ACC inline bool isValidLink(int32_t fedId, int32_t hs, int32_t el) const {
    const int32_t nHSLinksPerFED = nHSLinksPerSerenity / nFEDsPerSerenity;
    const int32_t fedInSerenity = (fedId - firstFedId) % nFEDsPerSerenity;
    return hs >= 0 && hs < nHSLinksPerSerenity && el >= 0 && el < nELinksPerHSLink &&
           hs / nHSLinksPerFED == fedInSerenity;
  }

  // Total number of entries in the mapping SoA.
  inline int32_t size() const {
    return (nFeds / nFEDsPerSerenity) * nHSLinksPerSerenity * nELinksPerHSLink * nChannelsPerELink;
  }
};

#endif
