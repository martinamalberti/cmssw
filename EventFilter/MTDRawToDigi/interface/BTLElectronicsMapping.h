#ifndef EventFilter_MTDRawToDigi_BTLElectronicsMapping_h
#define EventFilter_MTDRawToDigi_BTLElectronicsMapping_h

#include <array>
#include <utility>

#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "EventFilter/MTDRawToDigi/interface/BTLElectronicsSpecs.h"

/**
 * Helper functions for BTL electronics/DAQ mapping
 */

class BTLElectronicsMapping {
public:
  /** Default constructor -- invalid value */
  BTLElectronicsMapping(const BTLDetId::CrysLayout lay);

  // ------------------------------------------------------------------
  // SiPM channel / TOFHIR channel mapping helper functions
  // ------------------------------------------------------------------
  /** \brief BTL TOFHIR channel mapping with crystal BTLDetId
      Convention:
      SiPMside 0 == Minus Side
      SiPMside 1 == Plus Side
  */

  // -- Get SiPM Channel number from crystal
  int SiPMCh(uint32_t smodCopy, uint32_t crystal, uint32_t SiPMSide);
  int SiPMCh(BTLDetId det, uint32_t SiPMSide);
  int SiPMCh(uint32_t rawID, uint32_t SiPMSide);

  // -- Get TOFHIR Channel number from crystal
  int TOFHIRCh(uint32_t smodCopy, uint32_t crystal, uint32_t SiPMSide);
  int TOFHIRCh(BTLDetId det, uint32_t SiPMSide);
  int TOFHIRCh(uint32_t rawID, uint32_t SiPMSide);

  // -- Returns TOFHIR ASIC number in construction database
  int TOFHIRASIC(uint32_t dmodule, uint32_t smodCopy);
  int TOFHIRASIC(BTLDetId det);
  int TOFHIRASIC(uint32_t rawID);

  // ------------------------------------------------------------------
  // E-link mapping helper functions
  // ------------------------------------------------------------------

  // -- Get the e-link for a given SM
  int elinkFromSM(uint32_t dmodule, uint32_t smodCopy, int lpgbt_id = 0);
  int elink(BTLDetId det, int lpgbt_id = 0);
  int elink(uint32_t rawID, int lpgbt_id = 0);

  // ------------------------------------------------------------------
  // HS-link mapping helper functions
  // ------------------------------------------------------------------
  // -- Get the HS-link corresponding to a given RU/CC in a tray
  int opticalTxPosition(uint32_t tray, int optTxCh);
  int hslinkFromRU(uint32_t runit, uint32_t tray, int lpgbt_id = 0);
  int hslink(BTLDetId det, int lpgbt_id = 0);
  int hslink(uint32_t rawID, int lpgbt_id = 0);

  // ------------------------------------------------------------------
  // S-link mapping
  // ------------------------------------------------------------------
  /** one S-link corresponds to a group of 6 trays. One S-link = One FEDId **/
  // -- Get the S-link for a given tray
  int slinkFromTray(uint32_t tray, uint32_t zside);
  int slink(BTLDetId det);
  int slink(uint32_t rawID);

private:
};

#endif
