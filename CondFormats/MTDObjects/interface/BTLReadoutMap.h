#ifndef CondFormats_MTDObjects_BTLReadoutMap_h
#define CondFormats_MTDObjects_BTLReadoutMap_h

#include <cstdint>
#include <unordered_map>

#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "CondFormats/MTDObjects/interface/BTLElectronicsId.h"

// ------------------------------------------------------------
// Readout map: BTLDetId <-> BTLElectronicsId
// ------------------------------------------------------------
class BTLReadoutMap {
public:
  BTLReadoutMap();
  virtual ~BTLReadoutMap();

  // ----------------------------
  // Fill interface - inserts a new record in the readout map
  // ----------------------------
  void add(const BTLDetId& detId, const std::array<BTLElectronicsId, 2>& elecId);

  // ----------------------------
  // Forward lookup: DetId -> electronics
  // ----------------------------
  std::array<BTLElectronicsId, 2> getElectronicsId(const BTLDetId& detId) const;

  // ----------------------------
  // Reverse lookup: electronics -> DetId
  // ----------------------------
  BTLDetId getDetId(const BTLElectronicsId& elecId) const;

  // ----------------------------
  // Utilities
  // ----------------------------
  void clear();

  int size() const { return detToElec_.size(); };

private:
  // forward mapping
  std::unordered_map<uint32_t, std::array<BTLElectronicsId, 2>> detToElec_;

  // reverse mapping (packed electronics key -> detid)
  std::unordered_map<uint32_t, uint32_t> elecToDet_;
};

#endif
