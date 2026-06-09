#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CondFormats/MTDObjects/interface/BTLReadoutMap.h"

#include <stdexcept>
#include <ostream>

BTLReadoutMap::BTLReadoutMap() {}

BTLReadoutMap::~BTLReadoutMap() {}

// ------------------------------------------------------------
// Add full 2-channel mapping
// ------------------------------------------------------------
void BTLReadoutMap::add(const BTLDetId& detId, const std::array<BTLElectronicsId, 2>& elecIds) {
  const uint32_t detKey = detId.rawId();

  // ----------------------------
  // Forward: det -> 2 channels
  // ----------------------------
  if (detToElec_.find(detKey) != detToElec_.end()) {
    throw cms::Exception("BTLReadoutMap::add() - duplicate BTLDetId entry") << std::endl;
  }

  detToElec_[detKey] = elecIds;

  // ----------------------------
  // Reverse: each channel -> det
  // ----------------------------
  for (const auto& id : elecIds) {
    const uint32_t elecKey = id.rawId();
    auto it = elecToDet_.find(elecKey);
    if (it != elecToDet_.end()) {
      throw cms::Exception("BTLReadoutMap::add() - duplicate BTLElectronicsId entry") << std::endl;
    }

    elecToDet_[elecKey] = detKey;
  }
}

// ------------------------------------------------------------
// Forward lookup
// ------------------------------------------------------------
std::array<BTLElectronicsId, 2> BTLReadoutMap::getElectronicsId(const BTLDetId& detId) const {
  const auto it = detToElec_.find(detId.rawId());

  if (it == detToElec_.end()) {
    edm::LogWarning("BTLReadoutMap") << "BTLReadoutMap::getElectronicsId(): "
                                     << "******************  BTLDetId not found! ";  // log warning or exception?
    return {};
  }

  return it->second;
}

// ------------------------------------------------------------
// Reverse lookup
// ------------------------------------------------------------
BTLDetId BTLReadoutMap::getDetId(const BTLElectronicsId& elecId) const {
  const auto it = elecToDet_.find(elecId.rawId());

  if (it == elecToDet_.end()) {
    edm::LogWarning("BTLReadoutMap") << "BTLReadoutMap::getDetId(): "
                                     << "******************  BTLElectronicsId not found! ";  // log warning or exception?
    return BTLDetId();
  }

  return BTLDetId(it->second);
}

// ------------------------------------------------------------
// clear
// ------------------------------------------------------------
void BTLReadoutMap::clear() {
  detToElec_.clear();
  elecToDet_.clear();
}
