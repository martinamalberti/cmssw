#include "CondFormats/MTDObjects/interface/BTLDBChannelId.h"

#include <ostream>

// ------------------------------------------------------------
// Constructors
// ------------------------------------------------------------
BTLDBChannelId::BTLDBChannelId() : rawid_(0) {}

BTLDBChannelId::BTLDBChannelId(uint32_t rawid) : rawid_(rawid) {}

// ------------------------------------------------------------
// Comparison operators
// ------------------------------------------------------------

bool BTLDBChannelId::operator==(const BTLDBChannelId& other) const { return rawid_ == other.rawid_; }

bool BTLDBChannelId::operator!=(const BTLDBChannelId& other) const { return rawid_ != other.rawid_; }

// ------------------------------------------------------------
// Stream operator
// ------------------------------------------------------------

std::ostream& operator<<(std::ostream& os, const BTLDBChannelId& id) {
  os << "BTLDBChannelId: "
     << "rawId = " << static_cast<int>(id.rawId()) << ", BTLDetId = " << id.detId() << ", TOFHIR channel Id = " << id.channelId();

  return os;
}
