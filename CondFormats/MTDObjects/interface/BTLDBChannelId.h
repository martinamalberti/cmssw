#ifndef CondFormats_MTDObjects_BTLDBChannelId_h
#define CondFormats_MTDObjects_BTLDBChannelId_h

#include "CondFormats/Serialization/interface/Serializable.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/ForwardDetId/interface/MTDDetId.h"

#include <cstdint>
#include <iosfwd>

class BTLDBChannelId {
public:
  // Bit layout:
  //
  //  0 - 24 : bits 24:0 of the BTLDetId  // a meno che non ci serva "liberare" piu' bit ...
  // 25 - 29 : TOFHIR channel identifier [0-31]
  // 30 - 31 : unused
  //

  static constexpr uint32_t kBtlDetIdMask = 0x01FFFFFF;  // 25 bits
  static constexpr uint32_t kChannelMask  = 0x1F;        // 5 bits

  static constexpr unsigned kBtlDetIdShift = 0;
  static constexpr unsigned kChannelShift  = 25;

  static constexpr uint32_t kMTDPrefix = (DetId::Forward << 28) | (MTDDetId::FastTime << 25);
			   
  /** Default constructor **/
  BTLDBChannelId() ;

  /** Constructor from packed rawId **/
  explicit BTLDBChannelId(uint32_t rawid);

  /** Constructor from (BTLDetId, TOFHIR channel) **/
  BTLDBChannelId(BTLDetId detid, uint8_t chId) {
    rawid_ =
      (detid.rawId() & kBtlDetIdMask) |
      ((static_cast<uint32_t>(chId) & kChannelMask) << kChannelShift);
  }

  /// Accessors
  BTLDetId detId() const { return BTLDetId((rawid_ & kBtlDetIdMask) | kMTDPrefix); }

  int channelId() const {
    return (rawid_ >> kChannelShift) & kChannelMask;
  }

  uint32_t rawId() const { return rawid_; }

  bool operator==(const BTLDBChannelId&) const;
  bool operator!=(const BTLDBChannelId&) const;

private:
  uint32_t rawid_ = 0;

  COND_SERIALIZABLE;
};

#endif
