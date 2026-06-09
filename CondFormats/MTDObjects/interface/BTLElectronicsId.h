#ifndef CondFormats_MTDObjects_BTLElectronicsId_h
#define CondFormats_MTDObjects_BTLElectronicsId_h

#include <cstdint>
#include <iosfwd>

class BTLElectronicsId {
public:
  // Bit layout:
  //
  //  0 -  4 : channelId 0..31  (5 bits)
  //  5 -  9 : eLinkId   0..23  (5 bits)
  // 10 - 16 : hsLinkId  0..127 (7 bits)
  // 17 - 29 : fedId     ?      (13 bits)
  // 30 - 31 : free for now

  static constexpr uint32_t kChannelMask = 0x1F;  // 5 bits
  static constexpr uint32_t kELinkMask = 0x1F;    // 5 bits
  static constexpr uint32_t kHSLinkMask = 0x7F;   // 7 bits
  static constexpr uint32_t kFEDMask = 0x1FFF;    // 13 bits

  static constexpr unsigned kChannelShift = 0;
  static constexpr unsigned kELinkShift = 5;
  static constexpr unsigned kHSLinkShift = 10;
  static constexpr unsigned kFEDShift = 17;

  /** Default constructor **/
  BTLElectronicsId();

  /** Constructor from rawId **/
  explicit BTLElectronicsId(uint32_t rawid);

  /** Constructor from (FED id, HS-link id, e-link id, tofhir channel id) **/
  BTLElectronicsId(uint16_t fedId,  // sLinkId
                   uint8_t hsLinkId,
                   uint8_t eLinkId,
                   uint8_t channelId);

  /// Accessors
  int fedId() const;
  int hsLinkId() const;
  int eLinkId() const;
  int channelId() const;

  uint32_t rawId() const;

  bool operator==(const BTLElectronicsId&) const;
  bool operator!=(const BTLElectronicsId&) const;

private:
  uint32_t rawid_ = 0;

  //COND_SERIALIZABLE;
};

std::ostream& operator<<(std::ostream&, const BTLElectronicsId&);

#endif
