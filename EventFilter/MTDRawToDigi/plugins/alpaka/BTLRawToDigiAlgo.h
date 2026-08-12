#ifndef EventFilter_MTDRawToDigi_plugins_alpaka_BTLRawToDigiAlgo_h
#define EventFilter_MTDRawToDigi_plugins_alpaka_BTLRawToDigiAlgo_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

#include "DataFormats/FEDRawData/interface/RawDataBuffer.h"

#include "DataFormats/FTLDigiSoA/interface/alpaka/BTLDigiDeviceCollection.h"
#include "CondFormats/MTDObjects/interface/alpaka/BTLElectronicsToDetIdMappingDevice.h"

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

// ----------------------------------------------------------------
// Intermediate SoA produced by the decoding kernel.
// Each row corresponds to one 128-bit electronics channel payload.
// This is not yet a BTLDigi: two channel payloads belonging to the
// same ASIC will later be paired to build one BTLDigi entry.

// Keep local to this module rather than in DataFormats/: it's a purely as it is a pure auxiliary format.
// See also: https://github.com/cms-sw/cmssw/blob/master/RecoLocalTracker/SiStripClusterizer/plugins/alpaka/SiStripRawToClusterAlgo.h

GENERATE_SOA_LAYOUT(BTLChannelPayloadSoALayout,
                    SOA_COLUMN(uint16_t, bc0count),
                    SOA_COLUMN(bool, status),
                    SOA_COLUMN(uint32_t, bcCount),
                    SOA_COLUMN(uint8_t, chId),  // 0..31, as read from the payload
                    SOA_COLUMN(uint16_t, t1Coarse),
                    SOA_COLUMN(uint16_t, t2Coarse),
                    SOA_COLUMN(uint16_t, eoiCoarse),
                    SOA_COLUMN(uint16_t, charge),
                    SOA_COLUMN(uint16_t, t1Fine),
                    SOA_COLUMN(uint16_t, t2Fine),
                    SOA_COLUMN(uint16_t, idleTime),
                    SOA_COLUMN(uint8_t, prevTrigF),
                    SOA_COLUMN(uint8_t, tacId),

                    SOA_COLUMN(uint32_t, fedId),
                    SOA_COLUMN(uint8_t, hsLinkId),
                    SOA_COLUMN(uint8_t, eLinkId),

                    // Dense key encoding (fedId, hsLinkId, eLinkId).
                    // Used to group channels belonging to the same ASIC before pairing.
                    SOA_COLUMN(uint32_t, chipKey))

using BTLChannelPayloadSoA = BTLChannelPayloadSoALayout<>;
using BTLChannelPayloadHostCollection = PortableHostCollection<BTLChannelPayloadSoA>;

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using BTLChannelPayloadDeviceCollection = PortableCollection<BTLChannelPayloadSoA>;
}

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(BTLChannelPayloadDeviceCollection, BTLChannelPayloadHostCollection);
// ----------------------------------------------------------------


// ----------------------------------------------------------------
// *** Device pipeline: decode -> bucket by chip -> pair -> fill digis.   
//
namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using btldigi::BTLDigiDeviceCollection;

  // Dense per-chip/per-channel lookup: channelIndexByChip[chip * kChannelsPerChip + chId] holds the row
  // index (into BTLChannelPayloadSoA) of that channel, or -1 if that channel wasn't present in this event (sentinel).
  static constexpr int32_t kChannelsPerChip = 32;

  using BTLChannelIndexTable = cms::alpakatools::device_buffer<Device, int32_t[]>;

  //struct BTLChannelIndexTable {
  //  cms::alpakatools::device_buffer<Device, int32_t[]> channelIndexByChip;  // size nChips*kChannelsPerChip
  //  int32_t nChips = 0;
  //};
  

  class BTLRawToDigiAlgo {
  public:
    BTLChannelPayloadDeviceCollection decodeOnly(Queue& queue,
                                                 const uint64_t* rawWords_h,
                                                 const int32_t* channelFedId_h,
                                                 int32_t nChannels,
                                                 BTLElectronicsIndexer const& indexer) const;

    BTLChannelIndexTable buildChannelIndexTable(Queue& queue,
						BTLChannelPayloadDeviceCollection const& channelPayload_d,
						BTLElectronicsIndexer const& indexer) const;
    
    BTLDigiDeviceCollection pairChannels(Queue& queue,
					 BTLChannelPayloadDeviceCollection const& channelsPayload_d,
					 BTLChannelIndexTable const& table,
					 BTLElectronicsIndexer const& indexer,
					 BTLElectronicsToDetIdMappingDevice const& elecToDetId) const;
    
    BTLDigiDeviceCollection process(Queue& queue,
                                    const uint64_t* rawWords_h,
                                    const int32_t* channelFedId_h,
                                    int32_t nChannels,
                                    BTLElectronicsIndexer const& indexer,
                                    BTLElectronicsToDetIdMappingDevice const& elecToDetId) const;
												
  private:
    void debugChannelIndexTable(Queue& queue,
				BTLChannelIndexTable const& table,
				BTLChannelPayloadDeviceCollection const& channelsPayload_d,
				BTLElectronicsIndexer const& indexer) const;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
