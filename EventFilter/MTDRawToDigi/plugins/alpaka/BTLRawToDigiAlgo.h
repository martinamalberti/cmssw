#ifndef EventFilter_MTDRawToDigi_plugins_alpaka_BTLRawToDigiAlgo_h
#define EventFilter_MTDRawToDigi_plugins_alpaka_BTLRawToDigiAlgo_h

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

#include "DataFormats/FEDRawData/interface/RawDataBuffer.h"

#include "DataFormats/FTLDigiSoA/interface/alpaka/BTLDigiDeviceCollection.h"
#include "CondFormats/MTDObjects/interface/alpaka/BTLElectronicsToDetIdMappingDevice.h"
#include "EventFilter/MTDRawToDigi/interface/BTLElectronicsSpecs.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include <sstream>
#include <optional>

// ----------------------------------------------------------------
// Intermediate SoA produced by the decoding kernel.
// Each row corresponds to one 128-bit electronics channel payload.
// This is not yet a BTLDigi: two channel payloads belonging to the
// same ASIC will later be paired to build one BTLDigi entry.

// Keep local to this module rather than in DataFormats/: it's a pure auxiliary format.
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
// *** Device pipeline: decode per channel -> find segments corresponding to chip boundaries -> for each active chip: pairing and fill digis.
//
namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using btldigi::BTLDigiDeviceCollection;

  constexpr int kPlusSide = 1;
  constexpr int kMinusSide = 0;

  struct BTLChipSegments {
    int32_t* segmentStart = nullptr;  // segmentStart[k]/[k+1] define the segment k
    int32_t* nSegments = nullptr;     // numerber of segments found in the event
  };

  class BTLRawToDigiAlgo {
  public:
    BTLChannelPayloadDeviceCollection decodeChannels(Queue& queue,
                                                     const uint64_t* rawWords_h,
                                                     const int32_t* channelFedId_h,
                                                     int32_t nChannels,
                                                     BTLElectronicsIndexer const& indexer);

    BTLChipSegments findChipSegments(Queue& queue,
                                     BTLChannelPayloadDeviceCollection const& channelPayload_d,
                                     int32_t nChannels,
                                     BTLElectronicsIndexer const& indexer);

    void countDigis(Queue& queue,
                    BTLChannelPayloadDeviceCollection const& channelPayload_d,
                    BTLChipSegments const& segments,
                    BTLElectronicsToDetIdMappingDevice const& mapping,
                    BTLElectronicsIndexer const& indexer,
                    int32_t nActiveChips);

    // computes the per-segment digi offsets on device and enqueues the async copy of nDigisTot to the host
    void computeDigiOffsets(Queue& queue, int32_t nChips, cms::alpakatools::host_buffer<int32_t>& nDigisTot_h);

    BTLDigiDeviceCollection fillDigis(Queue& queue,
                                      BTLChannelPayloadDeviceCollection const& channelPayload_d,
                                      BTLChipSegments const& segments,
                                      BTLElectronicsToDetIdMappingDevice const& mapping,
                                      BTLElectronicsIndexer const& indexer,
                                      int32_t nActiveChips,
                                      int32_t nDigisTot);

    // ---- The processing is split in two parts, to be called from acquire() and produce() of a
    // SynchronizingEDProducer: the framework waits (non-blocking) for the queue between the two,
    // so that nDigisTot is available on the host without an explicit alpaka::wait.

    // STEP1 -> STEP2a.5: decode, find chip segments, pair and count digis, compute offsets.
    // Enqueues the async copy of the total number of digis into nDigisTot_h.
    void prepareDigis(Queue& queue,
                      const uint64_t* rawWords_h,
                      const int32_t* channelFedId_h,
                      int32_t nChannels,
                      BTLElectronicsToDetIdMappingDevice const& elecToDetId,
                      BTLElectronicsIndexer const& indexer,
                      cms::alpakatools::host_buffer<int32_t>& nDigisTot_h);

    // STEP2b: allocate the output collection with nDigisTot (read on the host after the
    // synchronization) and fill it. Must be called after prepareDigis().
    BTLDigiDeviceCollection produceDigis(Queue& queue,
                                         BTLElectronicsToDetIdMappingDevice const& elecToDetId,
                                         BTLElectronicsIndexer const& indexer,
                                         int32_t nDigisTot);

  private:
    // ---- intermediate products, kept from prepareDigis() (acquire) to produceDigis() (produce)
    std::optional<BTLChannelPayloadDeviceCollection> channelsPayload_d_;
    BTLChipSegments chipSegments_d_;

    // ----  std::optional: buffers are (re)allocated for each event with the caching allocator, and must not be
    // reused across events: produce() does not wait for the fill kernel, so the previous event may still be using them ---
    std::optional<cms::alpakatools::device_buffer<Device, int32_t[]>> segmentStart_d_;  // starting of chip segments
    std::optional<cms::alpakatools::device_buffer<Device, int32_t[]>>
        nSegments_d_;  // total number of chip segments ( = active chips)

    std::optional<cms::alpakatools::device_buffer<Device, int32_t[]>>
        pairPlusGlobalIndex_;  // global index in channelPayload_d of the plus-side channel
    std::optional<cms::alpakatools::device_buffer<Device, int32_t[]>>
        pairMinusGlobalIndex_;  // global index in channelPayload_d of the minus-side channel
    std::optional<cms::alpakatools::device_buffer<Device, int32_t[]>>
        digiCountPerSegment_;  // number of digis for each chip segment
    std::optional<cms::alpakatools::device_buffer<Device, int32_t>> nDigisTot_d_;  // total number of digis
    std::optional<cms::alpakatools::device_buffer<Device, int32_t[]>>
        segmentDigiOffset_;  // offset, needed to write digis in the same ordering as input raw data.

    //int32_t digiCountCapacity_ = 0;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
