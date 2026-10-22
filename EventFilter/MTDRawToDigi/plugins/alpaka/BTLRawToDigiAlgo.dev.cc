#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/prefixScan.h"

#include "BTLRawToDigiAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  ALPAKA_FN_ACC ALPAKA_FN_INLINE uint64_t extractBits(uint64_t lo, uint64_t hi, int pos, int width) {
    // width must be in [1, 64]
    const uint64_t mask = (width == 64) ? ~0ULL : ((1ULL << width) - 1ULL);

    if (pos < 64) {
      // Entire field is in the low word
      if (pos + width <= 64)
        return (lo >> pos) & mask;

      // Field spans the 64-bit boundary
      const uint64_t lowPart = lo >> pos;
      const uint64_t highPart = hi << (64 - pos);
      return (lowPart | highPart) & mask;
    }

    // Entire field is in the high word
    return (hi >> (pos - 64)) & mask;
  }

  //---------------------------------------------------------------------
  // KERNEL 1: one thread per raw channel payload.
  // Fills the device-side auxiliary SoA
  //---------------------------------------------------------------------
  class BTLDecodeChannelsKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  const uint64_t* rawWords,     // 2 x uint64 per channel
                                  const int32_t* channelFedId,  // one entry per channel, filled host-side
                                  int32_t nChannels,
                                  BTLElectronicsIndexer indexer,
                                  BTLChannelPayloadSoA::View out) const {
      for (int32_t i : cms::alpakatools::uniform_elements(acc, nChannels)) {
        const uint64_t lo = rawWords[2 * i];
        const uint64_t hi = rawWords[2 * i + 1];
        //***********
        // !!!!!!!!!!!!!!! FIXME: Update bit positions once the final BTL payload format is frozen. !!!!!!!!!!!!!!!!!!
        //***********
        out[i].bc0count() = static_cast<uint16_t>(extractBits(lo, hi, 118, 10));
        out[i].status() = extractBits(lo, hi, 117, 1) != 0;
        const uint8_t hslink = static_cast<uint8_t>(extractBits(lo, hi, 104, 6));
        const uint8_t elink = static_cast<uint8_t>(extractBits(lo, hi, 99, 5));
        out[i].hsLinkId() = hslink;
        out[i].eLinkId() = elink;
        out[i].chId() = static_cast<uint8_t>(extractBits(lo, hi, 94, 5));
        out[i].bcCount() = static_cast<uint32_t>(extractBits(lo, hi, 82, 12));
        out[i].t1Coarse() = static_cast<uint16_t>(extractBits(lo, hi, 67, 15));
        out[i].t2Coarse() = static_cast<uint16_t>(extractBits(lo, hi, 57, 10));
        out[i].eoiCoarse() = static_cast<uint16_t>(extractBits(lo, hi, 47, 10));
        out[i].charge() = static_cast<uint16_t>(extractBits(lo, hi, 37, 10));
        out[i].t1Fine() = static_cast<uint16_t>(extractBits(lo, hi, 27, 10));
        out[i].t2Fine() = static_cast<uint16_t>(extractBits(lo, hi, 17, 10));
        out[i].idleTime() = static_cast<uint16_t>(extractBits(lo, hi, 7, 10));
        out[i].prevTrigF() = static_cast<uint8_t>(extractBits(lo, hi, 3, 4));
        out[i].tacId() = static_cast<uint8_t>(extractBits(lo, hi, 0, 3));

        const int32_t fed = channelFedId[i];
        out[i].fedId() = fed;
        // Dense chip index.  Convert (fedId, hsLinkId, eLinkId) into a contiguous index in [0, nChips).
        // !!!!! NOTE: In the input, the channels belonging to the same chip are contiguous
        // (the chips themselves are NOT sorted by chipKey).
        // This STEP 1 preserves the input ordering because each input channel i
        // is decoded into output row i, so channels belonging to the same chip form a contiguous segment.
        // This is exploited in step 1.5 to detect boundaries between segments of data corresponding to different chips
        out[i].chipKey() =
            static_cast<uint32_t>(indexer.flatIndex(fed, hslink, elink, 0) / BTLElectronicsIndexer::nChannelsPerChip);

        // TEMPORARY PRINT JUST FOR DEBUGGING
        /*
	  printf("KERNEL i=%d fed=%u hs=%u elink=%u ch=%u chipKey=%u (0x%08x)\n",
	       i,
	       static_cast<unsigned>(out[i].fedId()),
	       static_cast<unsigned>(out[i].hsLinkId()),
	       static_cast<unsigned>(out[i].eLinkId()),
	       static_cast<unsigned>(out[i].chId()),
	       static_cast<unsigned>(out[i].chipKey()),
	       static_cast<unsigned>(out[i].chipKey())
	       );
	*/
      }
    }
  };

  //---------------------------------------------------------------------
  // STEP 1.5: finds the segments corresponding to each chip, in parallel.
  // In the BTLChannelPayloadSoA filled in STEP1 the channels of the same chip are contiguous, so a
  // segment starts at each channel whose chipKey differs from the previous one.
  // The segments are numbered in the order of the raw data (not sorted by chipKey):
  //   KERNEL 1.5a: one thread per channel, chipStartFlag[i] = 1 if channel i is the first of its chip
  //   inclusive prefix scan of chipStartFlag (multiBlockPrefixScan) -> chipStartRank
  //   KERNEL 1.5b: one thread per channel, the first channel of each chip writes the start of its
  //                segment, segmentStart[chipStartRank[i] - 1] = i; the last channel writes the end
  //                of the last segment and the number of segments.
  // About 332k channels in the worst case, typically ~30k.
  //-----------------------------------------------------------------------------------------------------
  class BTLFlagChipStartsKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  BTLChannelPayloadSoA::ConstView channelData,
                                  int32_t nChannels,
                                  uint32_t* chipStartFlag) const {
      for (int32_t i : cms::alpakatools::uniform_elements(acc, nChannels)) {
        chipStartFlag[i] = (i == 0 || channelData[i].chipKey() != channelData[i - 1].chipKey()) ? 1 : 0;
      }
    }
  };

  class BTLFindChipSegmentsKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  const uint32_t* chipStartFlag,
                                  const uint32_t* chipStartRank,  // inclusive prefix scan of chipStartFlag
                                  int32_t nChannels,
                                  int32_t* segmentStart,
                                  int32_t* nSegmentsOut) const {
      for (int32_t i : cms::alpakatools::uniform_elements(acc, nChannels)) {
        if (chipStartFlag[i])
          segmentStart[chipStartRank[i] - 1] = i;
        if (i == nChannels - 1) {
          const int32_t nSegments = chipStartRank[i];
          segmentStart[nSegments] = nChannels;  // end of last segment
          *nSegmentsOut = nSegments;
        }
      }
    }
  };

  //---------------------------------------------------------------------
  // KERNEL 2a: one thread per chip segment: pair plus/minus channels and count digis per segment.
  // For each segment, scans all channels, matches each plus-side channel with its partner minus-side
  // channel (if present), and stores the result as a list of (plus, minus) index pairs.
  // Complete pairs are counted once from the plus side; single-side hits are also kept.
  //
  // Output arrays (size nSegments * nPairsPerChip):
  //   pairPlusGlobalIndex[seg * nPairsPerChip + k]  = index in channelPayload_d of the plus-side channel of the k-th digi in segment seg (-1 if the digi has no plus-side hit)
  //   pairMinusGlobalIndex[seg * nPairsPerChip + k] = same for the minus-side channel (-1 if absent)
  //   digiCountPerSegment[seg]                      = number of digis found in segment seg
  //---------------------------------------------------------------------
  class BTLPairAndCountDigisKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  BTLChannelPayloadSoA::ConstView channelData,
                                  const int32_t* nSegments,
                                  const int32_t* segmentStart,
                                  BTLElectronicsToDetIdMappingSoA::ConstView mapping,
                                  BTLElectronicsIndexer indexer,
                                  int32_t nChips,
                                  int32_t* pairPlusGlobalIndex,  // size nSegmentsUpperBound * nPairsMax
                                  int32_t* pairMinusGlobalIndex,
                                  int32_t* digiCountPerSegment) const {
      for (int32_t segment : cms::alpakatools::uniform_elements(acc, nChips)) {
        if (segment >= *nSegments)
          continue;

        //      for (int32_t segment : cms::alpakatools::uniform_elements(acc, *nSegments)) {

        const int32_t start = segmentStart[segment];
        const int32_t end = segmentStart[segment + 1];

        const int32_t base = segment * BTLElectronicsIndexer::nPairsPerChip;

        int32_t nDigis = 0;  // local thread counter

        for (int32_t i = start; i < end; ++i) {
          const int32_t mappingIndex = indexer.flatIndex(
              channelData[i].fedId(), channelData[i].hsLinkId(), channelData[i].eLinkId(), channelData[i].chId());

          const auto& m = mapping[mappingIndex];

          if (!m.valid()) {
            printf("INVALID MAPPING: fed=%u hs=%u elink=%u ch=%u mappingIndex=%d\n",
                   static_cast<unsigned int>(channelData[i].fedId()),
                   static_cast<unsigned int>(channelData[i].hsLinkId()),
                   static_cast<unsigned int>(channelData[i].eLinkId()),
                   static_cast<unsigned int>(channelData[i].chId()),
                   mappingIndex);
            continue;
          }

          // Find the partner channel within this chip.
          int32_t partnerRow = -1;
          for (int32_t j = start; j < end; ++j) {
            if (j == i)
              continue;
            if (channelData[j].chId() == m.partnerChannelId()) {
              partnerRow = j;
              break;
            }
          }

          if (partnerRow >= 0) {
            // Complete pair: count it only from the plus side.
            if (m.side() == kPlusSide) {
              pairPlusGlobalIndex[base + nDigis] =
                  i;  // i = global index of this channel. Maybe better to change the name to pairPlusGlobalIndex?????
              pairMinusGlobalIndex[base + nDigis] = partnerRow;
              ++nDigis;
            }
          } else {
            // Single-side digi.
            if (m.side() == kPlusSide) {
              pairPlusGlobalIndex[base + nDigis] = i;
              pairMinusGlobalIndex[base + nDigis] = -1;
            } else {
              pairPlusGlobalIndex[base + nDigis] = -1;
              pairMinusGlobalIndex[base + nDigis] = i;
            }
            ++nDigis;
          }
        }

        digiCountPerSegment[segment] = nDigis;
      }
    }
  };

  //---------------------------------------------------------------------
  // KERNEL 2b: fill output BTLDigiSoA, one thread per chip segment
  // Reads the (plus, minus) index pairs built by BTLPairAndCountDigisKernel and
  // the per-segment digi offsets to write each digi at the correct position in
  // the output collection, preserving the segment ordering of the raw data.
  // No channel pairing is redone here.
  //---------------------------------------------------------------------
  class BTLFillDigisKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  BTLChannelPayloadSoA::ConstView channelData,
                                  const int32_t* nSegments,
                                  const int32_t* digiCountPerSegment,
                                  const int32_t* segmentDigiOffset,
                                  const int32_t* pairPlusGlobalIndex,
                                  const int32_t* pairMinusGlobalIndex,
                                  BTLElectronicsToDetIdMappingSoA::ConstView mapping,
                                  BTLElectronicsIndexer indexer,
                                  int32_t nChips,
                                  ::btldigi::BTLDigiSoA::View digisOut) const {
      for (int32_t segment : cms::alpakatools::uniform_elements(acc, nChips)) {
        if (segment >= *nSegments)
          continue;

        //for (int32_t segment : cms::alpakatools::uniform_elements(acc, *nSegments)) {
        const int32_t base = segment * BTLElectronicsIndexer::nPairsPerChip;
        const int32_t writeBase = segmentDigiOffset[segment];

        for (int32_t k = 0; k < digiCountPerSegment[segment]; ++k) {
          const int32_t idxPlus = pairPlusGlobalIndex[base + k];
          const int32_t idxMinus = pairMinusGlobalIndex[base + k];
          const int32_t validIdx = (idxMinus >= 0) ? idxMinus : idxPlus;

          const int32_t mappingIndex = indexer.flatIndex(channelData[validIdx].fedId(),
                                                         channelData[validIdx].hsLinkId(),
                                                         channelData[validIdx].eLinkId(),
                                                         channelData[validIdx].chId());
          const auto& m = mapping[mappingIndex];

          auto d = digisOut[writeBase + k];
          d.rawId() = m.rawId();

          const auto& cValid = channelData[validIdx];
          d.BC0count() = cValid.bc0count();
          d.status() = cValid.status();
          d.BCcount() = cValid.bcCount();

          // minus side - set to default values (0 for now) if not present
          if (idxMinus >= 0) {
            const auto& cm = channelData[idxMinus];
            d.chIDMinus() = cm.chId();
            d.T1coarseMinus() = cm.t1Coarse();
            d.T2coarseMinus() = cm.t2Coarse();
            d.EOIcoarseMinus() = cm.eoiCoarse();
            d.ChargeMinus() = cm.charge();
            d.T1fineMinus() = cm.t1Fine();
            d.T2fineMinus() = cm.t2Fine();
            d.IdleTimeMinus() = cm.idleTime();
            d.PrevTrigFMinus() = cm.prevTrigF();
            d.TACIDMinus() = cm.tacId();
          } else {
            d.chIDMinus() = 0;
            d.T1coarseMinus() = 0;
            d.T2coarseMinus() = 0;
            d.EOIcoarseMinus() = 0;
            d.ChargeMinus() = 0;
            d.T1fineMinus() = 0;
            d.T2fineMinus() = 0;
            d.IdleTimeMinus() = 0;
            d.PrevTrigFMinus() = 0;
            d.TACIDMinus() = 0;
          }
          // plus side - set to default values (0 for now) if not present
          if (idxPlus >= 0) {
            const auto& cp = channelData[idxPlus];
            d.chIDPlus() = cp.chId();
            d.T1coarsePlus() = cp.t1Coarse();
            d.T2coarsePlus() = cp.t2Coarse();
            d.EOIcoarsePlus() = cp.eoiCoarse();
            d.ChargePlus() = cp.charge();
            d.T1finePlus() = cp.t1Fine();
            d.T2finePlus() = cp.t2Fine();
            d.IdleTimePlus() = cp.idleTime();
            d.PrevTrigFPlus() = cp.prevTrigF();
            d.TACIDPlus() = cp.tacId();
          } else {
            d.chIDPlus() = 0;
            d.T1coarsePlus() = 0;
            d.T2coarsePlus() = 0;
            d.EOIcoarsePlus() = 0;
            d.ChargePlus() = 0;
            d.T1finePlus() = 0;
            d.T2finePlus() = 0;
            d.IdleTimePlus() = 0;
            d.PrevTrigFPlus() = 0;
            d.TACIDPlus() = 0;
          }
        }
      }
    }
  };

  //---------------------------------------------------------------------
  // BTLRawToDigiAlgo
  //---------------------------------------------------------------------
  BTLChannelPayloadDeviceCollection BTLRawToDigiAlgo::decodeChannels(Queue& queue,
                                                                     const uint64_t* rawWords_h,
                                                                     const int32_t* channelFedId_h,
                                                                     int32_t nChannels,
                                                                     BTLElectronicsIndexer const& indexer) {
#ifdef EDM_ML_DEBUG
    LogDebug("BTLRawToDigi") << "========== BTLRawToDigi [decodeChannels] START ========== \n"
                             << " decoding " << nChannels << " channels";
#endif

    // Copy host --> device
    auto rawWords_d = cms::alpakatools::make_device_buffer<uint64_t[]>(queue, 2 * nChannels);
    auto channelFedId_d = cms::alpakatools::make_device_buffer<int32_t[]>(queue, nChannels);
    alpaka::memcpy(queue, rawWords_d, cms::alpakatools::make_host_view(rawWords_h, 2 * nChannels));
    alpaka::memcpy(queue, channelFedId_d, cms::alpakatools::make_host_view(channelFedId_h, nChannels));

    BTLChannelPayloadDeviceCollection channelsPayload_d(queue, nChannels);

    auto workDiv =
        cms::alpakatools::make_workdiv<Acc1D>(cms::alpakatools::divide_up_by(uint32_t(nChannels), 256u), 256u);
    alpaka::exec<Acc1D>(queue,
                        workDiv,
                        BTLDecodeChannelsKernel{},
                        rawWords_d.data(),
                        channelFedId_d.data(),
                        nChannels,
                        indexer,
                        channelsPayload_d.view());

    return channelsPayload_d;
  }

  BTLChipSegments BTLRawToDigiAlgo::findChipSegments(Queue& queue,
                                                     BTLChannelPayloadDeviceCollection const& channelPayload_d,
                                                     int32_t nChannels,
                                                     BTLElectronicsIndexer const& indexer) {
    const int32_t nChips = indexer.size() / indexer.nChannelsPerELink;

#ifdef EDM_ML_DEBUG
    LogDebug("BTLRawToDigi") << "========== BTLRawToDigi [findChipSegments] START ==========";
    LogDebug("BTLRawToDigi") << "========== BTLRawToDigi [findChipSegments] : total number of chips = " << nChips;
#endif

    // allocated for each event (cheap with the caching allocator): the previous event may still be
    // using its buffers on another queue, since produce() does not wait for the fill kernel
    segmentStart_d_.emplace(cms::alpakatools::make_device_buffer<int32_t[]>(queue, nChips + 1));
    nSegments_d_.emplace(cms::alpakatools::make_device_buffer<int32_t[]>(queue, 1));

    // reset counters before launching the kernel
    alpaka::memset(queue, *nSegments_d_, 0);
    alpaka::memset(queue, *segmentStart_d_, 0);

    // temporary buffers, only used by the kernels enqueued here
    auto chipStartFlag_d = cms::alpakatools::make_device_buffer<uint32_t[]>(queue, nChannels);
    auto chipStartRank_d = cms::alpakatools::make_device_buffer<uint32_t[]>(queue, nChannels);

    // KERNEL 1.5a: flag the first channel of each chip
    auto workDiv =
        cms::alpakatools::make_workdiv<Acc1D>(cms::alpakatools::divide_up_by(uint32_t(nChannels), 256u), 256u);
    alpaka::exec<Acc1D>(
        queue, workDiv, BTLFlagChipStartsKernel{}, channelPayload_d.const_view(), nChannels, chipStartFlag_d.data());

    // inclusive prefix scan of the flags: chipStartRank[i] = number of segments starting at or before channel i
    const uint32_t threadsScan = 1024u;
    const uint32_t nBlocksScan = cms::alpakatools::divide_up_by(uint32_t(nChannels), threadsScan);
    const auto workDivScan = cms::alpakatools::make_workdiv<Acc1D>(nBlocksScan, threadsScan);
    auto blockCounter_d = cms::alpakatools::make_device_buffer<int32_t>(queue);
    alpaka::memset(queue, blockCounter_d, 0);
    cms::alpakatools::checkSharedMemoryPrefixScan<Acc1D>(uint32_t(nChannels), nBlocksScan, alpaka::getDev(queue));
    alpaka::exec<Acc1D>(queue,
                        workDivScan,
                        cms::alpakatools::multiBlockPrefixScan<uint32_t>(),
                        chipStartFlag_d.data(),
                        chipStartRank_d.data(),
                        uint32_t(nChannels),
                        int32_t(nBlocksScan),
                        blockCounter_d.data(),
                        alpaka::getPreferredWarpSize(alpaka::getDev(queue)));

    // KERNEL 1.5b: write the start of each segment, in the order of the raw data
    alpaka::exec<Acc1D>(queue,
                        workDiv,
                        BTLFindChipSegmentsKernel{},
                        chipStartFlag_d.data(),
                        chipStartRank_d.data(),
                        nChannels,
                        segmentStart_d_->data(),
                        nSegments_d_->data());

#ifdef EDM_ML_DEBUG
    LogDebug("BTLRawToDigi") << "========== BTLRawToDigi [findChipSegments] : kernel launched, returning "
                                "segmentStart_d_, nSegments_d_ ==========";
#endif
    return {segmentStart_d_->data(), nSegments_d_->data()};
  }

  void BTLRawToDigiAlgo::countDigis(Queue& queue,
                                    BTLChannelPayloadDeviceCollection const& channelPayload_d,
                                    BTLChipSegments const& segments,
                                    BTLElectronicsToDetIdMappingDevice const& mapping,
                                    BTLElectronicsIndexer const& indexer,
                                    //int32_t nActiveChips) {
                                    int32_t nChips) {
#ifdef EDM_ML_DEBUG
    LogDebug("BTLRawToDigi") << "========== BTLRawToDigi countDigis() START ==========";
#endif

    // -- buffers dimensioned on the total number of chips, allocated for each event (see findChipSegments)
    digiCountPerSegment_.emplace(cms::alpakatools::make_device_buffer<int32_t[]>(queue, nChips));
    pairPlusGlobalIndex_.emplace(
        cms::alpakatools::make_device_buffer<int32_t[]>(queue, nChips * BTLElectronicsIndexer::nPairsPerChip));
    pairMinusGlobalIndex_.emplace(
        cms::alpakatools::make_device_buffer<int32_t[]>(queue, nChips * BTLElectronicsIndexer::nPairsPerChip));
    segmentDigiOffset_.emplace(cms::alpakatools::make_device_buffer<int32_t[]>(queue, nChips + 1));

    // --- reset counters before launching the kernel
    alpaka::memset(queue, *digiCountPerSegment_, 0);
    alpaka::memset(queue, *pairPlusGlobalIndex_, 0xFF);   // -1 come int32_t
    alpaka::memset(queue, *pairMinusGlobalIndex_, 0xFF);  // -1 come int32_t

    auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(cms::alpakatools::divide_up_by(uint32_t(nChips), 256u), 256u);

    alpaka::exec<Acc1D>(queue,
                        workDiv,
                        BTLPairAndCountDigisKernel{},
                        channelPayload_d.const_view(),
                        segments.nSegments,
                        segments.segmentStart,
                        mapping.const_view(),
                        indexer,
                        nChips,
                        pairPlusGlobalIndex_->data(),
                        pairMinusGlobalIndex_->data(),
                        digiCountPerSegment_->data());

#ifdef EDM_ML_DEBUG
    LogDebug("BTLRawToDigi") << "========== BTLRawToDigi [countDigis] : kernel launched ==========";
#endif
  }

  void BTLRawToDigiAlgo::computeDigiOffsets(Queue& queue,
                                            int32_t nChips,
                                            cms::alpakatools::host_buffer<int32_t>& nDigisTot_h) {
#ifdef EDM_ML_DEBUG
    LogDebug("BTLRawToDigi") << "========== BTLRawToDigi computeDigiOffsets() START ==========";
#endif

    // STEP 2a.5: digi offsets per chip segment, as an exclusive prefix sum of digiCountPerSegment:
    //   segmentDigiOffset[s]      = index of the first digi of segment s in the output BTLDigiSoA
    //   segmentDigiOffset[nChips] = nDigisTot = total number of digis in the event
    // The inclusive prefix scan (multiBlockPrefixScan) is written shifted by one element,
    // into segmentDigiOffset + 1, so that segmentDigiOffset[0] = 0 and the result is the exclusive scan.
    // NOTE: digiCountPerSegment must be zeroed for segments >= nSegments
    // (memset before BTLPairAndCountDigisKernel), since the scan runs over all the nChips.
    alpaka::memset(queue, *segmentDigiOffset_, 0);

    const uint32_t threadsScan = 1024u;
    const uint32_t nBlocksScan = cms::alpakatools::divide_up_by(uint32_t(nChips), threadsScan);
    const auto workDivScan = cms::alpakatools::make_workdiv<Acc1D>(nBlocksScan, threadsScan);
    auto blockCounter_d = cms::alpakatools::make_device_buffer<int32_t>(queue);
    alpaka::memset(queue, blockCounter_d, 0);
    cms::alpakatools::checkSharedMemoryPrefixScan<Acc1D>(uint32_t(nChips), nBlocksScan, alpaka::getDev(queue));
    alpaka::exec<Acc1D>(queue,
                        workDivScan,
                        cms::alpakatools::multiBlockPrefixScan<int32_t>(),
                        digiCountPerSegment_->data(),
                        segmentDigiOffset_->data() + 1,
                        uint32_t(nChips),
                        int32_t(nBlocksScan),
                        blockCounter_d.data(),
                        alpaka::getPreferredWarpSize(alpaka::getDev(queue)));

    // only device -> host copy: nDigisTot (4 bytes, last element of segmentDigiOffset), needed to allocate
    // the output collection.
    // No wait here: the value is read in produce(), after the framework has synchronized the queue.
    alpaka::memcpy(queue, nDigisTot_h, cms::alpakatools::make_device_view(queue, segmentDigiOffset_->data()[nChips]));
  }

  BTLDigiDeviceCollection BTLRawToDigiAlgo::fillDigis(Queue& queue,
                                                      BTLChannelPayloadDeviceCollection const& channelPayload_d,
                                                      BTLChipSegments const& segments,
                                                      BTLElectronicsToDetIdMappingDevice const& mapping,
                                                      BTLElectronicsIndexer const& indexer,
                                                      int32_t nChips,
                                                      int32_t nDigisTot) {
#ifdef EDM_ML_DEBUG
    LogDebug("BTLRawToDigi") << "========== BTLRawToDigi fillDigis() START ==========";
    LogDebug("BTLRawToDigi") << " === BTLRawToDigi fillDigis(): total number of digis = " << nDigisTot;
#endif

    // protection against nDigisTot = 0 (should not happen, except if the mapping is wrong)
    if (nDigisTot == 0)
      return BTLDigiDeviceCollection(queue, 0);

    BTLDigiDeviceCollection digis_d(queue, nDigisTot);

    auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(cms::alpakatools::divide_up_by(uint32_t(nChips), 256u), 256u);
    alpaka::exec<Acc1D>(queue,
                        workDiv,
                        BTLFillDigisKernel{},
                        channelPayload_d.const_view(),
                        segments.nSegments,
                        digiCountPerSegment_->data(),
                        segmentDigiOffset_->data(),
                        pairPlusGlobalIndex_->data(),
                        pairMinusGlobalIndex_->data(),
                        mapping.const_view(),
                        indexer,
                        nChips,
                        digis_d.view());

    return digis_d;
  }

  void BTLRawToDigiAlgo::prepareDigis(Queue& queue,
                                      const uint64_t* rawWords_h,
                                      const int32_t* channelFedId_h,
                                      int32_t nChannels,
                                      BTLElectronicsToDetIdMappingDevice const& elecToDetIdMapping,
                                      BTLElectronicsIndexer const& indexer,
                                      cms::alpakatools::host_buffer<int32_t>& nDigisTot_h) {
#ifdef EDM_ML_DEBUG
    LogDebug("BTLRawToDigi") << "========== BTLRawToDigi prepareDigis() START ==========";
#endif

    // -- STEP1: Launch decode for each channel
    channelsPayload_d_.emplace(
        decodeChannels(queue, rawWords_h, channelFedId_h, nChannels, indexer));  // auxiliary soa on device
    auto const& channelsPayload_d = *channelsPayload_d_;

    // -- DEBUG STEP1: copy to host and print
#ifdef EDM_ML_DEBUG
    BTLChannelPayloadHostCollection channelsPayload_h(std::max(nChannels, 1));
    alpaka::memcpy(queue, channelsPayload_h.buffer(), channelsPayload_d.buffer());
    alpaka::wait(queue);
    auto soa = channelsPayload_h.view();

    std::ostringstream os;
    for (int i = 0; i < nChannels; ++i) {
      os << "i = " << i << " fedId = " << soa.fedId()[i] << " hsLinkId = " << static_cast<int>(soa.hsLinkId()[i])
         << " eLinkId = " << static_cast<int>(soa.eLinkId()[i]) << " chId = " << static_cast<int>(soa.chId()[i])
         << " chipKey = " << soa.chipKey()[i] << "\n";
    }
    LogDebug("BTLRawToDigi") << os.str();
#endif

    // -- STEP1.5: find segments corresponding to each active chip
    chipSegments_d_ = findChipSegments(queue, channelsPayload_d, nChannels, indexer);

    const int32_t nChips = indexer.size() / indexer.nChannelsPerELink;

    // -- DEBUG findChipSegments
#ifdef EDM_ML_DEBUG
    auto nSegments_h = cms::alpakatools::make_host_buffer<int32_t[]>(queue, 1);
    alpaka::memcpy(queue, nSegments_h, *nSegments_d_);
    alpaka::wait(queue);
    const int32_t nActiveChips = nSegments_h[0];

    LogDebug("BTLRawToDigi") << "Debugging findChipSegments...";
    LogDebug("BTLRawToDigi") << "Number of active chips = " << nActiveChips;

    auto segmentStart_h =
        cms::alpakatools::make_host_buffer<int32_t[]>(queue, nChips + 1);  // alloco con lunghezza nChips
    alpaka::memcpy(queue, segmentStart_h, *segmentStart_d_);
    alpaka::wait(queue);

    std::ostringstream os1;
    for (int32_t s = 0; s < nActiveChips; ++s) {  // loop only over nActiveChips elements
      os1 << "segment " << s << " : [" << segmentStart_h[s] << ", " << segmentStart_h[s + 1] << ")"
          << " size=" << segmentStart_h[s + 1] - segmentStart_h[s] << '\n';
    }
    LogDebug("BTLRawToDigi") << os1.str();
#endif

    // STEP 2a: pairing and count digis
    //
    // NOTE: nSegments is the number of active chips found in STEP 1.5,
    // while nChips is the maximum number of chips.
    //
    //
    // Two options:
    // 1) Use nChips and skip segments >= nSegments inside the kernel.
    // 2) Copy nSegments from device to host and use it to define the
    //    work division directly.

    countDigis(queue, channelsPayload_d, chipSegments_d_, elecToDetIdMapping, indexer, nChips);  // 1)
    //countDigis(queue, channelsPayload_d, chipSegments_d_, elecToDetIdMapping, indexer, nActiveChips); // 2)

    // -- DEBUG pair and count digis
#ifdef EDM_ML_DEBUG
    auto digiCount_h =
        cms::alpakatools::make_host_buffer<int32_t[]>(queue, alpaka::getExtentProduct(*digiCountPerSegment_));
    alpaka::memcpy(queue, digiCount_h, *digiCountPerSegment_);
    alpaka::wait(queue);

    std::ostringstream os2;
    for (int32_t s = 0; s < nActiveChips; ++s) {
      os2 << "segment " << s << "  number of digis in this segment : " << digiCount_h[s] << '\n';
    }
    LogDebug("BTLRawToDigi") << os2.str();
#endif

    // -- STEP2a.5: digi offsets per segment and total number of digis (async copy to nDigisTot_h)
    computeDigiOffsets(queue, nChips, nDigisTot_h);
  }

  BTLDigiDeviceCollection BTLRawToDigiAlgo::produceDigis(Queue& queue,
                                                         BTLElectronicsToDetIdMappingDevice const& elecToDetIdMapping,
                                                         BTLElectronicsIndexer const& indexer,
                                                         int32_t nDigisTot) {
#ifdef EDM_ML_DEBUG
    LogDebug("BTLRawToDigi") << "========== BTLRawToDigi produceDigis() START ==========";
#endif

    const int32_t nChips = indexer.size() / indexer.nChannelsPerELink;

    // -- STEP2b: fill digis
    auto digis_d =
        fillDigis(queue, *channelsPayload_d_, chipSegments_d_, elecToDetIdMapping, indexer, nChips, nDigisTot);

    // the intermediate SoA is not needed anymore: release it (the caching allocator frees it
    // only once the work already enqueued on the queue has completed)
    channelsPayload_d_.reset();

    // -- DEBUG FILL DIGIS
#ifdef EDM_ML_DEBUG
    ::btldigi::BTLDigiHostCollection digis_h(queue, std::max(nDigisTot, 1));
    alpaka::memcpy(queue, digis_h.buffer(), digis_d.buffer());
    alpaka::wait(queue);
    auto digis_h_view = digis_h.view();
    std::ostringstream os3;
    for (int32_t i = 0; i < nDigisTot; ++i) {
      os3 << "DIGI " << i << " rawId=" << digis_h_view[i].rawId() << " BC0=" << digis_h_view[i].BC0count()
          << " status=" << digis_h_view[i].status() << " BCcount=" << digis_h_view[i].BCcount()
          << " plus(ch=" << int(digis_h_view[i].chIDPlus()) << ", charge=" << digis_h_view[i].ChargePlus() << ")"
          << " minus(ch=" << int(digis_h_view[i].chIDMinus()) << ", charge=" << digis_h_view[i].ChargeMinus() << ")"
          << '\n';
    }
    LogDebug("BTLRawToDigi") << os3.str();
#endif

    // -- return the filled digi collection
    return digis_d;
  }
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
