#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

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
	// !!!!! NOTE: The input channels are ordered by FED, HSLink and ELink.
	// this STEP 1 preserves this ordering because each input channel i
	// is decoded into output row i. Therefore chipKey is non-decreasing in the output SoA, and
	// channels belonging to the same chip form a contiguous segment.
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
  // KERNEL 1.5:  finds segnments corresponding to each chip
  // The BTLChannelPayloadSoA filled in STEP1 is ordered by chip, so the chip boundaries
  // can be found with a simple sequential scan.
  // About 332k in the worst case, typically ~30k. To be understood if this step needs to be parallelized. For now it is not.
  //-----------------------------------------------------------------------------------------------------
  class BTLFindChipSegmentsKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
				  BTLChannelPayloadSoA::ConstView channelData,
				  int32_t* segmentStart,
				  int32_t* nSegmentsOut) const {
      const int32_t nChannels = channelData.metadata().size(); 
      // For now just one thread, to detect chip boundaries. To be checked if performance is acceptable.
      if (cms::alpakatools::once_per_grid(acc)) { // capire meglio cosa fa
	int32_t nSegments = 0;
	for (int32_t i = 0; i < nChannels; ++i) {
	  if (i == 0 || channelData[i].chipKey() != channelData[i - 1].chipKey()) {
	    segmentStart[nSegments] = i;
	    ++nSegments;
	  }
	}
	segmentStart[nSegments] = nChannels;  // end of last segment
	*nSegmentsOut = nSegments;
      }
    }
  };
  
  //---------------------------------------------------------------------
  // KERNEL 2a: one thread per chip: buil pairs and count digis       
  //---------------------------------------------------------------------
  class BTLPairAndCountDigisKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
				  BTLChannelPayloadSoA::ConstView channelData,
				  const int32_t* nSegments,
				  const int32_t* segmentStart,
				  BTLElectronicsToDetIdMappingSoA::ConstView mapping,
				  BTLElectronicsIndexer indexer,
				  //int32_t nChips,
				  int32_t* pairPlusLocalIndex,   // size nSegmentsUpperBound * nPairsMax
				  int32_t* pairMinusLocalIndex,
      				  int32_t* digiCountPerSegment) const {

    //    for (int32_t chip : cms::alpakatools::uniform_elements(acc, nChips)) {
    //if (chip >= *nSegments)
    //continue;
      
      for (int32_t segment : cms::alpakatools::uniform_elements(acc, *nSegments)) {
	
	const int32_t start = segmentStart[segment];
	const int32_t end = segmentStart[segment + 1];
	
	const int32_t base = segment * BTLElectronicsIndexer::nPairsPerChip;
	
	int32_t nDigis = 0; // local thread counter
	
	for (int32_t i = start; i < end; ++i) {
	  
	  const int32_t mappingIndex = indexer.flatIndex(channelData[i].fedId(),
							 channelData[i].hsLinkId(),
							 channelData[i].eLinkId(),
							 channelData[i].chId());
	  
	  const auto& m = mapping[mappingIndex];
	  
	  if (!m.valid())
	    continue;
	  
	  //const int32_t pairId = m.pairId(); // not needed anymore?
	  
	  
	  // Find the partner channel within this chip.
	  int32_t partnerRow = -1;
	  for (int32_t j = start; j < end; ++j) {
	    if (j == i) continue;
	    if (channelData[j].chId() == m.partnerChannelId()){
	      partnerRow = j;
	      break;
	    }
	  }
	  
	  if (partnerRow >= 0) {
	    // Complete pair: count it only from the plus side.
	    if (m.side() == kPlusSide) {
	      pairPlusLocalIndex[base + nDigis] = i; // i = global index of this channel
	      pairMinusLocalIndex[base + nDigis] = partnerRow;
	      ++nDigis;
	    }
	  } else {
	    // Single-side digi.
	    if (m.side() == kPlusSide) {
	      pairPlusLocalIndex[base + nDigis] = i;
	      pairMinusLocalIndex[base + nDigis] = -1;
	    } else {
	      pairPlusLocalIndex[base + nDigis] = -1;
	      pairMinusLocalIndex[base + nDigis] = i;
	    }
	    ++nDigis;
	  }
	}

	digiCountPerSegment[segment] = nDigis;
      }
    }
  };


  //---------------------------------------------------------------------
  // KERNEL 2b: fill digis, one thread per segment
  // Reuse indices found in BTLPairAndCountDigisKernel, no need of redoing the pairing
  //---------------------------------------------------------------------
  class BTLFillDigisKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
				  BTLChannelPayloadSoA::ConstView channelData,
				  const int32_t* nSegments,
				  const int32_t* digiCountPerSegment,
				  const int32_t* segmentDigiOffset,
				  const int32_t* pairPlusLocalIndex,
				  const int32_t* pairMinusLocalIndex,
				  BTLElectronicsToDetIdMappingSoA::ConstView mapping,
				  BTLElectronicsIndexer indexer,
				  ::btldigi::BTLDigiSoA::View digisOut) const {
      for (int32_t segment : cms::alpakatools::uniform_elements(acc, *nSegments)) {
	const int32_t base = segment * BTLElectronicsIndexer::nPairsPerChip;
	const int32_t writeBase = segmentDigiOffset[segment];
	
	for (int32_t k = 0; k < digiCountPerSegment[segment]; ++k) {
	  const int32_t idxPlus = pairPlusLocalIndex[base + k];
	  const int32_t idxMinus = pairMinusLocalIndex[base + k];
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
	  }
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
	  }
	}
      }
    }
  };
  
  
  //---------------------------------------------------------------------
  // Btlrawtodigialgo
  //---------------------------------------------------------------------
  BTLChannelPayloadDeviceCollection BTLRawToDigiAlgo::decodeOnly(Queue& queue,
                                                                 const uint64_t* rawWords_h,
                                                                 const int32_t* channelFedId_h,
                                                                 int32_t nChannels,
                                                                 BTLElectronicsIndexer const& indexer) {
    // Copy host --> device
    auto rawWords_d =
        cms::alpakatools::make_device_buffer<uint64_t[]>(queue, 2 * std::max(nChannels, 1));  // alloca memoria su gpu
    auto channelFedId_d = cms::alpakatools::make_device_buffer<int32_t[]>(queue, std::max(nChannels, 1));
    if (nChannels > 0) {
      alpaka::memcpy(queue, rawWords_d, cms::alpakatools::make_host_view(rawWords_h, 2 * nChannels));
      alpaka::memcpy(queue, channelFedId_d, cms::alpakatools::make_host_view(channelFedId_h, nChannels));
    }

    BTLChannelPayloadDeviceCollection channelsPayload_d(queue, std::max(nChannels, 1));
    if (nChannels > 0) {
      auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(cms::alpakatools::divide_up_by(uint32_t(nChannels), 256u), 256u);
      alpaka::exec<Acc1D>(queue,
                          workDiv,
                          BTLDecodeChannelsKernel{},
                          rawWords_d.data(),
                          channelFedId_d.data(),
                          nChannels,
                          indexer,
                          channelsPayload_d.view());
    }
    return channelsPayload_d;
  }
  


  BTLChipSegments BTLRawToDigiAlgo::findChipSegments(Queue& queue,
						     BTLChannelPayloadDeviceCollection const& channelPayload_d,
						     BTLElectronicsIndexer const& indexer) {

    std::cout << "[findChipSegments] START" << std::endl;

    const int32_t nChips = indexer.nFeds * indexer.nHsLinks * indexer.nELinks;

    std::cout << "[findChipSegments] nChips = "        << nChips << std::endl;
    
    if (!segmentStart_d_) {
      std::cout << "[findChipSegments] allocating segmentStart_d_"
              << std::endl;
      segmentStart_d_.emplace(cms::alpakatools::make_device_buffer<int32_t[]>(queue, nChips + 1));
      std::cout << "[findChipSegments] segmentStart_d_ allocated"
              << std::endl;
    }
    if (!nSegments_d_) {
      std::cout << "[findChipSegments] allocating nSegments_d_"
              << std::endl;
      nSegments_d_.emplace(cms::alpakatools::make_device_buffer<int32_t[]>(queue, 1));
      std::cout << "[findChipSegments] nSegments_d_ allocated"
              << std::endl;
    }
    
    std::cout << "[findChipSegments] creating workDiv"
            << std::endl;
    auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(1, 1);

    std::cout << "[findChipSegments] launching kernel"
	      << std::endl;
 
    alpaka::exec<Acc1D>(queue, workDiv, BTLFindChipSegmentsKernel{}, channelPayload_d.const_view(),
			segmentStart_d_->data(), nSegments_d_->data());

    std::cout << "[findChipSegments] kernel launched"
            << std::endl;
    return {segmentStart_d_->data(), nSegments_d_->data()};
  }
  

  
  void BTLRawToDigiAlgo::countDigis(Queue& queue,
				    BTLChannelPayloadDeviceCollection const& channelPayload_d,
				    BTLChipSegments const& segments,
				    BTLElectronicsToDetIdMappingDevice const& mapping,
				    BTLElectronicsIndexer const& indexer,
				    int32_t nActiveChips) {
    

    //buffer dimensionati su nActiveChips e riusati tra eventi,
    if (!digiCountPerSegment_ || digiCountCapacity_ < nActiveChips) {
      //digiCountPerSegment_.emplace(cms::alpakatools::make_device_buffer<int32_t[]>(queue, nChips)); // se voglio dimensionarlo nSegments (chip attivi) invece che nChips (nuemro totale di chip) dovrei conoscere nSegments lato host, pero' implica una sincronizzazione device-host.
      digiCountPerSegment_.emplace(cms::alpakatools::make_device_buffer<int32_t[]>(queue, nActiveChips)); // nActiveChips from host
      pairPlusLocalIndex_.emplace(cms::alpakatools::make_device_buffer<int32_t[]>(queue, nActiveChips * BTLElectronicsIndexer::nPairsPerChip));
      pairMinusLocalIndex_.emplace(cms::alpakatools::make_device_buffer<int32_t[]>(queue, nActiveChips * BTLElectronicsIndexer::nPairsPerChip));
      segmentDigiOffset_.emplace(cms::alpakatools::make_device_buffer<int32_t[]>(queue, nActiveChips + 1));
      totalDigis_.emplace(cms::alpakatools::make_device_buffer<int32_t[]>(queue, 1));
      digiCountCapacity_ = nActiveChips;
    }
    
    std::cout << "[countDigis]  digiCountPerSegment_ allocated" << std::endl; 
      
    
    //auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(nChips, 1); // se volessi solo nSegments work-items dorvei copiare da device a host nsegments, perche' workDiv vuole sapere quanti work-item lato host
    //  auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(nActiveChips, 1); // number of work-items = nActiveChips (workDiv vuole sapere quanti work-item lato host)
    auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(cms::alpakatools::divide_up_by(uint32_t(nActiveChips), 256u),256u);


    std::cout << "[countDigis] launching kernel"
              << std::endl;
     
    alpaka::exec<Acc1D>(queue,
			workDiv,
			BTLPairAndCountDigisKernel{},
			channelPayload_d.const_view(),
			segments.nSegments,
			segments.segmentStart,
			mapping.const_view(),
			indexer,
			pairPlusLocalIndex_->data(),
			pairMinusLocalIndex_->data(),
			digiCountPerSegment_->data());

    std::cout << "[countDigis] kernel launched"
              << std::endl;
  }
  

  BTLDigiDeviceCollection BTLRawToDigiAlgo::fillDigis(Queue& queue,
						      BTLChannelPayloadDeviceCollection const& channelPayload_d,
						      BTLChipSegments const& segments,
						      BTLElectronicsToDetIdMappingDevice const& mapping,
						      BTLElectronicsIndexer const& indexer,
						      int32_t nActiveChips) {

  // copy device to host: needed to correctly dimension the BTLDigiDeviceCollection 
  auto digiCount_h = cms::alpakatools::make_host_buffer<int32_t[]>(queue, std::max(nActiveChips, 1));
  alpaka::memcpy(queue, digiCount_h, *digiCountPerSegment_);
  alpaka::wait(queue);

  auto segmentDigiOffset_h = cms::alpakatools::make_host_buffer<int32_t[]>(queue, nActiveChips + 1);
  int32_t running = 0;
  for (int32_t s = 0; s < nActiveChips; ++s) {
    segmentDigiOffset_h[s] = running;
    running += digiCount_h[s];
  }
  segmentDigiOffset_h[nActiveChips] = running;
  const int32_t nDigisTot = running;
  
  std::cout << "Total number of digis = " << nDigisTot <<std::endl;
  
  // copy offset to device
  alpaka::memcpy(queue, *segmentDigiOffset_, segmentDigiOffset_h);
  
  BTLDigiDeviceCollection digis_d(queue, std::max(nDigisTot, 1));
  
  if (nDigisTot > 0) {
    auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(cms::alpakatools::divide_up_by(uint32_t(nActiveChips), 256u), 256u);
    alpaka::exec<Acc1D>(queue, workDiv, BTLFillDigisKernel{},
			channelPayload_d.const_view(),
			segments.nSegments,
			digiCountPerSegment_->data(),
			segmentDigiOffset_->data(),
			pairPlusLocalIndex_->data(),
			pairMinusLocalIndex_->data(),
			mapping.const_view(),
			indexer,
			digis_d.view());
  }
  
  return digis_d;
  }
  
  
  BTLDigiDeviceCollection BTLRawToDigiAlgo::process(Queue& queue,
                                                    const uint64_t* rawWords_h,
                                                    const int32_t* channelFedId_h,
                                                    int32_t nChannels,
                                                    BTLElectronicsToDetIdMappingDevice const& elecToDetIdMapping,
						    BTLElectronicsIndexer const& indexer) {
    // -- STEP1: Launch decode for each channel
    auto channelsPayload_d = decodeOnly(queue, rawWords_h, channelFedId_h, nChannels, indexer); // auxiliary soa on device

    // -- DEBUG STEP1: copy to host and print
    BTLChannelPayloadHostCollection channelsPayload_h(std::max(nChannels, 1 ));
    alpaka::memcpy(queue, channelsPayload_h.buffer(), channelsPayload_d.buffer());
    alpaka::wait(queue);
    auto soa = channelsPayload_h.view();   
    
    for (int i = 0; i < nChannels; ++i) {
      std::cout << "i = " << i
		<< " fedId = " << soa.fedId()[i]
		<< " hsLinkId = " << static_cast<int>(soa.hsLinkId()[i])
		<< " eLinkId = " << static_cast<int>(soa.eLinkId()[i])
		<< " chId = " << static_cast<int>(soa.chId()[i])
		<< " chipKey = " << soa.chipKey()[i]
		<< std::endl;
    }
    

    // -- STEP1.5: find segments corresponding to each active chip
    auto chipSegments_d = findChipSegments(queue, channelsPayload_d, indexer); 

    // -- DEBUG findChipSegments
    auto nSegments_h = cms::alpakatools::make_host_buffer<int32_t[]>(queue, 1);
    alpaka::memcpy(queue, nSegments_h, *nSegments_d_);
    alpaka::wait(queue);
    const int32_t nActiveChips = nSegments_h[0];

    std::cout << "Number of active chips = " << nActiveChips << std::endl;

    const int32_t nChips = indexer.nFeds * indexer.nHsLinks * indexer.nELinks;
    auto segmentStart_h = cms::alpakatools::make_host_buffer<int32_t[]>(queue, nChips + 1);// alloco con lunghezza nChips
    alpaka::memcpy(queue, segmentStart_h, *segmentStart_d_);
    alpaka::wait(queue);

    std::cout << "Debugging findChipSegments..." <<std::endl; 
    
    for (int32_t s = 0; s < nActiveChips; ++s) { // solo gli elementi corrispondenti al numero di segmenti (chip attivi)
      std::cout << "segment " << s
		<< " : [" << segmentStart_h[s]
		<< ", " << segmentStart_h[s + 1]
		<< ")"
		<< " size=" << segmentStart_h[s + 1] - segmentStart_h[s]
		<< '\n';
    }

    std::cout << "Number of active chips = " << nActiveChips << std::endl;
    
 
    // STEP 2a: count digis
    //
    // NOTE: nSegments is the number of active chips found in STEP 1.5,
    // while nChips is the maximum number of chips.
    //
    //
    // Two options:
    // 1) Use nChips and skip segments >= nSegments inside the kernel.
    // 2) Copy nSegments from device to host and use it to define the
    //    work division directly.
    
    //const int32_t nChips = indexer.nFeds * indexer.nHsLinks * indexer.nELinks;
    //countDigis(queue, channelsPayload_d, chipSegments_d, nChips); // qui uso nChips. Another possibility is to copy chipSegments_d from device to host and use the host-side nSegments
    countDigis(queue, channelsPayload_d, chipSegments_d, elecToDetIdMapping, indexer, nActiveChips);

    
    // -- DEBUG pair and count digis  
    auto digiCount_h = cms::alpakatools::make_host_buffer<int32_t[]>(queue, nActiveChips);
    alpaka::memcpy(queue, digiCount_h, *digiCountPerSegment_);
    alpaka::wait(queue);
    
    int32_t nDigisTot = 0;
    for (int32_t s = 0; s < nActiveChips; ++s) {
      std::cout << "segment " << s
		<< "  number of digis in this segment : " << digiCount_h[s]
		<<std::endl;
      
      nDigisTot+=digiCount_h[s];
    }


    
    
    // -- STEP2b: fill digis
    auto digis_d = fillDigis(queue, channelsPayload_d, chipSegments_d, elecToDetIdMapping, indexer, nActiveChips);


    // -- DEBUGGING FILL DIGIS
    ::btldigi::BTLDigiHostCollection digis_h(queue, std::max(nDigisTot, 1));
    alpaka::memcpy(queue, digis_h.buffer(), digis_d.buffer());
    alpaka::wait(queue);
    auto digis_h_view = digis_h.view();

    /*for (int32_t i = 0; i < nDigisTot; ++i) {
      LogDebug("BTLRawToDigi")
	<< "DIGI " << i
	<< " rawId=" << digis_h_view[i].rawId()
	<< " BC0=" << digis_h_view[i].BC0count()
	<< " status=" << digis_h_view[i].status()
	<< " BCcount=" << digis_h_view[i].BCcount()
	<< " plus(ch=" << int(digis_h_view[i].chIDPlus())
	<< ", charge=" << digis_h_view[i].ChargePlus()
	<< ")"
	<< " minus(ch=" << int(digis_h_view[i].chIDMinus())
	<< ", charge=" << digis_h_view[i].ChargeMinus()
	<< ")";
	}*/


    for (int32_t i = 0; i < nDigisTot; ++i) {
      std::cout 
	<< "DIGI " << i
	<< " rawId=" << digis_h_view[i].rawId()
	<< " BC0=" << digis_h_view[i].BC0count()
	<< " status=" << digis_h_view[i].status()
	<< " BCcount=" << digis_h_view[i].BCcount()
	<< " plus(ch=" << int(digis_h_view[i].chIDPlus())
	<< ", charge=" << digis_h_view[i].ChargePlus()
	<< ")"
	<< " minus(ch=" << int(digis_h_view[i].chIDMinus())
	<< ", charge=" << digis_h_view[i].ChargeMinus()
	<< ")"
	<<std::endl;
    }

    

  
    
    return digis_d;
    
    
  }}
// namespace ALPAKA_ACCELERATOR_NAMESPACE
