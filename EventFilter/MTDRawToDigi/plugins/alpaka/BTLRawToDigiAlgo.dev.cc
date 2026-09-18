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
            static_cast<uint32_t>(indexer.flatIndex(fed, hslink, elink, 0) / BTLElectronicsIndexer::nPairs);

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
				  BTLChannelPayloadSoA::ConstView channels,
				  int32_t* segmentStart,
				  int32_t* nSegmentsOut) const {
      const int32_t nChannels = channels.metadata().size(); 
      // For now just one thread, to detect chip boundaries. To be checked if performance is acceptable.
      if (cms::alpakatools::once_per_grid(acc)) { // capire meglio cosa fa
	int32_t nSegments = 0;
	for (int32_t i = 0; i < nChannels; ++i) {
	  if (i == 0 || channels[i].chipKey() != channels[i - 1].chipKey()) {
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
  class BTLCountDigisKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
				  BTLChannelPayloadSoA::ConstView channels,
				  const int32_t* nSegments,
				  const int32_t* segmentStart,
				  int32_t nChips,
				  //int32_t* pairFirstIndex,
				  //int32_t* pairSecondIndex,
				  int32_t* digiCount) const {
      
      for (int32_t segment : cms::alpakatools::uniform_elements(acc, nChips)) { //uniform_elements() distribuisce gli indici degli elementi tra le diverse esecuzioni parallele del kernel (i work-item).
	
	if (segment >= *nSegments) continue;
	
	const int32_t start = segmentStart[segment];
	const int32_t end = segmentStart[segment + 1];
	
	int32_t nDigis = 0;
	
	for (int32_t i = start; i < end; ++i) {
	  const uint8_t chId = channels[i].chId();
	  
	  const int8_t partnerChId = BTLElectronicsSpecs::kTofhirChannelPartner[chId];
	  
	  // Process each electronics pair only once.
	  if (static_cast<int32_t>(chId) > partnerChId)
	    continue;
	  
	  //int32_t partnerRow = -1;
	  
	  for (int32_t j = start; j < end; ++j) {
            if (j == i) continue;
            if (channels[j].chId() == partnerChId) {
	      //int32_t partnerRow = -1;
	      break;
	    }
          }
	  
	  
	  //const int32_t pairIndex = pairBase + nPairs;
	  //pairFirstIndex[pairIndex] = i;
	  //pairSecondIndex[pairIndex] = partnerRow;
	  
	  // A digi is produced even if the partner channel is missing.
	  ++nDigis;
	}
	
	digiCount[segment] = nDigis;
      }
    }
  };


  //---------------------------------------------------------------------
  // KERNEL 2b: one thread per chip: buil pairs and fill digis    
  //
  // TODO

  
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
    }
    return channelsPayload_d;
  }



  BTLChipSegments BTLRawToDigiAlgo::findChipSegments(Queue& queue,
						     BTLChannelPayloadDeviceCollection const& channelPayload_d,
						     BTLElectronicsIndexer const& indexer) {
  
    const int32_t nChips = indexer.nFeds * indexer.nHsLinks * indexer.nELinks;
  
    if (!segmentStart_->data()) {
      segmentStart_.emplace(cms::alpakatools::make_device_buffer<int32_t[]>(queue, nChips + 1));
    }
    if (!nSegments_->data()) {
      nSegments_.emplace(cms::alpakatools::make_device_buffer<int32_t[]>(queue, 1));
    }
    
    auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(1, 1);
    alpaka::exec<Acc1D>(queue, workDiv, BTLFindChipSegmentsKernel{}, channelPayload_d.const_view(),
			segmentStart_->data(), nSegments_->data());
    
    return {segmentStart_->data(), nSegments_->data()};
  }
  

  
  void BTLRawToDigiAlgo::countDigis(Queue& queue,
				    BTLChannelPayloadDeviceCollection const& channelPayload_d,
				    BTLChipSegments const& segments,
				    int32_t nChips) {
    
    if (!digiCount_) {
      digiCount_.emplace(cms::alpakatools::make_device_buffer<int32_t[]>(queue, nChips)); // se voglio dimensionarlo nSegments invece che nChips dovrei conoscere nSegments lato host, pero' implica una sincronizzazione device-host.
    }
    
    auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(nChips, 1); // se volessi solo nSegments work-items dorvei copiare da device a host nsegments, perche' workDiv vuole sapere quanti work-item lato host
    
    alpaka::exec<Acc1D>(queue,
			workDiv,
			BTLCountDigisKernel{},
			channelPayload_d.const_view(),
			segments.nSegments,
			segments.segmentStart,
			nChips,
			digiCount_->data());
  }
  
  
  BTLDigiDeviceCollection BTLRawToDigiAlgo::process(Queue& queue,
                                                    const uint64_t* rawWords_h,
                                                    const int32_t* channelFedId_h,
                                                    int32_t nChannels,
                                                    BTLElectronicsIndexer const& indexer,
                                                    BTLElectronicsToDetIdMappingDevice const& elecToDetId) {
    // -- STEP1: Launch decode for each channel
    auto channelsPayload_d = decodeOnly(queue, rawWords_h, channelFedId_h, nChannels, indexer);

    // -- DEBUG STEP1: copy to host
    BTLChannelPayloadHostCollection channelsPayload_h( std::max(nChannels, 1 ));
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
    auto chipSegments = findChipSegments(queue, channelsPayload_d, indexer); 


    // STEP 2a: count digis
    //
    // NOTE: nSegments is the number of active chips found in STEP 1.5,
    // while nChips is the maximum number of chips.
    //
    // Here I use nChips to define the work division, because nChips is
    // known on the host. In principle, I would like to launch only
    // nSegments work-items.
    //
    // Two options:
    // 1) Use nChips and skip segments >= nSegments inside the kernel.
    // 2) Copy nSegments from device to host and use it to define the
    //    work division directly.
    const int32_t nChips = indexer.nFeds * indexer.nHsLinks * indexer.nELinks;
    countDigis(queue, channelsPayload_d, chipSegments, nChips); // qui uso nChips. Another possibility is to copy chipSegments from device to host and use the host-side nSegments. 
    
    // STEP2b: pair and fill digis
    // TO DO Fill digis
    
    
    // -- STEP2: Launch pairing and fill digis
    //auto digis_d = pairChannels(queue, channelsPayload_d, table, indexer, elecToDetId);
    BTLDigiDeviceCollection digis_d{queue, 0}; /// temp
    return digis_d;

  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
