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
	// This is needed in STEP 2 which runs one work item per chip, so its `chip` index is exactly this chipKey. 
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
  // KERNEL 1.5:  build the dense per-chip channel lookup. Each thread
  // writes its own channel
  //---------------------------------------------------------------------
  class BTLBuildChannelIndexTableKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
				  BTLChannelPayloadSoA::ConstView channels,
				  int32_t nChannels,
				  int32_t* channelIndexByChip) const {
      for (int32_t i : cms::alpakatools::uniform_elements(acc, nChannels)) {
        const uint32_t chip = channels[i].chipKey();
        const uint8_t chId = channels[i].chId();
        channelIndexByChip[chip * kChannelsPerChip + chId] = i;
      }
    }
  };
  
  
  //---------------------------------------------------------------------
  // KERNEL 2: one thread per chip.  
  // For each crystal belonging to the chip:
  //   1. get the electronics channel IDs for minus/plus from the
  //      electronics -> DetId mapping;
  //   2. look up the corresponding payload rows through
  //      channelIndexByChip;
  //   3. determine which sides are present.
  //
  // channelIndexByChip[chip * 32 + chId] contains the row index in
  // BTLChannelPayloadSoA, or -1 if that channel is absent.
  //---------------------------------------------------------------------
  class BTLPairChannelsKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
				  BTLChannelPayloadSoA::ConstView channels,
				  int32_t nChips,
				  const int32_t* channelIndexByChip,
				  BTLElectronicsToDetIdMappingSoA::ConstView elecToDetId,
				  BTLElectronicsIndexer indexer,
				  ::btldigi::BTLDigiSoA::View digisOut) const {

      for (int32_t chip : cms::alpakatools::uniform_elements(acc, nChips)) {

	// Dense channel lookup for this chip.
	const int32_t channelBase = chip * kChannelsPerChip;

	const int32_t eLinkId = chip % indexer.nELinks;
	const int32_t tmp = chip / indexer.nELinks;

	const int32_t hsLinkId = tmp % indexer.nHsLinks;
	const int32_t fed = tmp / indexer.nHsLinks;

	const int32_t fedId = indexer.firstFedId + fed;
	const int32_t hsLink = indexer.hsLinkOffset + hsLinkId;

	for (int32_t pairIdx = 0; pairIdx < BTLElectronicsIndexer::nPairs; ++pairIdx) { // at most 16 crystals

	  const int32_t mapIndex = indexer.flatIndex(fedId, hsLink, eLinkId, pairIdx);
	  const auto& mapping = elecToDetId[mapIndex];

	  if (!mapping.valid())
	    continue;

	  const uint8_t minusCh = mapping.minusChannelId();
	  const uint8_t plusCh = mapping.plusChannelId();

	  const int32_t minusIndex = channelIndexByChip[channelBase + minusCh];

	  const int32_t plusIndex = channelIndexByChip[channelBase + plusCh];
	  
	  const bool hasMinus = minusIndex >= 0;
	  const bool hasPlus = plusIndex >= 0;

	  // now fill the digis

	  if (hasMinus && hasPlus){
	    
	  }
	  
	}
  
      }	
    }
  };
  

  void BTLRawToDigiAlgo::debugChannelIndexTable(Queue& queue,
						BTLChannelIndexTable const& table,
						BTLChannelPayloadDeviceCollection const& channelsPayload_d,
						BTLElectronicsIndexer const& indexer) const {
    
    const int32_t nChannels = channelsPayload_d->metadata().size();
    const int32_t nChips = indexer.nFeds * indexer.nHsLinks * indexer.nELinks;
    const int32_t tableSize = nChips * kChannelsPerChip;
    
    // --------------------------------------------------------------
    // Copy lookup table device -> host
    // --------------------------------------------------------------
    auto table_h = cms::alpakatools::make_host_buffer<int32_t[]>(queue, std::max(tableSize, 1));
    alpaka::memcpy(queue, table_h, table);
    alpaka::wait(queue);
    
    // --------------------------------------------------------------
    // Copy payload device -> host
    // --------------------------------------------------------------
    BTLChannelPayloadHostCollection channelsPayload_h(std::max(nChannels, 1));
    alpaka::memcpy(queue, channelsPayload_h.buffer(), channelsPayload_d.buffer());
    alpaka::wait(queue);
    auto channels = channelsPayload_h.view();
    
    // --------------------------------------------------------------
    // Print table: one line per occupied chip
    // Example:  chip 42: [ 15,  -1,  17,  -1, ... ]
    // The value is the payload index i. -1 if channel not present
    // --------------------------------------------------------------
    std::cout << "\n===== BTL channel index table =====\n";
    
    for (int32_t chip = 0; chip < nChips; ++chip) {
      bool chipHasChannels = false;
      for (int32_t ch = 0; ch < kChannelsPerChip; ++ch) {
	const int32_t i = table_h[chip * kChannelsPerChip + ch];
	
	if (i >= 0) {
	  chipHasChannels = true;
	  break;
	}
      }
      
      if (!chipHasChannels)
	continue;
      
      std::cout << "chip " << chip << ": [";
      
      for (int32_t ch = 0; ch < kChannelsPerChip; ++ch) {
	const int32_t i =
          table_h[chip * kChannelsPerChip + ch];
	
	if (ch != 0)
	  std::cout << ", ";
	
	std::cout << std::setw(4) << i;
      }
      
      std::cout << " ]\n";
    }
    
    // --------------------------------------------------------------
    // Validate - for every payload i: table[chipKey(i) * 32 + chId(i)] == i
    // --------------------------------------------------------------
    bool tableOK = true;
    
    for (int32_t i = 0; i < nChannels; ++i) {
      const uint32_t chip = channels.chipKey()[i];
      const uint8_t chId = channels.chId()[i];
      const int32_t tableIndex = table_h[chip * kChannelsPerChip + chId];
      
      if (tableIndex != i) {
	
	std::cout << "ERROR: table mismatch:"
		  << " i = " << i
		  << " chip = " << chip
		  << " chId = " << static_cast<int>(chId)
		  << " tableIndex = " << tableIndex
		  << '\n';
	
	tableOK = false;
      }
    }
    
    if (tableOK) {
      std::cout << "===== BTL channel index table: OK =====\n";
    } else {
      std::cout << "===== BTL channel index table: FAILED =====\n";
    }
  }
  
  //---------------------------------------------------------------------
  // BTLRawToDigiAlgo
  //---------------------------------------------------------------------
  BTLChannelPayloadDeviceCollection BTLRawToDigiAlgo::decodeOnly(Queue& queue,
                                                                 const uint64_t* rawWords_h,
                                                                 const int32_t* channelFedId_h,
                                                                 int32_t nChannels,
                                                                 BTLElectronicsIndexer const& indexer) const {
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


  BTLChannelIndexTable BTLRawToDigiAlgo::buildChannelIndexTable(Queue& queue,
								BTLChannelPayloadDeviceCollection const& channelPayload_d,
								BTLElectronicsIndexer const& indexer) const {
    const int32_t nChannels = channelPayload_d->metadata().size();
    const int32_t nChips = indexer.nFeds * indexer.nHsLinks * indexer.nELinks;
    
    BTLChannelIndexTable table = cms::alpakatools::make_device_buffer<int32_t[]>(queue, std::max(nChips * kChannelsPerChip, 1));

    // Initialize all entries to -1:
    // -1 means that this channel was not present in the event.
    alpaka::memset(queue, table, 0xFF);

    if (nChannels > 0) {
      auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(cms::alpakatools::divide_up_by(uint32_t(nChannels), 256u), 256u);
      alpaka::exec<Acc1D>(queue,
			  workDiv,
			  BTLBuildChannelIndexTableKernel{},
			  channelPayload_d.const_view(),
			  nChannels,
			  table.data());
    }
    
    return table;

  }


  BTLDigiDeviceCollection BTLRawToDigiAlgo::pairChannels(Queue& queue,
							 BTLChannelPayloadDeviceCollection const& channelsPayload_d,
							 BTLChannelIndexTable const& table,
							 BTLElectronicsIndexer const& indexer,
							 BTLElectronicsToDetIdMappingDevice const& elecToDetId) const {
    
    const int32_t nChips = indexer.nFeds * indexer.nHsLinks * indexer.nELinks;
    const int32_t maxDigis = nChips * 32;// temporaneamente mettiamo max, poi puo'
    //const int32_t maxDigis = nChannels; 
    
    BTLDigiDeviceCollection digis_d(queue, std::max(maxDigis, 1));
    
    if (nChips > 0) {
      auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(cms::alpakatools::divide_up_by(uint32_t(nChips), 256u), 256u);
      alpaka::exec<Acc1D>(queue,
			  workDiv,
			  BTLPairChannelsKernel{},
			  channelsPayload_d.const_view(),
			  nChips,
			  table.data(),
			  elecToDetId.const_view(),
			  indexer,
			  digis_d.view());
    }

    return digis_d;
  }

  
  BTLDigiDeviceCollection BTLRawToDigiAlgo::process(Queue& queue,
                                                    const uint64_t* rawWords_h,
                                                    const int32_t* channelFedId_h,
                                                    int32_t nChannels,
                                                    BTLElectronicsIndexer const& indexer,
                                                    BTLElectronicsToDetIdMappingDevice const& elecToDetId) const {
    // -- STEP1: Launch decode for each channel
    auto channelsPayload_d = decodeOnly(queue, rawWords_h, channelFedId_h, nChannels, indexer);

    // -- DEBUG STEP1
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
    

    // -- STEP1.5: Build the dense per-chip channel lookup
    auto table = buildChannelIndexTable(queue, channelsPayload_d, indexer);

    // DEBUG STEP1.5:
    debugChannelIndexTable(queue, table, channelsPayload_d, indexer);
    
    // -- STEP2: Launch pairing and fill digis
    auto digis_d = pairChannels(queue, channelsPayload_d, table, indexer, elecToDetId);
    return digis_d;
    
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
