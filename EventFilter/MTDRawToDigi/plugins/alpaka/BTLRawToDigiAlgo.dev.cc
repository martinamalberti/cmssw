#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/OneToManyAssoc.h"

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
	out[i].chipKey() = static_cast<uint32_t>(indexer.flatIndex(fed, hslink, elink, 0) / BTLElectronicsIndexer::nPairs); // flatIndex DA CONTROLLARE
      }
    }
  };
    


  //---------------------------------------------------------------------
  // KERNEL 1.5: fill the chip -> [channel indices] association.
  //---------------------------------------------------------------------
    
  // DA IMPLEMENTARE


    
  //---------------------------------------------------------------------
  // KERNEL 2: one block per occupied chip bucket. Threads 0..count-1 load
  // their channel into shared memory, then each thread checks whether it
  // is the "minus" side of its pair and, if so, looks for the partner
  // chId within the (<=32-element) shared block and fills one digi.
  //---------------------------------------------------------------------
  
  /*
    class BTLPairChannelsKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
				  BTLChannelPayloadSoA::ConstView channels,
				  ChipAssoc::ConstView chipBuckets,  // bucket -> channel indices
				  BTLElectronicsToDetIdSoA::ConstView elecToDetId,
				  BTLElectronicsIndexer indexer,
				  BTLDigiSoA::View digisOut) const {      
      // pairing and filling digis
      // ...
      // ...
    };
  };
  */

  //---------------------------------------------------------------------
  // BTLRawToDigiAlgo
  //---------------------------------------------------------------------
  BTLChannelPayloadDeviceCollection BTLRawToDigiAlgo::decodeOnly(Queue& queue,
								 const uint64_t* rawWords_h,
								 const int32_t* channelFedId_h,
								 int32_t nChannels,
								 BTLElectronicsIndexer const& indexer) const {
    
    // Copy host --> device
    auto rawWords_d = cms::alpakatools::make_device_buffer<uint64_t[]>(queue, 2 * std::max(nChannels, 1)); // alloca memoria su gpu
    auto channelFedId_d = cms::alpakatools::make_device_buffer<int32_t[]>(queue, std::max(nChannels, 1));
    if (nChannels > 0) {
      alpaka::memcpy(queue, rawWords_d, cms::alpakatools::make_host_view(rawWords_h, 2 * nChannels));
      alpaka::memcpy(queue, channelFedId_d, cms::alpakatools::make_host_view(channelFedId_h, nChannels));
    }
    
    BTLChannelPayloadDeviceCollection channelsPayload_d(queue, std::max(nChannels, 1));
    if (nChannels > 0) {
      auto workDiv = cms::alpakatools::make_workdiv<Acc1D>(cms::alpakatools::divide_up_by(uint32_t(nChannels), 256u), 256u);
      alpaka::exec<Acc1D>(queue, workDiv, BTLDecodeChannelsKernel{}, rawWords_d.data(), channelFedId_d.data(),
			  nChannels, indexer, channelsPayload_d.view());
    }
    return channelsPayload_d;
  }
  

  BTLDigiDeviceCollection BTLRawToDigiAlgo::process(Queue& queue,
						    const uint64_t* rawWords_h,
						    const int32_t* channelFedId_h,
						    int32_t nChannels,
						    BTLElectronicsIndexer const& indexer,
						    BTLElectronicsToDetIdMappingDevice const& elecToDetId) const {
    
    // launch decode for each channel
    auto channels_d = decodeOnly(queue, rawWords_h, channelFedId_h, nChannels, indexer);
    
    // launch bucket by chip
    //alpaka::exec(..., BTLBuildChipAssocKernel{}, ...);
    
    // launch pairing and fill digis
    //alpaka::exec(..., BTLPairChannelsKernel{}, ...);


    // temporaneo 
    return BTLDigiDeviceCollection(queue, 0);

    //return digis_d;
  }
  
  
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

