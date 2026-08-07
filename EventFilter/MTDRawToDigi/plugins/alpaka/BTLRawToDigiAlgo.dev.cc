#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
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
      
      // !!!!!!!!!!!!!!! FIXME: Update bit positions once the final BTL payload format is frozen. !!!!!!!!!!!!!!!!!!
      out[i].bc0count() = static_cast<uint16_t>(extractBits(118, 10));
      out[i].status() = extractBits(117, 1) != 0;
      const int16_t hslink = static_cast<int16_t>(extractBits(104, 6));
      const int16_t elink = static_cast<int16_t>(extractBits(99, 5));
      out[i].hsLinkId() = hslink;
      out[i].eLinkId() = elink;
      out[i].chId() = static_cast<uint8_t>(extractBits(94, 5));
      out[i].bcCount() = static_cast<uint32_t>(extractBits(82, 12));
      out[i].t1Coarse() = static_cast<uint16_t>(extractBits(67, 15));
      out[i].t2Coarse() = static_cast<uint16_t>(extractBits(57, 10));
      out[i].eoiCoarse() = static_cast<uint16_t>(extractBits(47, 10));
      out[i].charge() = static_cast<uint16_t>(extractBits(37, 10));
      out[i].t1Fine() = static_cast<uint16_t>(extractBits(27, 10));
      out[i].t2Fine() = static_cast<uint16_t>(extractBits(17, 10));
      out[i].idleTime() = static_cast<uint16_t>(extractBits(7, 10));
      out[i].prevTrigF() = static_cast<uint8_t>(extractBits(3, 4));
      out[i].tacId() = static_cast<uint8_t>(extractBits(0, 3));
      
      const int32_t fed = channelFedId[i];
      out[i].fedId() = fed;
      out[i].chipKey() = static_cast<uint32_t>(indexer.flatIndex(fed, hslink, elink, 0) / BTLElectronicsIndexer::nPairs); // flatIndex DA IMPLEMENTARE
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
  // chId within the (<=32-element) shared block and emits one digi.
  //---------------------------------------------------------------------
  
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
    
    
  BTLDigiDeviceCollection BTLRawToDigiAlgo::process(Queue& queue,
                                    const uint64_t* rawWords_h,
                                    const int32_t* channelFedId_h,
                                    int32_t nChannels,
                                    BTLElectronicsIndexer const& indexer,
                                    BTLElectronicsToDetIdDeviceCollection const& elecToDetId) {
    
    // Move the data to the device
    
    // launch decode
    alpaka::exec(..., BTLDecodeChannelsKernel{}, ...);
    
    // launch bucket
    alpaka::exec(..., BTLBuildChipAssocKernel{}, ...);
    
    // launch pairing and fill digis
    alpaka::exec(..., BTLPairChannelsKernel{}, ...);
    
    return digis;
  }
  
  
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
