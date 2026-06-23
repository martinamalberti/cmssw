#include <alpaka/alpaka.hpp>

#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

#include "DataFormats/FTLRecHitSoA/interface/BTLBaseRecHitSoA.h"
#include "DataFormats/FTLRecHitSoA/interface/BTLRecHitSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "RecoLocalFastTime/Records/interface/MTDTimeCalibRecord.h"
#include "RecoLocalFastTime/FTLCommonAlgos/interface/MTDTimeCalib.h"

#include "BTLRecHitSoAProducerAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::btlrechit {

  using namespace ::btlrechit;
  ALPAKA_FN_ACC float timeResolutionInNs(float amp) { return 0.0593858 * pow(amp, -1.02826) + 0.0156719; }

  ALPAKA_FN_ACC float getTimeCalib() { return 0.25; }

  ALPAKA_FN_ACC float timeWalkCorr(float amp) {
    // taken from SLHCUpgradeSimulations/Configuration/python/aging.py
    // for 1000 fb-1 scenario
    return 2.40863 * pow(amp,-0.583148) + 0.0545682;
  }

  class BTLBaseToRecoKernel {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  BTLBaseRecHitSoA::ConstView input,
                                  BTLRecHitSoA::View output,
                                  const double c_LYSO_,
                                  const double thresholdToKeep_,
                                  const double calibration_,
                                  const double npeSaturationCorr0_,
                                  const double npeSaturationCorr1_,
                                  const double npePerGeV_) const {  // when condformat for calib ready, add also tdc and qdc in inputs
      // make a strided loop over the kernel grid, covering up to "size" elements

      for (int32_t i : cms::alpakatools::uniform_elements(acc, input.metadata().size())) {
        auto entry = input[i];
        float time1R = entry.time1R();
        float time1L = entry.time1L();
        float time2R = entry.time2R();
        float time2L = entry.time2L();
        float ampR = entry.ampR();
        float ampL = entry.ampL();

        // Apply time and energy corrections
        //   apply amp walk corrections
        auto corrR = timeWalkCorr(ampR);
        auto corrL = timeWalkCorr(ampL);
        time1R = time1R - corrR;
        time1L = time1L - corrL;
        time2R = time2R - corrR;
        time2L = time2L - corrL;

        //   correction for SiPM saturation (just invert the function used to model this effect in BTLElectronicsSim)
        float dR = npeSaturationCorr1_ * npeSaturationCorr1_ + 4. * npeSaturationCorr0_ * ampR;
        ampR = (-npeSaturationCorr1_ + sqrt(dR)) / (2. * (npeSaturationCorr0_));
        ampR /= npePerGeV_;
        float dL = npeSaturationCorr1_ * npeSaturationCorr1_ + 4. * npeSaturationCorr0_ * ampL;
        ampL = (-npeSaturationCorr1_ + sqrt(dL)) / (2. * (npeSaturationCorr0_));
        ampL /= npePerGeV_;

        float time1 = 0;
        float time2 = 0;
        float position = -1.;
        float position_error = -1.;
        float time_error = 0;
        float energy = 0;
        uint8_t flag = 0;

        //!!!!!!! time error calculation to be added
        //!!!!!!! position error calculation to be added

        // -- if you have both sipm info and they are not saturated
        if (entry.flagsR() == 0x1 && entry.flagsL() == 0x1) {
          time1 = 0.5f * (time1L + time1R);
          time2 = 0.5f * (time2L + time2R);  // to be discussed
          position = 0.5f * c_LYSO_ * (time1L - time1R);
          position_error = 0.6;  // as in the std btl uncalibrated hit producer
          energy = (ampR + ampL) / 2.;
          flag |= 0x3;

        }
        // --- If only one SiPM has good not saturated signal
        else if (entry.flagsL() == 0x1 && (time1R == 0x3 || time1R == 0)) {
          time1 = time1L;
          time2 = time2L;
          energy = ampL;
          flag |= (0x1 << 1);
        }

        else if (entry.flagsR() == 0x1 && (entry.flagsL() == 0x3 || entry.flagsR() == 0)) {
          time1 = time1R;
          time2 = time2R;
          energy = ampR;
          flag |= 0x1;
        }

        // energy calibration
        energy *= calibration_;

        // --- Time calibration: for the time being just removes a time offset in BTL
        time1 -= getTimeCalib();

        time_error = timeResolutionInNs(energy);

        // Now fill flags
        // good is 1--> 2 channels && over threshold, bad is 0
        if (energy > thresholdToKeep_ && flag == 0x3) {
          flag = 1;
        } else {
          flag = 0;
        }

#ifdef EDM_ML_DEBUG

        printf("RecHit SoA with raw id %i \n", entry.detId().rawId());
        printf(
            "Time 1  L,R (%f, %f) and average, error (%f, %f) \n", time1L, time1R, time1, time_error);
        printf("Time 2  L,R (%f, %f) and average %f \n", time2L, time2R, time2);
        printf("Energy  L,R (%f, %f) and average %f \n", ampL, ampR, energy);
        printf("Position and error (%f, %f) \n", position, position_error);

#endif

        // fill the rechit
        output[i] = {entry.detId(),
                     entry.row(),  //dummy
                     time1,
                     time2,
                     energy,
                     position,
                     time_error,
                     position_error,  // dummy
                     flag};
      }
    }
  };

  void BTLRecHitSoAProducerAlgo::fromBaseToReco(Queue& queue,
                                                BTLBaseRecHitSoA::ConstView const& input,
                                                BTLRecHitSoA::View& output,
                                                const double c_LYSO_,
                                                const double thresholdToKeep_,
                                                const double calibration_,
                                                const double npeSaturationCorr0_,
                                                const double npeSaturationCorr1_,
                                                const double npePerGeV_) {
    // Use 64 items per group.
    // This value is arbitrary, but it's a reasonable starting point.
    uint32_t items = 64;

    // Use as many groups as needed to cover the whole problem.
    // If this value is too large, a smaller number of blocks can give better performance.
    uint32_t groups = cms::alpakatools::divide_up_by(input.metadata().size(), items);

    auto grid = cms::alpakatools::make_workdiv<Acc1D>(groups, items);
    alpaka::exec<Acc1D>(queue,
                        grid,
                        BTLBaseToRecoKernel{},
                        input,
                        output,
                        c_LYSO_,
                        thresholdToKeep_,
                        calibration_,
                        npeSaturationCorr0_,
                        npeSaturationCorr1_,
                        npePerGeV_);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::btlrechit
