#include <cstdio>

#include <alpaka/alpaka.hpp>

#include "DataFormats/FTLRecHitSoA/interface/alpaka/BTLBaseRecHitDeviceCollection.h"
#include "DataFormats/FTLDigiSoA/interface/alpaka/BTLDigiDeviceCollection.h"

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "DataFormats/ForwardDetId/interface/BTLDetId.h"

#include "BTLBaseRecHitSoAProducerAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::btlrechit {

  using namespace ::btlrechit;

  ALPAKA_FN_ACC uint8_t rowFromId(uint32_t rawId) {  // NB working only with new geometry
    int crys = ((rawId >> BTLDetId::kBTLCrystalOffset) & BTLDetId::kBTLCrystalMask);
    uint8_t row = crys % BTLDetId::kCrystalsPerModuleV2;
    return row;
  }

  ALPAKA_FN_ACC float TcoarseTfineToTime(
      uint32_t rawId, uint8_t chID, uint8_t TACID, uint16_t tcoarse, uint16_t tfine, bool isT1) {
    // tdc calibration parameters
    // (to be modified: these parameters are evaluated by channel and stored in parquet files)
    static constexpr float a0 = 57.244545;
    static constexpr float a1 = 511.27832;
    static constexpr float a2 = -7.8838577;
    static constexpr float t0 = -0.048343264;

    float const qT = (-a1 + sqrt(a1 * a1 - 4.0 * (a0 - float(tfine)) * a2)) / (2.0 * a2);
    float const time = tcoarse - qT - t0;
    return time;
  }

  ALPAKA_FN_ACC uint32_t
  QfineToADC(uint32_t rawId, uint8_t chID, uint8_t TACID, uint16_t qfine, float time1, uint16_t timeEndQ) {
    // qdc calibration parameters
    // (to be modified: these parameters are evaluated by channel and stored in parquet files)
    static constexpr float p0 = 49.542229;
    static constexpr float p1 = -0.323424;
    static constexpr float p2 = 0.062578;
    static constexpr float p3 = -0.002484;
    static constexpr float p4 = 0.0;
    static constexpr float p5 = 0.0;
    static constexpr float p6 = 0.0;
    static constexpr float p7 = 0.0;
    static constexpr float p8 = 0.0;
    static constexpr float p9 = 0.0;
    float const ti = float(timeEndQ) - time1;

    uint32_t pedestal = (  // check the type
        p0 + p1 * ti + p2 * ti * ti + p3 * ti * ti * ti + p4 * ti * ti * ti * ti + p5 * ti * ti * ti * ti * ti +
        p6 * ti * ti * ti * ti * ti * ti + p7 * ti * ti * ti * ti * ti * ti * ti +
        p8 * ti * ti * ti * ti * ti * ti * ti * ti + p9 * ti * ti * ti * ti * ti * ti * ti * ti * ti);

    const uint32_t adc = qfine - pedestal;

    return adc;
  }

  class BTLdigiToBaseKernel {
  public:
    ALPAKA_FN_ACC void operator()(
        Acc1D const& acc,
        ::btldigi::BTLDigiSoA::ConstView input,
        BTLBaseRecHitSoA::View output,
        const double npeToADC0_,
        const double npeToADC1_) const {

      static constexpr uint32_t adcBitSaturation_ = 1023;
      static constexpr float tclock = 6.25;
      // make a strided loop over the kernel grid, covering up to "size" elements
      for (int32_t i : cms::alpakatools::uniform_elements(acc, input.metadata().size())) {
        auto entry = input[i];
        // here you should call your functions to apply TDC and QDC, reading timecoarse, fine, ... etc from digi

        // for the times at first and second th, still in clock units
        // atm tdc and qdc calibs are fixed to dummy values for each channel, hence rawId, ch, and the bool to select branch 1 or 2 are not used.
        auto time1Plus =
            TcoarseTfineToTime(entry.rawId(), entry.chIDPlus(), entry.TACIDPlus(), entry.T1coarsePlus(), entry.T1finePlus(), true);
        auto time1Minus =
            TcoarseTfineToTime(entry.rawId(), entry.chIDMinus(), entry.TACIDMinus(), entry.T1coarseMinus(), entry.T1fineMinus(), true);

        auto time2Plus =
            TcoarseTfineToTime(entry.rawId(), entry.chIDPlus(), entry.TACIDPlus(), entry.T2coarsePlus(), entry.T2finePlus(), false);
        auto time2Minus =
            TcoarseTfineToTime(entry.rawId(), entry.chIDMinus(), entry.TACIDMinus(), entry.T2coarseMinus(), entry.T2fineMinus(), false);

        // from qfine to energy in adc, NB you need to pass calibrated time
        auto ampMinus =
            QfineToADC(entry.rawId(), entry.chIDMinus(), entry.TACIDMinus(), entry.ChargeMinus(), time1Minus, entry.EOIcoarseMinus());
        auto ampPlus =
            QfineToADC(entry.rawId(), entry.chIDPlus(), entry.TACIDPlus(), entry.ChargePlus(), time1Plus, entry.EOIcoarsePlus());

        uint8_t row = rowFromId(entry.rawId());

        // flags for the usability of the channel uint_8: atm 2 bit are used
        //  first bit is channel has signal (1) or not (0)
        //  second bit channel was saturated (1) or not (0)
        uint8_t flagsMinus = 0;
        uint8_t flagsPlus = 0;

        if (ampMinus > 0)
          flagsMinus |= 0x1;
        if (ampMinus == adcBitSaturation_)
          flagsMinus |= (0x1 << 1);
        if (ampPlus > 0)
          flagsPlus |= 0x1;
        if (ampPlus == adcBitSaturation_)
          flagsPlus |= (0x1 << 1);

        // detId from rawId
        DetId detId(entry.rawId());

        // convert from clock units to ps
        time1Plus *= tclock;
        time1Minus *= tclock;
        time2Plus *= tclock;
        time2Minus *= tclock;

        // converting the energy from ADC to energy
        auto energyPlus = float((float(ampPlus) - npeToADC0_) / npeToADC1_);
        auto energyMinus = float((float(ampMinus) - npeToADC0_) / npeToADC1_);

#ifdef EDM_ML_DEBUG
        printf("Base recHit SoA with raw id %i \n", entry.rawId());
        printf("Time 1 before corrections -,+ (%f, %f) - ", time1Minus, time1Plus);
        printf("Time 2 before corrections -,+ (%f, %f) - ", time2Minus, time2Plus);

        printf("Amplidute in ADC -,+ (%i, %i) - ", ampMinus, ampPlus);
        printf("converting to energy -,+ (%f, %f) --> ", npeToADC0_, invADCPerMeV_);
        printf("Energy in MeV -,+ (%f, %f) \n", energyMinus, energyPlus);
#endif

        // fill the base rechit
        output[i] = {
            detId,
            row,
            time1Plus,  // in ns
            time2Plus,
            energyPlus,  // energy
            entry.IdleTimePlus(),
            flagsPlus,
            time1Minus,  // in ns
            time2Minus,
            energyMinus,  // energy
            entry.IdleTimeMinus(),
            flagsMinus,

        };
      }
    }
  };

  void BTLBaseRecHitSoAProducerAlgo::fromDigiToBase(Queue& queue,
                                                    ::btldigi::BTLDigiSoA::ConstView const& input,
                                                    BTLBaseRecHitSoA::View& output,
                                                    const double npeToADC0_,
                                                    const double npeToADC1_) {
    //,
    //Table const& tdc,
    //Table const& qdc) {
    // Use 64 items per group.
    // This value is arbitrary, but it's a reasonable starting point.
    uint32_t items = 64;

    // Use as many groups as needed to cover the whole problem.
    // If this value is too large, a smaller number of blocks can give better performance.
    uint32_t groups = cms::alpakatools::divide_up_by(input.metadata().size(), items);

    auto grid = cms::alpakatools::make_workdiv<Acc1D>(groups, items);
    alpaka::exec<Acc1D>(queue,
                        grid,
                        BTLdigiToBaseKernel{},
                        input,
                        output,
                        npeToADC0_,
                        npeToADC1_);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::btlrechit
