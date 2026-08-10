#ifndef CondFormats_MTDObjects_alpaka_BTLElectronicsToDetIdMappingDevice_h
#define CondFormats_MTDObjects_alpaka_BTLElectronicsToDetIdMappingDevice_h

#include "CondFormats/MTDObjects/interface/BTLElectronicsToDetIdMappingHost.h"
#include "CondFormats/MTDObjects/interface/BTLElectronicsToDetIdMappingSoA.h"
#include "DataFormats/Portable/interface/PortableCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  using BTLElectronicsToDetIdDeviceCollection = PortableCollection<BTLElectronicsToDetIdSoA>;
}

#endif
