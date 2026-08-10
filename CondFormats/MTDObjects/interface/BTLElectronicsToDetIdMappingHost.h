#ifndef CondFormats_MTDObjects_interface_BTLElectronicsToDetIdMappingHost_h
#define CondFormats_MTDObjects_interface_BTLElectronicsToDetIdMappingHost_h

#include "CondFormats/EcalObjects/interface/BTLElectronicsToDetIdMappingSoA.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

using BTLElectronicsToDetIdMappingHost = PortableHostCollection<BTLElectronicsToDetIdMappingSoA>;

#endif
