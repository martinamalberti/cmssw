#include <vector>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ModuleFactory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

#include "CondFormats/MTDObjects/interface/BTLReadoutMap.h"
#include "CondFormats/DataRecord/interface/BTLReadoutMapRcd.h"
#include "CondFormats/MTDObjects/interface/alpaka/BTLElectronicsToDetIdMappingDevice.h"
#include "EventFilter/MTDRawToDigi/interface/BTLElectronicsSpecs.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // Builds the dense electronics(fedId,hsLinkId,eLinkId,pairIdx) -> rawId
  // Table used by BTLPairChannelsKernel (STEP2). 
  class BTLElectronicsToDetIdMappingESProducer : public ESProducer {
  public:
    explicit BTLElectronicsToDetIdMappingESProducer(const edm::ParameterSet& iConfig) : ESProducer(iConfig) {
      auto cc = setWhatProduced(this);
      readoutMapToken_ = cc.consumes();
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      descriptions.addWithDefaultLabel(desc);
    }

    std::unique_ptr<BTLElectronicsToDetIdMappingHost> produce(const BTLReadoutMapRcd& iRecord) {
      const BTLReadoutMap& readoutMap = iRecord.get(readoutMapToken_);

      BTLElectronicsIndexer indexer;
      indexer.firstFedId = BTLElectronicsSpecs::kFirstFEDId;
      indexer.nFeds = BTLElectronicsSpecs::kNumberOfFEDs;
      indexer.nHsLinks = BTLElectronicsSpecs::kNumberOfHsLinks;  // check real name
      indexer.nELinks = BTLElectronicsSpecs::kNumberOfELinks;    // check real name

      auto product = std::make_unique<BTLElectronicsToDetIdMappingHost>(cms::alpakatools::host(), indexer.size());

      // -- defaults
      for (int32_t i = 0; i < indexer.size(); ++i) {
        product->view()[i].valid() = false;
        product->view()[i].rawId() = 0;
      }

      // -- Fill from the crystals actually present in the readout map.
      int32_t nFilled = 0, nSkippedInconsistent = 0;
      for (const auto& detId : readoutMap.getListOfDetIds()) {  
        BTLElectronicsIdPair elecIds = readoutMap.getElectronicsId(detId);

	// -- check consistent links
        const bool consistentLinks = elecIds.minus.fedId() == elecIds.plus.fedId() &&
                                     elecIds.minus.hsLinkId() == elecIds.plus.hsLinkId() &&
        	                     elecIds.minus.eLinkId() == elecIds.plus.eLinkId();

        if (!consistentLinks) {
          ++nSkippedInconsistent;
          edm::LogWarning("BTLElectronicsToDetIdMappingESProducer")
              << "crystal " << std::hex << detId.rawId() << std::dec << ": electronics mapping inconsistent"
              << " (minus/plus fed/hs/e-link or pair index mismatch) - excluded from device lookup table.";
          continue;
        }

        //const int32_t flat = indexer.flatIndex(elecIds.minus.fedId(), elecIds.minus.hsLinkId(), elecIds.minus.eLinkId(), btlPairIdx(chIdMinus));
	const int32_t idx = indexer.flatIndex(elecIds.minus.fedId(), elecIds.minus.hsLinkId(), elecIds.minus.eLinkId(), detId.crystal());
	product->view()[idx].rawId() = detId.rawId();
	product->view()[idx].valid() = true;
        ++nFilled;
      }

      LogDebug("BTLElectronicsToDetIdMappingESProducer")
          << "filled " << nFilled << " / " << indexer.size() << " table entries, skipped " << nSkippedInconsistent
          << " inconsistent crystals.";

      return product;
    }

  private:
    edm::ESGetToken<BTLReadoutMap, BTLReadoutMapRcd> readoutMapToken_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_EVENTSETUP_ALPAKA_MODULE(BTLElectronicsToDetIdMappingESProducer);
