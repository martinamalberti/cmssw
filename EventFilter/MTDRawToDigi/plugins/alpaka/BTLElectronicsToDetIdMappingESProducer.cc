#include <vector>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ModuleFactory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

#include "CondFormats/MTDObjects/interface/BTLReadoutMap.h"
#include "CondFormats/DataRecord/interface/BTLReadoutMapRcd.h"
#include "CondFormats/MTDObjects/interface/alpaka/BTLElectronicsToDetIdMappingDevice.h"
#include "CondFormats/MTDObjects/interface/BTLElectronicsSpecs.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // Builds the dense electronics(fedId,hsLinkId,eLinkId,channelId) -> rawId
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
      indexer.nFEDsPerSerenity = BTLElectronicsSpecs::nFEDsPerSerenity;
      indexer.nHSLinksPerSerenity = BTLElectronicsSpecs::nHSLinksPerSerenity;
      indexer.nELinksPerHSLink = BTLElectronicsSpecs::nELinksPerHSLink;
      indexer.nChannelsPerELink = BTLElectronicsSpecs::nChannelsPerELink;

      auto product = std::make_unique<BTLElectronicsToDetIdMappingHost>(cms::alpakatools::host(), indexer.size());

      // -- defaults
      for (int32_t i = 0; i < indexer.size(); ++i) {
        product->view()[i].valid() = false;
        product->view()[i].rawId() = 0;
        product->view()[i].partnerChannelId() = 0;
        product->view()[i].side() = 0;
      }

      // -- Fill from the crystals actually present in the readout map.
      int32_t nFilled = 0;
      for (const auto& detId : readoutMap.getListOfDetIds()) {
        BTLElectronicsIdPair elecIds = readoutMap.getElectronicsId(detId);

        // -- check consistent links
        const bool consistentLinks = elecIds.minus.fedId() == elecIds.plus.fedId() &&
                                     elecIds.minus.hsLinkId() == elecIds.plus.hsLinkId() &&
                                     elecIds.minus.eLinkId() == elecIds.plus.eLinkId();

        if (!consistentLinks) {
          throw cms::Exception("BTLElectronicsToDetIdMappingESProducer")
              << "crystal " << std::hex << detId.rawId() << std::dec << ": minus and plus sides on different links"
              << " (fed " << elecIds.minus.fedId() << " vs " << elecIds.plus.fedId() << ", hs "
              << elecIds.minus.hsLinkId() << " vs " << elecIds.plus.hsLinkId() << ", el " << elecIds.minus.eLinkId()
              << " vs " << elecIds.plus.eLinkId() << ")";
        }

        // -- isValidLink() is meant to check the ids read from the payload: a readout map that does not respect
        // the layout assumed by the indexer must not be used
        if (!indexer.isValidLink(elecIds.minus.fedId(), elecIds.minus.hsLinkId(), elecIds.minus.eLinkId())) {
          throw cms::Exception("BTLElectronicsToDetIdMappingESProducer")
              << "crystal " << std::hex << detId.rawId() << std::dec << ": fed=" << elecIds.minus.fedId()
              << " hs=" << elecIds.minus.hsLinkId() << " el=" << elecIds.minus.eLinkId()
              << " not compatible with the layout of BTLElectronicsIndexer";
        }

        // -- the unpacker pairs a channel with the channel kTofhirChannelPartner[chId] of the same chip, assuming
        // that it is the other side of the same crystal (kTofhirChannelPartner is an involution, so checking
        // the plus -> minus direction is enough)
        const int32_t chPlus = elecIds.plus.channelId();
        const int32_t chMinus = elecIds.minus.channelId();
        const int32_t nChannels = static_cast<int32_t>(BTLElectronicsSpecs::kTofhirChannelPartner.size());
        if (chPlus < 0 || chPlus >= nChannels || chMinus < 0 || chMinus >= nChannels ||
            BTLElectronicsSpecs::kTofhirChannelPartner[chPlus] != chMinus) {
          throw cms::Exception("BTLElectronicsToDetIdMappingESProducer")
              << "crystal " << std::hex << detId.rawId() << std::dec << ": plus channel " << chPlus
              << " and minus channel " << chMinus << " are not partners in kTofhirChannelPartner";
        }

        const int32_t idxplus = indexer.flatIndex(
            elecIds.plus.fedId(), elecIds.plus.hsLinkId(), elecIds.plus.eLinkId(), elecIds.plus.channelId());
        const int32_t idxminus = indexer.flatIndex(
            elecIds.minus.fedId(), elecIds.minus.hsLinkId(), elecIds.minus.eLinkId(), elecIds.minus.channelId());

        // -- each electronics channel must belong to a single crystal
        for (const int32_t idx : {idxplus, idxminus}) {
          if (product->view()[idx].valid()) {
            throw cms::Exception("BTLElectronicsToDetIdMappingESProducer")
                << "crystal " << std::hex << detId.rawId() << std::dec << ": electronics channel (fed="
                << elecIds.minus.fedId() << " hs=" << elecIds.minus.hsLinkId() << " el=" << elecIds.minus.eLinkId()
                << " ch=" << (idx == idxplus ? chPlus : chMinus) << ") already assigned to crystal " << std::hex << product->view()[idx].rawId() << std::dec;
          }
        }

        product->view()[idxplus].rawId() = detId.rawId();
        product->view()[idxplus].valid() = true;
        product->view()[idxplus].partnerChannelId() =
            BTLElectronicsSpecs::kTofhirChannelPartner[elecIds.plus.channelId()];
        product->view()[idxplus].side() = 1;                                        // plus side
        ++nFilled;

        product->view()[idxminus].rawId() = detId.rawId();
        product->view()[idxminus].valid() = true;
        product->view()[idxminus].partnerChannelId() =
            BTLElectronicsSpecs::kTofhirChannelPartner[elecIds.minus.channelId()];
        product->view()[idxminus].side() = 0;                                        // minus side
        ++nFilled;
      }

      LogDebug("BTLElectronicsToDetIdMappingESProducer")
          << "filled " << nFilled << " / " << indexer.size() << " table entries.";

      return product;
    }

  private:
    edm::ESGetToken<BTLReadoutMap, BTLReadoutMapRcd> readoutMapToken_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_EVENTSETUP_ALPAKA_MODULE(BTLElectronicsToDetIdMappingESProducer);
