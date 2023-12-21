// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "MtdSimLayerClusterToTPAssociatorByUniqueIdImpl.h"

#include "FWCore/Utilities/interface/EDGetToken.h"

//
// class decleration
//

class MtdSimLayerClusterToTPAssociatorByUniqueIdProducer : public edm::global::EDProducer<> {
public:
  explicit MtdSimLayerClusterToTPAssociatorByUniqueIdProducer (const edm::ParameterSet &);
  ~MtdSimLayerClusterToTPAssociatorByUniqueIdProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  
private:
  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;
};


MtdSimLayerClusterToTPAssociatorByUniqueIdProducer::MtdSimLayerClusterToTPAssociatorByUniqueIdProducer(const edm::ParameterSet &ps)
{

  // Register the product
  produces<reco::MtdSimLayerClusterToTPAssociator>();

}



MtdSimLayerClusterToTPAssociatorByUniqueIdProducer::~MtdSimLayerClusterToTPAssociatorByUniqueIdProducer() {}



void MtdSimLayerClusterToTPAssociatorByUniqueIdProducer::produce(edm::StreamID,
								 edm::Event &iEvent,
								 const edm::EventSetup &es) const {

  auto impl = std::make_unique<MtdSimLayerClusterToTPAssociatorByUniqueIdImpl>(iEvent.productGetter());
  auto toPut = std::make_unique<reco::MtdSimLayerClusterToTPAssociator>(std::move(impl));
  iEvent.put(std::move(toPut));
  
}


void MtdSimLayerClusterToTPAssociatorByUniqueIdProducer::fillDescriptions(edm::ConfigurationDescriptions &cfg) {
  edm::ParameterSetDescription desc;
  
  cfg.add("mtdSimLayerClusterToTPAssociatorByUniqueId", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MtdSimLayerClusterToTPAssociatorByUniqueIdProducer);
