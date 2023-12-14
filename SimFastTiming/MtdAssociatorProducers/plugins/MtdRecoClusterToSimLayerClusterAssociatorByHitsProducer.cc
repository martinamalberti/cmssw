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

#include "MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl.h"

#include "FWCore/Utilities/interface/EDGetToken.h"

//
// class decleration
//

class MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer : public edm::global::EDProducer<> {
public:
  explicit MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer (const edm::ParameterSet &);
  ~MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  
private:
  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;
  const double energyCut_;
  const double timeCut_;
};


MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer::MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer(const edm::ParameterSet &ps)
  : energyCut_(ps.getParameter<double>("energyCut")), 
    timeCut_(ps.getParameter<double>("timeCut"))  {

  // Register the product
  produces<reco::MtdRecoClusterToSimLayerClusterAssociator>();

}



MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer::~MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer() {}



void MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer::produce(edm::StreamID,
								 edm::Event &iEvent,
								 const edm::EventSetup &es) const {

  auto impl = std::make_unique<MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl>(iEvent.productGetter(), energyCut_, timeCut_);
  auto toPut = std::make_unique<reco::MtdRecoClusterToSimLayerClusterAssociator>(std::move(impl));
  iEvent.put(std::move(toPut));
  
}


void MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer::fillDescriptions(edm::ConfigurationDescriptions &cfg) {
  edm::ParameterSetDescription desc;

  desc.add<double>("energyCut", 9999.);
  desc.add<double>("timeCut", 9999.);
  
  cfg.add("mtdClusterAssociatorByHits", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer);
