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
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "Geometry/MTDNumberingBuilder/interface/MTDTopology.h"
#include "MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl.h"

//
// Class declaration
//

class MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer : public edm::global::EDProducer<> {
public:
  explicit MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer (const edm::ParameterSet &);
  ~MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  
private:
  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;
  //edm::EDGetTokenT<FTLRecHitCollection> btlRecHitsToken_;
  //edm::EDGetTokenT<FTLRecHitCollection> etlRecHitsToken_;
  //edm::ESGetToken<MTDTopology, MTDTopologyRcd> topoToken_;
  const double energyCut_;
  const double timeCut_;
};


MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer::MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer(const edm::ParameterSet &ps)
  : energyCut_(ps.getParameter<double>("energyCut")), 
    timeCut_(ps.getParameter<double>("timeCut"))  {

  //btlRecHitsToken_ = consumes<FTLRecHitCollection>(ps.getParameter<edm::InputTag>("btlRecHitsTag"));
  //etlRecHitsToken_ = consumes<FTLRecHitCollection>(ps.getParameter<edm::InputTag>("etlRecHitsTag"));
  //topoToken_ = esConsumes<MTDTopology, MTDTopologyRcd>();
 
  // Register the product
  produces<reco::MtdRecoClusterToSimLayerClusterAssociator>();

}



MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer::~MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer() {}



void MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer::produce(edm::StreamID,
								      edm::Event &iEvent,
								      const edm::EventSetup &es) const {

  /*
    edm::Handle<FTLRecHitCollection> btlRecHitsH;
  iEvent.getByToken(btlRecHitsToken_, btlRecHitsH);

  edm::Handle<FTLRecHitCollection> etlRecHitsH;
  iEvent.getByToken(etlRecHitsToken_, etlRecHitsH);

  auto topologyHandle = es.getTransientHandle(topoToken_);
  const MTDTopology* topology = topologyHandle.product();
  
  auto impl = std::make_unique<MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl>(iEvent.productGetter(), btlRecHitsH, etlRecHitsH, topology, energyCut_, timeCut_);
  */
  auto impl = std::make_unique<MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl>(iEvent.productGetter(), energyCut_, timeCut_);
  auto toPut = std::make_unique<reco::MtdRecoClusterToSimLayerClusterAssociator>(std::move(impl));
  iEvent.put(std::move(toPut));
  
}


void MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer::fillDescriptions(edm::ConfigurationDescriptions &cfg) {
  edm::ParameterSetDescription desc;
  //desc.add<edm::InputTag>("btlRecHitsTag", edm::InputTag("mtdClusters:FTLBarrel"));
  //desc.add<edm::InputTag>("etlRecHitsTag", edm::InputTag("mtdClusters:FTLEndcap"));
  desc.add<double>("energyCut", 9999.);
  desc.add<double>("timeCut", 9999.);
  
  //cfg.add("mtdClusterAssociatorByHits", desc);
  cfg.add("mtdRecoClusterToSimLayerClusterAssociatorByHits", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MtdRecoClusterToSimLayerClusterAssociatorByHitsProducer);
