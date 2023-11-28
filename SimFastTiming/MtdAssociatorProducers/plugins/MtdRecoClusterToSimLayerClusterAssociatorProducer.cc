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

#include "SimDataFormats/Associations/interface/MtdRecoClusterToSimLayerClusterAssociator.h"

#include "FWCore/Utilities/interface/EDGetToken.h"

#include "MtdRecoClusterToSimLayerClusterAssociatorImpl.h"

//
// class decleration
//

class MtdRecoClusterToSimLayerClusterAssociatorProducer : public edm::global::EDProducer<> {
public:
  explicit MtdRecoClusterToSimLayerClusterAssociatorProducer(const edm::ParameterSet &);
  ~MtdRecoClusterToSimLayerClusterAssociatorProducer() override;

private:
  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;
  
  edm::EDGetTokenT<FTLClusterCollection> btlRecoClustersToken_;
  edm::EDGetTokenT<FTLClusterCollection> etlRecoClustersToken_;
  edm::EDGetTokenT<FTLRecHitCollection> btlRecHitsToken_;
  edm::EDGetTokenT<FTLRecHitCollection> etlRecHitsToken_;
  edm::EDGetTokenT<MtdSimLayerClusterCollection> simClustersToken_;

  edm::EDGetTokenT<reco::MtdRecoClusterToSimLayerClusterAssociator> associatorToken_;
};

MtdRecoClusterToSimLayerClusterAssociatorProducer::MtdRecoClusterToSimLayerClusterAssociatorProducer(const edm::ParameterSet &pset) {
  produces<reco::SimToRecoCollectionMtd>();
  produces<reco::RecoToSimCollectionMtd>();

  btlRecoClustersToken_ = consumes<FTLClusterCollection>(pset.getParameter<edm::InputTag>("btlRecoClustersTag"));
  etlRecoClustersToken_ = consumes<FTLClusterCollection>(pset.getParameter<edm::InputTag>("etlRecoClustersTag"));
  btlRecHitsToken_ = consumes<FTLRecHitCollection>(pset.getParameter<edm::InputTag>("btlRecHitsTag"));
  etlRecHitsToken_ = consumes<FTLRecHitCollection>(pset.getParameter<edm::InputTag>("etlRecHitsTag"));
  simClustersToken_  = consumes<MtdSimLayerClusterCollection>(pset.getParameter<edm::InputTag>("mtdSimClustersTag"));
  associatorToken_ = consumes<reco::MtdRecoClusterToSimLayerClusterAssociator>(pset.getParameter<edm::InputTag>("associator"));
}

MtdRecoClusterToSimLayerClusterAssociatorProducer::~MtdRecoClusterToSimLayerClusterAssociatorProducer() {}

//
// member functions
//

// ------------ method called to produce the data  ------------
void MtdRecoClusterToSimLayerClusterAssociatorProducer::produce(edm::StreamID, edm::Event &iEvent, const edm::EventSetup &iSetup) const {
  using namespace edm;

  edm::Handle<reco::MtdRecoClusterToSimLayerClusterAssociator> theAssociator;
  iEvent.getByToken(associatorToken_, theAssociator);

  edm::Handle<FTLClusterCollection> btlRecoClusters;
  iEvent.getByToken(btlRecoClustersToken_, btlRecoClusters);

  edm::Handle<FTLClusterCollection> etlRecoClusters;
  iEvent.getByToken(etlRecoClustersToken_, etlRecoClusters);

  edm::Handle<FTLRecHitCollection> btlRecHits;
  iEvent.getByToken(btlRecHitsToken_, btlRecHits);

  edm::Handle<FTLRecHitCollection> etlRecHits;
  iEvent.getByToken(etlRecHitsToken_, etlRecHits);

  edm::Handle<MtdSimLayerClusterCollection> simClusters;
  iEvent.getByToken(simClustersToken_, simClusters);


  // associate reco clus to sim layer clus
  reco::RecoToSimCollectionMtd recoToSimColl = theAssociator->associateRecoToSim(btlRecoClusters, etlRecoClusters, simClusters, btlRecHits, etlRecHits);
  reco::SimToRecoCollectionMtd simToRecoColl = theAssociator->associateSimToReco(btlRecoClusters, etlRecoClusters, simClusters, btlRecHits, etlRecHits);

  auto r2s = std::make_unique<reco::RecoToSimCollectionMtd>(recoToSimColl);
  auto s2r = std::make_unique<reco::SimToRecoCollectionMtd>(simToRecoColl);

  iEvent.put(std::move(r2s));
  iEvent.put(std::move(s2r));
}

// define this as a plug-in
DEFINE_FWK_MODULE(MtdRecoClusterToSimLayerClusterAssociatorProducer);
