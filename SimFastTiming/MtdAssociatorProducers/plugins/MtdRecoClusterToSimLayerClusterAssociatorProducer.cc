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

//
// class decleration
//

class MtdRecoClusterToSimLayerClusterAssociatorProducer : public edm::global::EDProducer<> {
public:
  explicit MtdRecoClusterToSimLayerClusterAssociatorProducer(const edm::ParameterSet &);
  ~MtdRecoClusterToSimLayerClusterAssociatorProducer() override;

private:
  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;
  
  edm::EDGetTokenT<FTLClusterCollection> recoClusCollectionToken_;
  edm::EDGetTokenT<MtdSimLayerClusterCollection> simClusCollectionToken_;
  edm::EDGetTokenT<reco::MtdRecoClusterToSimLayerClusterAssociator> associatorToken_;
};

MtdRecoClusterToSimLayerClusterAssociatorProducer::MtdRecoClusterToSimLayerClusterAssociatorProducer(const edm::ParameterSet &pset) {
  produces<reco::SimToRecoCollection>();
  produces<reco::RecoToSimCollection>();

  recoClusCollectionToken_ = consumes<FTLClusterCollection>(pset.getParameter<edm::InputTag>("mtdRecoClustersTag"));
  simClusCollectionToken_  = consumes<SimClusterCollection>(pset.getParameter<edm::InputTag>("mtdSimClustersTag"));
  
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

  Handle<FTLClusterCollection> recoClusCollection;
  iEvent.getByToken(recoClusCollectionToken_, recoClusCollection);

  Handle<MtdSimLayerClusterCollection> simClusCollection;
  iEvent.getByToken(simClusCollectionToken_, simClusCollection);

  // associate reco clus to sim layer clus
  reco::RecoToSimCollection recoToSimColl = theAssociator->associateRecoToSim(recoClusCollection, simClusCollection);
  reco::SimToRecoCollection simToRecoColl = theAssociator->associateSimToReco(recoClusCollection, simClusCollection);

  auto r2s = std::make_unique<reco::RecoToSimCollection>(recoToSimColl);
  auto s2r = std::make_unique<reco::SimToRecoCollection>(simToRecoColl);

  iEvent.put(std::move(r2s));
  iEvent.put(std::move(s2r));
}

// define this as a plug-in
DEFINE_FWK_MODULE(MtdRecoClusterToSimLayerClusterAssociatorProducer);
