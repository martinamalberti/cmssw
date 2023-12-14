// -*- C++ -*-
//
// Package:    Validation/MtdValidation
// Class:      TestClusterAssociation
//
/**\class TestClusterAssociation TestClusterAssociation.cc Validation/MtdValidation/plugins/TestClusterAssociation.cc

 Description: BTL RECO hits and clusters validation

 Implementation:
     [Notes on implementation]
*/

#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/ValidHandle.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "DataFormats/TrackerRecHit2D/interface/MTDTrackingRecHit.h"

#include "SimDataFormats/Associations/interface/MtdRecoClusterToSimLayerClusterAssociationMap.h"

using namespace std;

class TestClusterAssociation : public edm::one::EDAnalyzer<edm::one::WatchRuns,
       							   edm::one::SharedResources> {
public:
  explicit TestClusterAssociation(const edm::ParameterSet&);
  ~TestClusterAssociation() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void beginRun(const edm::Run&, const edm::EventSetup&) override{};
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override{};
  void endJob() override ;
  

  // ------------ member data ------------
  edm::EDGetTokenT<FTLClusterCollection> btlRecCluToken_;
  //edm::EDGetTokenT<FTLClusterCollection> etlRecCluToken_;
  edm::EDGetTokenT<MTDTrackingDetSetVector> mtdTrackingHitToken_;
  edm::EDGetTokenT<MtdRecoClusterToSimLayerClusterAssociationMap> clusterAssociationMapToken_;

};

// ------------ constructor and destructor --------------
TestClusterAssociation::TestClusterAssociation(const edm::ParameterSet& iConfig) {
  btlRecCluToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("bltRecoClusTag"));
  //etlRecCluToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("eltRecoClusTag"));
  clusterAssociationMapToken_ = consumes<MtdRecoClusterToSimLayerClusterAssociationMap>(iConfig.getParameter<edm::InputTag>("clusterAssociationMapTag"));
  mtdTrackingHitToken_ = consumes<MTDTrackingDetSetVector>(iConfig.getParameter<edm::InputTag>("trkHitTag"));
}

TestClusterAssociation::~TestClusterAssociation() {}



// ------------ method called once each job just before starting event loop  ------------
void TestClusterAssociation::beginJob() {

}


// ------------ method called for each event  ------------
void TestClusterAssociation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  
  auto btlRecCluHandle = makeValid(iEvent.getHandle(btlRecCluToken_));
  //auto etlRecCluHandle = makeValid(iEvent.getHandle(etlRecCluToken_));
  auto mtdTrkHitHandle = makeValid(iEvent.getHandle(mtdTrackingHitToken_));
  auto clusterAssociationMap = iEvent.get(clusterAssociationMapToken_);
  
  // --- Loop over the BTL RECO clusters ---
  for (const auto& DetSetClu : *btlRecCluHandle) {
    for (const auto& cluster : DetSetClu) {

      float recoClusEnergy = cluster.energy();
      float recoClusTime = cluster.time();

      // get the reco cluster ref
      edm::Ref<edmNew::DetSetVector<FTLCluster>, FTLCluster> clusterRef = edmNew::makeRefTo(btlRecCluHandle, &cluster);

      // get the corresponding sim cluster ref
      auto it = std::find_if( clusterAssociationMap.begin(), clusterAssociationMap.end(),
			      [&](const std::pair<FTLClusterRef, MtdSimLayerClusterRef>& p) { return p.first == clusterRef; });
      
      auto simClusterRef = (*it).second;
      
      float simClusEnergy = (*simClusterRef).simLCEnergy();
      float simClusTime = (*simClusterRef).simLCTime();

      

      std::cout << "reco cluster energy = " << recoClusEnergy << "    sim cluster energy = " << simClusEnergy <<std::endl;
      std::cout << "reco cluster time = " << recoClusTime << "    sim cluster time = " << simClusTime <<std::endl;
      
    }  // Cluster loop

  }  // DetSetClu loop

}


void TestClusterAssociation::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TestClusterAssociation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("recCluTag", edm::InputTag("mtdClusters", "FTLBarrel"));
  desc.add<edm::InputTag>("trkHitTag", edm::InputTag("mtdTrackingRecHits"));
  descriptions.add("TestClusterAssociation", desc);
}

DEFINE_FWK_MODULE(TestClusterAssociation);
