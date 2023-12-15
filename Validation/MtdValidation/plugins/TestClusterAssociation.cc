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
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "DataFormats/TrackerRecHit2D/interface/MTDTrackingRecHit.h"
#include "DataFormats/ForwardDetId/interface/MTDDetId.h"

#include "SimDataFormats/CaloAnalysis/interface/MtdSimLayerCluster.h"

#include "SimDataFormats/Associations/interface/MtdRecoClusterToSimLayerClusterAssociationMap.h"
#include "SimDataFormats/Associations/interface/MtdSimLayerClusterToRecoClusterAssociationMap.h"

#include "TTree.h"
#include "TH1F.h"

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
  edm::EDGetTokenT<FTLClusterCollection> btlRecoClusToken_;
  edm::EDGetTokenT<FTLClusterCollection> etlRecoClusToken_;
  edm::EDGetTokenT<MTDTrackingDetSetVector> mtdTrackingHitToken_;
  edm::EDGetTokenT<std::vector<MtdSimLayerCluster>> mtdSimLayerClustersToken_;
  edm::EDGetTokenT<MtdRecoClusterToSimLayerClusterAssociationMap> r2sAssociationMapToken_;
  edm::EDGetTokenT<MtdSimLayerClusterToRecoClusterAssociationMap> s2rAssociationMapToken_;


  // --  Output histos
  edm::Service<TFileService> fs_;
  TH1F *h_energyResol_BTL;
  TH1F *h_energyResol_ETL;

  TH1F *h_timeResol_BTL; 
  TH1F *h_timeResol_ETL;

  TH1F *h_nSimClusPerRecoClus_BTL;
  TH1F *h_nSimClusPerRecoClus_ETL;

  TH1F *h_nRecoClusPerSimClus_BTL;
  TH1F *h_nRecoClusPerSimClus_ETL;

  TH1F *h_etaSimClus;
  TH1F *h_etaSimClus_assocToReco;
    
};


  
// ------------ constructor and destructor --------------
TestClusterAssociation::TestClusterAssociation(const edm::ParameterSet& iConfig) {
  btlRecoClusToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("btlRecoClustersTag"));
  etlRecoClusToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("etlRecoClustersTag"));
  r2sAssociationMapToken_ = consumes<MtdRecoClusterToSimLayerClusterAssociationMap>(iConfig.getParameter<edm::InputTag>("r2sAssociationMapTag"));
  s2rAssociationMapToken_ = consumes<MtdSimLayerClusterToRecoClusterAssociationMap>(iConfig.getParameter<edm::InputTag>("s2rAssociationMapTag"));
  mtdTrackingHitToken_ = consumes<MTDTrackingDetSetVector>(iConfig.getParameter<edm::InputTag>("trkHitTag"));
  mtdSimLayerClustersToken_ = consumes<std::vector<MtdSimLayerCluster> >(iConfig.getParameter<edm::InputTag>("mtdSimLayerClustersTag"));


  // -- book histograms
  h_energyResol_BTL =  fs_->make<TH1F>("h_energyResol_BTL","h_energyResol_BTL", 200, -1, 1);
  h_energyResol_ETL = fs_->make<TH1F>("h_energyResol_ETL","h_energyResol_ETL", 200, -1, 1);

  h_timeResol_BTL = fs_->make<TH1F>("h_timeResol_BTL","h_timeResol_BTL", 200, -1, 1);
  h_timeResol_ETL = fs_->make<TH1F>("h_timeResol_ETL","h_timeResol_ETL", 200, -1, 1);

  h_nSimClusPerRecoClus_BTL = fs_->make<TH1F>("h_nSimClusPerRecoClus_BTL","h_nSimClusPerRecoClus_BTL", 5, -0.0, 4.5);
  h_nSimClusPerRecoClus_ETL = fs_->make<TH1F>("h_nSimClusPerRecoClus_ETL","h_nSimClusPerRecoClus_ETL", 5, -0.0, 4.5);

  h_nRecoClusPerSimClus_BTL = fs_->make<TH1F>("h_nRecoClusPerSimClus_BTL","h_nRecoClusPerSimClus_BTL", 5, -0.0, 4.5);
  h_nRecoClusPerSimClus_ETL = fs_->make<TH1F>("h_nRecoClusPerSimClus_ETL","h_nRecoClusPerSimClus_ETL", 5, -0.0, 4.5);

  h_etaSimClus = fs_->make<TH1F>("h_etaSimClus","h_etaSimClus", 100, -3.0, 3.0);

  h_etaSimClus_assocToReco = fs_->make<TH1F>("h_etaSimClus_assocToReco","h_etaSimClus_assocToReco", 100, -3.0, 3.0);
}

TestClusterAssociation::~TestClusterAssociation() {}



// ------------ method called once each job just before starting event loop  ------------
void TestClusterAssociation::beginJob() {

  
}


// ------------ method called for each event  ------------
void TestClusterAssociation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // -- get the collections
  auto btlRecoClusH = iEvent.getHandle(btlRecoClusToken_);
  auto etlRecoClusH = iEvent.getHandle(etlRecoClusToken_);
  auto mtdTrkHitH = iEvent.getHandle(mtdTrackingHitToken_);
  auto r2sAssociationMap = iEvent.get(r2sAssociationMapToken_);
  auto s2rAssociationMap = iEvent.get(s2rAssociationMapToken_);

  auto mtdSimLayerClustersH = iEvent.getHandle(mtdSimLayerClustersToken_);
  std::vector<MtdSimLayerCluster> mtdSimLayerClusters =  *mtdSimLayerClustersH;
      

  std::array<edm::Handle<FTLClusterCollection>, 2> inputRecoClusH{{btlRecoClusH, etlRecoClusH}};

  // --- Loop over the RECO clusters ---
  for (auto const& recoClusH : inputRecoClusH) {
    for (const auto& DetSetClu : *recoClusH) {
      for (const auto& cluster : DetSetClu) {

	MTDDetId mtdDetId = cluster.id();	
	float recoClusEnergy = cluster.energy();
	float recoClusTime = cluster.time();

	// get the reco cluster ref
	edm::Ref<edmNew::DetSetVector<FTLCluster>, FTLCluster> clusterRef = edmNew::makeRefTo(recoClusH, &cluster);

	// get the corresponding sim cluster ref
	auto it = std::find_if( r2sAssociationMap.begin(), r2sAssociationMap.end(),
				[&](const std::pair<FTLClusterRef, std::vector<MtdSimLayerClusterRef>>& p) { return p.first == clusterRef; });

	std::vector<MtdSimLayerClusterRef> simClustersRefs = (*it).second;
	
	std::cout << "Number of associated MtdSimLayerclusters: " << simClustersRefs.size() <<std::endl;

	if (mtdDetId.mtdSubDetector() == MTDDetId::BTL)	h_nSimClusPerRecoClus_BTL->Fill(simClustersRefs.size());
	if (mtdDetId.mtdSubDetector() == MTDDetId::ETL)	h_nSimClusPerRecoClus_ETL->Fill(simClustersRefs.size());
	
	for (unsigned int i = 0; i < simClustersRefs.size(); i++){
	  auto simClusterRef = simClustersRefs[i];
	  
	  float simClusEnergy = (*simClusterRef).simLCEnergy();
	  float simClusTime = (*simClusterRef).simLCTime();
	  
	  std::cout << "reco cluster energy = " << recoClusEnergy << "    sim cluster energy = " << simClusEnergy <<std::endl;
	  std::cout << "reco cluster time = " << recoClusTime << "    sim cluster time = " << simClusTime <<std::endl;

	  // -- BTL
	  if (mtdDetId.mtdSubDetector() == MTDDetId::BTL){
	    h_energyResol_BTL->Fill( (0.001*recoClusEnergy-simClusEnergy)/simClusEnergy);
	    h_timeResol_BTL->Fill(recoClusTime-simClusTime);
	  }
	  
	  // -- ETL
	  if (mtdDetId.mtdSubDetector() == MTDDetId::ETL){
	    h_energyResol_ETL->Fill( (0.001*recoClusEnergy-simClusEnergy)/simClusEnergy);
	    h_timeResol_ETL->Fill(recoClusTime-simClusTime);
	  }
	  
	}
	  
      }  // reco cluster loop

    }  // DetSetClu loop
  }
  

  
  // --- Loop over the SIM clusters ---    
  
  // -- Loop over MtdSimLayerClusters
  edm::Ref<MtdSimLayerClusterCollection>::key_type simClusIndex = 0;
  for (auto simClus : mtdSimLayerClusters){

    // -- get the sim cluster ref
    edm::Ref<MtdSimLayerClusterCollection> simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(mtdSimLayerClustersH, simClusIndex);

    // -- get the corresponding reco cluster ref
    auto it = std::find_if( s2rAssociationMap.begin(), s2rAssociationMap.end(),
			    [&](const std::pair<MtdSimLayerClusterRef, std::vector<FTLClusterRef>>& p) { return p.first == simClusterRef; });

    std::vector<FTLClusterRef> recoClustersRefs = (*it).second;
	
    std::cout << "Number of associated reco clusters: " << recoClustersRefs.size() <<std::endl;
      
    if (std::abs(simClus.eta())<1.5 ) h_nRecoClusPerSimClus_BTL->Fill(recoClustersRefs.size());
    if (std::abs(simClus.eta())>1.5 ) h_nRecoClusPerSimClus_ETL->Fill(recoClustersRefs.size());

    h_etaSimClus->Fill(simClus.eta());
    if (recoClustersRefs.size()>0) h_etaSimClus_assocToReco->Fill(simClus.eta());
    
    
    simClusIndex++;

  }



}




void TestClusterAssociation::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TestClusterAssociation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("btlRecoClustersTag", edm::InputTag("mtdClusters", "FTLBarrel"));
  desc.add<edm::InputTag>("etlRecoClustersTag", edm::InputTag("mtdClusters", "FTLEndcap"));
  desc.add<edm::InputTag>("r2sAssociationMapTag", edm::InputTag("mtdRecoClusterToSimLayerClusterAssociation"));
  desc.add<edm::InputTag>("s2rAssociationMapTag", edm::InputTag("mtdRecoClusterToSimLayerClusterAssociation"));
  desc.add<edm::InputTag>("trkHitTag", edm::InputTag("mtdTrackingRecHits"));
  desc.add<edm::InputTag>("mtdSimLayerClustersTag", edm::InputTag("mix","MergedMtdTruthLC"));
  descriptions.add("TestClusterAssociation", desc);
}

DEFINE_FWK_MODULE(TestClusterAssociation);
