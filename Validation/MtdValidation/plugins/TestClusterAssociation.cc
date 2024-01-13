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
#include "SimDataFormats/Associations/interface/MtdSimLayerClusterToTPAssociatorBaseImpl.h"

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
  edm::EDGetTokenT<reco::SimToTPCollectionMtd> sim2TPAssociationMapToken_;
  edm::EDGetTokenT<reco::TPToSimCollectionMtd> tp2SimAssociationMapToken_;

  // --  Output histos
  edm::Service<TFileService> fs_;
  TH1F *h_energy_BTL_all;
  TH1F *h_energy_ETL_all;

  TH1F *h_energy_BTL;
  TH1F *h_energy_ETL;

  TH1F *h_energy_TPmatched_BTL;
  TH1F *h_energy_TPmatched_ETL;
  
  TH1F *h_time_BTL_all;
  TH1F *h_time_ETL_all;

  TH1F *h_time_BTL;
  TH1F *h_time_ETL;

  TH1F *h_time_TPmatched_BTL;
  TH1F *h_time_TPmatched_ETL;
  
  TH1F *h_energyResol_BTL;
  TH1F *h_energyResol_ETL;

  TH1F *h_timeResol_BTL; 
  TH1F *h_timeResol_ETL;

  TH1F *h_nSimClusPerRecoClus_BTL;
  TH1F *h_nSimClusPerRecoClus_ETL;

  TH1F *h_nRecoClusPerSimClus_BTL;
  TH1F *h_nRecoClusPerSimClus_ETL;

  TH1F *h_energySimClus_BTL;
  TH1F *h_energySimClus_ETL;

  TH1F *h_timeSimClus_BTL;
  TH1F *h_timeSimClus_ETL;

  TH1F *h_etaSimClus;
  TH1F *h_etaSimClus_assocToReco;

  int nRecoClus = 0;
  int nSimClus = 0;
  
};


  
// ------------ constructor and destructor --------------
TestClusterAssociation::TestClusterAssociation(const edm::ParameterSet& iConfig) {
  btlRecoClusToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("btlRecoClustersTag"));
  etlRecoClusToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("etlRecoClustersTag"));
  r2sAssociationMapToken_ = consumes<MtdRecoClusterToSimLayerClusterAssociationMap>(iConfig.getParameter<edm::InputTag>("r2sAssociationMapTag"));
  s2rAssociationMapToken_ = consumes<MtdSimLayerClusterToRecoClusterAssociationMap>(iConfig.getParameter<edm::InputTag>("s2rAssociationMapTag"));
  sim2TPAssociationMapToken_ = consumes<reco::SimToTPCollectionMtd>(iConfig.getParameter<edm::InputTag>("sim2TPAssociationMapTag"));
  tp2SimAssociationMapToken_ = consumes<reco::TPToSimCollectionMtd>(iConfig.getParameter<edm::InputTag>("tp2SimAssociationMapTag"));
  mtdTrackingHitToken_ = consumes<MTDTrackingDetSetVector>(iConfig.getParameter<edm::InputTag>("trkHitTag"));
  mtdSimLayerClustersToken_ = consumes<std::vector<MtdSimLayerCluster> >(iConfig.getParameter<edm::InputTag>("mtdSimLayerClustersTag"));

  // -- book histograms
  h_energy_BTL_all = fs_->make<TH1F>("h_energy_BTL_all","h_energy_BTL_all", 200, 0, 20);
  h_energy_ETL_all = fs_->make<TH1F>("h_energy_ETL_all","h_energy_ETL_all", 200, 0, 20);
  
  h_energy_BTL = fs_->make<TH1F>("h_energy_BTL","h_energy_BTL", 200, 0, 20);
  h_energy_ETL = fs_->make<TH1F>("h_energy_ETL","h_energy_ETL", 200, 0, 20);

  h_energy_TPmatched_BTL = fs_->make<TH1F>("h_energy_TPmatched_BTL","h_energy_TPmatched_BTL", 200, 0, 20);
  h_energy_TPmatched_ETL = fs_->make<TH1F>("h_energy_TPmatched_ETL","h_energy_TPmatched_ETL", 200, 0, 20);
  
  h_time_BTL_all = fs_->make<TH1F>("h_time_BTL_all","h_time_BTL_all", 100, 0, 25);
  h_time_ETL_all = fs_->make<TH1F>("h_time_ETL_all","h_time_ETL_all", 100, 0, 25);
  
  h_time_BTL = fs_->make<TH1F>("h_time_BTL","h_time_BTL", 100, 0, 25);
  h_time_ETL = fs_->make<TH1F>("h_time_ETL","h_time_ETL", 100, 0, 25);

  h_time_TPmatched_BTL = fs_->make<TH1F>("h_time_TPmatched_BTL","h_time_TPmatched_BTL", 100, 0, 25);
  h_time_TPmatched_ETL = fs_->make<TH1F>("h_time_TPmatched_ETL","h_time_TPmatched_ETL", 100, 0, 25);

  h_energyResol_BTL =  fs_->make<TH1F>("h_energyResol_BTL","h_energyResol_BTL", 800, 0, 4);
  h_energyResol_ETL = fs_->make<TH1F>("h_energyResol_ETL","h_energyResol_ETL", 800, 0, 4);

  h_timeResol_BTL = fs_->make<TH1F>("h_timeResol_BTL","h_timeResol_BTL", 200, -1, 1);
  h_timeResol_ETL = fs_->make<TH1F>("h_timeResol_ETL","h_timeResol_ETL", 200, -1, 1);

  h_nSimClusPerRecoClus_BTL = fs_->make<TH1F>("h_nSimClusPerRecoClus_BTL","h_nSimClusPerRecoClus_BTL", 5, -0.5, 4.5);
  h_nSimClusPerRecoClus_ETL = fs_->make<TH1F>("h_nSimClusPerRecoClus_ETL","h_nSimClusPerRecoClus_ETL", 5, -0.5, 4.5);

  h_nRecoClusPerSimClus_BTL = fs_->make<TH1F>("h_nRecoClusPerSimClus_BTL","h_nRecoClusPerSimClus_BTL", 5, -0.5, 4.5);
  h_nRecoClusPerSimClus_ETL = fs_->make<TH1F>("h_nRecoClusPerSimClus_ETL","h_nRecoClusPerSimClus_ETL", 5, -0.5, 4.5);

  h_energySimClus_BTL = fs_->make<TH1F>("h_energySimClus_BTL","h_energySimClus_BTL", 200, 0., 20.0);
  h_energySimClus_ETL = fs_->make<TH1F>("h_energySimClus_ETL","h_energySimClus_ETL", 200, 0., 20.0);

  h_timeSimClus_BTL = fs_->make<TH1F>("h_timeSimClus_BTL","h_timeSimClus_BTL", 100, 0., 25.0);
  h_timeSimClus_ETL = fs_->make<TH1F>("h_timeSimClus_ETL","h_timeSimClus_ETL", 100, 0., 25.0);
  
  h_etaSimClus = fs_->make<TH1F>("h_etaSimClus","h_etaSimClus", 150, -3.0, 3.0);
  h_etaSimClus_assocToReco = fs_->make<TH1F>("h_etaSimClus_assocToReco","h_etaSimClus_assocToReco", 150, -3.0, 3.0);
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
  auto sim2TPAssociationMap = iEvent.get(sim2TPAssociationMapToken_);
  auto tp2SimAssociationMap = iEvent.get(tp2SimAssociationMapToken_);

  auto mtdSimLayerClustersH = iEvent.getHandle(mtdSimLayerClustersToken_);
  std::vector<MtdSimLayerCluster> mtdSimLayerClusters =  *mtdSimLayerClustersH;
      
  std::array<edm::Handle<FTLClusterCollection>, 2> inputRecoClusH{{btlRecoClusH, etlRecoClusH}};
  
  std::cout << "recoToSim map size: " << r2sAssociationMap.size() << "           simToReco map size: " << s2rAssociationMap.size() << std::endl;
	
	
  // --- Loop over the RECO clusters ---
  for (auto const& recoClusH : inputRecoClusH) {
    for (const auto& DetSetClu : *recoClusH) {
      for (const auto& cluster : DetSetClu) {

	MTDDetId mtdDetId = cluster.id();	
	float recoClusEnergy = cluster.energy();
	float recoClusTime = cluster.time();

	//std::cout << " ----- NEW RECO CLUSTER ----  recoClusEnergy  = " << recoClusEnergy << "    recoClusTime = " << recoClusTime <<std::endl;

	if (mtdDetId.mtdSubDetector() == MTDDetId::BTL){
            h_energy_BTL_all->Fill(recoClusEnergy);
            h_time_BTL_all->Fill(recoClusTime);
	}

	if (mtdDetId.mtdSubDetector() == MTDDetId::ETL){
            h_energy_ETL_all->Fill(recoClusEnergy);
            h_time_ETL_all->Fill(recoClusTime);
	}
	
	// get the reco cluster ref
	edm::Ref<edmNew::DetSetVector<FTLCluster>, FTLCluster> clusterRef = edmNew::makeRefTo(recoClusH, &cluster);

	// get the corresponding sim clusters ref
	auto it = std::find_if( r2sAssociationMap.begin(), r2sAssociationMap.end(),
				[&](const std::pair<FTLClusterRef, std::vector<MtdSimLayerClusterRef>>& p) { return p.first == clusterRef; });

	if ( it == r2sAssociationMap.end()) {
	  std::cout << " reco cluster not found "<<std::endl;
	  continue;
	}
	
	std::vector<MtdSimLayerClusterRef> simClustersRefs = (*it).second;
	
	//std::cout << "Number of associated MtdSimLayerclusters: " << simClustersRefs.size() <<std::endl;
	  
	if (mtdDetId.mtdSubDetector() == MTDDetId::BTL)	h_nSimClusPerRecoClus_BTL->Fill(simClustersRefs.size());
	if (mtdDetId.mtdSubDetector() == MTDDetId::ETL)	h_nSimClusPerRecoClus_ETL->Fill(simClustersRefs.size());
	
	for (unsigned int i = 0; i < simClustersRefs.size(); i++){
	  auto simClusterRef = simClustersRefs[i];
	  
	  float simClusEnergy = (*simClusterRef).simLCEnergy();
	  float simClusTime = (*simClusterRef).simLCTime();
	  
	  //std::cout << "reco cluster energy = " << recoClusEnergy << "    sim cluster energy = " << simClusEnergy <<std::endl;
	  //std::cout << "reco cluster time   = " << recoClusTime << "    sim cluster time = " << simClusTime <<std::endl;

	  // -- BTL
	  if (mtdDetId.mtdSubDetector() == MTDDetId::BTL){
            h_time_BTL->Fill(recoClusTime);
	    h_energy_BTL->Fill(recoClusEnergy);
	    h_energyResol_BTL->Fill( 0.001*recoClusEnergy/simClusEnergy);
	    h_timeResol_BTL->Fill(recoClusTime-simClusTime);
	  }
	  
	  // -- ETL
	  if (mtdDetId.mtdSubDetector() == MTDDetId::ETL){
	    h_time_ETL->Fill(recoClusTime);
	    h_energy_ETL->Fill(recoClusEnergy);
	    h_energyResol_ETL->Fill( 0.001*recoClusEnergy/simClusEnergy);
	    h_timeResol_ETL->Fill(recoClusTime-simClusTime);
	  }

	  // -- check matching with tracking particles
	  auto found = sim2TPAssociationMap.find(simClusterRef);
	  //if (found == sim2TPAssociationMap.end()) cout<< "Sim cluster not matched to any TP"<<std::endl;
	  if (found != sim2TPAssociationMap.end()) {
	    if (mtdDetId.mtdSubDetector() == MTDDetId::BTL){
	      h_time_TPmatched_BTL->Fill(recoClusTime);
	      h_energy_TPmatched_BTL->Fill(recoClusEnergy);
	    }
	    if (mtdDetId.mtdSubDetector() == MTDDetId::ETL){
	      h_time_TPmatched_ETL->Fill(recoClusTime);
	      h_energy_TPmatched_ETL->Fill(recoClusEnergy);
	    }
	    //  cout << "Sim layer cluster matched to a tracking particle"<<std::endl;
	    //  for (const auto& tpRef : found->val) {
	    //    cout <<  (*tpRef).g4Tracks()[0].trackId() << std::endl;
	    //  }
	  }
	  
	} // end loop over simClustersrefs
	  
      }  // reco cluster loop

    }  // DetSetClu loop
  }
  

  
  // -- Loop over MtdSimLayerClusters
  edm::Ref<MtdSimLayerClusterCollection>::key_type simClusIndex = 0;
  for (auto simClus : mtdSimLayerClusters){

    // -- get the sim cluster ref
    edm::Ref<MtdSimLayerClusterCollection> simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(mtdSimLayerClustersH, simClusIndex);

    // -- get the corresponding reco cluster ref
    auto it = std::find_if( s2rAssociationMap.begin(), s2rAssociationMap.end(),
			    [&](const std::pair<MtdSimLayerClusterRef, std::vector<FTLClusterRef>>& p) { return p.first == simClusterRef; });

    std::vector<FTLClusterRef> recoClustersRefs = (*it).second;
	
    //std::cout << "Number of associated reco clusters: " << recoClustersRefs.size() <<std::endl;

    if (std::abs(simClus.eta())<1.5 ) h_nRecoClusPerSimClus_BTL->Fill(recoClustersRefs.size());
    if (std::abs(simClus.eta())>1.5 ) h_nRecoClusPerSimClus_ETL->Fill(recoClustersRefs.size());

    h_etaSimClus->Fill(simClus.eta());
    if (recoClustersRefs.size()>0) h_etaSimClus_assocToReco->Fill(simClus.eta());
    
    
    if (std::abs(simClus.eta())<1.5 ) {
	h_energySimClus_BTL ->Fill(simClus.simLCEnergy()*1000);
	if (simClus.simLCEnergy()*1000 > 1.)	h_timeSimClus_BTL ->Fill(simClus.simLCTime());
      }
    if (std::abs(simClus.eta())>1.5 ) {
      h_energySimClus_ETL ->Fill(simClus.simLCEnergy()*1000); // GeV -> MeV
      h_timeSimClus_ETL ->Fill(simClus.simLCTime());
    }
    
    simClusIndex++;

  }


  // loop over tracking particles Ref
  

}




void TestClusterAssociation::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TestClusterAssociation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("btlRecoClustersTag", edm::InputTag("mtdClusters", "FTLBarrel"));
  desc.add<edm::InputTag>("etlRecoClustersTag", edm::InputTag("mtdClusters", "FTLEndcap"));
  desc.add<edm::InputTag>("r2sAssociationMapTag", edm::InputTag("mtdRecoClusterToSimLayerClusterAssociation"));
  desc.add<edm::InputTag>("s2rAssociationMapTag", edm::InputTag("mtdRecoClusterToSimLayerClusterAssociation"));
  desc.add<edm::InputTag>("sim2TPAssociationMapTag", edm::InputTag("mtdSimLayerClusterToTPAssociation"));
  desc.add<edm::InputTag>("tp2SimAssociationMapTag", edm::InputTag("mtdSimLayerClusterToTPAssociation"));
  desc.add<edm::InputTag>("trkHitTag", edm::InputTag("mtdTrackingRecHits"));
  desc.add<edm::InputTag>("mtdSimLayerClustersTag", edm::InputTag("mix","MergedMtdTruthLC"));
  descriptions.add("TestClusterAssociation", desc);
}

DEFINE_FWK_MODULE(TestClusterAssociation);
