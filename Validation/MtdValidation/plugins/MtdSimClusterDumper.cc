#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/ValidHandle.h"
#include "DataFormats/Math/interface/GeantUnits.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"

#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimLayerCluster.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/Records/interface/MTDTopologyRcd.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDNumberingBuilder/interface/MTDTopology.h"

#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

#include "Geometry/MTDCommonData/interface/MTDTopologyMode.h"

#include "DataFormats/ForwardDetId/interface/MTDDetId.h"

#include "TTree.h"

using namespace std;

struct MTDSimCluster {
  float energy;
  int nClusHits;
  vector<float> clusHitEnergy;
  vector<int> simHitDetId;
  vector<float> simHitEnergy;
  vector<unsigned int> simHitTrackId;
  vector<unsigned int> simHitOriginalTrackId;
  vector<unsigned int> simHitOffsetTrackId;
};




class MtdSimClusterDumper : public edm::one::EDAnalyzer<edm::one::WatchRuns,
							edm::one::SharedResources> {
public:
  explicit MtdSimClusterDumper(const edm::ParameterSet&);
  ~MtdSimClusterDumper() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void beginRun(const edm::Run&, const edm::EventSetup&) override{};
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override{};
  void endJob() override ;
    
  void clearTreeVariables();


  // ------------ member data ------------
  const float hitMinEnergy_;

  edm::EDGetTokenT<CrossingFrame<PSimHit> > btlSimHitsToken_;
  edm::EDGetTokenT<std::vector<MtdSimLayerCluster>> mtdSimLayerClustersToken_;
    
  // --  Output tree
  TTree* tree_; 
  MTDSimCluster *mtdSimClusterInfo;  

};



// ------------ constructor and destructor --------------
MtdSimClusterDumper::MtdSimClusterDumper(const edm::ParameterSet& iConfig):
  hitMinEnergy_(iConfig.getParameter<double>("hitMinEnergy"))
{
  btlSimHitsToken_ = consumes<CrossingFrame<PSimHit> >(iConfig.getParameter<edm::InputTag>("btlSimHitsTag"));
  mtdSimLayerClustersToken_ = consumes<std::vector<MtdSimLayerCluster> >(iConfig.getParameter<edm::InputTag>("mtdSimLayerClustersTag"));
}

MtdSimClusterDumper::~MtdSimClusterDumper() {}

// ------------ method called once each job just before starting event loop  ------------
void MtdSimClusterDumper::beginJob() {

  mtdSimClusterInfo = new MTDSimCluster;
  
  edm::Service<TFileService> fs_;
  tree_ = fs_->make<TTree>("clusters", "simClusters tree");

  tree_ -> Branch( "energy",                &mtdSimClusterInfo->energy);
  tree_ -> Branch( "nClusHits",             &mtdSimClusterInfo->nClusHits);
  tree_ -> Branch( "clusHitEnergy",         &mtdSimClusterInfo->clusHitEnergy);
  tree_ -> Branch( "simHitDetId",           &mtdSimClusterInfo->simHitDetId);
  tree_ -> Branch( "simHitEnergy",          &mtdSimClusterInfo->simHitEnergy);
  tree_ -> Branch( "simHitTrackId",         &mtdSimClusterInfo->simHitTrackId);
  tree_ -> Branch( "simHitOriginalTrackId", &mtdSimClusterInfo->simHitOriginalTrackId);
  tree_ -> Branch( "simHitOffsetTrackId",   &mtdSimClusterInfo->simHitOffsetTrackId);
  
}

// ------------ method called for each event  ------------
void MtdSimClusterDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

 
  using namespace edm;
  using namespace std;
  
  // -- get simHits and simClusters
  auto btlSimHitsHandle = makeValid(iEvent.getHandle(btlSimHitsToken_));
  MixCollection<PSimHit> btlSimHits(btlSimHitsHandle.product());

  auto mtdSimLayerClustersHandle = makeValid(iEvent.getHandle(mtdSimLayerClustersToken_));
  std::vector<MtdSimLayerCluster> mtdSimLayerClusters =  *mtdSimLayerClustersHandle;

  std::unordered_map<uint32_t, std::vector<int> > map_btlSimHitTrackId;
  std::unordered_map<uint32_t, std::vector<int> > map_btlSimHitOrigTrackId;
  std::unordered_map<uint32_t, std::vector<int> > map_btlSimHitOffsetTrackId;
  std::unordered_map<uint32_t, std::vector<float> > map_btlSimHitEnergy;
    
  // -- Loop over BTL sim hits
  for (auto const& simHit : btlSimHits) {

    // --- Use only hits compatible with the in-time bunch-crossing
    if (simHit.tof() < 0 || simHit.tof() > 25.)
      continue;

    // -- min energy cut?
    if (simHit.energyLoss() < hitMinEnergy_)
      continue;

    // -- for each cell, save the simHit trackId and energyLoss  
    DetId id = simHit.detUnitId();
    map_btlSimHitTrackId[id].push_back(simHit.trackId()); 
    map_btlSimHitOrigTrackId[id].push_back(simHit.originalTrackId()); 
    map_btlSimHitOffsetTrackId[id].push_back(simHit.offsetTrackId()); 
    map_btlSimHitEnergy[id].push_back(simHit.energyLoss()); 
  }

  
  // -- Loop over MtdSimLayerClusters
  for (auto simClus : mtdSimLayerClusters){

    // -- initialize output tree
    clearTreeVariables();
    
    // -- get detIds for this cluster
    std::vector<std::pair<uint32_t, std::pair<uint8_t, uint8_t>>> simClusterDetIds = simClus.detIds_and_rows();

    std::vector<std::pair<uint64_t, float>> simClusterHitsAndEnergies = simClus.hits_and_energies();

    std::cout << "Sim cluster energy : " << simClus.simLCEnergy() <<std::endl;
    std::cout << "Number of simHits in simLayerCluster :" << simClusterDetIds.size() << std::endl;
    for (auto it: simClusterHitsAndEnergies){
      std::cout << " Sim cluster hit " << (uint32_t)it.first << "   hit energy : " << it.second << std::endl;
      //mtdSimClusterInfo->clusHitEnergy.push_back( it.second );
    }
    
    mtdSimClusterInfo->energy = simClus.simLCEnergy();        
    mtdSimClusterInfo->nClusHits = simClusterDetIds.size();


    bool isInBTL = true;
    
    // -- loop over detIds of this cluster
    for (auto simClusterDetId : simClusterDetIds){
      DetId thisDetId = simClusterDetId.first;

      // -- check if is BTL
      MTDDetId mtdId = MTDDetId(thisDetId);
      if (mtdId.mtdSubDetector() != MTDDetId::BTL) {
	isInBTL = false;
	continue;
      }

      std::cout << "  Sim hit DetID : " << simClusterDetId.first
        	<< "  row, col : " << (uint32_t)((simClusterDetId.second).first) << "," << (uint32_t)((simClusterDetId.second).second)
      	        << std::endl;

      for (auto it: simClusterHitsAndEnergies){
	std::cout << (uint32_t)it.first << "   " << (uint32_t)(simClusterDetId.second).second <<std::endl;
	if ((uint32_t)it.first == (uint32_t)(simClusterDetId.second).second) mtdSimClusterInfo->clusHitEnergy.push_back( it.second );
    }

          
      mtdSimClusterInfo->simHitDetId.push_back((int)thisDetId);
	
      for (auto const ene : map_btlSimHitEnergy[thisDetId]){
	 std::cout <<  "  --> simHit energy: " <<  ene << std::endl;
	 mtdSimClusterInfo->simHitEnergy.push_back(ene);
      }
	  
      for (auto const tkId : map_btlSimHitTrackId[thisDetId]){
	std::cout <<  "  --> simHit trackId: " <<  tkId << std::endl;
	mtdSimClusterInfo->simHitTrackId.push_back(tkId);
      }

      for (auto const tkId : map_btlSimHitOrigTrackId[thisDetId]){
	std::cout <<  "  --> simHit original trackId: " <<  tkId << std::endl;
	mtdSimClusterInfo->simHitOriginalTrackId.push_back(tkId);
      }

      for (auto const tkId : map_btlSimHitOffsetTrackId[thisDetId]){
	std::cout <<  "  --> simHit offset trackId: " <<  tkId << std::endl;
	mtdSimClusterInfo->simHitOffsetTrackId.push_back(tkId);
      }

    }

    // -- fill the tree for each cluster
    if (isInBTL){
      tree_->Fill();
    }

    std::cout << "  " <<std::endl;
    
  }//---end loop over sim clusters


 }

void MtdSimClusterDumper::clearTreeVariables() {

  mtdSimClusterInfo->energy    = -999.;
  mtdSimClusterInfo->nClusHits = 0;
  mtdSimClusterInfo->clusHitEnergy.clear();
  mtdSimClusterInfo->simHitDetId.clear();
  mtdSimClusterInfo->simHitEnergy.clear();
  mtdSimClusterInfo->simHitTrackId.clear();
  mtdSimClusterInfo->simHitOriginalTrackId.clear();
  mtdSimClusterInfo->simHitOffsetTrackId.clear();
  
}

void MtdSimClusterDumper::endJob() {}

void MtdSimClusterDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<double>("hitMinEnergy", 1.);  // [MeV]
  desc.add<edm::InputTag>("btlSimHitsTag", edm::InputTag("mix:g4SimHitsFastTimerHitsBarrel"));
  desc.add<edm::InputTag>("mtdSimLayerClustersTag", edm::InputTag("mix:MergedMtdTruthLC"));
    
}

DEFINE_FWK_MODULE(MtdSimClusterDumper);
