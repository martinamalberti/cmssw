//
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl.h"


using namespace reco;
using namespace std;


/* Constructor */

MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl::MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl(edm::EDProductGetter const& productGetter,
													 const edm::Handle<FTLRecHitCollection>& btlRecHitsH,
													 const edm::Handle<FTLRecHitCollection>& etlRecHitsH,
													 double energyCut,
													 double timeCut)
: productGetter_(&productGetter), btlRecHitsH_(btlRecHitsH), etlRecHitsH_(etlRecHitsH), energyCut_(energyCut), timeCut_(timeCut){}


//
//---member functions
//

reco::RecoToSimCollectionMtd MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl::associateRecoToSim(
      const edm::Handle<FTLClusterCollection>& btlRecoClusH, const edm::Handle<FTLClusterCollection>& etlRecoClusH,
      const edm::Handle<MtdSimLayerClusterCollection>& simClusH) const {  
												     
  RecoToSimCollectionMtd outputCollection;

  // -- get the collections
  const auto& simClusters  = *simClusH.product();
    
  std::array<edm::Handle<FTLClusterCollection>, 2> inputRecoClusH{{btlRecoClusH, etlRecoClusH}};
  std::array<edm::Handle<FTLRecHitCollection>, 2> inputRecHitH{{btlRecHitsH_, etlRecHitsH_}};
 
  for (auto const& recoClusH : inputRecoClusH) {
    // -- loop over detSetVec
    for (const auto& detSet : *recoClusH) { 
      // -- loop over reco clusters
      for (const auto& recoClus : detSet) {

	MTDDetId clusId = recoClus.id();
	
	LogDebug("MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl") << "Reco cluster : " << clusId;
	
	std::vector<uint64_t> recoClusHitIds;

	// -- loop over hits in the reco cluster and find their ids
	for (int ihit = 0; ihit < recoClus.size(); ++ihit) {
	  int hit_row = recoClus.minHitRow() + recoClus.hitOffset()[ihit * 2];
	  int hit_col = recoClus.minHitCol() + recoClus.hitOffset()[ihit * 2 + 1];
	  
	  for (const auto& recHitsH : inputRecHitH) {
	    for (auto recHit : *recHitsH) {
	      MTDDetId hitId(recHit.id().rawId());
		    
	      // -- check the hit position
	      if (hitId.mtdSide() != clusId.mtdSide() || hitId.mtdRR() != clusId.mtdRR() || recHit.row() != hit_row || recHit.column() != hit_col)
		continue;
	    
	      // -- check the hit energy and time
	      if (recHit.energy() != recoClus.hitENERGY()[ihit] || recHit.time() != recoClus.hitTIME()[ihit])
		continue;
	      
	      // -- Get an unique id: for BTL the detId is unique (one for each crystal), for ETL the detId is not enough
	      //    also row and column are needed. An unique number is created from detId, row, col
	      uint64_t uniqueId = static_cast<uint64_t>(hitId.rawId()) << 32;
	      uniqueId |= recHit.row() << 16;
	      uniqueId |= recHit.column();
	      recoClusHitIds.push_back(uniqueId);
	      LogDebug("MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl") << " ======= recHit raw Id : " << hitId.rawId()
									      << "  row : " << recHit.row() << "   column: " <<  recHit.column()
									      << " recHit uniqueId : " << uniqueId;

	    }
	  }
	}// end loop over rec hits in reco cluster
	
	// -- loop over sim clusters and if this reco clus shares some hits
	edm::Ref<MtdSimLayerClusterCollection>::key_type simClusIndex = 0;
	std::vector<MtdSimLayerClusterRef> simClusterRefs;
	simClusterRefs.clear();

	for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){
	  auto simClus = *simClusIt;

	  std::vector<std::pair<uint64_t, float>> hitsAndFrac = simClus.hits_and_fractions();
	  std::vector<uint64_t> simClusHitIds(hitsAndFrac.size());
	  std::transform(hitsAndFrac.begin(), hitsAndFrac.end(), simClusHitIds.begin(), [](const std::pair<uint64_t, float>& pair) {
											  return pair.first;});

	  // -- Get shared hits
	  std::vector<uint64_t> sharedHitIds;
	  std::set_intersection(recoClusHitIds.begin(), recoClusHitIds.end(), simClusHitIds.begin(), simClusHitIds.end(), std::back_inserter(sharedHitIds));

	  float dE = recoClus.energy()*0.001/simClus.simLCEnergy(); // reco cluster energy is in MeV
	  float dtSig = std::abs((recoClus.time()-simClus.simLCTime())/recoClus.timeError());

	  // -- If the sim and reco clusters have common hits, fill the std:vector of sim clusters refs 
	  if (!sharedHitIds.empty() && dE < energyCut_ && dtSig < timeCut_){ // at least one hit in common + requirement on energy and time compatibility  	

	    edm::Ref<MtdSimLayerClusterCollection> simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIndex); 
	    simClusterRefs.push_back(simClusterRef);

	    LogDebug("MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl") << "RecoToSim --> Found " << sharedHitIds.size() << " shared hits";
	    LogDebug("MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl") << "E_recoClus = " << recoClus.energy() << "   E_simClus = " << simClus.simLCEnergy() << "   E_recoClus/E_simClus = " << dE;
	    LogDebug("MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl") << "(t_recoClus-t_simClus)/sigma_t = " << dtSig;
	    
	  }
	  
	  simClusIndex++;
	  
	}// -- end loop over sim clus
	
	// -- Now fill the output collection
	edm::Ref<edmNew::DetSetVector<FTLCluster>, FTLCluster> recoClusterRef = edmNew::makeRefTo(recoClusH, &recoClus); 
	outputCollection.emplace_back(recoClusterRef, simClusterRefs);
	  
      }// -- end loop over reco clus
    }// -- end loop over detsetclus
  }	  
    
  return outputCollection;

}



reco::SimToRecoCollectionMtd MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl::associateSimToReco(
      const edm::Handle<FTLClusterCollection>& btlRecoClusH, const edm::Handle<FTLClusterCollection>& etlRecoClusH,
      const edm::Handle<MtdSimLayerClusterCollection>& simClusH) const {

  SimToRecoCollectionMtd outputCollection;
  
  // -- get the collections
  const auto& simClusters  = *simClusH.product();
  
  std::array<edm::Handle<FTLClusterCollection>, 2> inputH{{btlRecoClusH, etlRecoClusH}};
  std::array<edm::Handle<FTLRecHitCollection>, 2> inputRecHitH{{btlRecHitsH_, etlRecHitsH_}};
  
  // -- loop over MtdSimLayerClusters
  edm::Ref<MtdSimLayerClusterCollection>::key_type simClusIndex = 0;
  
  for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){
    auto simClus = *simClusIt;
    
    std::vector<std::pair<uint64_t, float>> hitsAndFrac = simClus.hits_and_fractions();
    std::vector<uint64_t> simClusHitIds(hitsAndFrac.size());
    std::transform(hitsAndFrac.begin(), hitsAndFrac.end(), simClusHitIds.begin(), [](const std::pair<uint64_t, float>& pair) {
										    return pair.first;});

    std::vector<FTLClusterRef> recoClusterRefs;
    recoClusterRefs.clear();
    
    // -- loop over reco clusters
    for (auto const& recoClusH : inputH) {
      // -- loop over detSetVec
      for (const auto& detSet : *recoClusH) {
	// -- loop over reco clusters
	for (const auto& recoClus : detSet) {
 	  MTDDetId clusId = recoClus.id();

	  std::vector<uint64_t> recoClusHitIds;
 
	  // -- loop over hits in the reco cluster and find their ids
	  for (int ihit = 0; ihit < recoClus.size(); ++ihit) {
	    int hit_row = recoClus.minHitRow() + recoClus.hitOffset()[ihit * 2];
	    int hit_col = recoClus.minHitCol() + recoClus.hitOffset()[ihit * 2 + 1];
	    
	    for (const auto& recHitsH : inputRecHitH) {
	      for (auto recHit : *recHitsH) {
		MTDDetId hitId(recHit.id().rawId());
		
		// -- check the hit position
		if (hitId.mtdSide() != clusId.mtdSide() || hitId.mtdRR() != clusId.mtdRR() || recHit.row() != hit_row || recHit.column() != hit_col)
		  continue;
		
		// -- check the hit energy and time
		if (recHit.energy() != recoClus.hitENERGY()[ihit] || recHit.time() != recoClus.hitTIME()[ihit])
		  continue;
		
		// -- Get an unique id: for BTL the detId is unique (one for each crystal), for ETL the detId is not enough
		//    also row and column are needed. An unique number is created from detId, row, col
		uint64_t uniqueId = static_cast<uint64_t>(hitId.rawId()) << 32;
		uniqueId |= recHit.row() << 16;
		uniqueId |= recHit.column();
		recoClusHitIds.push_back(uniqueId);
	      }
	    }
	  } // -- end loop over hits in the reco clus

	  // -- Get shared hits
          std::vector<uint64_t> sharedHitIds;
          std::set_intersection(simClusHitIds.begin(), simClusHitIds.end(), recoClusHitIds.begin(), recoClusHitIds.end(), std::back_inserter(sharedHitIds));
	  
          float dE = recoClus.energy()*0.001/simClus.simLCEnergy(); // reco cluster energy is in MeV
          float dtSig = std::abs((recoClus.time()-simClus.simLCTime())/recoClus.timeError());
	  
          // -- If the sim and reco clusters have common hits, fill the std:vector of reco clusters refs
          if (!sharedHitIds.empty() && dE < energyCut_ && dtSig < timeCut_){
	  
	    // Create a persistent edm::Ref to the cluster
	    edm::Ref<edmNew::DetSetVector<FTLCluster>, FTLCluster> recoClusterRef = edmNew::makeRefTo(recoClusH, &recoClus);
	    recoClusterRefs.push_back(recoClusterRef);

            LogDebug("MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl") << "SimToReco --> Found " << sharedHitIds.size() << " shared hits";
	    LogDebug("MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl") << "E_recoClus = " << recoClus.energy() << "   E_simClus = " << simClus.simLCEnergy() << "   E_recoClus/E_simClus = " << dE;
	    LogDebug("MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl") << "(t_recoClus-t_simClus)/sigma_t = " << dtSig;
	  }
	  
	}// end loop ove reco clus
      }// end loop over detsets
    }

    // -- Now fill the output collection
    edm::Ref<MtdSimLayerClusterCollection> simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIndex); 
    outputCollection.emplace_back(simClusterRef, recoClusterRefs);
    
    simClusIndex++;
    
  } // -- end loop over sim clusters
  
  return outputCollection;
}
