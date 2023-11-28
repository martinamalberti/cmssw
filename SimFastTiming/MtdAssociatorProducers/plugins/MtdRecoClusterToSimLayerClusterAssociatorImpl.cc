//
//
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "MtdRecoClusterToSimLayerClusterAssociatorImpl.h"



using namespace reco;
using namespace std;


/* Constructor */

MtdRecoClusterToSimLayerClusterAssociatorImpl::MtdRecoClusterToSimLayerClusterAssociatorImpl(edm::EDProductGetter const& productGetter)
  : productGetter_(&productGetter){}


//
//---member functions
//

reco::RecoToSimCollectionMtd MtdRecoClusterToSimLayerClusterAssociatorImpl::associateRecoToSim(
      const edm::Handle<FTLClusterCollection>& btlRecoClusH, const edm::Handle<FTLClusterCollection>& etlRecoClusH,
      const edm::Handle<MtdSimLayerClusterCollection>& simClusH,
      const edm::Handle<FTLRecHitCollection>& btlRecHitsH, const edm::Handle<FTLRecHitCollection>& etlRecHitsH) const {

  RecoToSimCollectionMtd outputCollection(productGetter_);

  // -- get the collections
  //const auto& btlRecHits   = *btlRecHitsH.product();
  //const auto& etlRecHits   = *etlRecHitsH.product();
  const auto& simClusters  = *simClusH.product();
    
  std::array<edm::Handle<FTLClusterCollection>, 2> inputRecoClusH{{btlRecoClusH, etlRecoClusH}};
  std::array<edm::Handle<FTLRecHitCollection>, 2> inputRecHitH{{btlRecHitsH, etlRecHitsH}};
 
  for (auto const& recoClusH : inputRecoClusH) {

    const auto& detSetVecs = *recoClusH.product();

    // -- loop over detSetVec
    for (auto detSetVecIt = detSetVecs.begin(); detSetVecIt != detSetVecs.end(); detSetVecIt++) {
      
      auto detSetVec = *detSetVecIt;

      // -- loop over reco clusters
      for (const auto& recoClus : detSetVec) {

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
	}// end loop over rec hits in reco cluster

	// -- loop over sim clusters and if this reco clus shares some hits
	edm::Ref<MtdSimLayerClusterCollection>::key_type simClusIndex = 0;
	int nSharedHits = 0;
	float dE = 0;
	float dt = 0;
	for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){
	  auto simClus = *simClusIt;
	  std::vector<std::pair<uint64_t, float>> hitsAndFrac = simClus.hits_and_fractions();
	  std::vector<uint64_t> simClusHitIds(hitsAndFrac.size());
	  std::transform(hitsAndFrac.begin(), hitsAndFrac.end(), simClusHitIds.begin(), [](const std::pair<int, float>& pair) {
											  return pair.first;});
	  // -- get shared hits
	  std::vector<uint64_t> sharedHitIds;
	  std::set_intersection(recoClusHitIds.begin(), recoClusHitIds.end(), simClusHitIds.begin(), simClusHitIds.end(), std::back_inserter(sharedHitIds));
	  nSharedHits = sharedHitIds.size();

	  // -- may add some requirement on energy and/or time compatibility between the sim cluster and the reco cluster
	  dE = recoClus.energy()/simClus.simLCEnergy();
	  dt = (recoClus.time()-simClus.simLCTime())/recoClus.time();

	  simClusIndex++;

	  if (!sharedHitIds.empty() ){ 	
	    break; 
	  }
	}// -- end loop over sim clus

	
	// -- if they share at least one hit fill the output collection
	if ( nSharedHits > 0 ){ // at least one hit in common
	  //if ( nSharedHits > 0 && dE < 10 && dt <10){ // at least one hit in common + compat.
	  edm::Ref<MtdSimLayerClusterCollection> simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIndex); // this works.
	  	  
	  // Create a persistent edm::Ref to the cluster
	  // --> voglio : edm::Ref<edmNew::DetSetVector<FTLCluster>, edmNew::DetSet<FTLCluster> ?????
	  edm::Ref<edmNew::DetSetVector<FTLCluster>, FTLCluster> recoClusterRef = edmNew::makeRefTo(recoClusH, &recoClus); 
	  outputCollection.insert(recoClusterRef, std::make_pair(simClusterRef, float(nSharedHits)));
	  
	}
      }// -- end loop over reco clus
    }// -- end loop over detsetclus
  }	  
    
  return outputCollection;

}



reco::SimToRecoCollectionMtd MtdRecoClusterToSimLayerClusterAssociatorImpl::associateSimToReco(
      const edm::Handle<FTLClusterCollection>& btlRecoClusH, const edm::Handle<FTLClusterCollection>& etlRecoClusH,
      const edm::Handle<MtdSimLayerClusterCollection>& simClusH,
      const edm::Handle<FTLRecHitCollection>& btlRecHitsH, const edm::Handle<FTLRecHitCollection>& etlRecHitsH) const {

  SimToRecoCollectionMtd outputCollection(productGetter_);
  
  // -- get the collections
  const auto& simClusters  = *simClusH.product();
  
  std::array<edm::Handle<FTLClusterCollection>, 2> inputH{{btlRecoClusH, etlRecoClusH}};
  std::array<edm::Handle<FTLRecHitCollection>, 2> inputRecHitH{{btlRecHitsH, etlRecHitsH}};
  
  // -- loop over MtdSimLayerClusters
  edm::Ref<MtdSimLayerClusterCollection>::key_type simClusIndex = 0;
  
  for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){

    auto simClus = *simClusIt;
    
    std::vector<std::pair<uint64_t, float>> hitsAndFrac = simClus.hits_and_fractions();
    std::vector<uint64_t> simClusHitIds(hitsAndFrac.size());
    std::transform(hitsAndFrac.begin(), hitsAndFrac.end(), simClusHitIds.begin(), [](const std::pair<int, float>& pair) {
										    return pair.first;});

    // -- loop over reco clusters
    for (auto const& recoClusH : inputH) {

      const auto& detSetVecs = *recoClusH.product();

      // -- loop over detSetVec
      for (auto detSetVecIt = detSetVecs.begin(); detSetVecIt != detSetVecs.end(); detSetVecIt++) {
      
	auto detSetVec = *detSetVecIt;
    
	// -- loop over reco clusters
	for (const auto& recoClus : detSetVec) {

	  MTDDetId clusId = recoClus.id();

	  int nSharedHits = 0;

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
		if (std::find(simClusHitIds.begin(), simClusHitIds.end(), uniqueId) != simClusHitIds.end())
		  nSharedHits++;
	      }
	    }
	  } // -- end loop over hits in the reco clus

	  // -- may add some requirement on energy and/or time compatibility between the sim cluster and the reco cluster
          //float dE = recoClus.energy()/simClus.simLCEnergy();
          //float dt = (recoClus.time()-simClus.simLCTime())/recoClus.time();
	  
	  // -- if the reco clus shares at least one hit with the sim clus fill the output collection
	   if ( nSharedHits > 0 ){
	     //if ( nSharedHits > 0 && dE < 10 && dt <10){ // at least one hit in common + compat.
	    edm::Ref<MtdSimLayerClusterCollection> simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIndex); // OK
	    
	    // Create a persistent edm::Ref to the cluster
	    edm::Ref<edmNew::DetSetVector<FTLCluster>, FTLCluster> recoClusterRef = edmNew::makeRefTo(recoClusH, &recoClus);
	    outputCollection.insert(simClusterRef, std::make_pair(recoClusterRef, float(nSharedHits)));
	  }
      
	  
	  ///	  
	}
      }
    }

    simClusIndex++;
    
  } // -- end loop over sim clusters
  
  return outputCollection;
}
