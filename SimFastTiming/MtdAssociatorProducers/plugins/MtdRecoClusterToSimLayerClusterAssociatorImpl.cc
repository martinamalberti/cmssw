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

  // do stuff to fill the output collection

  // -- get the collections
  const auto& btlRecHits   = *btlRecHitsH.product();
  const auto& etlRecHits   = *etlRecHitsH.product();
  const auto& simClusters  = *simClusH.product();
    
  std::array<edm::Handle<FTLClusterCollection>, 2> inputH{{btlRecoClusH, etlRecoClusH}};

  for (auto const& recoClusH : inputH) {
    
    const auto& detSetVecs = *recoClusH.product();

    // -- loop over detSetVec
    for (auto detSetVecIt = detSetVecs.begin(); detSetVecIt != detSetVecs.end(); detSetVecIt++) {
      
      auto detSetVec = *detSetVecIt;

      // -- loop over reco clusters
      for (const auto& recoClus : detSetVec) {

	BTLDetId clusId = recoClus.id(); // DetId dei reco clus e' il sensor module
	MTDDetId mtdId(clusId);

	// === Clusters in BTL   -- testing only BTL for now.
	if (mtdId.mtdSubDetector() != MTDDetId::BTL)  continue;

	std::vector<uint64_t> recoClusHitIds;

	// -- loop over hits in the reco cluster and find their ids
	for (int ihit = 0; ihit < recoClus.size(); ++ihit) {
	  int hit_row = recoClus.minHitRow() + recoClus.hitOffset()[ihit * 2];
	  int hit_col = recoClus.minHitCol() + recoClus.hitOffset()[ihit * 2 + 1];
	  
	  for (auto recHit : btlRecHits) {
	    BTLDetId hitId(recHit.id().rawId());
	    
	    // -- check the hit position
	    if (hitId.mtdSide() != clusId.mtdSide() || hitId.mtdRR() != clusId.mtdRR() || recHit.row() != hit_row || recHit.column() != hit_col)
	      continue;
	    
	    // -- check the hit energy and time
	    if (recHit.energy() != recoClus.hitENERGY()[ihit] || recHit.time() != recoClus.hitTIME()[ihit])
	      continue;
	    
	    recoClusHitIds.push_back(hitId);
	  }
	} // end loop over hits in reco cluster

	// -- loop over sim clusters and if this reco clus shares some hits
	edm::Ref<MtdSimLayerClusterCollection>::key_type simClusIndex = 0;
	float quality = 0;
	int nSharedHits = 0;
	//for (const auto& simClus  : simClusters){
	for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){
	  auto simClus = *simClusIt;
	  simClusIndex++;
	  std::vector<std::pair<uint64_t, float>> hitsAndFrac = simClus.hits_and_fractions();
	  std::vector<uint64_t> simClusHitIds(hitsAndFrac.size());
	  std::transform(hitsAndFrac.begin(), hitsAndFrac.end(), simClusHitIds.begin(), [](const std::pair<int, float>& pair) {
											  return pair.first;});
	  std::vector<uint64_t> sharedHitIds;
	  std::set_intersection(recoClusHitIds.begin(), recoClusHitIds.end(), simClusHitIds.begin(), simClusHitIds.end(), std::back_inserter(sharedHitIds));
	  if (!sharedHitIds.empty()){ 	// NB : may add some requirement on energy and/or time compatibility between the sim cluster and the reco cluster
	    nSharedHits = sharedHitIds.size();
	    quality = sharedHitIds.size()/recoClusHitIds.size();
	    break; 
	  }
	}

	// -- if they share at least one hit fill the output collection
	if ( nSharedHits > 0 ){ // at least one hit in common
	  edm::Ref<MtdSimLayerClusterCollection> simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIndex); // OK
	  
	  // Create a persistent edm::Ref to the cluster
	  // --> voglio : edm::Ref<edmNew::DetSetVector<FTLCluster>, edmNew::DetSet<FTLCluster>
	  edm::Ref<edmNew::DetSetVector<FTLCluster>, FTLCluster> recoClusterRef = edmNew::makeRefTo(recoClusH, &recoClus);
	  outputCollection.insert(recoClusterRef, std::make_pair(simClusterRef, quality));
	  
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

  // do stuff to fill the output collection
  


  return outputCollection;
}
