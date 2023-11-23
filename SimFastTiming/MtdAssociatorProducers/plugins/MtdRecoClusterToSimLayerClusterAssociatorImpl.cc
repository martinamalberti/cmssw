//
//
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "MtdRecoClusterToSimLayerClusterAssociatorImpl.h"

using namespace reco;
using namespace std;


/* Constructor */

MtdRecoClusterToSimLayerClusterAssociatorImpl::MtdRecoClusterToSimLayerClusterAssociatorImpl(edm::EDProductGetter const& productGetter)
  : productGetter_(&productGetter){}


//
//---member functions
//

reco:RecoToSimCollectionMtd MtdRecoClusterToSimLayerClusterAssociatorImpl::associateRecoToSim(
      const edm::Handle<FTLClusterCollection>& rCCH, const edm::Handle<MtdSimLayerClusterCollection>& sCCH,
      const edm::Handle<FTLRecHitCollection>& btlRecHitsHandle, const edm::Handle<FTLRecHitCollection>& etlRecHitsHandle) const {

  RecoToSimCollectionMtd outputCollection(productGetter_);

  // do stuff to fill the output collection

  // -- get collections
  const auto& recoClusters = *rCCH.product();
  const auto& simClusters  = *sCCH.product();
  const auto& btlRecHits   = *btlRecHitsHandle.product();
  const auto& etlRecHits   = *etlRecHitsHandle.product();

  // -- loop over reco clusters
  for (auto recoClus : recoClusters){

    std::vector<uint64_t> recoClusHitIds;

    BTLDetId clusId = cluster.id();
    MTDDetId mtdDetId(clusId);

    // === Clusters in BTL 
    if (mtdId.mtdSubDetector() != MTDDetId::BTL)  continue;
 
    // -- loop over hits in the reco cluster and find their ids
    for (int ihit = 0; ihit < cluster.size(); ++ihit) {
      int hit_row = cluster.minHitRow() + cluster.hitOffset()[ihit * 2];
      int hit_col = cluster.minHitCol() + cluster.hitOffset()[ihit * 2 + 1];

      for (auto recHit : btlRecHits) {
	BTLDetId hitId(recHit.id().rawId());
	
	// -- check the hit position
	if (hitId.mtdSide() != clusId.mtdSide() || hitId.mtdRR() != clusId.mtdRR() || recHit.row() != hit_row || recHit.column() != hit_col)
	  continue;
 
        // -- check the hit energy and time
	if (recHit.energy() != cluster.hitENERGY()[ihit] || recHit.time() != cluster.hitTIME()[ihit])
	  continue;

	recoClusHitIds.push_back(hitId);
      }
    }

    // -- loop over sim clusters and if this reco clus shares some hits
    int simClusIndex = -1;
    for (auto simClus : simClusters){
      simClusIndex++;
      std::vector<std::pair<uint64_t, float>> hitsAndFrac = simClus.hits_and_fractions();
      //const auto& hitsAndFrac = simClus.hits_and_fractions();
      std::vector<uint64_t> simClusHitIds(hitsAndFrac.size());
      std::transform(hitsAndFrac.begin(), hitsAndFrac.end(), simClusIds.begin(), [](const std::pair<int, float>& pair) {
          return pair.first;
      });
      std::vector<uint64_t> sharedHitIds;
      std::set_intersection(recoClusHitIds.begin(), recoClusHitIds.end(), simClusHitIds.begin(), simClusHitIds.end(), std::back_inserter(sharedHitIds));
      if (sharedHitIds.empty()) continue;
    }

    // -- if they share some hits fill the output collection
    if ( simClusIndex != -1){
      outputCollection.insert(recoClus, simClus[simClusIndex]);
    }
    
  }// -- end loop over reco clusters
  
  return outputCollection;
}



reco:SimToRecoCollectionMtd MtdRecoClusterToSimLayerClusterAssociatorImpl::associateSimToReco(
  const edm::Handle<FTLlusterCollection>& rCCH, const edm::Handle<MtdSimLayerClusterCollection>& sCCH) const {

  SimToRecoCollectionMtd outputCollection(productGetter_);

  // do stuff to fill the output collection
  


  return outputCollection;
}
