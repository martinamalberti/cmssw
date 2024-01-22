#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"

using namespace reco;
using namespace std;


/* Constructor */

MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl::MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl(edm::EDProductGetter const& productGetter,
													 double energyCut,
													 double timeCut)
: productGetter_(&productGetter),
  energyCut_(energyCut),
  timeCut_(timeCut){}


//
//---member functions
//

reco::RecoToSimCollectionMtd MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl::associateRecoToSim(
      const edm::Handle<FTLClusterCollection>& btlRecoClusH, const edm::Handle<FTLClusterCollection>& etlRecoClusH,
      const edm::Handle<MtdSimLayerClusterCollection>& simClusH) const {  
												     
  RecoToSimCollectionMtd outputCollection;
    
  // -- get the collections
  std::array<edm::Handle<FTLClusterCollection>, 2> inputRecoClusH{{btlRecoClusH, etlRecoClusH}};

  const auto& simClusters  = *simClusH.product();
  
  // -- create temporary map  DetId, SimClusterRef (praticamente ... il DetSetVector dei poveri)
  std::map< uint32_t, std::vector<MtdSimLayerClusterRef>> simClusIdsMap;
  for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){

    auto simClus = *simClusIt;

    edm::Ref<MtdSimLayerClusterCollection> simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIt-simClusters.begin());

    std::vector<std::pair<uint32_t, std::pair<uint8_t, uint8_t>>> detIdsAndRows = simClus.detIds_and_rows();
    std::vector<uint32_t> detIds(detIdsAndRows.size());
    std::transform(detIdsAndRows.begin(), detIdsAndRows.end(), detIds.begin(), [](const std::pair<uint32_t, std::pair<uint8_t, uint8_t>>& pair) {
                                                                                 return pair.first;} );
    
    simClusIdsMap[detIds[0]].push_back(simClusterRef);                                                                                                                                                                  
    LogDebug("MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl")<< "Sim cluster  detId = " << detIds[0] << std::endl;   
    
  }


  
  for (auto const& recoClusH : inputRecoClusH) {
    // -- loop over detSetVec
    for (const auto& detSet : *recoClusH) { 
      // -- loop over reco clusters
      for (const auto& recoClus : detSet) {

	MTDDetId clusId = recoClus.id();
	
	LogDebug("MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl") << "Reco cluster : " << clusId;
	
	std::vector<uint64_t> recoClusHitIds;

	// -- loop over hits in the reco cluster and find their unique ids
	for (int ihit = 0; ihit < recoClus.size(); ++ihit) {
	  int hit_row = recoClus.minHitRow() + recoClus.hitOffset()[ihit * 2];
	  int hit_col = recoClus.minHitCol() + recoClus.hitOffset()[ihit * 2 + 1];

	  // -- Get an unique id from sensor module detId , row, column
	  uint64_t uniqueId = static_cast<uint64_t>(clusId.rawId()) << 32;
	  uniqueId |= hit_row << 16;
	  uniqueId |= hit_col;
	  recoClusHitIds.push_back(uniqueId);

	  LogDebug("MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl") << " ======= cluster raw Id : " << clusId.rawId()
									  << "  row : " << hit_row << "   column: " <<  hit_col
									  << " recHit uniqueId : " << uniqueId;
	}
	

	// -- loop over sim clusters and check if this reco clus shares some hits
        std::vector<MtdSimLayerClusterRef> simClusterRefs;
		
	for  ( auto simClusterRef : simClusIdsMap[clusId.rawId()] ){
	  auto simClus = *simClusterRef;

	  std::vector<std::pair<uint64_t, float>> hitsAndFrac = simClus.hits_and_fractions();
	  std::vector<uint64_t> simClusHitIds(hitsAndFrac.size());
	  std::transform(hitsAndFrac.begin(), hitsAndFrac.end(), simClusHitIds.begin(), [](const std::pair<uint64_t, float>& pair) {
											  return pair.first;});
	  
	  // -- Get shared hits
	  std::vector<uint64_t> sharedHitIds;
	  std::set_intersection(recoClusHitIds.begin(), recoClusHitIds.end(), simClusHitIds.begin(), simClusHitIds.end(), std::back_inserter(sharedHitIds));

	  if (sharedHitIds.empty()) continue;
	  	  
	  float dE = recoClus.energy()*0.001/simClus.simLCEnergy(); // reco cluster energy is in MeV!
	  float dtSig = std::abs((recoClus.time()-simClus.simLCTime())/recoClus.timeError());

	  // -- If the sim and reco clusters have common hits, fill the std:vector of sim clusters refs 
	  if (!sharedHitIds.empty() && dE < energyCut_ && dtSig < timeCut_){ // at least one hit in common + requirement on energy and time compatibility  	
	  
	    simClusterRefs.push_back(simClusterRef);

	    LogDebug("MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl") << "RecoToSim --> Found " << sharedHitIds.size() << " shared hits";
	    LogDebug("MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl") << "E_recoClus = " << recoClus.energy() << "   E_simClus = " << simClus.simLCEnergy() << "   E_recoClus/E_simClus = " << dE;
	    LogDebug("MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl") << "(t_recoClus-t_simClus)/sigma_t = " << dtSig;
	    
	  }
	} // -- end loop over sim clus refs
	
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
  
  // -- loop over MtdSimLayerClusters
  for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){
    auto simClus = *simClusIt;

    // - get the uniqueIds of the hits in the sim layer cluster
    std::vector<std::pair<uint64_t, float>> hitsAndFrac = simClus.hits_and_fractions();
    std::vector<uint64_t> simClusHitIds(hitsAndFrac.size());
    std::transform(hitsAndFrac.begin(), hitsAndFrac.end(), simClusHitIds.begin(), [](const std::pair<uint64_t, float>& pair) {return pair.first;} );

    std::vector<std::pair<uint32_t, std::pair<uint8_t, uint8_t>>> detIdsAndRows = simClus.detIds_and_rows();
    std::vector<uint32_t> detIds(detIdsAndRows.size());
    std::transform(detIdsAndRows.begin(), detIdsAndRows.end(), detIds.begin(), [](const std::pair<uint32_t, std::pair<uint8_t, uint8_t>>& pair) {return pair.first;});

    DetId simClusId  = detIds[0];
    
    std::vector<FTLClusterRef> recoClusterRefs;
    
    // -- loop over reco clusters
    for (auto const& recoClusH : inputH) {
      // -- loop over detSetVec
      for (const auto& detSet : *recoClusH) {

	if (detSet.id() != simClusId.rawId())  continue;
	
	// -- loop over reco clusters
	for (const auto& recoClus : detSet) {
 	  MTDDetId clusId = recoClus.id();
	  
	  std::vector<uint64_t> recoClusHitIds;
	  
	  // -- loop over hits in the reco cluster and find their ids
	  for (int ihit = 0; ihit < recoClus.size(); ++ihit) {
	    int hit_row = recoClus.minHitRow() + recoClus.hitOffset()[ihit * 2];
	    int hit_col = recoClus.minHitCol() + recoClus.hitOffset()[ihit * 2 + 1];

	  // -- Get an unique id from sensor module detId , row, column  
	  uint64_t uniqueId = static_cast<uint64_t>(clusId.rawId()) << 32;
	  uniqueId |= hit_row << 16;
	  uniqueId |= hit_col;
	  recoClusHitIds.push_back(uniqueId);

	  LogDebug("MtdRecoClusterToSimLayerClusterAssociatorByHitsImpl") << " ======= cluster raw Id : " << clusId.rawId()
									  << "  row : " << hit_row << "   column: " <<  hit_col
									  << " recHit uniqueId : " << uniqueId;

	  }
	  
	  // -- Get shared hits
          std::vector<uint64_t> sharedHitIds;
          std::set_intersection(simClusHitIds.begin(), simClusHitIds.end(), recoClusHitIds.begin(), recoClusHitIds.end(), std::back_inserter(sharedHitIds));
	  if (sharedHitIds.empty()) continue;
	  
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
    edm::Ref<MtdSimLayerClusterCollection> simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIt-simClusters.begin()); 
    outputCollection.emplace_back(simClusterRef, recoClusterRefs);
      
  } // -- end loop over sim clusters
  
  return outputCollection;

}
