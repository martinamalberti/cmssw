//
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "MtdSimLayerClusterToTPAssociatorByTrackIdImpl.h"

#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

#include "SimFastTiming/MtdAssociatorProducers/interface/MtdAssociatorTools.h"

using namespace reco;
using namespace std;


/* Constructor */

/*MtdSimLayerClusterToTPAssociatorByTrackIdImpl::MtdSimLayerClusterToTPAssociatorByTrackIdImpl(edm::EDProductGetter const& productGetter,
  const edm::Handle<CrossingFrame<PSimHit>> &btlSimHitsH,
  const edm::Handle<CrossingFrame<PSimHit>> &etlSimHitsH,
  mtd::MTDGeomUtil &geomTools)
  : productGetter_(&productGetter), btlSimHitsH_(btlSimHitsH), etlSimHitsH_(etlSimHitsH), geomTools_(geomTools) {}
*/


MtdSimLayerClusterToTPAssociatorByTrackIdImpl::MtdSimLayerClusterToTPAssociatorByTrackIdImpl(edm::EDProductGetter const& productGetter)
  : productGetter_(&productGetter) {}


//
//---member functions
//

reco::SimToTPCollectionMtd MtdSimLayerClusterToTPAssociatorByTrackIdImpl::associateSimToTP(
      const edm::Handle<MtdSimLayerClusterCollection>& simClusH,
      const edm::Handle<TrackingParticleCollection>& trackingParticleH) const {

  
  SimToTPCollectionMtd outputCollection(productGetter_);

  // -- get the collections
  const auto& simClusters  = *simClusH.product();
  const auto& trackingParticles = *trackingParticleH.product();

  // -- Loop over tracking particles and build a temporary map of trackId, eventId  --> tpRef
  //std::map<int, std::pair<unsigned int, uint32_t>> tpIdMap;
  std::map< std::pair<unsigned int, uint32_t>, TrackingParticleRef> tpIdMap;
  size_t tpIndex(0);
  for (auto tpIt = trackingParticles.begin(); tpIt != trackingParticles.end(); tpIt++) {
    //std::cout << tpIndex << std::endl;
    auto tp = *tpIt;
    unsigned int tpTrackId = tp.g4Tracks()[0].trackId();
    EncodedEventId tpEventId = tp.eventId();
    TrackingParticleRef tpRef = edm::Ref<TrackingParticleCollection>(trackingParticleH, tpIndex);
    //tpIdMap[tpIndex] = std::make_pair(tpTrackId, tpEventId.rawId());
    tpIdMap[std::make_pair(tpTrackId, tpEventId.rawId())] = tpRef ;
    tpIndex++;
  }

  // -- loop over sim clusters and get the trackId, eventId
  LogDebug("MtdSimLayerClusterToTPAssociator") << " Found " << simClusters.size() << " sim layer clusters in the event";
  
  size_t simClusIndex(0);
  for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){
    auto simClus = *simClusIt;
    MtdSimLayerClusterRef simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIndex);
    unsigned int simClusTrackId = simClus.g4Tracks()[0].trackId();
    EncodedEventId simClusEventId = simClus.eventId();

    // -- Check the trackId offset of the sim hits and keep only clusters with "direct" hits (offset == 0)    
    /*    !!!!!  need to implement offset method in MtdSimLayerClusters 
     */
    
    // -- Look for the associated tracking particle
    /*
    std::pair uniqueId = std::make_pair(simClusTrackId, simClusEventId.rawId());
    auto it = std::find_if( tpIdMap.begin(), tpIdMap.end(),
			    [&](const auto& m) {return m.second == uniqueId; });
			    
    if ( it == tpIdMap.end()) continue;
			    
    size_t tpIndex = (*it).first;
    TrackingParticleRef tpRef = edm::Ref<TrackingParticleCollection>(trackingParticleH, tpIndex);
    */

    std::pair uniqueId = std::make_pair(simClusTrackId, simClusEventId.rawId());
    auto it = tpIdMap.find(uniqueId);

    if ( it != tpIdMap.end() ) {
    
	TrackingParticleRef tpRef = tpIdMap[uniqueId];
    
	LogDebug("MtdSimLayerClusterToTPAssociator::associateSimToTP") << "MtdSimLayerCluster: index = " << simClusIndex << "   simClusTrackId = " << simClusTrackId << " simClusEventId = " << simClusEventId.rawId()
								       << " simClusEta = "<< simClus.eta() << " simClusPhi = " << simClus.phi() << "  simClusTime = " << simClus.simLCTime() <<  "  simClusEnergy = " << simClus.simLCEnergy() << std::endl;                                            
	LogDebug("MtdSimLayerClusterToTPAssociator::associateSimToTP") << "  --> Found associated tracking particle:  index = " << tpIndex << "    tpTrackId = " << (*tpRef).g4Tracks()[0].trackId() << " tpEventId = " << (*tpRef).eventId().rawId() << std::endl;

	outputCollection.insert(simClusterRef, tpRef);  

    }
   
    
    /* size_t tpIndex(0);
    for (auto tpIt = trackingParticles.begin(); tpIt != trackingParticles.end(); tpIt++) {
      auto tp = *tpIt;
      unsigned int tpTrackId = tp.g4Tracks()[0].trackId();
      EncodedEventId tpEventId = tp.eventId();

      // -- fill the output collection
      if ( tpTrackId == simClusTrackId && tpEventId.rawId() == simClusEventId.rawId() ){
	LogDebug("MtdSimLayerClusterToTPAssociator") << "MtdSimLayerCluster: index = " << simClusIndex << "   simClusTrackId = " << simClusTrackId << " simClusEventId = " << simClusEventId.rawId()
						     <<  " simClusEta = "<< simClus.eta() << " simClusPhi = " << simClus.phi() << "  simClusTime = " << simClus.simLCTime() << std::endl;
	LogDebug("MtdSimLayerClusterToTPAssociator") << "  --> Found associated tracking particle:  index = " << tpIndex << "    tpTrackId = " << tpTrackId << " tpEventId = " << tpEventId.rawId() << std::endl;

	TrackingParticleRef tpRef = edm::Ref<TrackingParticleCollection>(trackingParticleH, tpIndex);
    	outputCollection.insert(simClusterRef, tpRef);
	break;
      }

      tpIndex++;
      
    }// -- end loop over tracking particles
    */
    
    simClusIndex++;

  }// -- end loop over sim clus
      
  return outputCollection;
  
}



reco::TPToSimCollectionMtd MtdSimLayerClusterToTPAssociatorByTrackIdImpl::associateTPToSim(
      const edm::Handle<MtdSimLayerClusterCollection>& simClusH,
      const edm::Handle<TrackingParticleCollection>& trackingParticleH) const {

  TPToSimCollectionMtd outputCollection(productGetter_);

  // -- get the collections
  const auto& simClusters  = *simClusH.product();
  const auto& trackingParticles = *trackingParticleH.product();

  // -- Loop over MtdSimLayerClusters and build a temporary map of trackId, eventId --> simClusterRef
  //std::map<int, std::pair<unsigned int, uint32_t>> simClusIdMap;
  std::map<std::pair<unsigned int, uint32_t>, MtdSimLayerClusterRef> simClusIdMap;
  size_t simClusIndex(0);
  for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){
    auto simClus = *simClusIt;
    unsigned int simClusTrackId = simClus.g4Tracks()[0].trackId();
    EncodedEventId simClusEventId = simClus.eventId();
    MtdSimLayerClusterRef simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIndex);
    //simClusIdMap[simClusIndex] = std::make_pair(simClusTrackId, simClusEventId.rawId());
    simClusIdMap[std::make_pair(simClusTrackId, simClusEventId.rawId())] = simClusterRef;
    simClusIndex++;
  }
  
  
  // -- Loop over the tracking particles
  size_t tpIndex(0);
  for (auto tpIt = trackingParticles.begin(); tpIt != trackingParticles.end(); tpIt++) {
    auto tp = *tpIt;
    TrackingParticleRef tpRef = edm::Ref<TrackingParticleCollection>(trackingParticleH, tpIndex);
    unsigned int tpTrackId = tp.g4Tracks()[0].trackId();
    EncodedEventId tpEventId = tp.eventId();

    // -- Look for the associated mtd sim layer clus
    /*
    auto it = std::find_if( simClusIdMap.begin(), simClusIdMap.end(),
                            [&](const auto& m) {return m.second == std::make_pair(tpTrackId, tpEventId.rawId()); });

    if ( it == simClusIdMap.end()) continue;

    size_t simClusIndex = (*it).first;
    */
        
    std::pair uniqueId = std::make_pair(tpTrackId, tpEventId.rawId());
    auto it = simClusIdMap.find(uniqueId);

    if ( it != simClusIdMap.end() ) {
      MtdSimLayerClusterRef simClusterRef = simClusIdMap[uniqueId];
      
      LogDebug("MtdSimLayerClusterToTPAssociator") << "Tracking particle:  index = " << tpIndex << "  tpTrackId = " << tpTrackId << "  tpEventId = " << tpEventId.rawId();
      LogDebug("MtdSimLayerClusterToTPAssociator") << " --> Found associated MtdSimLayerCluster: index = " << simClusIndex << "   simClusTrackId = " << (*simClusterRef).g4Tracks()[0].trackId() << " simClusEventId = " << (*simClusterRef).eventId().rawId() <<  " simClusEta = "<< (*simClusterRef).eta() << " simClusPhi = " << (*simClusterRef).phi() << "  simClusTime = " << (*simClusterRef).simLCTime() <<  "  simClusEnergy = " << (*simClusterRef).simLCEnergy() << std::endl;
      
      // -- Check the trackId offset of the sim hits and keep only clusters with "direct" hits (offset == 0)    
      /*    !!!!!  need to implement offset method in MtdSimLayerClusters
       */
      
      outputCollection.insert(tpRef, simClusterRef);
      
    }
    
    // -- loop over MtdSimLayerClusters
    /*
      size_t simClusIndex(0);
      for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){
      auto simClus = *simClusIt;
      MtdSimLayerClusterRef simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIndex);
      unsigned int simClusTrackId = simClus.g4Tracks()[0].trackId();
      EncodedEventId simClusEventId = simClus.eventId();
      
      // -- Fill the output collection
      if ( tpTrackId == simClusTrackId && tpEventId.rawId() == simClusEventId.rawId() ){
        LogDebug("MtdSimLayerClusterToTPAssociator") << "Tracking particle:  index = " << tpIndex << "  tpTrackId = " << tpTrackId << "  tpEventId = " << tpEventId.rawId();
	LogDebug("MtdSimLayerClusterToTPAssociator") << "--> Found associated MtdSimLayerCluster: index = " << simClusIndex << "   simClusTrackId = " << simClusTrackId << " simClusEventId = " << simClusEventId.rawId() <<  " simClusEta = "<< simClus.eta();
	outputCollection.insert(tpRef, simClusterRef);
	break;
      }
      
      simClusIndex++;
     
    } // -- end loop over sim clusters
     */

    tpIndex++;
    
  }
  
  return outputCollection;
  
}
