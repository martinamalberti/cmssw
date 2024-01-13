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
  std::map< std::pair<unsigned int, uint32_t>, TrackingParticleRef> tpIdMap;
  size_t tpIndex(0);
  for (auto tpIt = trackingParticles.begin(); tpIt != trackingParticles.end(); tpIt++) {
    auto tp = *tpIt;
    unsigned int tpTrackId = tp.g4Tracks()[0].trackId();
    EncodedEventId tpEventId = tp.eventId();
    TrackingParticleRef tpRef = edm::Ref<TrackingParticleCollection>(trackingParticleH, tpIndex);
    tpIdMap[std::make_pair(tpTrackId, tpEventId.rawId())] = tpRef ;
    tpIndex++;
  }

  // -- loop over sim clusters and get the trackId, eventId

  LogDebug("MtdSimLayerClusterToTPAssociator") << " Found " << simClusters.size() << " MtdSimLayerClusters in the event";

  size_t simClusIndex(0);
  for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){
    auto simClus = *simClusIt;
    MtdSimLayerClusterRef simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIndex);
    unsigned int simClusTrackId = simClus.g4Tracks()[0].trackId();
    EncodedEventId simClusEventId = simClus.eventId();

    // -- Check the trackId offset of the sim hits and keep only clusters with "direct" hits (offset == 0)    
    /*    !!!!!  need to implement offset method in MtdSimLayerClusters 
     */

    std::pair uniqueId = std::make_pair(simClusTrackId, simClusEventId.rawId());
    auto it = tpIdMap.find(uniqueId);
    if ( it == tpIdMap.end() ) continue;
    
    TrackingParticleRef tpRef = tpIdMap[uniqueId];
    
    LogDebug("MtdSimLayerClusterToTPAssociator::associateSimToTP") << "MtdSimLayerCluster: index = " << simClusIndex << "   simClusTrackId = " << simClusTrackId << " simClusEventId = " << simClusEventId.rawId() << " simClusEta = "<< simClus.eta() << " simClusPhi = " << simClus.phi() << "  simClusTime = " << simClus.simLCTime() <<  "  simClusEnergy = " << simClus.simLCEnergy() << std::endl;                                            
    LogDebug("MtdSimLayerClusterToTPAssociator::associateSimToTP") << "  --> Found associated tracking particle:  index = " << tpIndex << "    tpTrackId = " << (*tpRef).g4Tracks()[0].trackId() << " tpEventId = " << (*tpRef).eventId().rawId() << std::endl;
    
    outputCollection.insert(simClusterRef, tpRef);  
    
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
  std::map<std::pair<unsigned int, uint32_t>, MtdSimLayerClusterRef> simClusIdMap;
  size_t simClusIndex(0);
  for (auto simClusIt = simClusters.begin(); simClusIt != simClusters.end(); simClusIt++){
    auto simClus = *simClusIt;
    unsigned int simClusTrackId = simClus.g4Tracks()[0].trackId();
    EncodedEventId simClusEventId = simClus.eventId();
    MtdSimLayerClusterRef simClusterRef = edm::Ref<MtdSimLayerClusterCollection>(simClusH, simClusIndex);
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

    std::pair uniqueId = std::make_pair(tpTrackId, tpEventId.rawId());
    auto it = simClusIdMap.find(uniqueId);

    if ( it == simClusIdMap.end() ) continue;
    
    MtdSimLayerClusterRef simClusterRef = simClusIdMap[uniqueId];
    
    LogDebug("MtdSimLayerClusterToTPAssociator") << "Tracking particle:  index = " << tpIndex << "  tpTrackId = " << tpTrackId << "  tpEventId = " << tpEventId.rawId();
    LogDebug("MtdSimLayerClusterToTPAssociator") << " --> Found associated MtdSimLayerCluster: index = " << simClusIndex << "   simClusTrackId = " << (*simClusterRef).g4Tracks()[0].trackId() << " simClusEventId = " << (*simClusterRef).eventId().rawId() <<  " simClusEta = "<< (*simClusterRef).eta() << " simClusPhi = " << (*simClusterRef).phi() << "  simClusTime = " << (*simClusterRef).simLCTime() <<  "  simClusEnergy = " << (*simClusterRef).simLCEnergy() << std::endl;
    
    // -- Check the trackId offset of the sim hits and keep only clusters with "direct" hits (offset == 0)    
    /*    !!!!!  need to implement offset method in MtdSimLayerClusters
     */
    
    outputCollection.insert(tpRef, simClusterRef);

    tpIndex++;
  }
  
  return outputCollection;
  
}
