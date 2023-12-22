// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "Geometry/MTDCommonData/interface/MTDTopologyMode.h"
#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeomUtil.h"
#include "Geometry/MTDNumberingBuilder/interface/MTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

#include "MtdSimLayerClusterToTPAssociatorByTrackIdImpl.h"


//
// Class declaration
//

class MtdSimLayerClusterToTPAssociatorByTrackIdProducer : public edm::global::EDProducer<> {
public:
  explicit MtdSimLayerClusterToTPAssociatorByTrackIdProducer (const edm::ParameterSet &);
  ~MtdSimLayerClusterToTPAssociatorByTrackIdProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  
private:
  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;
  
  edm::EDGetTokenT<CrossingFrame<PSimHit> > btlSimHitsToken_;
  edm::EDGetTokenT<CrossingFrame<PSimHit> > etlSimHitsToken_;

  edm::ESGetToken<MTDGeometry, MTDDigiGeometryRecord> mtdgeomToken_;
  edm::ESGetToken<MTDTopology, MTDTopologyRcd> mtdtopoToken_;

};


MtdSimLayerClusterToTPAssociatorByTrackIdProducer::MtdSimLayerClusterToTPAssociatorByTrackIdProducer(const edm::ParameterSet &pset)
{
  btlSimHitsToken_ = consumes<CrossingFrame<PSimHit> >(pset.getParameter<edm::InputTag>("btlSimHitsTag"));
  etlSimHitsToken_ = consumes<CrossingFrame<PSimHit> >(pset.getParameter<edm::InputTag>("etlSimHitsTag"));
  mtdgeomToken_ = esConsumes<MTDGeometry, MTDDigiGeometryRecord>();
  mtdtopoToken_ = esConsumes<MTDTopology, MTDTopologyRcd>();

  // Register the product
  produces<reco::MtdSimLayerClusterToTPAssociator>();
    
}


MtdSimLayerClusterToTPAssociatorByTrackIdProducer::~MtdSimLayerClusterToTPAssociatorByTrackIdProducer() {}


void MtdSimLayerClusterToTPAssociatorByTrackIdProducer::produce(edm::StreamID,
								 edm::Event &iEvent,
								 const edm::EventSetup &es) const {
  
  using namespace angle_units::operators;

  auto geometryHandle = es.getTransientHandle(mtdgeomToken_);
  const MTDGeometry* geom = geometryHandle.product();
  
  auto topologyHandle = es.getTransientHandle(mtdtopoToken_);
  const MTDTopology* topology = topologyHandle.product();
  
  mtd::MTDGeomUtil geomTools_;
  geomTools_.setGeometry(geom);
  geomTools_.setTopology(topology);
  
  edm::Handle<CrossingFrame<PSimHit>> btlSimHitsH;
  iEvent.getByToken(btlSimHitsToken_, btlSimHitsH);
  
  edm::Handle<CrossingFrame<PSimHit>> etlSimHitsH;
  iEvent.getByToken(etlSimHitsToken_, etlSimHitsH);

  auto impl = std::make_unique<MtdSimLayerClusterToTPAssociatorByTrackIdImpl>(iEvent.productGetter(), btlSimHitsH, etlSimHitsH, geomTools_);
  auto toPut = std::make_unique<reco::MtdSimLayerClusterToTPAssociator>(std::move(impl));
  iEvent.put(std::move(toPut));
  
}


void MtdSimLayerClusterToTPAssociatorByTrackIdProducer::fillDescriptions(edm::ConfigurationDescriptions &cfg) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("btlSimHitsTag", edm::InputTag("mix:g4SimHitsFastTimerHitsBarrel"));
  desc.add<edm::InputTag>("etlSimHitsTag", edm::InputTag("mix:g4SimHitsFastTimerHitsEndcap"));

  cfg.add("mtdSimLayerClusterToTPAssociatorByTrackId", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MtdSimLayerClusterToTPAssociatorByTrackIdProducer);
