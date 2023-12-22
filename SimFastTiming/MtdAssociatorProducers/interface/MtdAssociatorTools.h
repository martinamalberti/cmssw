#ifndef SimFastTiming_MtdAssociatorProducers_interface_MtdAssociatorTools_h
#define SimFastTiming_MtdAssociatorProducers_interface_MtdAssociatorTools_h

#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "Geometry/MTDCommonData/interface/MTDTopologyMode.h"
#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeomUtil.h"
#include "Geometry/MTDNumberingBuilder/interface/MTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

#include <map>

void fillSimHitTkIdMap(edm::Handle<CrossingFrame<PSimHit>> btlSimHitsH_,
		       edm::Handle<CrossingFrame<PSimHit>> etlSimHitsH_,
		       mtd::MTDGeomUtil geomTools_,
		       std::unordered_map<uint64_t, std::vector<int>>& tkIdMap,
		       std::unordered_map<uint64_t, std::vector<int>>& offsetMap) {

  using namespace angle_units::operators;

  std::array<edm::Handle<CrossingFrame<PSimHit>>, 2> inputSimHitsH{{btlSimHitsH_, etlSimHitsH_}};

  for (auto const& simHitsH : inputSimHitsH) {
    MixCollection<PSimHit> simHits(simHitsH.product());

    // -- loop over sim hits
    for (const auto& simHit : simHits) {
      DetId id = simHit.detUnitId();

      // --- Use only hits compatible with the in-time bunch-crossing
      if (simHit.tof() < 0 || simHit.tof() > 25.) continue;

      // -- Get an unique id: for BTL the detId is unique (one for each crystal), for ETL the detId is not enough
      //    also row and column are needed. An unique number is created from detId, row, col
      // -  Get row and column
      const auto &position = simHit.localPosition();
      LocalPoint simscaled(convertMmToCm(position.x()), convertMmToCm(position.y()), convertMmToCm(position.z()));
      std::pair<uint8_t, uint8_t> pixel = geomTools_.pixelInModule(id, simscaled);
      // -  create the unique id
      uint64_t uniqueId = static_cast<uint64_t>(id.rawId()) << 32;
      uniqueId |= pixel.first << 16;
      uniqueId |= pixel.second;
      //std::cout << id.rawId() << "   uniqueID = " << uniqueId <<std::endl;
      tkIdMap[uniqueId].push_back(simHit.trackId());
      offsetMap[uniqueId].push_back(simHit.offsetTrackId());
    }
  }
}

#endif
