#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/Common/interface/View.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include <vector>
#include <iostream>

template <typename T>
class PositionFromDetIDTableProducer : public edm::stream::EDProducer<> {
public:
  PositionFromDetIDTableProducer(edm::ParameterSet const& params)
      : name_(params.getParameter<std::string>("name")),
        doc_(params.getParameter<std::string>("doc")),
        src_(consumes<T>(params.getParameter<edm::InputTag>("src"))),
        cut_(params.getParameter<std::string>("cut"), true) {
    produces<nanoaod::FlatTable>();
  }

  ~PositionFromDetIDTableProducer() override {}

  void beginRun(const edm::Run&, const edm::EventSetup& iSetup) {
    // TODO: check that the geometry exists
    iSetup.get<CaloGeometryRecord>().get(caloGeom_);
  }

  GlobalPoint positionFromHit(const CaloRecHit& hit) {
    return positionFromDetId(hit.detid());
  }

  GlobalPoint positionFromDetId(DetId id) {

    hgcal::RecHitTools rhtools_;
    rhtools_.setGeometry(*caloGeom_);
    return rhtools_.getPosition(id);

  }
  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    edm::Handle<T> objs;
    iEvent.getByToken(src_, objs);

    std::vector<float> xvals;
    std::vector<float> yvals;
    std::vector<float> zvals;
    for (const auto& obj : *objs) {
      if (cut_(obj)) {
        auto position = positionFromHit(obj);
        xvals.emplace_back(position.x());
        yvals.emplace_back(position.y());
        zvals.emplace_back(position.z());
      }
    }

    auto tab = std::make_unique<nanoaod::FlatTable>(xvals.size(), name_, false, true);
    tab->addColumn<float>("x", xvals, "x position");
    tab->addColumn<float>("y", yvals, "y position");
    tab->addColumn<float>("z", zvals, "z position");

    iEvent.put(std::move(tab));
  }

protected:
  const std::string name_, doc_;
  const edm::EDGetTokenT<T> src_;
  const StringCutObjectSelector<typename T::value_type> cut_;
  edm::ESHandle<CaloGeometry> caloGeom_;
  //  edm::ESHandle<GlobalTrackingGeometry> trackGeom_;

};

#include "FWCore/Framework/interface/MakerMacros.h"
typedef PositionFromDetIDTableProducer<HGCRecHitCollection> HGCRecHitPositionFromDetIDTableProducer;
DEFINE_FWK_MODULE(HGCRecHitPositionFromDetIDTableProducer);
