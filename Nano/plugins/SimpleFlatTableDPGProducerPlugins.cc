#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
typedef SimpleFlatTableProducer<CaloParticle> SimpleCaloParticleFlatTableProducer;

#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
typedef SimpleFlatTableProducer<CaloRecHit> SimpleCaloRecHitFlatTableProducer;

#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
typedef SimpleFlatTableProducer<HFRecHit> SimpleHFRecHitFlatTableProducer;

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
typedef SimpleFlatTableProducer<reco::CaloCluster> SimpleCaloClusterFlatTableProducer;

#include "DataFormats/HGCalReco/interface/Trackster.h"
typedef SimpleFlatTableProducer<ticl::Trackster> SimpleTracksterFlatTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(SimpleCaloParticleFlatTableProducer);
DEFINE_FWK_MODULE(SimpleHFRecHitFlatTableProducer);
DEFINE_FWK_MODULE(SimpleCaloRecHitFlatTableProducer);
DEFINE_FWK_MODULE(SimpleCaloClusterFlatTableProducer);
DEFINE_FWK_MODULE(SimpleTracksterFlatTableProducer);
