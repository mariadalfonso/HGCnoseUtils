#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
typedef SimpleFlatTableProducer<CaloParticle> SimpleCaloParticleFlatTableProducer;

#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
typedef SimpleFlatTableProducer<CaloRecHit> SimpleCaloRecHitFlatTableProducer;

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
typedef SimpleFlatTableProducer<reco::CaloCluster> SimpleCaloClusterFlatTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(SimpleCaloParticleFlatTableProducer);
DEFINE_FWK_MODULE(SimpleCaloRecHitFlatTableProducer);
DEFINE_FWK_MODULE(SimpleCaloClusterFlatTableProducer);
