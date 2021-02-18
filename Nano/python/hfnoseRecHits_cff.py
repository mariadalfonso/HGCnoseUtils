import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var,P3Vars

hfnoseRecHitsTable = cms.EDProducer("SimpleCaloRecHitFlatTableProducer",
    src = cms.InputTag("HGCalRecHit:HGCHFNoseRecHits"),
    cut = cms.string(""), 
    name = cms.string("RecHitHFNose"),
    doc  = cms.string("RecHits in HFNose"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        detId = Var('detid().rawId()', 'int', precision=-1, doc='detId'),
        energy = Var('energy', 'float', precision=14, doc='energy'),
        time = Var('time', 'float', precision=14, doc='hit time'),
    )
)

hfnoseRecHitsPositionTable = cms.EDProducer("HGCRecHitPositionFromDetIDTableProducer",
    src = hfnoseRecHitsTable.src,
    cut = hfnoseRecHitsTable.cut, 
    name = hfnoseRecHitsTable.name,
    doc  = hfnoseRecHitsTable.doc,
)

hfnoseRecHitsSequence = cms.Sequence(
    hfnoseRecHitsTable
    +hfnoseRecHitsPositionTable
)

