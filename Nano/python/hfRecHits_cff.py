import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var,P3Vars

hfRecHitsTable = cms.EDProducer("SimpleHFRecHitFlatTableProducer",
    src = cms.InputTag("hfreco"),
    cut = cms.string(""), 
    name = cms.string("RecHitHF"),
    doc  = cms.string("RecHits in HF"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        detId = Var('detid().rawId()', 'int', precision=-1, doc='detId'),
        energy = Var('energy', 'float', precision=14, doc='energy'),
        time = Var('time', 'float', precision=14, doc='hit time'),
        ieta = Var('id().ieta()', 'int', precision=-1, doc='ieta'),
        iphi = Var('id().iphi()', 'int', precision=-1, doc='iphi'),
        depth = Var('id().depth()', 'int', precision=-1, doc='depth'),
    )
)

hfRecHitsTask = cms.Task(
    hfRecHitsTable
)

