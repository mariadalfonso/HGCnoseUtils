import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var,P3Vars


HFNoseTracksterEMTable = cms.EDProducer("SimpleTracksterFlatTableProducer",
    src = cms.InputTag("ticlTrackstersHFNoseEM"),
    cut = cms.string(""), 
    name = cms.string("TICLTrackstersHFNoseEM"),
    doc  = cms.string("TICLTrackstersEM in HFNose"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        raw_energy = Var('raw_energy', 'float', precision=14, doc='raw energy'),
        raw_em_energy = Var('raw_em_energy', 'float', precision=14, doc='raw energy'),
        time = Var('time', 'float', precision=14, doc='time'),
        time_error = Var('timeError', 'float', precision=14, doc='time error'),
    )
)


HFNoseTracksterHADTable = cms.EDProducer("SimpleTracksterFlatTableProducer",
    src = cms.InputTag("ticlTrackstersHFNoseHAD"),
    cut = cms.string(""), 
    name = cms.string("TICLTrackstersHFNoseHAD"),
    doc  = cms.string("TICLTrackstersHAD in HFNose"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        raw_energy = Var('raw_energy', 'float', precision=14, doc='raw energy'),
        raw_em_energy = Var('raw_em_energy', 'float', precision=14, doc='raw energy'),
        time = Var('time', 'float', precision=14, doc='time'),
        time_error = Var('timeError', 'float', precision=14, doc='time error'),
    )
)



HFNoseTracksterTrkTable = cms.EDProducer("SimpleTracksterFlatTableProducer",
    src = cms.InputTag("ticlTrackstersHFNoseTrk"),
    cut = cms.string(""), 
    name = cms.string("TICLTrackstersHFNoseTrk"),
    doc  = cms.string("TICLTrackstersTrk in HFNose"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        raw_energy = Var('raw_energy', 'float', precision=14, doc='raw energy'),
        raw_em_energy = Var('raw_em_energy', 'float', precision=14, doc='raw energy'),
        time = Var('time', 'float', precision=14, doc='time'),
        time_error = Var('timeError', 'float', precision=14, doc='time error'),
    )
)

hfnoseTrackstersTask = cms.Task( HFNoseTracksterEMTable, HFNoseTracksterHADTable, HFNoseTracksterTrkTable)
