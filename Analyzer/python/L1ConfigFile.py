import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C10_cff import Phase2C10

process = cms.Process("Demo",Phase2C10)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet (
    wantSummary = cms.untracked.bool(False),
    numberOfThreads = cms.untracked.uint32(8),
    numberOfStreams = cms.untracked.uint32(8)
)

process.load('Configuration.Geometry.GeometryExtended2026D60Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D60_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#if hasattr(process,'MessageLogger'):
#process.MessageLogger.categories.append('HGCalGeom')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

#        'file:$CMSSW_BASE/src/test_L1.root'
#        'file:/eos/user/d/dalfonso/HGCnose/11_2_calTow/gamma_pt20_D60_Scale058_L1_TOWER.root'
        'file:/eos/user/d/dalfonso/HGCnose/11_2_calTow/pion_pt20_D60_Scale058_L1_TOWER.root'

    )
)

process.comparisonPlots = cms.EDAnalyzer("L1Analyzer",
#                                         pdgId = cms.untracked.int32(22)
                                         pdgId = cms.untracked.int32(211)
                                     )

#process.TFileService = cms.Service('TFileService', fileName = cms.string('test_Photon20_L1.root'))
process.TFileService = cms.Service('TFileService', fileName = cms.string('test_Pion20_L1.root'))

process.p = cms.Path(process.comparisonPlots)
