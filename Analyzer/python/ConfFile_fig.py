import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C10_cff import Phase2C10

process = cms.Process('DEMO',Phase2C10)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet (
    wantSummary = cms.untracked.bool(False),
    numberOfThreads = cms.untracked.uint32(8),
    numberOfStreams = cms.untracked.uint32(8)
)

process.load('Configuration.Geometry.GeometryExtended2026D47Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D47_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100


#if hasattr(process,'MessageLogger'):
#process.MessageLogger.categories.append('HGCalGeom')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     'file:/eos/user/d/dalfonso/HGCnose/11_1_pre7/step3_D47_PionPt20.root'
    )
)

process.comparisonPlots = cms.EDAnalyzer("GenAnalyzer")
#process.comparisonPlots = cms.EDAnalyzer("GenAnalyzer",
#                                         GenParticleTag = cms.InputTag("genParticles"),
#                                         GenJetTag = cms.InputTag("genJets")
#)


process.TFileService = cms.Service('TFileService', fileName = cms.string('plots_recHits_D47_PionPt20.root'))

process.p = cms.Path(process.comparisonPlots)
