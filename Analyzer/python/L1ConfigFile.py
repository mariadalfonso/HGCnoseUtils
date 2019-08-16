import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C6_cff import Phase2C6

process = cms.Process("Demo",Phase2C6)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet (
    wantSummary = cms.untracked.bool(False),
    numberOfThreads = cms.untracked.uint32(8),
    numberOfStreams = cms.untracked.uint32(8)
)

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.load("Geometry.ForwardCommonData.hfnoseXML_cfi")
process.load("Geometry.ForwardCommonData.hfnoseParametersInitialization_cfi")
process.load("Geometry.ForwardCommonData.hfnoseNumberingInitialization_cfi")
process.load("Geometry.CaloEventSetup.HFNoseTopology_cfi")
process.load("Geometry.ForwardGeometry.HFNoseGeometryESProducer_cfi")

#if hasattr(process,'MessageLogger'):
#process.MessageLogger.categories.append('HGCalGeom')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

        'file:/tmp/dalfonso/test.root'

    )
)

process.comparisonPlots = cms.EDAnalyzer("L1Analyzer")
process.TFileService = cms.Service('TFileService', fileName = cms.string('test_Photon100_L1.root'))

process.p = cms.Path(process.comparisonPlots)
