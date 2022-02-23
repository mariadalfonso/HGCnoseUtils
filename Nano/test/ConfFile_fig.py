import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C10_cff import Phase2C10

process = cms.Process('DEMO',Phase2C10)

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#-----------------------------------------
# INPUT / OUTPUT
#-----------------------------------------

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        "file:/eos/user/d/dalfonso/HGCnose/11_3_ticl/step3_PionPt20_D60.root"
        "file:$CMSSW_BASE/src/step3_ttbar10ev.root"
    )
)

process.output = cms.OutputModule("NanoAODOutputModule",
     compressionAlgorithm = cms.untracked.string('LZMA'),
     compressionLevel = cms.untracked.int32(9),
     dataset = cms.untracked.PSet(
         dataTier = cms.untracked.string('NANOAODSIM'),
         filterName = cms.untracked.string('')
     ),
     fileName = cms.untracked.string('nanoHFNose.root'),
     outputCommands = cms.untracked.vstring(
         'drop *',
         "keep nanoaodFlatTable_*Table_*_*",     # event data
     )
 )

#-----------------------------------------
# HGCstuff
#-----------------------------------------

####
process.load("HGCnoseUtils.Nano.nanoHFNose_cff")

process.hgc = cms.Path(process.nanoHFNoseTask)

process.finalize = cms.EndPath(process.output)

process.schedule = cms.Schedule(
    process.hgc
    ,process.finalize
)
