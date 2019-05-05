import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000


process.load("Geometry.ForwardCommonData.hfnoseXML_cfi")
process.load("Geometry.ForwardCommonData.hfnoseParametersInitialization_cfi")
process.load("Geometry.ForwardCommonData.hfnoseNumberingInitialization_cfi")
process.load("Geometry.CaloEventSetup.HFNoseTopology_cfi")
process.load("Geometry.ForwardGeometry.HFNoseGeometryESProducer_cfi")

#if hasattr(process,'MessageLogger'):
#process.MessageLogger.categories.append('HGCalGeom')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/eos/user/d/dalfonso/HGCnose/SinglePtNoSmearing/K0L_Pt3_Eta35_Step3.root'
#        'file:/eos/user/d/dalfonso/HGCnose/SinglePtNoSmearing/Gamma_Pt3_Eta35_Step3.root'

#        'file:/eos/user/d/dalfonso/HGCnose/SinglePt/K0L_Pt5_Eta35_Step3.root'
#        'file:/eos/user/d/dalfonso/HGCnose/SinglePt/K0L_Pt10_Eta35_Step3.root'
#        'file:/eos/user/d/dalfonso/HGCnose/SinglePt/K0L_Pt20_Eta35_Step3.root'

#        'file:/eos/user/d/dalfonso/VBS/WWjj_SS_ewk/NOSE/WWjj_SS_ewk_RECO_HGCNose_1to10.root',
#        'file:/eos/user/d/dalfonso/VBS/WWjj_SS_ewk/NOSE/WWjj_SS_ewk_RECO_HGCNose_0.root'
#        'file:/eos/user/d/dalfonso/VBS/WWjj_SS_ewk/HF/WWjj_SS_ewk_RECO_phase2.root',
#        'file:/eos/user/d/dalfonso/VBS/WWjj_SS_ewk/HF/WWjj_SS_ewk_RECO_phase2_0.root'

#    'file:/eos/user/d/dalfonso/HGCnose/SingleGammaE/NOSE/Muon_E10_Step3.root'
#    'file:/eos/user/d/dalfonso/HGCnose/SingleGammaE/NOSE/gamma_E10_step3.root'
#    'file:/eos/user/d/dalfonso/HGCnose/SingleGammaE/NOSE/Gamma_E100_Step3.root'
#    'file:/eos/user/d/dalfonso/HGCnose/SingleGammaE/HF/Gamma_E10_Step3_phase2.root'
#    'file:/eos/user/d/dalfonso/HGCnose/SingleGammaE/HF/Gamma_E60_Step3_phase2.root'
#    'file:/eos/user/d/dalfonso/HGCnose/SingleGammaE/gamma_E10_step3_ev1000.root'
#    'file:/eos/user/d/dalfonso/HGCnose/SinglePionE/NOSE/Pion_E10_Step3.root'

    )
)

process.comparisonPlots = cms.EDAnalyzer("GenAnalyzer")
#process.comparisonPlots = cms.EDAnalyzer("GenAnalyzer",
#                                         GenParticleTag = cms.InputTag("genParticles"),
#                                         GenJetTag = cms.InputTag("genJets")
#)

#process.TFileService = cms.Service('TFileService', fileName = cms.string('test_Muon.root') )
#process.TFileService = cms.Service('TFileService', fileName = cms.string('SCAN/test_Pion_Nose_E10_W12.root') )   
#process.TFileService = cms.Service('TFileService', fileName = cms.string('RESOLUTION/test_Gamma_HF_E60.root') )
#process.TFileService = cms.Service('TFileService', fileName = cms.string('RESOLUTION/test_Pion_HF_E60.root') )
#process.TFileService = cms.Service('TFileService', fileName = cms.string('RESOLUTION/test_Pion_Nose_E10.root') )
#process.TFileService = cms.Service('TFileService', fileName = cms.string('TIME2/test_gamma_Pt3.root') )  
process.TFileService = cms.Service('TFileService', fileName = cms.string('TIME2/test_K0L_Pt3.root') )  

process.p = cms.Path(process.comparisonPlots)
