# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: TTbar_14TeV_TuneCUETP8M1_cfi --conditions auto:phase2_realistic -n 10 --era Phase2C6 --eventcontent FEVTDEBUG --relval 9000,100 -s GEN,SIM --datatier GEN-SIM --beamspot HLLHC14TeV --geometry Extended2023D31 --fileout file:step1.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('SIM',eras.Phase2C6)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D31Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D31_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
##process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.load('IOMC.EventVertexGenerators.VtxSmearedFlat_cfi')

process.VtxSmeared.MaxZ = cms.double(0.)
process.VtxSmeared.MinZ = cms.double(0.)
process.VtxSmeared.MaxT = cms.double(0.)
process.VtxSmeared.MinT = cms.double(0.)


process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(100)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 10

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(8),
    numberOfStreams = cms.untracked.uint32(8),
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('TTbar_14TeV_TuneCUETP8M1_cfi nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:/tmp/dalfonso/QCD_GEN_SIM.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')



process.generator = cms.EDFilter("Pythia8GeneratorFilter",
    PythiaParameters = cms.PSet(
        parameterSets = cms.vstring(
            'pythia8CommonSettings', 
            'pythia8CUEP8M1Settings', 
            'processParameters'
        ),
        processParameters = cms.vstring(
            'HardQCD:all = off',
            'HardQCD:gg2gg = on',
            'HardQCD:qq2qq = on',
#            'HardQCD:all = on',
            'PhaseSpace:pTHatMin = 15.',
#            'PhaseSpace:pTHatMax = .'
#                'MSEL=0          !  for user specification of sub-processes',
#                'MSUB(11)=1      !  qq->qq   ON, if one needs quark jets',
#                'MSUB(68)=1      !  gg->gg   ON, if one needs gluon jets',
#                'CKIN(3)=100.    !  Pt low cut but also the Et jet required',
#                'CKIN(13)=3.49    !  3.49',
#                'CKIN(14)=3.51    !  3.51',
#                'CKIN(15)=-3.51   ! -3.51',
#HardQCD                'CKIN(16)=-3.49   ! -3.49'
#            'HardQCD:all = on', 
#            'PhaseSpace:pTHatMin = 30.'
        ),
        pythia8CUEP8M1Settings = cms.vstring(
            'Tune:pp 14', 
            'Tune:ee 7', 
            'MultipartonInteractions:pT0Ref=2.4024', 
            'MultipartonInteractions:ecmPow=0.25208', 
            'MultipartonInteractions:expPow=1.6'
        ),
        pythia8CommonSettings = cms.vstring(
            'Tune:preferLHAPDF = 2', 
            'Main:timesAllowErrors = 10000', 
            'Check:epTolErr = 0.01', 
            'Beams:setProductionScalesFromLHEF = off', 
            'SLHA:keepSM = on', 
            'SLHA:minMassSM = 1000.', 
            'ParticleDecays:limitTau0 = on', 
            'ParticleDecays:tau0Max = 10', 
            'ParticleDecays:allowPhotonRadiation = on'
        )
    ),
    comEnergy = cms.double(13000.0),
    crossSection = cms.untracked.double(57200000.0),
    filterEfficiency = cms.untracked.double(1.0),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    pythiaPylistVerbosity = cms.untracked.int32(1)
)


#process.generator = cms.EDProducer("FlatRandomEGunProducer",
#    PGunParameters = cms.PSet(
#            PartID = cms.vint32(211),
#        MinPhi = cms.double(-3.14159265359),
#        MaxPhi = cms.double(3.14159265359), ## in radians
#        MinEta = cms.double(3.499),
#        MaxEta = cms.double(3.501),
#        MinE = cms.double(99.99), # in GeV
#        MaxE = cms.double(100.01)
##        MinE = cms.double(299.99), # in GeV
##        MaxE = cms.double(300.01)
#    ),
#    Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts
#    AddAntiParticle = cms.bool(True),
#)

# Path and EndPath definitions                                                                                                                                                  

process.load('GeneratorInterface.GenFilters.VBFGenJetFilter_cfi')

process.vbfGenJetFilterA.oppositeHemisphere = cms.untracked.bool(True)
process.vbfGenJetFilterA.minPt = cms.untracked.double(10)
process.vbfGenJetFilterA.minInvMass = cms.untracked.double(100)
process.vbfGenJetFilterA.minDeltaPhi = cms.untracked.double(  -9999.), # Minimum dijet delta phi
process.vbfGenJetFilterA.maxDeltaPhi = cms.untracked.double(  9999.), # Maximum dijet delta phi

#vbfGenJetFilterA = cms.EDFilter("VBFGenJetFilter",
#  inputTag_GenJetCollection = cms.untracked.InputTag('ak4GenJetsNoNu'),
#  oppositeHemisphere = cms.untracked.bool  ( False), # Require j1_eta*j2_eta<0
#  minPt              = cms.untracked.double(    40), # Minimum dijet jet_pt
#  minEta             = cms.untracked.double(  -4.8), # Minimum dijet jet_eta
#  maxEta             = cms.untracked.double(   4.8), # Maximum dijet jet_eta
#  minInvMass         = cms.untracked.double( 1000.), # Minimum dijet invariant mass
#  maxInvMass         = cms.untracked.double(99999.), # Maximum dijet invariant mass
#  minDeltaPhi        = cms.untracked.double(  -1.0), # Minimum dijet delta phi
#  maxDeltaPhi        = cms.untracked.double(  2.15), # Maximum dijet delta phi
#  minDeltaEta        = cms.untracked.double(   3.0), # Minimum dijet delta eta
#  maxDeltaEta        = cms.untracked.double(99999.)  # Maximum dijet delta eta
#)
                                                         
process.generation_step = cms.Path(process.pgen*process.vbfGenJetFilterA)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)


process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.FEVTDEBUGoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path).insert(0, process.ProductionFilterSequence)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
