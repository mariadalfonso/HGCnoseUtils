from __future__ import print_function
import FWCore.ParameterSet.Config as cms

###
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import *
from PhysicsTools.NanoAOD.genVertex_cff import *
###
from HGCnoseUtils.Nano.hfnoseRecHits_cff import *
from HGCnoseUtils.Nano.caloParticles_cff import *
from HGCnoseUtils.Nano.layerClusters_cff import *

nanoMetadata = cms.EDProducer("UniqueStringProducer",
    strings = cms.PSet(
        tag = cms.string("untagged"),
    )
)

genParticleTable.src = "genParticles"
genParticleTable.variables = cms.PSet(genParticleTable.variables,
    charge = CandVars.charge)


#from Validation.Configuration.hgcalSimValid_cff import _hfnose_hgcalAssociatorsTask
#hfnoseAssociatorsSequence = cms.Sequence(_hfnose_hgcalAssociatorsTask)

nanoHFNoseTask = cms.Task(nanoMetadata,genVertexTablesTask,genParticleTablesTask
                          ,caloParticleTable
                          #                                  +hfnoseAssociatorsSequence
                          ,hfnoseRecHitsTask
                          ,hfnoseLayerClusterTask
)
