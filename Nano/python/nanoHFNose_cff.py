from __future__ import print_function
import FWCore.ParameterSet.Config as cms

###
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import genParticleTable
from PhysicsTools.NanoAOD.genVertex_cff import *
###
from hfnoseRecHits_cff import *
from caloParticles_cff import *
from layerClusters_cff import *

nanoMetadata = cms.EDProducer("UniqueStringProducer",
    strings = cms.PSet(
        tag = cms.string("untagged"),
    )
)

genParticleTable.src = "genParticles"
genParticleTable.variables = cms.PSet(genParticleTable.variables,
    charge = CandVars.charge)

nanoHFNoseSequence = cms.Sequence(nanoMetadata+genVertexTables+genParticleTable
                                  +caloParticleTable
                                  +hfnoseRecHitsSequence
                                  +hfnoseLayerClusterSequence
)
