import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoHRT.ak8_cff import setupCustomizedAK8
from PhysicsTools.NanoHRT.ca15_cff import setupCA15
from PhysicsTools.NanoHRT.hotvr_cff import setupHOTVR


def nanoHRT_customizeCommon(process, runOnMC):
    setupCustomizedAK8(process, runOnMC=runOnMC)
    setupCA15(process, runOnMC=runOnMC)
    setupHOTVR(process, runOnMC=runOnMC)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process


def nanoHRT_customizeData(process):
    process = nanoHRT_customizeCommon(process, False)
    return process


def nanoHRT_customizeMC(process):
    process = nanoHRT_customizeCommon(process, True)
    return process
