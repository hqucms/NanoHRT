import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoHRT.ak8_cff import setupCustomizedAK8
from PhysicsTools.NanoHRT.ca15_cff import setupCA15
from PhysicsTools.NanoHRT.hotvr_cff import setupHOTVR


def nanoHRT_customizeCommon(process, runOnMC):
    setupCustomizedAK4(process, runOnMC=runOnMC)
    setupCustomizedAK8(process, runOnMC=runOnMC)
    setupCA15(process, runOnMC=runOnMC)
    setupHOTVR(process, runOnMC=runOnMC)
    # update MET w/ JEC
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    runMetCorAndUncFromMiniAOD(process, isData=not runOnMC)
    return process


def nanoHRT_customizeData(process):
    process = nanoHRT_customizeCommon(process, False)
    process.NANOAODoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process


def nanoHRT_customizeData_METMuEGClean(process):
    process = nanoHRT_customizeCommon(process, False)

    from PhysicsTools.PatUtils.tools.corMETFromMuonAndEG import corMETFromMuonAndEG
    corMETFromMuonAndEG(process,
                        pfCandCollection="",  # not needed
                        electronCollection="slimmedElectronsBeforeGSFix",
                        photonCollection="slimmedPhotonsBeforeGSFix",
                        corElectronCollection="slimmedElectrons",
                        corPhotonCollection="slimmedPhotons",
                        allMETEGCorrected=True,
                        muCorrection=False,
                        eGCorrection=True,
                        runOnMiniAOD=True,
                        postfix="MuEGClean"
                        )
    process.slimmedMETsMuEGClean = process.slimmedMETs.clone()
    process.slimmedMETsMuEGClean.src = cms.InputTag("patPFMetT1MuEGClean")
    process.slimmedMETsMuEGClean.rawVariation = cms.InputTag("patPFMetRawMuEGClean")
    process.slimmedMETsMuEGClean.t1Uncertainties = cms.InputTag("patPFMetT1%sMuEGClean")
    del process.slimmedMETsMuEGClean.caloMET
    process.metTable.src = cms.InputTag('slimmedMETsMuEGClean')

    process.NANOAODoutput.fakeNameForCrab = cms.untracked.bool(True)  # hack for crab publication
    return process


def nanoHRT_customizeMC(process):
    process = nanoHRT_customizeCommon(process, True)
    process.NANOAODSIMoutput.fakeNameForCrab = cms.untracked.bool(True)  # needed for crab publication
    return process
