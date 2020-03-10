import FWCore.ParameterSet.Config as cms
from PhysicsTools.Pancakes.ak8_cff import setupCustomizedAK8, addCustomizedAK8PF
from PhysicsTools.NanoAOD.common_cff import Var


def nanoHRT_customizeSV(process):
    process.vertexTable.dlenMin = -1
    process.vertexTable.dlenSigMin = -1
    process.svCandidateTable.variables.ntracks = Var("numberOfDaughters()", int, doc="number of tracks")


def nanoHRT_customizeCommon(process, runOnMC):
    nanoHRT_customizeSV(process)
    setupCustomizedAK8(process, runOnMC=runOnMC)
    addCustomizedAK8PF(process)

    # fix genParticles: keep first gen decay product for all top/W/Z/H
    process.finalGenParticles.select.append('keep+ (abs(pdgId) == 6 || abs(pdgId) == 23 || abs(pdgId) == 24 || abs(pdgId) == 25)')
    # update MET w/ JEC
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    runMetCorAndUncFromMiniAOD(process, isData=not runOnMC)
    # remove regular fat jets
    process.jetTables.remove(process.fatJetTable)
    process.jetTables.remove(process.subJetTable)
    # don't drop bacon loose taus, but keep the other (era-dependent?) selection as well
    baconcut = "pt > 18 && tauID('decayModeFinding') && (tauID('byCombinedIsolationDeltaBetaCorrRaw3Hits') < 5)"
    process.finalTaus.cut = cms.string("(%s) || (%s)" % (baconcut, process.finalTaus.cut.value()))
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

def nanoHRT_customizeMINLOnnlops(process):
    wnames = [
        "nnlops-11-1",
        "nnlops-11-2",
        "nnlops-11-3",
        "nnlops-11-4",
        "nnlops-11-5",
        "nnlops-11-6",
        "nnlops-11-7",
        "nnlops-11-8",
        "nnlops-11-9",
        "nnlops-22-1",
        "nnlops-22-2",
        "nnlops-22-3",
        "nnlops-22-4",
        "nnlops-22-5",
        "nnlops-22-6",
        "nnlops-22-7",
        "nnlops-22-8",
        "nnlops-22-9",
        "nnlops-0505-1",
        "nnlops-0505-2",
        "nnlops-0505-3",
        "nnlops-0505-4",
        "nnlops-0505-5",
        "nnlops-0505-6",
        "nnlops-0505-7",
        "nnlops-0505-8",
        "nnlops-0505-9",
    ]
    process.genWeightsTable.namedWeightIDs = cms.vstring(wnames)
    process.genWeightsTable.namedWeightLabels = cms.vstring(wnames)
    return process
