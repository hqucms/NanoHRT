import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process('NANO', eras.Run2_2016, eras.run2_miniAOD_80XLegacy)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Services_cff')

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1000))
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
'/store/mc/RunIISummer16MiniAODv2/ZprimeToTT_M-3000_W-30_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/D6D620EF-73BE-E611-8BFB-B499BAA67780.root',
]

process.load("PhysicsTools.NanoAOD.nano_cff")

# customized for HRT
runOnMC = True

if runOnMC:
    process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v2', '')  # use v3 for new JEC
    process.nanoPath = cms.Path(process.nanoSequenceMC)
    from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeMC
    process = nanoAOD_customizeMC(process)
else:
    # for data:
    process.GlobalTag.globaltag = GlobalTag(process.GlobalTag, '94X_dataRun2_v4', '')  # use v10 for new JEC
    process.nanoPath = cms.Path(process.nanoSequence)
    from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeData
    process = nanoAOD_customizeData(process)

from PhysicsTools.NanoHRT.ak8_cff import setupCustomizedAK8
from PhysicsTools.NanoHRT.ca15_cff import setupCA15
from PhysicsTools.NanoHRT.hotvr_cff import setupHOTVR
setupCustomizedAK8(process, runOnMC=runOnMC, path='nanoPath')
setupCA15(process, runOnMC=runOnMC, path='nanoPath')
setupHOTVR(process, runOnMC=runOnMC, path='nanoPath')

from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
patTask = getPatAlgosToolsTask(process)

process.out = cms.OutputModule("NanoAODOutputModule",
    fileName = cms.untracked.string('nano.root'),
    outputCommands = process.NanoAODEDMEventContent.outputCommands,
    compressionLevel=cms.untracked.int32(9),
    compressionAlgorithm=cms.untracked.string("LZMA"),
)
process.end = cms.EndPath(process.out, patTask)
