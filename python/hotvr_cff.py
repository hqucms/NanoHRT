import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

# ---------------------------------------------------------


def setupHOTVR(process, runOnMC=False, path=None):

    # HOTVR
    process.hotvrPuppi = cms.EDProducer('HOTVRProducer',
        src=cms.InputTag("puppi")
    )

    from  PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import patJetCorrFactors
    process.jetCorrFactorsHOTVRSubjets = patJetCorrFactors.clone(
        src=cms.InputTag('hotvrPuppi', 'RecoSubJets'),
        levels=cms.vstring('L2Relative', 'L3Absolute', 'L2L3Residual'),
        payload=cms.string('AK4PFPuppi'),
        primaryVertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
    )
    from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import _patJets
    process.updatedHOTVRSubjets = _patJets.clone(
        jetSource=cms.InputTag('hotvrPuppi', 'RecoSubJets'),
        jetCorrFactorsSource=cms.VInputTag(cms.InputTag("jetCorrFactorsHOTVRSubjets")),
        addBTagInfo=cms.bool(False),
        addAssociatedTracks=cms.bool(False),
        addJetCharge=cms.bool(False),
        addGenPartonMatch=cms.bool(False),
        embedGenPartonMatch=cms.bool(False),
        addGenJetMatch=cms.bool(False),
        embedGenJetMatch=cms.bool(False),
        getJetMCFlavour=cms.bool(False),
        addJetFlavourInfo=cms.bool(False),
    )

    process.finalHOTVR = cms.EDProducer('HOTVRUpdater',
        src=cms.InputTag("hotvrPuppi"),
        subjets=cms.InputTag("updatedHOTVRSubjets"),
    )

    process.hotvrTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src=cms.InputTag("finalHOTVR"),
        name=cms.string("HOTVRPuppi"),
        cut=cms.string(""),
        doc=cms.string("HOTVR Puppi jets"),
        singleton=cms.bool(False),  # the number of entries is variable
        extension=cms.bool(False),  # this is the main table for the jets
        variables=cms.PSet(P4Vars,
            tau1=Var("userFloat('tau1')", float, doc="Nsubjettiness (1 axis)", precision=10),
            tau2=Var("userFloat('tau2')", float, doc="Nsubjettiness (2 axis)", precision=10),
            tau3=Var("userFloat('tau3')", float, doc="Nsubjettiness (3 axis)", precision=10),
            fpt=Var("userFloat('fpt')", float, doc="pT(sj1)/pT(jet)", precision=10),
            mmin=Var("userFloat('mmin')", float, doc="min mass of subjet pairs", precision=10),
            nsubjets=Var("?nSubjetCollections()>0?subjets().size():0", int, doc="number of subjets"),
            subJetIdx1=Var("?nSubjetCollections()>0 && subjets().size()>0?subjets()[0].key():-1", int,
                 doc="index of first subjet"),
            subJetIdx2=Var("?nSubjetCollections()>0 && subjets().size()>1?subjets()[1].key():-1", int,
                 doc="index of second subjet"),
            subJetIdx3=Var("?nSubjetCollections()>0 && subjets().size()>2?subjets()[2].key():-1", int,
                 doc="index of third subjet"),
        )
    )
    process.hotvrTable.variables.pt.precision = 10

    process.hotvrSubJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src=cms.InputTag("updatedHOTVRSubjets"),
        cut=cms.string(""),
        name=cms.string("HOTVRPuppiSubJet"),
        doc=cms.string("HOTVR Puppi subjets"),
        singleton=cms.bool(False),  # the number of entries is variable
        extension=cms.bool(False),  # this is the main table for the jets
        variables=cms.PSet(P4Vars,
            rawFactor=Var("1.-jecFactor('Uncorrected')", float, doc="1 - Factor to get back to raw pT", precision=6),
            area=Var("jetArea()", float, doc="jet catchment area, for JECs", precision=10),
        )
    )
    process.hotvrSubJetTable.variables.pt.precision = 10

    process.hotvrTask = cms.Task(
        process.hotvrPuppi,
        process.jetCorrFactorsHOTVRSubjets,
        process.updatedHOTVRSubjets,
        process.finalHOTVR,
        process.hotvrTable,
        process.hotvrSubJetTable
        )

    if path is None:
        process.schedule.associate(process.hotvrTask)
    else:
        getattr(process, path).associate(process.hotvrTask)

# ---------------------------------------------------------
