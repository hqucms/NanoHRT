import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *
from Configuration.Eras.Modifier_run2_miniAOD_80XLegacy_cff import run2_miniAOD_80XLegacy

# ---------------------------------------------------------


def setupCustomizedAK4(process, runOnMC=False, path=None):

    # Mu subtraction
    process.MuSubProducer = cms.EDProducer('MuSubProducer',
        src=cms.InputTag('slimmedMuons'),
        pfcs=cms.InputTag('packedPFCandidates'),
        vtxs=cms.InputTag('offlineSlimmedPrimaryVertices'),
        ptmin=cms.double(55.0)
    )


    bTagDiscriminators = [
        'pfCombinedInclusiveSecondaryVertexV2BJetTags',
        'pfDeepCSVJetTags:probb',
        'pfDeepCSVJetTags:probbb'
    ]
    JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']

    from PhysicsTools.NanoHRT.jetToolbox_cff import jetToolbox
    jetToolbox(process, 'ak4', 'dummySeq', 'out', associateTask=False,
               PUMethod='CHS', JETCorrPayload='AK4PFchs', JETCorrLevels=JETCorrLevels,
	       postFix='MuSub',	
	       newPFCollection=True, nameNewPFCollection="MuSubProducer",
               Cut='pt > 20.0 && abs(rapidity()) < 2.4',
               bTagDiscriminators=bTagDiscriminators)

    srcJets = cms.InputTag('selectedPatJetsAK4PFCHSMuSub')



    # jetID
    process.looseJetIdCustomAK4 = cms.EDProducer("PatJetIDValueMapProducer",
              filterParams=cms.PSet(
                version = cms.string('WINTER16'),
                quality = cms.string('LOOSE'),
              ),
              src=srcJets
    )

    process.tightJetIdCustomAK4 = cms.EDProducer("PatJetIDValueMapProducer",
              filterParams=cms.PSet(
                version=cms.string('WINTER17'),
                quality = cms.string('TIGHT'),
              ),
              src=srcJets
    )
    run2_miniAOD_80XLegacy.toModify(process.tightJetIdCustomAK4.filterParams, version="WINTER16")

    process.tightJetIdLepVetoCustomAK4 = cms.EDProducer("PatJetIDValueMapProducer",
              filterParams=cms.PSet(
                version = cms.string('WINTER17'),
                quality = cms.string('TIGHTLEPVETO'),
              ),
              src=srcJets
    )


    process.customAK4WithUserData = cms.EDProducer("PATJetUserDataEmbedder",
        src=srcJets,
        userFloats=cms.PSet(),
        userInts=cms.PSet(
           tightId=cms.InputTag("tightJetIdCustomAK4"),
           tightIdLepVeto=cms.InputTag("tightJetIdLepVetoCustomAK4"),
        ),
    )


    run2_miniAOD_80XLegacy.toModify(process.customAK4WithUserData.userInts,
        looseId=cms.InputTag("looseJetIdCustomAK4"),
        tightIdLepVeto=None,
    )

    process.customAK4Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src=cms.InputTag("customAK4WithUserData"),
        name=cms.string("CustomAK4CHS"),
        cut=cms.string(""),
        doc=cms.string("lepton subtracted AK4"),
        singleton=cms.bool(False),  # the number of entries is variable
        extension=cms.bool(False),  # this is the main table for the jets
        variables=cms.PSet(P4Vars,
            jetId=Var("userInt('tightId')*2+4*userInt('tightIdLepVeto')", int, doc="Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto"),
            area=Var("jetArea()", float, doc="jet catchment area, for JECs", precision=10),
            rawFactor=Var("1.-jecFactor('Uncorrected')", float, doc="1 - Factor to get back to raw pT", precision=6),
            btagCSVV2=Var("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')", float, doc=" pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)", precision=10),
            btagDeepB=Var("bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb')", float, doc="DeepCSV b+bb tag discriminator", precision=10)
        )
    )
    run2_miniAOD_80XLegacy.toModify(process.customAK4Table.variables, jetId=Var("userInt('tightId')*2+userInt('looseId')", int, doc="Jet ID flags bit1 is loose, bit2 is tight"))
    process.customAK4Table.variables.pt.precision = 10

    process.customizedAK4Task = cms.Task(
        process.MuSubProducer,
        process.tightJetIdCustomAK4,
        process.tightJetIdLepVetoCustomAK4,
        process.customAK4WithUserData,
        process.customAK4Table
        )


    _customizedAK4Task_80X = process.customizedAK4Task.copy()
    _customizedAK4Task_80X.replace(process.tightJetIdLepVetoCustomAK4, process.looseJetIdCustomAK4)
    run2_miniAOD_80XLegacy.toReplaceWith(process.customizedAK4Task, _customizedAK4Task_80X)

    if path is None:
        process.schedule.associate(process.customizedAK4Task)
    else:
        getattr(process, path).associate(process.customizedAK4Task)
