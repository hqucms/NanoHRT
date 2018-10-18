import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *
from Configuration.Eras.Modifier_run2_miniAOD_80XLegacy_cff import run2_miniAOD_80XLegacy

# ---------------------------------------------------------


def setupCA15(process, runOnMC=False, path=None):
    # recluster Puppi jets, add N-Subjettiness and ECF
    bTagDiscriminators = [
        'pfCombinedInclusiveSecondaryVertexV2BJetTags',
    ]
    subjetBTagDiscriminators = [
        'pfCombinedInclusiveSecondaryVertexV2BJetTags',
        'pfDeepCSVJetTags:probb',
        'pfDeepCSVJetTags:probbb',
    ]
    JETCorrLevels = ['L2Relative', 'L3Absolute', 'L2L3Residual']

    from PhysicsTools.NanoHRT.jetToolbox_cff import jetToolbox
    jetToolbox(process, 'ca15', 'dummySeq', 'out', associateTask=False,
               PUMethod='Puppi', JETCorrPayload='AK8PFPuppi', JETCorrLevels=JETCorrLevels,
               Cut='pt > 170.0 && abs(rapidity()) < 2.4',
               miniAOD=True, runOnMC=runOnMC,
               addNsub=False, maxTau=3, addEnergyCorrFunc=False,
               addSoftDrop=True, addSoftDropSubjets=True, subJETCorrPayload='AK4PFPuppi', subJETCorrLevels=JETCorrLevels,
               bTagDiscriminators=bTagDiscriminators, subjetBTagDiscriminators=subjetBTagDiscriminators)

    if runOnMC:
        process.ca15GenJetsNoNu.jetPtMin = 100
        process.ca15GenJetsNoNuSoftDrop.jetPtMin = 100

    # ECFTopTag
    process.ecfTopTagCA15Puppi = cms.EDProducer('ECFTopTagsProducer',
        src=cms.InputTag('packedPatJetsCA15PFPuppiSoftDrop'),
        bdt_path=cms.untracked.FileInPath('PhysicsTools/NanoHRT/data/ECFTopTag/top_ecfbdt_v8_BDT.weights.xml'),
        bdt_name=cms.untracked.string('v0'),
    )

    # src
    srcJets = cms.InputTag('ecfTopTagCA15Puppi')

    # jetID
    process.looseJetIdCA15Puppi = cms.EDProducer("PatJetIDValueMapProducer",
              filterParams=cms.PSet(
                version=cms.string('WINTER16'),
                quality=cms.string('LOOSE'),
              ),
              src=srcJets
    )

    process.tightJetIdCA15Puppi = cms.EDProducer("PatJetIDValueMapProducer",
              filterParams=cms.PSet(
                version=cms.string('WINTER17'),
                quality=cms.string('TIGHT'),
              ),
              src=srcJets
    )
    run2_miniAOD_80XLegacy.toModify(process.tightJetIdCA15Puppi.filterParams, version="WINTER16")

    process.tightJetIdLepVetoCA15Puppi = cms.EDProducer("PatJetIDValueMapProducer",
              filterParams=cms.PSet(
                version=cms.string('WINTER17'),
                quality=cms.string('TIGHTLEPVETO'),
              ),
              src=srcJets
    )

    process.ca15WithUserData = cms.EDProducer("PATJetUserDataEmbedder",
        src=srcJets,
        userFloats=cms.PSet(),
        userInts=cms.PSet(
           tightId=cms.InputTag("tightJetIdCA15Puppi"),
           tightIdLepVeto=cms.InputTag("tightJetIdLepVetoCA15Puppi"),
        ),
    )
    run2_miniAOD_80XLegacy.toModify(process.ca15WithUserData.userInts,
        looseId=cms.InputTag("looseJetIdCA15Puppi"),
        tightIdLepVeto=None,
    )

    process.ca15Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src=cms.InputTag("ca15WithUserData"),
        name=cms.string("CA15Puppi"),
        cut=cms.string(""),
        doc=cms.string("ca15 puppi jets"),
        singleton=cms.bool(False),  # the number of entries is variable
        extension=cms.bool(False),  # this is the main table for the jets
        variables=cms.PSet(P4Vars,
            jetId=Var("userInt('tightId')*2+4*userInt('tightIdLepVeto')", int, doc="Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto"),
            area=Var("jetArea()", float, doc="jet catchment area, for JECs", precision=10),
            rawFactor=Var("1.-jecFactor('Uncorrected')", float, doc="1 - Factor to get back to raw pT", precision=6),
            msoftdrop=Var("groomedMass()", float, doc="Corrected soft drop mass with PUPPI", precision=10),
            btagCSVV2=Var("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')", float, doc=" pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)", precision=10),
            subJetIdx1=Var("?nSubjetCollections()>0 && subjets().size()>0?subjets()[0].key():-1", int,
                 doc="index of first subjet"),
            subJetIdx2=Var("?nSubjetCollections()>0 && subjets().size()>1?subjets()[1].key():-1", int,
                 doc="index of second subjet"),
            ecf0=Var("userFloat('ecf_0')", float, doc="ecfN_1_2_20/pow(ecfN_1_2_10,2.00)", precision=10),
            httFRec=Var("userFloat('httFRec')", float, doc="HTT frec", precision=10),
            tau32sd=Var("userFloat('tau32sd')", float, doc="soft drop tau32", precision=10),
            ecfTopTagBDT=Var("userFloat('ecfTopTagBDT')", float, doc="ECF top BDT", precision=10),
        )
    )
    run2_miniAOD_80XLegacy.toModify(process.ca15Table.variables, jetId=Var("userInt('tightId')*2+userInt('looseId')", int, doc="Jet ID flags bit1 is loose, bit2 is tight"))
    process.ca15Table.variables.pt.precision = 10

    process.ca15SubJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src=cms.InputTag("selectedPatJetsCA15PFPuppiSoftDropPacked", "SubJets"),
        cut=cms.string(""),
        name=cms.string("CA15PuppiSubJet"),
        doc=cms.string("ca15 puppi subjets for HRT"),
        singleton=cms.bool(False),  # the number of entries is variable
        extension=cms.bool(False),  # this is the main table for the jets
        variables=cms.PSet(P4Vars,
            btagDeepB=Var("bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb')", float, doc="DeepCSV b+bb tag discriminator", precision=10),
            btagCSVV2=Var("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')", float, doc=" pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)", precision=10),
            rawFactor=Var("1.-jecFactor('Uncorrected')", float, doc="1 - Factor to get back to raw pT", precision=6),
            area=Var("jetArea()", float, doc="jet catchment area, for JECs", precision=10),
        )
    )
    process.ca15SubJetTable.variables.pt.precision = 10

    process.ca15Task = cms.Task(
        process.ecfTopTagCA15Puppi,
        process.tightJetIdCA15Puppi,
        process.tightJetIdLepVetoCA15Puppi,
        process.ca15WithUserData,
        process.ca15Table,
        process.ca15SubJetTable
        )

    if runOnMC:
        process.genJetCA15Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
            src=cms.InputTag("ca15GenJetsNoNu"),
            cut=cms.string("pt > 100."),
            name=cms.string("GenJetCA15"),
            doc=cms.string("CA15 GenJets made with visible genparticles"),
            singleton=cms.bool(False),  # the number of entries is variable
            extension=cms.bool(False),  # this is the main table for the genjets
            variables=cms.PSet(P4Vars,
            )
        )
        process.genJetCA15Table.variables.pt.precision = 10

        process.genSubJetCA15Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
            src=cms.InputTag("ca15GenJetsNoNuSoftDrop", "SubJets"),
            cut=cms.string(""),
            name=cms.string("GenSubJetCA15"),
            doc=cms.string("CA15 Gen-SubJets made with visible genparticles"),
            singleton=cms.bool(False),  # the number of entries is variable
            extension=cms.bool(False),  # this is the main table for the genjets
            variables=cms.PSet(P4Vars,
            )
        )
        process.genSubJetCA15Table.variables.pt.precision = 10

        process.ca15Task.add(process.genJetCA15Table)
        process.ca15Task.add(process.genSubJetCA15Table)

    _ca15Task_80X = process.ca15Task.copy()
    _ca15Task_80X.replace(process.tightJetIdLepVetoCA15Puppi, process.looseJetIdCA15Puppi)
    run2_miniAOD_80XLegacy.toReplaceWith(process.ca15Task, _ca15Task_80X)

    if path is None:
        process.schedule.associate(process.ca15Task)
    else:
        getattr(process, path).associate(process.ca15Task)

# ---------------------------------------------------------
