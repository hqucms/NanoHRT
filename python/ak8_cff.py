import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *
from Configuration.Eras.Modifier_run2_miniAOD_80XLegacy_cff import run2_miniAOD_80XLegacy

# ---------------------------------------------------------


def setupCustomizedAK8(process, runOnMC=False, path=None):
    # recluster Puppi jets, add N-Subjettiness and ECF
    bTagDiscriminators = [
        'pfCombinedInclusiveSecondaryVertexV2BJetTags',
        'pfBoostedDoubleSecondaryVertexAK8BJetTags',
    ]
    subjetBTagDiscriminators = [
        'pfCombinedInclusiveSecondaryVertexV2BJetTags',
        'pfDeepCSVJetTags:probb',
        'pfDeepCSVJetTags:probbb',
    ]
    JETCorrLevels = ['L2Relative', 'L3Absolute', 'L2L3Residual']

    from PhysicsTools.NanoHRT.jetToolbox_cff import jetToolbox
    jetToolbox(process, 'ak8', 'dummySeq', 'out', associateTask=False,
               PUMethod='Puppi', JETCorrPayload='AK8PFPuppi', JETCorrLevels=JETCorrLevels,
               Cut='pt > 170.0 && abs(rapidity()) < 2.4',
               miniAOD=True, runOnMC=runOnMC,
               addNsub=True, maxTau=4, addEnergyCorrFunc=True,
               addSoftDrop=True, addSoftDropSubjets=True, subJETCorrPayload='AK4PFPuppi', subJETCorrLevels=JETCorrLevels,
               bTagDiscriminators=bTagDiscriminators, subjetBTagDiscriminators=subjetBTagDiscriminators)

    # DeepAK8
    process.deepBoostedJetsAK8Puppi = cms.EDProducer('DeepBoostedJetProducer',
        src=cms.InputTag('packedPatJetsAK8PFPuppiSoftDrop'),
        hasPuppiWeightedDaughters=cms.bool(True),
        jet_radius=cms.untracked.double(0.8),
        nominal_nn_path=cms.untracked.string('NNKit/data/ak8/full'),
        decorrelated_nn_path=cms.untracked.string('NNKit/data/ak8/decorrelated'),
    )

    # BEST
    process.boostedEventShapeJetsAK8Puppi = cms.EDProducer('BESTProducer',
        src=cms.InputTag('deepBoostedJetsAK8Puppi'),
        config_path=cms.untracked.FileInPath('PhysicsTools/NanoHRT/data/BEST/config.txt'),
        dnn_path=cms.untracked.FileInPath('PhysicsTools/NanoHRT/data/BEST/BEST_6bin_PUPPI.json'),
    )

    # src
    srcJets = cms.InputTag('boostedEventShapeJetsAK8Puppi')

    # jetID
    process.looseJetIdCustomAK8 = cms.EDProducer("PatJetIDValueMapProducer",
              filterParams=cms.PSet(
                version = cms.string('WINTER16'),
                quality = cms.string('LOOSE'),
              ),
              src=srcJets
    )

    process.tightJetIdCustomAK8 = cms.EDProducer("PatJetIDValueMapProducer",
              filterParams=cms.PSet(
                version=cms.string('WINTER17'),
                quality = cms.string('TIGHT'),
              ),
              src=srcJets
    )
    run2_miniAOD_80XLegacy.toModify(process.tightJetIdCustomAK8.filterParams, version="WINTER16")

    process.tightJetIdLepVetoCustomAK8 = cms.EDProducer("PatJetIDValueMapProducer",
              filterParams=cms.PSet(
                version = cms.string('WINTER17'),
                quality = cms.string('TIGHTLEPVETO'),
              ),
              src=srcJets
    )

    process.customAK8WithUserData = cms.EDProducer("PATJetUserDataEmbedder",
        src=srcJets,
        userFloats=cms.PSet(),
        userInts=cms.PSet(
           tightId=cms.InputTag("tightJetIdCustomAK8"),
           tightIdLepVeto=cms.InputTag("tightJetIdLepVetoCustomAK8"),
        ),
    )
    run2_miniAOD_80XLegacy.toModify(process.customAK8WithUserData.userInts,
        looseId=cms.InputTag("looseJetIdCustomAK8"),
        tightIdLepVeto=None,
    )

    process.customAK8Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src=cms.InputTag("customAK8WithUserData"),
        name=cms.string("CustomAK8Puppi"),
        cut=cms.string(""),
        doc=cms.string("customized ak8 puppi jets for HRT"),
        singleton=cms.bool(False),  # the number of entries is variable
        extension=cms.bool(False),  # this is the main table for the jets
        variables=cms.PSet(P4Vars,
            jetId=Var("userInt('tightId')*2+4*userInt('tightIdLepVeto')", int, doc="Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto"),
            area=Var("jetArea()", float, doc="jet catchment area, for JECs", precision=10),
            rawFactor=Var("1.-jecFactor('Uncorrected')", float, doc="1 - Factor to get back to raw pT", precision=6),
            tau1=Var("userFloat('NjettinessAK8Puppi:tau1')", float, doc="Nsubjettiness (1 axis)", precision=10),
            tau2=Var("userFloat('NjettinessAK8Puppi:tau2')", float, doc="Nsubjettiness (2 axis)", precision=10),
            tau3=Var("userFloat('NjettinessAK8Puppi:tau3')", float, doc="Nsubjettiness (3 axis)", precision=10),
            tau4=Var("userFloat('NjettinessAK8Puppi:tau4')", float, doc="Nsubjettiness (4 axis)", precision=10),
            n2b1=Var("userFloat('ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN2')", float, doc="N2 with beta=1", precision=10),
            n3b1=Var("userFloat('ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN3')", float, doc="N3 with beta=1", precision=10),
            msoftdrop=Var("groomedMass()", float, doc="Corrected soft drop mass with PUPPI", precision=10),
            btagCSVV2=Var("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')", float, doc=" pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)", precision=10),
            btagHbb=Var("bDiscriminator('pfBoostedDoubleSecondaryVertexAK8BJetTags')", float, doc="Higgs to BB tagger discriminator", precision=10),
            subJetIdx1=Var("?nSubjetCollections()>0 && subjets().size()>0?subjets()[0].key():-1", int,
                 doc="index of first subjet"),
            subJetIdx2=Var("?nSubjetCollections()>0 && subjets().size()>1?subjets()[1].key():-1", int,
                 doc="index of second subjet"),
            # DeepAK8: nominal
            nnTbcq=Var("userFloat('DeepBoostedJet:nn_Top_bcq')", float, doc="DeepAK8 score Top_bcq", precision=-1),
            nnTbqq=Var("userFloat('DeepBoostedJet:nn_Top_bqq')", float, doc="DeepAK8 score Top_bqq", precision=-1),
            nnTbc=Var("userFloat('DeepBoostedJet:nn_Top_bc')", float, doc="DeepAK8 score Top_bc", precision=-1),
            nnTbq=Var("userFloat('DeepBoostedJet:nn_Top_bq')", float, doc="DeepAK8 score Top_bq", precision=-1),
            nnWcq=Var("userFloat('DeepBoostedJet:nn_W_cq')", float, doc="DeepAK8 score W_cq", precision=-1),
            nnWqq=Var("userFloat('DeepBoostedJet:nn_W_qq')", float, doc="DeepAK8 score W_qq", precision=-1),
            nnZbb=Var("userFloat('DeepBoostedJet:nn_Z_bb')", float, doc="DeepAK8 score Z_bb", precision=-1),
            nnZcc=Var("userFloat('DeepBoostedJet:nn_Z_cc')", float, doc="DeepAK8 score Z_cc", precision=-1),
            nnZqq=Var("userFloat('DeepBoostedJet:nn_Z_qq')", float, doc="DeepAK8 score Z_qq", precision=-1),
            nnHbb=Var("userFloat('DeepBoostedJet:nn_H_bb')", float, doc="DeepAK8 score nnHbb", precision=-1),
            nnHcc=Var("userFloat('DeepBoostedJet:nn_H_cc')", float, doc="DeepAK8 score nnHcc", precision=-1),
            nnHqqqq=Var("userFloat('DeepBoostedJet:nn_H_qqqq')", float, doc="DeepAK8 score nnHqqqq", precision=-1),
            nnQCDbb=Var("userFloat('DeepBoostedJet:nn_QCD_bb')", float, doc="DeepAK8 score QCD_bb", precision=-1),
            nnQCDcc=Var("userFloat('DeepBoostedJet:nn_QCD_cc')", float, doc="DeepAK8 score QCD_cc", precision=-1),
            nnQCDb=Var("userFloat('DeepBoostedJet:nn_QCD_b')", float, doc="DeepAK8 score QCD_b", precision=-1),
            nnQCDc=Var("userFloat('DeepBoostedJet:nn_QCD_c')", float, doc="DeepAK8 score QCD_c", precision=-1),
            nnQCDothers=Var("userFloat('DeepBoostedJet:nn_QCD_others')", float, doc="DeepAK8 score QCD_others", precision=-1),
            # DeepAK8: decorrelated
            nnMDTbcq=Var("userFloat('DeepBoostedJet:decorr_nn_Top_bcq')", float, doc="Mass-decorrelated DeepAK8 score Top_bcq", precision=-1),
            nnMDTbqq=Var("userFloat('DeepBoostedJet:decorr_nn_Top_bqq')", float, doc="Mass-decorrelated DeepAK8 score Top_bqq", precision=-1),
            nnMDTbc=Var("userFloat('DeepBoostedJet:decorr_nn_Top_bc')", float, doc="Mass-decorrelated DeepAK8 score Top_bc", precision=-1),
            nnMDTbq=Var("userFloat('DeepBoostedJet:decorr_nn_Top_bq')", float, doc="Mass-decorrelated DeepAK8 score Top_bq", precision=-1),
            nnMDWcq=Var("userFloat('DeepBoostedJet:decorr_nn_W_cq')", float, doc="Mass-decorrelated DeepAK8 score W_cq", precision=-1),
            nnMDWqq=Var("userFloat('DeepBoostedJet:decorr_nn_W_qq')", float, doc="Mass-decorrelated DeepAK8 score W_qq", precision=-1),
            nnMDZbb=Var("userFloat('DeepBoostedJet:decorr_nn_Z_bb')", float, doc="Mass-decorrelated DeepAK8 score Z_bb", precision=-1),
            nnMDZcc=Var("userFloat('DeepBoostedJet:decorr_nn_Z_cc')", float, doc="Mass-decorrelated DeepAK8 score Z_cc", precision=-1),
            nnMDZqq=Var("userFloat('DeepBoostedJet:decorr_nn_Z_qq')", float, doc="Mass-decorrelated DeepAK8 score Z_qq", precision=-1),
            nnMDHbb=Var("userFloat('DeepBoostedJet:decorr_nn_H_bb')", float, doc="Mass-decorrelated DeepAK8 score nnHbb", precision=-1),
            nnMDHcc=Var("userFloat('DeepBoostedJet:decorr_nn_H_cc')", float, doc="Mass-decorrelated DeepAK8 score nnHcc", precision=-1),
            nnMDHqqqq=Var("userFloat('DeepBoostedJet:decorr_nn_H_qqqq')", float, doc="Mass-decorrelated DeepAK8 score nnHqqqq", precision=-1),
            nnMDQCDbb=Var("userFloat('DeepBoostedJet:decorr_nn_QCD_bb')", float, doc="Mass-decorrelated DeepAK8 score QCD_bb", precision=-1),
            nnMDQCDcc=Var("userFloat('DeepBoostedJet:decorr_nn_QCD_cc')", float, doc="Mass-decorrelated DeepAK8 score QCD_cc", precision=-1),
            nnMDQCDb=Var("userFloat('DeepBoostedJet:decorr_nn_QCD_b')", float, doc="Mass-decorrelated DeepAK8 score QCD_b", precision=-1),
            nnMDQCDc=Var("userFloat('DeepBoostedJet:decorr_nn_QCD_c')", float, doc="Mass-decorrelated DeepAK8 score QCD_c", precision=-1),
            nnMDQCDothers=Var("userFloat('DeepBoostedJet:decorr_nn_QCD_others')", float, doc="Mass-decorrelated DeepAK8 score QCD_others", precision=-1),
            # BEST Tagger
            bestT=Var("userFloat('BEST:dnn_top')", float, doc="Boosted Event Shape Tagger score Top", precision=10),
            bestW=Var("userFloat('BEST:dnn_w')", float, doc="Boosted Event Shape Tagger score W", precision=10),
            bestZ=Var("userFloat('BEST:dnn_z')", float, doc="Boosted Event Shape Tagger score Z", precision=10),
            bestH=Var("userFloat('BEST:dnn_higgs')", float, doc="Boosted Event Shape Tagger score Higgs", precision=10),
            bestQCD=Var("userFloat('BEST:dnn_qcd')", float, doc="Boosted Event Shape Tagger score QCD", precision=10),
            bestB=Var("userFloat('BEST:dnn_b')", float, doc="Boosted Event Shape Tagger score B", precision=10),
        )
    )
    run2_miniAOD_80XLegacy.toModify(process.customAK8Table.variables, jetId=Var("userInt('tightId')*2+userInt('looseId')", int, doc="Jet ID flags bit1 is loose, bit2 is tight"))
    process.customAK8Table.variables.pt.precision = 10

    process.customAK8SubJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src=cms.InputTag("selectedPatJetsAK8PFPuppiSoftDropPacked", "SubJets"),
        cut=cms.string(""),
        name=cms.string("CustomAK8PuppiSubJet"),
        doc=cms.string("customized ak8 puppi subjets for HRT"),
        singleton=cms.bool(False),  # the number of entries is variable
        extension=cms.bool(False),  # this is the main table for the jets
        variables=cms.PSet(P4Vars,
            btagDeepB=Var("bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb')", float, doc="DeepCSV b+bb tag discriminator", precision=10),
            btagCSVV2=Var("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')", float, doc=" pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)", precision=10),
            rawFactor=Var("1.-jecFactor('Uncorrected')", float, doc="1 - Factor to get back to raw pT", precision=6),
            area=Var("jetArea()", float, doc="jet catchment area, for JECs", precision=10),
        )
    )
    process.customAK8SubJetTable.variables.pt.precision = 10

    process.customizedAK8Task = cms.Task(
        process.deepBoostedJetsAK8Puppi,
        process.boostedEventShapeJetsAK8Puppi,
        process.tightJetIdCustomAK8,
        process.tightJetIdLepVetoCustomAK8,
        process.customAK8WithUserData,
        process.customAK8Table,
        process.customAK8SubJetTable
        )

    _customizedAK8Task_80X = process.customizedAK8Task.copy()
    _customizedAK8Task_80X.replace(process.tightJetIdLepVetoCustomAK8, process.looseJetIdCustomAK8)
    run2_miniAOD_80XLegacy.toReplaceWith(process.customizedAK8Task, _customizedAK8Task_80X)

    if path is None:
        process.schedule.associate(process.customizedAK8Task)
    else:
        getattr(process, path).associate(process.customizedAK8Task)
# ---------------------------------------------------------
