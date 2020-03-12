import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *
from Configuration.Eras.Modifier_run2_miniAOD_80XLegacy_cff import run2_miniAOD_80XLegacy
from Configuration.Eras.Modifier_run2_nanoAOD_94X2016_cff import run2_nanoAOD_94X2016
from Configuration.Eras.Modifier_run2_nanoAOD_102Xv1_cff import run2_nanoAOD_102Xv1 

# ---------------------------------------------------------


def setupCustomizedAK8(process, runOnMC=False, path=None):
    # recluster Puppi jets, add N-Subjettiness and ECF
    bTagDiscriminators = [
        'pfCombinedInclusiveSecondaryVertexV2BJetTags',
        'pfBoostedDoubleSecondaryVertexAK8BJetTags',
        'pfMassIndependentDeepDoubleBvLJetTags:probHbb',
        'pfMassIndependentDeepDoubleCvLJetTags:probHcc',
        'pfMassIndependentDeepDoubleCvBJetTags:probHcc',
    ]
    subjetBTagDiscriminators = [
        'pfCombinedInclusiveSecondaryVertexV2BJetTags',
        'pfDeepCSVJetTags:probb',
        'pfDeepCSVJetTags:probbb',
    ]
    JETCorrLevels = ['L2Relative', 'L3Absolute', 'L2L3Residual']

    from PhysicsTools.Pancakes.jetToolbox_cff import jetToolbox
    jetToolbox(process, 'ak8', 'dummySeq', 'out', associateTask=False,
               PUMethod='Puppi', JETCorrPayload='AK8PFPuppi', JETCorrLevels=JETCorrLevels,
               Cut='pt > 170.0 && abs(rapidity()) < 2.4',
               miniAOD=True, runOnMC=runOnMC,
               addNsub=True, maxTau=4, addEnergyCorrFunc=True,
	       GetSubjetMCFlavour=False,GetJetMCFlavour=False,
               addSoftDrop=True, addSoftDropSubjets=True, subJETCorrPayload='AK4PFPuppi', subJETCorrLevels=JETCorrLevels,
               bTagDiscriminators=bTagDiscriminators, subjetBTagDiscriminators=subjetBTagDiscriminators)

    if runOnMC:
        process.ak8GenJetsNoNu.jetPtMin = 100
        process.ak8GenJetsNoNuSoftDrop.jetPtMin = 100

    # DeepAK8
    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    from RecoBTag.MXNet.pfDeepBoostedJet_cff import _pfDeepBoostedJetTagsProbs, _pfMassDecorrelatedDeepBoostedJetTagsProbs

    Bdiscs = ['pfDeepFlavourJetTags:probb', 'pfDeepFlavourJetTags:probbb', 'pfDeepFlavourJetTags:probuds', 'pfDeepFlavourJetTags:probg' , 'pfDeepFlavourJetTags:problepb', 'pfDeepFlavourJetTags:probc','pfCombinedInclusiveSecondaryVertexV2BJetTags']

    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

    jetToolbox( process, 'ak8', 'ak8JetSubs', 'out', associateTask=False, 
		updateCollection='packedPatJetsAK8PFPuppiSoftDrop', JETCorrPayload='AK8PFPuppi', JETCorrLevels=JETCorrLevels,
                Cut='pt > 170.0 && abs(rapidity()) < 2.4',
                miniAOD=True, runOnMC=runOnMC,bTagDiscriminators=bTagDiscriminators + _pfDeepBoostedJetTagsProbs + _pfMassDecorrelatedDeepBoostedJetTagsProbs,
		updateCollectionSubjets='selectedPatJetsAK8PFPuppiSoftDropPacked:SubJets', subjetBTagDiscriminators=Bdiscs, 
		subJETCorrPayload='AK4PFPuppi', subJETCorrLevels=JETCorrLevels,postFix='AK8WithPuppiDaughters')

    # src
    srcJets = cms.InputTag('selectedUpdatedPatJetsAK8PFPuppiAK8WithPuppiDaughters')

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
                version=cms.string('WINTER17PUPPI'),
                quality = cms.string('TIGHT'),
              ),
              src=srcJets
    )
    run2_miniAOD_80XLegacy.toModify(process.tightJetIdCustomAK8.filterParams, version="WINTER16")

    process.tightJetIdLepVetoCustomAK8 = cms.EDProducer("PatJetIDValueMapProducer",
              filterParams=cms.PSet(
                version = cms.string('WINTER17PUPPI'),
                quality = cms.string('TIGHTLEPVETO'),
              ),
              src=srcJets
    )

    for modifier in run2_miniAOD_80XLegacy, run2_nanoAOD_94X2016:
    	modifier.toModify( process.tightJetIdCustomAK8.filterParams, version = "WINTER16" )
    	modifier.toModify( process.tightJetIdLepVetoCustomAK8.filterParams, version = "WINTER16" )
    run2_nanoAOD_102Xv1.toModify( process.tightJetIdCustomAK8.filterParams, version = "SUMMER18PUPPI" )
    run2_nanoAOD_102Xv1.toModify( process.tightJetIdLepVetoCustomAK8.filterParams, version = "SUMMER18PUPPI" )

    process.lepInJetVars = cms.EDProducer("LepInJetProducer",
                                          src = srcJets,
                                          srcEle = cms.InputTag("finalElectrons"),
                                          srcMu = cms.InputTag("finalMuons")
                                          )

    jetToolbox(process, 'ak8', 'leptonSubtractedJetSeq', 'out', associateTask=False,
               newPFCollection=True, nameNewPFCollection='lepInJetVars:pfCandsNoLep',	
               PUMethod='Puppi', JETCorrPayload='AK8PFPuppi', JETCorrLevels=JETCorrLevels,
               Cut='pt > 170.0 && abs(rapidity()) < 2.4',
               miniAOD=True, runOnMC=runOnMC,
               addNsub=True, maxTau=4, addEnergyCorrFunc=True,
	       GetSubjetMCFlavour=False,GetJetMCFlavour=False,
               addSoftDrop=True, addSoftDropSubjets=True, subJETCorrPayload='AK4PFPuppi', subJETCorrLevels=JETCorrLevels,
               bTagDiscriminators=bTagDiscriminators, subjetBTagDiscriminators=subjetBTagDiscriminators,
               postFix='LeptonSubtracted',
               )

    leptonSubtractedToEmbed = {
        'tau1': 'userFloat("NjettinessAK8PuppiLeptonSubtracted:tau1")',
        'tau2': 'userFloat("NjettinessAK8PuppiLeptonSubtracted:tau2")',
        'tau3': 'userFloat("NjettinessAK8PuppiLeptonSubtracted:tau3")',
        'tau4': 'userFloat("NjettinessAK8PuppiLeptonSubtracted:tau4")',
        'n2b1': 'userFloat("ak8PFJetsPuppiLeptonSubtractedSoftDropValueMap:nb1AK8PuppiLeptonSubtractedSoftDropN2")',
        'n3b1': 'userFloat("ak8PFJetsPuppiLeptonSubtractedSoftDropValueMap:nb1AK8PuppiLeptonSubtractedSoftDropN3")',
        'msoftdrop': 'groomedMass()',
        'msoftdropraw': 'userFloat("ak8PFJetsPuppiLeptonSubtractedSoftDropMass")',
        'subJet1btagDeepB': '?nSubjetCollections()>0 && subjets().size()>0?subjets()[0].get().bDiscriminator("pfDeepCSVJetTags:probb")+subjets()[0].get().bDiscriminator("pfDeepCSVJetTags:probbb"):-2',
        'subJet2btagDeepB': '?nSubjetCollections()>0 && subjets().size()>1?subjets()[1].get().bDiscriminator("pfDeepCSVJetTags:probb")+subjets()[1].get().bDiscriminator("pfDeepCSVJetTags:probbb"):-2',
        'pt': 'pt()',
    }
    process.matchLeptonSubtracted = cms.EDProducer("PatJetDeltaRValueMapProducer",
                                                   src=srcJets,
                                                   matched=cms.InputTag('packedPatJetsAK8PFPuppiLeptonSubtractedSoftDrop'),
                                                   distMax=cms.double(0.8),
                                                   values=cms.vstring(list(leptonSubtractedToEmbed.values())),
                                                   valueLabels=cms.vstring(list(leptonSubtractedToEmbed.keys())),
                                                   )

    process.customAK8WithUserData = cms.EDProducer("PATJetUserDataEmbedder",
                                                   src=srcJets,
                                                   userFloats = cms.PSet(lsf3 = cms.InputTag("lepInJetVars:lsf3"),
                                                                         dRLep = cms.InputTag("lepInJetVars:dRLep"),
                                                                         **{('leptonSubtracted:' + k): cms.InputTag("matchLeptonSubtracted:"+k) for k in leptonSubtractedToEmbed.keys()}
                                                                         ),
                                                   userInts = cms.PSet(tightId = cms.InputTag("tightJetIdCustomAK8"),
                                                                       tightIdLepVeto = cms.InputTag("tightJetIdLepVetoCustomAK8"),
                                                                       muonIdx3SJ = cms.InputTag("lepInJetVars:muIdx3SJ"),
                                                                       electronIdx3SJ = cms.InputTag("lepInJetVars:eleIdx3SJ"),
                                                                       idLep = cms.InputTag("lepInJetVars:idLep"),
                                                                       ),
    )
    run2_miniAOD_80XLegacy.toModify(process.customAK8WithUserData.userInts,
        looseId=cms.InputTag("looseJetIdCustomAK8"),
        tightIdLepVeto=None,
    )

    process.customAK8Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src=cms.InputTag("customAK8WithUserData"),
        name=cms.string("FatJet"),
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
            rawmsoftdrop=Var("userFloat('ak8PFJetsPuppiSoftDropMass')", float, doc="Raw soft drop mass with PUPPI", precision=10),
            #btagDeepB = Var("bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb')",float,doc="DeepCSV b+bb tag discriminator",precision=10),
            #btagCSVV2=Var("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')", float, doc=" pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)", precision=10),
            btagHbb=Var("bDiscriminator('pfBoostedDoubleSecondaryVertexAK8BJetTags')", float, doc="old Higgs to BB tagger discriminator", precision=10),
            btagDDBvL = Var("bDiscriminator('pfMassIndependentDeepDoubleBvLJetTags:probHbb')",float,doc="DeepDoubleX (mass-decorrelated) discriminator for H(Z)->bb vs QCD",precision=10),
            btagDDCvL = Var("bDiscriminator('pfMassIndependentDeepDoubleCvLJetTags:probHcc')",float,doc="DeepDoubleX (mass-decorrelated) discriminator for H(Z)->cc vs QCD",precision=10),
            btagDDCvB = Var("bDiscriminator('pfMassIndependentDeepDoubleCvBJetTags:probHcc')",float,doc="DeepDoubleX (mass-decorrelated) discriminator for H(Z)->cc vs H(Z)->bb",precision=10),
            subJetIdx1=Var("?nSubjetCollections()>0 && subjets().size()>0?subjets()[0].key():-1", int,
                 doc="index of first subjet"),
            subJetIdx2=Var("?nSubjetCollections()>0 && subjets().size()>1?subjets()[1].key():-1", int,
                 doc="index of second subjet"),
            nPFConstituents=Var("numberOfDaughters()", int, doc="Number of PF daughter constituents"),
            nBHadrons=Var("jetFlavourInfo().getbHadrons().size()", int, doc="number of b-hadrons"),
            nCHadrons=Var("jetFlavourInfo().getcHadrons().size()", int, doc="number of c-hadrons"),
            lsf3 = Var("userFloat('lsf3')",float, doc="LSF (3 subjets)",precision=10),
            dRLep = Var("userFloat('dRLep')", float, doc="dR(lep,jet)",precision=10),
            muonIdx3SJ = Var("userInt('muonIdx3SJ')",int, doc="index of muon matched (3 subjets)"),
            electronIdx3SJ = Var("userInt('electronIdx3SJ')",int,doc="index of electron matched (3 subjets)"),
            idLep =  Var("userInt('idLep')", int, doc="id of pf particle matched to that reco electron or muon"),
            LStau1 = Var("userFloat('leptonSubtracted:tau1')", float, doc="Nsubjettiness, after subtracting highest pt lepton from jet", precision=10),
            LStau2 = Var("userFloat('leptonSubtracted:tau2')", float, doc="Nsubjettiness, after subtracting highest pt lepton from jet", precision=10),
            LStau3 = Var("userFloat('leptonSubtracted:tau3')", float, doc="Nsubjettiness, after subtracting highest pt lepton from jet", precision=10),
            LStau4 = Var("userFloat('leptonSubtracted:tau4')", float, doc="Nsubjettiness, after subtracting highest pt lepton from jet", precision=10),
            LSn2b1 = Var("userFloat('leptonSubtracted:n2b1')", float, doc="N2 with beta=1, after subtracting highest pt lepton from jet", precision=10),
            LSn3b1 = Var("userFloat('leptonSubtracted:n3b1')", float, doc="N3 with beta=1, after subtracting highest pt lepton from jet", precision=10),
            LSmsoftdrop = Var("userFloat('leptonSubtracted:msoftdrop')", float, doc="Softdrop mass with PUPPI, after subtracting highest pt lepton from jet", precision=10),
            LSrawmsoftdrop = Var("userFloat('leptonSubtracted:msoftdropraw')", float, doc="Softdrop mass with PUPPI, after subtracting highest pt lepton from jet", precision=10),
            LSsubJet1btagDeepB = Var("userFloat('leptonSubtracted:subJet1btagDeepB')", float, doc="Subjet DeepCSV b-tag score (b+bb), after subtracting highest pt lepton from jet", precision=10),
            LSsubJet2btagDeepB = Var("userFloat('leptonSubtracted:subJet2btagDeepB')", float, doc="Subjet DeepCSV b-tag score (b+bb), after subtracting highest pt lepton from jet", precision=10),
            LSpt = Var("userFloat('leptonSubtracted:pt')", float, doc="Jet pt after sub", precision=10),
        )
    )

    #run2_miniAOD_80XLegacy.toModify( fatJetTable.variables, n2b1 = None)
    #run2_miniAOD_80XLegacy.toModify( fatJetTable.variables, n3b1 = None)
    run2_miniAOD_80XLegacy.toModify(process.customAK8Table.variables, jetId=Var("userInt('tightId')*2+userInt('looseId')", int, doc="Jet ID flags bit1 is loose, bit2 is tight"))
    process.customAK8Table.variables.pt.precision = 10

    _pfDeepBoostedJetTagsProbs_new = ['pfDeepBoostedJetTags:probWcq','pfDeepBoostedJetTags:probWqq', 'pfDeepBoostedJetTags:probZbb', 'pfDeepBoostedJetTags:probZcc', 'pfDeepBoostedJetTags:probZqq', 'pfDeepBoostedJetTags:probHbb','pfDeepBoostedJetTags:probHcc', 'pfDeepBoostedJetTags:probHqqqq', 'pfDeepBoostedJetTags:probQCDbb', 'pfDeepBoostedJetTags:probQCDcc']

    _pfMassDecorrelatedDeepBoostedJetTagsProbs_new = ['pfMassDecorrelatedDeepBoostedJetTags:probWcq', 'pfMassDecorrelatedDeepBoostedJetTags:probWqq', 'pfMassDecorrelatedDeepBoostedJetTags:probZbb', 'pfMassDecorrelatedDeepBoostedJetTags:probZcc', 'pfMassDecorrelatedDeepBoostedJetTags:probZqq', 'pfMassDecorrelatedDeepBoostedJetTags:probHbb', 'pfMassDecorrelatedDeepBoostedJetTags:probHcc', 'pfMassDecorrelatedDeepBoostedJetTags:probHqqqq', 'pfMassDecorrelatedDeepBoostedJetTags:probQCDbb', 'pfMassDecorrelatedDeepBoostedJetTags:probQCDcc']

    # add DeepAK8 scores: nominal
    for prob in _pfDeepBoostedJetTagsProbs_new:
        name = prob.split(':')[1].replace('prob', 'deepTag')
        setattr(process.customAK8Table.variables, name, Var("bDiscriminator('%s')" % prob, float, doc=prob, precision=-1))

    # add DeepAK8 scores: mass decorrelated
    for prob in _pfMassDecorrelatedDeepBoostedJetTagsProbs_new:
        name = prob.split(':')[1].replace('prob', 'deepTagMD')
        setattr(process.customAK8Table.variables, name, Var("bDiscriminator('%s')" % prob, float, doc=prob, precision=-1))

    process.customAK8SubJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src=cms.InputTag("selectedPatJetsAK8PFPuppiSoftDropPacked", "SubJets"),
        cut=cms.string(""),
        name=cms.string("SubJet"),
        doc=cms.string("customized ak8 puppi subjets for HRT"),
        singleton=cms.bool(False),  # the number of entries is variable
        extension=cms.bool(False),  # this is the main table for the jets
        variables=cms.PSet(P4Vars,
            btagDeepB=Var("bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb')", float, doc="DeepCSV b+bb tag discriminator", precision=10),
            btagCSVV2=Var("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')", float, doc=" pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)", precision=10),
            rawFactor=Var("1.-jecFactor('Uncorrected')", float, doc="1 - Factor to get back to raw pT", precision=6),
            area=Var("jetArea()", float, doc="jet catchment area, for JECs", precision=10),
            nBHadrons=Var("jetFlavourInfo().getbHadrons().size()", int, doc="number of b-hadrons"),
            nCHadrons=Var("jetFlavourInfo().getcHadrons().size()", int, doc="number of c-hadrons"),
        )
    )
    process.customAK8SubJetTable.variables.pt.precision = 10

    process.customizedAK8Task = cms.Task(
        process.tightJetIdCustomAK8,
        process.tightJetIdLepVetoCustomAK8,
        process.lepInJetVars,
        process.matchLeptonSubtracted,
        process.customAK8WithUserData,
        process.customAK8Table,
        process.customAK8SubJetTable,
        )

    if runOnMC:
        process.customGenJetAK8Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
            src=cms.InputTag("ak8GenJetsNoNu"),
            cut=cms.string("pt > 100."),
            name=cms.string("CustomGenJetAK8"),
            doc=cms.string("AK8 GenJets made with visible genparticles"),
            singleton=cms.bool(False),  # the number of entries is variable
            extension=cms.bool(False),  # this is the main table for the genjets
            variables=cms.PSet(P4Vars,
                               #msoftdrop = Var("groomedMass('SoftDropPuppi')",float, doc="Corrected soft drop mass with PUPPI",precision=10), #wishlist
                               )
        )
        process.customGenJetAK8Table.variables.pt.precision = 10

        process.customGenSubJetAK8Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
            src=cms.InputTag("ak8GenJetsNoNuSoftDrop", "SubJets"),
            cut=cms.string(""),
            name=cms.string("CustomGenSubJetAK8"),
            doc=cms.string("AK8 Gen-SubJets made with visible genparticles"),
            singleton=cms.bool(False),  # the number of entries is variable
            extension=cms.bool(False),  # this is the main table for the genjets
            variables=cms.PSet(P4Vars,
            )
        )
        process.customGenSubJetAK8Table.variables.pt.precision = 10

        process.customGenJetAK8Constituents = cms.EDProducer("GenJetPackedConstituentPtrSelector",
                                                             #src = cms.InputTag("slimmedGenJetsAK8"),
                                                             src = cms.InputTag("ak8GenJetsNoNu"),
                                                             cut = cms.string("pt > 100.0")
                                                             )

        process.customGenJetAK8ParticleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                                              src = cms.InputTag("customGenJetAK8Constituents", "constituents"),
                                                              cut = cms.string(""), #we should not filter after pruning
                                                              name= cms.string("GenPartAK8"),
                                                              doc = cms.string("interesting gen particles from AK8 jets"),
                                                              singleton = cms.bool(False), # the number of entries is variable
                                                              extension = cms.bool(False), # this is the main table for the AK8 constituents
                                                              variables = cms.PSet(CandVars,
                                                                                   )
                                                              )

        #process.customizedAK8Task.add(process.customGenJetAK8Table)
        #process.customizedAK8Task.add(process.customGenSubJetAK8Table)

        #process.customizedAK8Task.add(process.customGenJetAK8Constituents)
        #process.customizedAK8Task.add(process.customGenJetAK8ParticleTable)

    _customizedAK8Task_80X = process.customizedAK8Task.copy()
    _customizedAK8Task_80X.replace(process.tightJetIdLepVetoCustomAK8, process.looseJetIdCustomAK8)
    run2_miniAOD_80XLegacy.toReplaceWith(process.customizedAK8Task, _customizedAK8Task_80X)

    # Bacon 15 config, deprecated
    # process.pfInclusiveSecondaryVertexFinderAK8TagInfosAK8PFPuppiAK8WithPuppiDaughters.extSVDeltaRToJet = cms.double(0.3)
    # process.pfInclusiveSecondaryVertexFinderAK8TagInfosAK8PFPuppiAK8WithPuppiDaughters.trackSelection.jetDeltaRMax = cms.double(0.3)
    # process.pfInclusiveSecondaryVertexFinderAK8TagInfosAK8PFPuppiAK8WithPuppiDaughters.vertexCuts.maxDeltaRToJetAxis = cms.double(0.4)
    # process.pfImpactParameterAK8TagInfosAK8PFPuppiAK8WithPuppiDaughters.computeProbabilities = cms.bool(True)

    if path is None:
        process.schedule.associate(process.customizedAK8Task)
    else:
        getattr(process, path).associate(process.customizedAK8Task)


def addCustomizedAK8PF(process):
    if not hasattr(process, 'customizedAK8Task'):
        raise RuntimeError("Call setupCustomizedAK8 first")

    process.customAK8ConstituentsTable = cms.EDProducer("JetConstituentTableProducer",
                                                   #src = cms.InputTag("updatedJetsAK8"),
                                                   src = cms.InputTag("customAK8WithUserData"),
                                                   cut = cms.string("pt()>170"),
                                                   name = cms.string("FatJetPFCands"),
                                                   )

    process.customAK8ConstituentsExtTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                                                        src = cms.InputTag("customAK8ConstituentsTable"),
                                                        cut = cms.string(""), #we should not filter after pruning
                                                        name = cms.string("FatJetPFCands"),
                                                        doc = cms.string("interesting particles from AK8 jets"),
                                                        singleton = cms.bool(False), # the number of entries is variable
                                                        extension = cms.bool(True), # this is the extension table for the AK8 constituents
                                                        variables = cms.PSet(CandVars,
                                                                             puppiWeight = Var("puppiWeight()", float, doc="Puppi weight",precision=10),
                                                                             puppiWeightNoLep = Var("puppiWeightNoLep()", float, doc="Puppi weight removing leptons",precision=10),
                                                                             vtxChi2 = Var("?hasTrackDetails()?vertexChi2():-1", float, doc="vertex chi2", precision=10),
                                                                             trkChi2 = Var("?hasTrackDetails()?pseudoTrack().normalizedChi2():-1", float, doc="normalized trk chi2", precision=10),
                                                                             dz=Var("dz()", float, doc="pf dz", precision=10),
                                                                             dzErr = Var("?hasTrackDetails()?dzError():-1", float, doc="pf dz err", precision=10),
                                                                             d0=Var("dxy()", float, doc="pf d0", precision=10),
                                                                             d0Err = Var("?hasTrackDetails()?dxyError():-1", float, doc="pf d0 err", precision=10),
                                                                             pvAssocQuality = Var("pvAssociationQuality()", int, doc="primary vertex association quality"),
                                                                             lostInnerHits = Var("lostInnerHits()", int, doc="lost inner hits"),
                                                                             trkQuality = Var("?hasTrackDetails()?pseudoTrack().qualityMask():0", int, doc="track quality mask"),
                                                                             )
                                                        )

    process.customizedAK8Task.add(process.customAK8ConstituentsTable)
    process.customizedAK8Task.add(process.customAK8ConstituentsExtTable)
    return process
