import FWCore.ParameterSet.Config as cms

pfDeepBoostedJetPreprocessParams = cms.PSet(
    input_names = cms.vstring('pfcand'),
    pfcand = cms.PSet(
        input_shape=cms.vuint32(1, 9, 100, 1),
        var_infos = cms.PSet(
            pfcand_abseta=cms.PSet(
                median=cms.double(0.599505603313),
                upper=cms.double(1.21494185925)
            ),
            pfcand_deltaR = cms.PSet(
                median = cms.double(0.22575956583),
                upper = cms.double(0.488191870451)
            ),
            pfcand_drsubjet1 = cms.PSet(
                median = cms.double(0.231124095619),
                upper = cms.double(0.549522156715)
            ),
            pfcand_drsubjet2 = cms.PSet(
                median = cms.double(0.263272643089),
                upper = cms.double(0.605471189022)
            ),
            pfcand_erel_log = cms.PSet(
                median = cms.double(-5.38983869553),
                upper = cms.double(-3.53490426064)
            ),
            pfcand_etarel = cms.PSet(
                median = cms.double(-0.0054658302106),
                upper = cms.double(0.174858552814)
            ),
            pfcand_phirel = cms.PSet(
                median = cms.double(-5.10289683007e-05),
                upper = cms.double(0.215602903366)
            ),
            pfcand_pt_log = cms.PSet(
                median = cms.double(1.09469842911),
                upper = cms.double(3.02194809914)
            ),
            pfcand_ptrel_log = cms.PSet(
                median = cms.double(-5.38205528259),
                upper = cms.double(-3.52304198265)
            )
        ),
        var_length = cms.uint32(100),
        var_names = cms.vstring(
            'pfcand_pt_log', 
            'pfcand_ptrel_log', 
            'pfcand_erel_log', 
            'pfcand_phirel', 
            'pfcand_etarel', 
            'pfcand_deltaR', 
            'pfcand_abseta',
            'pfcand_drsubjet1', 
            'pfcand_drsubjet2'
        )
    )
)
