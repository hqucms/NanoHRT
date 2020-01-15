cmsDriver.py test_nanoHRT_mc -n 1000 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM \
  --conditions 102X_mc2017_realistic_v7 --step NANO --nThreads 4 --era Run2_2017,run2_nanoAOD_94XMiniAODv2 \
  --customise PhysicsTools/NanoHRT/nanoHRT_cff.nanoHRT_customizeMC \
  --customise PhysicsTools/NanoHRT/ak8_cff.addCustomizedAK8PF \
  --filein "dbs:/ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" \
  --fileout file:nano_mc.root \
  --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring

# lumis
# 2016: https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt
# 2017: https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt
# 2018: https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt
# nevts per job
#"2016_sig":20000,"2017_sig":20000,"2018_sig":20000,
#"2016_ttbar":200000,"2017_ttbar":300000,"2018_ttbar":300000,
#"2016_qcd":50000,"2017_qcd":50000,"2018_qcd":50000,
#"2016_mc":200000,"2017_mc":300000,"2018_mc":300000,
#"2016_data":200000,"2017_data":200000,"2018_data":200000,
python crab.py -p test_nanoHRT_mc_NANO.py -o /store/user/ncsmith/pancakes/02/ -t pancakes-02 -i mc_2017.txt \
  --num-cores 4 --max-memory 6000 --send-external -s EventAwareLumiBased -n 40000 --work-area crab_projects_mc_2017 \
  --dryrun
