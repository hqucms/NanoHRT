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
#python crab.py -p RunIIAutumn18NanoAODv5_pancakes01_mc_cfg.py -o /store/group/lpchbb/cmantill/pancakes/01/ -t pancakes-01 -i signal_2018.txt --num-cores 2 --send-external -s EventAwareLumiBased -n 50000 --work-area crab_projects_mc_2018 --max-memory 8000
python crab.py -p RunIIAutumn18NanoAODv5_pancakes01_data_cfg.py  -o /store/group/lpchbb/cmantill/pancakes/01/ -t pancakes-01 -i data_2018.txt --num-cores 2 --send-external -s EventAwareLumiBased -n 200000 --work-area crab_projects_data_2018 --max-memory 8000 --json https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt