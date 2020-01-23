cmsDriver.py pancakes_2017_data_UL -n -1 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 106X_dataRun2_v20 --step NANO --nThreads 2 \
    --era Run2_2017,run2_nanoAOD_106Xv1 \
    --customise PhysicsTools/Pancakes/nanoHRT_cff.nanoHRT_customizeData \
    --customise PhysicsTools/Pancakes/ak8_cff.addCustomizedAK8PF \
    --filein file:miniAOD.root --fileout file:nano_data_2017.root \
    --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring

python crab.py -p pancakes_2017_data_UL_NANO.py -o /store/group/lpcbacon/pancakes/02/2017/ -t pancakes-02 -i data_2017_UL.txt  --num-cores 2 --send-external \
    -s EventAwareLumiBased \
    -n 200000 --work-area crab_projects_data_2017 --max-memory 5000 \
    --json https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt 
#--dryrun

# cmsDriver.py pancakes_2018_data_UL -n -1 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 106X_dataRun2_v20 --step NANO --nThreads 2 \
#     --era Run2_2018,run2_nanoAOD_106Xv1 \
#     --customise PhysicsTools/Pancakes/nanoHRT_cff.nanoHRT_customizeData \
#     --customise PhysicsTools/Pancakes/ak8_cff.addCustomizedAK8PF \
#     --filein file:miniAOD.root --fileout file:nano_data_2018.root \
#     --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring

# python crab.py -p pancakes_2018_data_UL_NANO.py -o /store/group/lpcbacon/pancakes/02/2018/ -t pancakes-02 -i data_2018_UL.txt  --num-cores 2 --send-external \
#     -s EventAwareLumiBased \
#     -n 200000 --work-area crab_projects_data_2018 --max-memory 6000 \
#     --json https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt --dryrun

# cmsDriver.py pancakes_2016_data_UL -n -1 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 106X_dataRun2_v20 --step NANO --nThreads 2 \
#     --era Run2_2016,run2_nanoAOD_106Xv1 \
#     --customise PhysicsTools/Pancakes/nanoHRT_cff.nanoHRT_customizeData \
#     --customise PhysicsTools/Pancakes/ak8_cff.addCustomizedAK8PF \
#     --filein file:miniAOD.root --fileout file:nano_data_2016.root \
#     --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring

# python crab.py -p pancakes_2016_data_UL_NANO.py -o /store/group/lpcbacon/pancakes/02/2016/ -t pancakes-02 -i data_2016_UL.txt  --num-cores 2 --send-external \
#     -s EventAwareLumiBased \
#     -n 200000 --work-area crab_projects_data_2016 --max-memory 6000 \
#     --json https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt --dryrun
