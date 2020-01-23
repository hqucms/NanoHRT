cmsDriver.py pancakes_2017_mc_UL -n 100 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 106X_mc2017_realistic_v6 --step NANO --nThreads 4 \
    --era Run2_2017 --customise PhysicsTools/Pancakes/nanoHRT_cff.nanoHRT_customizeMC \
    --customise PhysicsTools/Pancakes/ak8_cff.addCustomizedAK8PF \
    --filein file:miniAOD.root --fileout file:nano_mc_2017.root \
    --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring

python crab.py -p pancakes_2017_mc_UL_NANO.py -o /store/group/lpcbacon/pancakes/02/2017/ -t pancakes-02 -i mc_2017_UL.txt --num-cores 4  -s EventAwareLumiBased \
    -n 100000 --work-area crab_projects_mc_2017 --max-memory 5000 --dryrun

cmsDriver.py pancakes_2017_mc -n 100 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_mc2017_realistic_v7 --step NANO --nThreads 4 \
    --era Run2_2017,run2_nanoAOD_94XMiniAODv2 --customise PhysicsTools/Pancakes/nanoHRT_cff.nanoHRT_customizeMC \
    --customise PhysicsTools/Pancakes/ak8_cff.addCustomizedAK8PF \
    --filein file:miniAOD.root --fileout file:nano_mc_2017.root \
    --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring

python crab.py -p pancakes_2017_mc_NANO.py -o /store/group/lpcbacon/pancakes/02/2017/ -t pancakes-02 -i mc_2017.txt \
    --num-cores 4 --max-memory 6000 --send-external -s EventAwareLumiBased -n 50000 --work-area crab_projects_mc_2017 --dryrun
