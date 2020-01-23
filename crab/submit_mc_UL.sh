cmsDriver.py pancakes_2017_mc_UL -n 100 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 106X_mc2017_realistic_v6 --step NANO --nThreads 4 \
    --era Run2_2017 --customise PhysicsTools/Pancakes/nanoHRT_cff.nanoHRT_customizeMC --fileout file:nano_mc_2017.root \
    --customise PhysicsTools/Pancakes/ak8_cff.addCustomizedAK8PF \
    --filein file:miniAOD.root --fileout file:nano_mc_2017.root \
    --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring

python crab.py -p pancakes_2017_mc_UL_NANO.py -o /store/group/lpcbacon/pancakes/02/2017/ -t pancakes-02 -i mc_2017_UL.txt  --num-cores 4  -s EventAwareLumiBased \
    -n 200000 --work-area crab_projects_mc_2017_UL --max-memory 6000 --dryrun
