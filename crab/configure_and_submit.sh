cmsDriver.py test_nanoHRT_mc -n 1000 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM \
  --conditions 102X_mc2017_realistic_v7 --step NANO --nThreads 4 --era Run2_2017,run2_nanoAOD_94XMiniAODv2 \
  --customise PhysicsTools/NanoHRT/nanoHRT_cff.nanoHRT_customizeMC --filein file:step-1.root --fileout file:nano_mc.root \
  --no_exec
python crab.py -p test_nanoHRT_mc_NANO.py -o /store/user/ncsmith/pancakes/01/ -t pancakes-01 -i mc_2017.txt \
  --num-cores 4 --send-external -s EventAwareLumiBased -n 50000 --work-area crab_projects_mc_2017
