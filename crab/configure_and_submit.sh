cmsDriver.py test_nanoHRT_mc -n 1000 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM \
  --conditions 102X_mc2017_realistic_v7 --step NANO --nThreads 4 --era Run2_2017,run2_nanoAOD_94XMiniAODv2 \
  --customise PhysicsTools/NanoHRT/nanoHRT_cff.nanoHRT_customizeMC --fileout file:nano_mc.root \
  --filein dummy.root # "dbs:/QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/MINIAODSIM" \
  --no_exec # --customise Configuration/DataProcessing/Utils.addMonitoring

python crab.py -p test_nanoHRT_mc_NANO.py -o /store/user/ncsmith/pancakes/01/ -t pancakes-01 -i mc_2017.txt \
  --num-cores 4 --max-memory 8000 --send-external -s EventAwareLumiBased -n 100000 --work-area crab_projects_mc_2017 # --dryrun
