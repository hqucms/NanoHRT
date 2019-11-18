# 2016 sample is junk anyway (too small, wrong PDF)
# cmsDriver.py pancakes_2016_mc -n 100 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM \
#   --conditions 102X_mcRun2_asymptotic_v7 --step NANO --nThreads 4 --era Run2_2016,run2_nanoAOD_94X2016 \
#   --customise PhysicsTools/NanoHRT/nanoHRT_cff.nanoHRT_customizeMC --fileout file:nano_mc_2016.root \
#   --filein "dbs:/GluGluHToBB_M-125_13TeV_powheg_MINLO_NNLOPS_pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM" \
#   --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring
# 
# python crab.py -p pancakes_2016_mc_NANO.py -o /store/user/ncsmith/pancakes/02/ -t pancakes-02 -i minlo2016.txt \
#   --num-cores 4 --max-memory 6000 --send-external -s EventAwareLumiBased -n 40000 --work-area crab_projects_mc_2016 \
#   --dryrun

cmsDriver.py pancakes_2017_mc -n 100 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM \
  --conditions 102X_mc2017_realistic_v7 --step NANO --nThreads 4 --era Run2_2017,run2_nanoAOD_94XMiniAODv2 \
  --customise PhysicsTools/NanoHRT/nanoHRT_cff.nanoHRT_customizeMC --fileout file:nano_mc_2017.root \
  --filein "dbs:/GluGluHToBB_M-125_13TeV_powheg_MINLO_NNLOPS_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM" \
  --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring \
  --customise PhysicsTools/NanoHRT/nanoHRT_cff.nanoHRT_customizeMINLOnnlops

python crab.py -p pancakes_2017_mc_NANO.py -o /store/user/ncsmith/pancakes/02/ -t pancakes-02 -i minlo2017.txt \
  --num-cores 4 --max-memory 6000 --send-external -s EventAwareLumiBased -n 40000 --work-area crab_projects_mc_2017 \

cmsDriver.py pancakes_2018_mc -n 100 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM \
  --conditions 102X_upgrade2018_realistic_v19 --step NANO --nThreads 4 --era Run2_2018,run2_nanoAOD_102Xv1 \
  --customise PhysicsTools/NanoHRT/nanoHRT_cff.nanoHRT_customizeMC --fileout file:nano_mc_2018.root \
  --filein "dbs:/GluGluHToBB_M-125_13TeV_powheg_MINLO_NNLOPS_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM" \
  --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring \
  --customise PhysicsTools/NanoHRT/nanoHRT_cff.nanoHRT_customizeMINLOnnlops

python crab.py -p pancakes_2018_mc_NANO.py -o /store/user/ncsmith/pancakes/02/ -t pancakes-02 -i minlo2018.txt \
  --num-cores 4 --max-memory 6000 --send-external -s EventAwareLumiBased -n 40000 --work-area crab_projects_mc_2018 \
