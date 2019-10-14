# DASZLE Pancakes

### Set up CMSSW

```bash
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git-cms-addpkg PhysicsTools/NanoAOD 
git-cms-addpkg PhysicsTools/SelectorUtils 
```

### Get customized NanoAOD producer

```bash
git clone https://github.com/DAZSLE/NanoHRT.git PhysicsTools/NanoHRT
```

### Compile

```bash
scram b -j 16
```

### Test

```bash
mkdir PhysicsTools/NanoHRT/test
cd PhysicsTools/NanoHRT/test
```

MC:

2018: 102X_upgrade2018_realistic_v19

ex. file: /store/mc/RunIIAutumn18MiniAOD/BulkGravTohhTohWWhbb_narrow_M-2300_TuneCP2_13TeV-madgraph_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/260000/24516D8A-B053-1A42-A10F-07EA8D96FE6C.root
```bash
cmsDriver.py test_nanoHRT_mc --filein /store/mc/RunIIAutumn18MiniAOD/BulkGravTohhTohWWhbb_narrow_M-2300_TuneCP2_13TeV-madgraph_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/260000/24516D8A-B053-1A42-A10F-07EA8D96FE6C.root --fileout file:RunIIAutumn18NanoAODv5_BulkGravTohhTohWWhbb.root --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v19 --step NANO --nThreads 2 --era Run2_2018,run2_nanoAOD_102Xv1 --python_filename RunIIAutumn18NanoAODv5_pancakes01_mc_cfg.py --customise PhysicsTools/NanoHRT/nanoHRT_cff.nanoHRT_customizeMC -n 10 --no_exec 
```

Data:

2018: 102X_dataRun2_v11

ex. file: /store/data/Run2018B/JetHT/MINIAOD/17Sep2018-v1/60000/FE3C69F0-A0BC-8941-92AF-B0DA1A6270BF.root
```bash
cmsDriver.py test_nanoHRT_data --filein /store/data/Run2018B/JetHT/MINIAOD/17Sep2018-v1/60000/FE3C69F0-A0BC-8941-92AF-B0DA1A6270BF.root --fileout file:RunIIAutumn18NanoAODv5_JetHTRun2018B.root --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v11 --eventcontent NANOAOD --era Run2_2018,run2_nanoAOD_102Xv1 --customise PhysicsTools/NanoHRT/nanoHRT_cff.nanoHRT_customizeData --python_filename=RunIIAutumn18NanoAODv5_pancakes01_data_cfg.py -s NANO --no_exec 
```

### Production

**Step 0**: switch to the crab production directory and set up grid proxy, CRAB environment, etc.

```bash
mkdir $CMSSW_BASE/src/PhysicsTools/NanoHRT/crab
cd $CMSSW_BASE/src/PhysicsTools/NanoHRT/crab
# set up grid proxy
voms-proxy-init -rfc -voms cms --valid 168:00
# set up CRAB env (must be done after cmsenv)
source /cvmfs/cms.cern.ch/crab3/crab.sh
```

**Step 1**: generate the python config file with `cmsDriver.py` with the following commands - these are included for now:

MC (102X, MiniAODv2):

```bash
cmsDriver.py mc -n -1 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v19 --step NANO --nThreads 2 --era Run2_2018,run2_nanoAOD_102Xv1 --customise PhysicsTools/NanoHRT/nanoHRT_cff.nanoHRT_customizeMC --filein file:miniAOD.root --fileout file:RunIIAutumn18NanoAODv5.root --python_filename RunIIAutumn18NanoAODv5_pancakes01_mc_cfg.py --no_exec
```

Data (`2018` PromptReco-17Jul2018):

```bash
cmsDriver.py data -n -1 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v11 --step NANO --nThreads 2 --era Run2_2018,run2_nanoAOD_102Xv1 --customise PhysicsTools/NanoHRT/nanoHRT_cff.nanoHRT_customizeData --filein file:miniAOD.root --fileout file:RunIIAutumn18NanoAODv5.root --python_filename RunIIAutumn18NanoAODv5_pancakes01_data_cfg.py --no_exec
```

**Step 2**: use the `crab.py` script to submit the CRAB jobs (use --dry-run to test), e.g.: 

version = 01

For MC (2018):
`python crab.py -p RunIIAutumn18NanoAODv5_pancakes01_mc_cfg.py -o /store/group/[group]/[username]/pancakes/[version]/ -t pancakes-[version] -i signal_2018.txt --num-cores 2 --send-external -s EventAwareLumiBased -n 50000 --work-area crab_projects_mc_2018 --max-memory 5000 --dryrun`

For data (2018):

`python crab.py -p RunIIAutumn18NanoAODv5_pancakes01_data_cfg.py  -o /store/group/[group]/[username]/pancakes/[version] -t pancakes-[version] -i data_2018.txt --num-cores 2 --send-external -s EventAwareLumiBased -n 200000 --work-area crab_projects_data_2018 --max-memory 5000 --json https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt --dryrun`

A JSON file can be applied for data samples with the `-j` options. By default, we use the golden JSON for 2016:

```
https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
```

These command will perform a "dryrun" to print out the CRAB configuration files. Please check everything is correct (e.g., the output path, version number, requested number of cores, etc.) before submitting the actual jobs. To actually submit the jobs to CRAB, just remove the `--dryrun` option at the end.

To check status
`https://monit-grafana.cern.ch/d/cmsTMDetail/cms-task-monitoring-task-view`