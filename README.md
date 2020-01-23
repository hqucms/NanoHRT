# DASZLE Pancakes

### Set up CMSSW

```bash
setenv SCRAM_ARCH slc7_amd64_gcc700
cmsrel CMSSW_10_6_5
cd CMSSW_10_6_5/src
cmsenv
git cms-addpkg PhysicsTools/NanoAOD 
git cms-addpkg PhysicsTools/SelectorUtils 
```

### Get customized NanoAOD producer

```bash
git clone https://github.com/DAZSLE/Pancakes.git PhysicsTools/Pancakes
```

### Compile

```bash
scram b -j 16
```

### Test

```bash
mkdir PhysicsTools/Pancakes/test
cd PhysicsTools/Pancakes/test
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

**Step 1**: generate the python config file with `cmsDriver.py`.
**Step 2**: use the `crab.py` script to submit the CRAB jobs (use --dry-run to test)

These commands are condensed in e.g. `submit_data_UL.sh` or `submit_mc_UL.sh`. They will perform a "dryrun" to print out the CRAB configuration files. Please check everything is correct (e.g., the output path, version number, requested number of cores, etc.) before submitting the actual jobs. To actually submit the jobs to CRAB, just remove the `--dryrun` option at the end.

To check status
`https://monit-grafana.cern.ch/d/cmsTMDetail/cms-task-monitoring-task-view`
