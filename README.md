# NanoHRT

### Set up CMSSW

```bash
cmsrel CMSSW_10_2_14
cd CMSSW_10_2_14/src
cmsenv
```

### Apply latest changes from offical NanoAOD repo

```bash
# pull updates from official NanoAOD repo
git cms-merge-topic -u cms-nanoAOD:master-102X
```

### Get customized NanoAOD producers for HeavyResTagging

```bash
git clone https://github.com/hqucms/NanoHRT.git PhysicsTools/NanoHRT -b prod/102X
```

### Compile

```bash
scram b -j16
```

### Test

```bash
mkdir PhysicsTools/NanoHRT/test
cd PhysicsTools/NanoHRT/test
```

MC (2017, 94X, MiniAODv2):

```bash
cmsDriver.py test_nanoHRT_mc2017 -n 1000 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_mc2017_realistic_v6 --step NANO --nThreads 4 --era Run2_2017,run2_nanoAOD_94XMiniAODv2 --customise PhysicsTools/NanoHRT/nanoTuples_cff.nanoTuples_customizeMC --customise_commands 'process.particleLevelSequence.remove(process.genParticles2HepMCHiggsVtx);process.particleLevelSequence.remove(process.rivetProducerHTXS);process.particleLevelTables.remove(process.HTXSCategoryTable)' --filein /store/mc/RunIIFall17MiniAODv2/ZprimeToTT_M3000_W30_TuneCP2_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/20FF99D9-702A-E911-B801-0025904CFB86.root --fileout file:nano_mc2017.root --customise_commands "process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))" >& test_mc2017.log &

less +F test_mc2017.log
```

Data (2017, 94X, MiniAODv2):

```bash
cmsDriver.py test_nanoHRT_data2017 -n 1000 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v8 --step NANO --nThreads 4 --era Run2_2017,run2_nanoAOD_94XMiniAODv2 --customise PhysicsTools/NanoHRT/nanoTuples_cff.nanoTuples_customizeData --filein /store/data/Run2017D/JetHT/MINIAOD/31Mar2018-v1/60000/1EEE02D3-E539-E811-9859-0025905A6066.root --fileout file:nano_data2017.root --customise_commands "process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))" >& test_data2017.log &

less +F test_data2017.log
```


MC (2018, 102X):

```bash
cmsDriver.py test_nanoHRT_mc2018 -n 1000 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v18 --step NANO --nThreads 4 --era Run2_2018,run2_nanoAOD_102Xv1 --customise PhysicsTools/NanoHRT/nanoTuples_cff.nanoTuples_customizeMC --customise_commands 'process.particleLevelSequence.remove(process.genParticles2HepMCHiggsVtx);process.particleLevelSequence.remove(process.rivetProducerHTXS);process.particleLevelTables.remove(process.HTXSCategoryTable)' --filein /store/mc/RunIIAutumn18MiniAOD/ZprimeToTT_M3000_W30_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/90000/1CFAC15C-895C-CD44-BC86-58EE90CBF456.root --fileout file:nano_mc2018.root --customise_commands "process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))" >& test_mc2018.log &

less +F test_mc2018.log
```

Data (2018, 102X):

```bash
cmsDriver.py test_nanoHRT_data2018 -n 1000 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_Sep2018ABC_v2 --step NANO --nThreads 4 --era Run2_2018 --customise PhysicsTools/NanoHRT/nanoTuples_cff.nanoTuples_customizeData --filein /store/data/Run2018C/JetHT/MINIAOD/17Sep2018-v1/80000/DDC38B74-3A1C-BF4B-9B01-11A3A6A4078A.root --fileout file:nano_data2018.root --customise_commands "process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))" >& test_data2018.log &

less +F test_data2018.log
```


### Production

**Step 0**: switch to the crab production directory and set up grid proxy, CRAB environment, etc.

```bash
cd $CMSSW_BASE/src/PhysicsTools/NanoHRT/crab
# set up grid proxy
voms-proxy-init -rfc -voms cms --valid 168:00
# set up CRAB env (must be done after cmsenv)
source /cvmfs/cms.cern.ch/crab3/crab.sh
```

**Step 1**: generate the python config file with `cmsDriver.py` with the following commands:

MC (2017, 94X, MiniAODv2):

```bash
cmsDriver.py test_nanoHRT_mc2017 -n -1 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_mc2017_realistic_v6 --step NANO --nThreads 4 --era Run2_2017,run2_nanoAOD_94XMiniAODv2 --customise PhysicsTools/NanoHRT/nanoTuples_cff.nanoTuples_customizeMC --customise_commands 'process.particleLevelSequence.remove(process.genParticles2HepMCHiggsVtx);process.particleLevelSequence.remove(process.rivetProducerHTXS);process.particleLevelTables.remove(process.HTXSCategoryTable)' --filein file:step-1.root --fileout file:nano.root --no_exec
```

Data (2017, 94X, MiniAODv2):

```bash
cmsDriver.py data_2017 -n -1 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_v8 --step NANO --nThreads 4 --era Run2_2017,run2_nanoAOD_94XMiniAODv2 --customise PhysicsTools/NanoHRT/nanoTuples_cff.nanoTuples_customizeData --filein file:step-1.root --fileout file:nano.root --no_exec
```


MC (2018, 102X):

```bash
cmsDriver.py mc_2018 -n -1 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 102X_upgrade2018_realistic_v18 --step NANO --nThreads 4 --era Run2_2018,run2_nanoAOD_102Xv1 --customise PhysicsTools/NanoHRT/nanoTuples_cff.nanoTuples_customizeMC --customise_commands 'process.particleLevelSequence.remove(process.genParticles2HepMCHiggsVtx);process.particleLevelSequence.remove(process.rivetProducerHTXS);process.particleLevelTables.remove(process.HTXSCategoryTable)' --filein file:step-1.root --fileout file:nano.root --no_exec
```

Data (2018, 102X):

[Note] GT with JEC V8 not ready yet. For now use separate GTs for 2018ABC and 2018D. 

For 2018ABC:

```bash
cmsDriver.py data_2018ABC -n -1 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_Sep2018ABC_v2 --step NANO --nThreads 4 --era Run2_2018 --customise PhysicsTools/NanoHRT/nanoTuples_cff.nanoTuples_customizeData --filein file:step-1.root --fileout file:nano.root --no_exec
```

For 2018D:

```bash
cmsDriver.py data_2018D -n -1 --data --eventcontent NANOAOD --datatier NANOAOD --conditions 102X_dataRun2_Prompt_v13 --step NANO --nThreads 4 --era Run2_2018 --customise PhysicsTools/NanoHRT/nanoTuples_cff.nanoTuples_customizeData --filein file:step-1.root --fileout file:nano.root --no_exec
```

**Step 2**: use the `crab.py` script to submit the CRAB jobs:

For MC:

`python crab.py -p mc_[yyyy]_NANO.py -o /store/group/lpcjme/noreplica/NanoHRT/mc/[version] -t NanoHRT-[version] -i mc_[ABC].txt --num-cores 4 --send-external -s EventAwareLumiBased -n 50000 --work-area crab_projects_mc_[ABC] --dryrun`

For data:

`python crab.py -p data_[yyyy]_NANO.py -j [JSON] -o /store/group/lpcjme/noreplica/NanoHRT/data/[version] -t NanoHRT-[version] -i data.txt --num-cores 4 --send-external -s EventAwareLumiBased -n 50000 --work-area crab_projects_data --dryrun`

A JSON file can be applied for data samples with the `-j` options.

Golden JSON, 2016:

```
https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
```

Golden JSON, 2017:

```
https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt
```

Golden JSON, 2018:

```
https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt
```

These command will perform a "dryrun" to print out the CRAB configuration files. Please check everything is correct (e.g., the output path, version number, requested number of cores, etc.) before submitting the actual jobs. To actually submit the jobs to CRAB, just remove the `--dryrun` option at the end.

**Step 3**: check job status

The status of the CRAB jobs can be checked with:

```bash
./crab.py --status --work-area crab_projects_[ABC]
```

Note that this will also resubmit failed jobs automatically.

The crab dashboard can also be used to get a quick overview of the job status:
`https://dashb-cms-job.cern.ch/dashboard/templates/task-analysis`

More options of this `crab.py` script can be found with:

```bash
./crab.py -h
```
