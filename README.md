# NanoHRT

### Set up CMSSW

```bash
cmsrel CMSSW_9_4_10
cd CMSSW_9_4_10/src
cmsenv
```

### Set up DeepAK8 recipe

```bash
# MXNet
scram setup /cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw/CMSSW_10_3_0_pre4/config/toolbox/slc6_amd64_gcc700/tools/selected/mxnet-predict.xml

git clone ssh://git@gitlab.cern.ch:7999/DeepAK8/NNKit.git -b for94X-reclustered-jets
```

### Get customized NanoAOD producers for HeavyResTagging

```bash
git clone https://github.com/hqucms/NanoHRT.git PhysicsTools/NanoHRT
```

### Compile

```bash
scram b -j16
```

### Test

```bash
cd PhysicsTools/NanoHRT/test
cmsRun nanoHRT_cfg.py
```

### Production

MC:

```bash
cmsDriver.py mc -n -1 --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --conditions 94X_mcRun2_asymptotic_v2 --step NANO --nThreads 4 --era Run2_2016,run2_miniAOD_80XLegacy --customise PhysicsTools/NanoHRT/nanoHRT_cff.nanoHRT_customizeMC --filein file:step-1.root --fileout file:nano.root --no_exec
```