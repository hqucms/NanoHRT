# NanoHRT

### Set up CMSSW

```bash
cmsrel CMSSW_9_4_10
cd CMSSW_9_4_10/src
cmsenv
```

### Use cutomized fastjet-contrib for HOTVR

```bash
scram setup /afs/cern.ch/user/h/hqu/public/tools/fastjet-contrib/slc6_amd64_gcc630/external/fastjet-contrib/1.033HOTVR/fastjet-contrib.xml
```

### Set up MXNet and DeepAK8 recipe

```bash
# MXNet
scram setup /afs/cern.ch/user/h/hqu/public/tools/mxnet/slc6_amd64_gcc630/external/mxnet-predict/1.2.1.mod3/mxnet-predict.xml

# NNKit
git clone ssh://git@gitlab.cern.ch:7999/DeepAK8/NNKit.git -b for94X
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