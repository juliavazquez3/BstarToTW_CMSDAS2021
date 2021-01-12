# BstarToTW_CMSDAS2021

## Getting started (in bash shell)

### Setup CMSSW environment:
```
ssh -XY USERNAME@cmslpc-sl7.fnal.gov
cexport $SCRAM_ARCH=slc7_amd64_gcc820 
cd nobackup/
mkdir b2g_exercise/
cd b2g_exercise/
cmsrel CMSSW_11_0_1
cd CMSSW_11_0_1/src
cmsenv
```

### In the `src` directory, clone repo:
```
git clone https://github.com/cmantill/BstarToTW_CMSDAS2021.git
```
OR fork the code onto your personal project space and set the upstream:
```
git clone https://github.com/<GitHubUsername>/BstarToTW_CMSDAS2021.git
cd BstarToTW_CMSDAS2021
git remote add upstream https://github.com/cmantill/BstarToTW_CMSDAS2021.git
git remote -v
```

### In the `src` directory, create environment
```
git clone https://github.com/lcorcodilos/TIMBER.git
python -m virtualenv timber-env
source timber-env/bin/activate
cd TIMBER
source setup.sh
cd ..
python -c 'import TIMBER.Analyzer'
```

## Starting up once environment is set:

Once you have an environment:
```
cd CMSSW_11_0_1/src/
cmsenv
source timber-env/bin/activate
cd TIMBER/
git fetch --all
git checkout master
python setup.py install
cd ../BstarToTW_CMSDAS2020
git fetch --all
git pull origin master
```
