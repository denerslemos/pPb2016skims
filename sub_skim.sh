#!/bin/bash

echo "Setup CMSSW (ROOT version)"
cd /afs/cern.ch/work/d/ddesouza/UIC/pPbMultAna/CMSSW_13_0_5/src/
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/d/ddesouza/UIC/pPbMultAna/CMSSW_13_0_5/src/JetSkims/pPb2016skims
echo "Submit skim jobs at "
echo PWD: $PWD

./pPbSkim $1 $2 $3 $4 $5
