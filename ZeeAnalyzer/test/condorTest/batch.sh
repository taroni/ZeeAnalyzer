#! /bin/bash

sbox=$PWD
cd /afs/cern.ch/user/t/taroni/work/private/newDeadCh102X/src/ZeeAnalyzer/ZeeAnalyzer/test/condorTest
eval `scramv1 runtime -sh`
export XRD_NETWORKSTACK=IPv4
export X509_USER_PROXY=x509up_u29820
cd $sbox

cmsRun /afs/cern.ch/user/t/taroni/work/private/newDeadCh102X/src/ZeeAnalyzer/ZeeAnalyzer/test/condorTest/electronTree_Zee_official.py $@
