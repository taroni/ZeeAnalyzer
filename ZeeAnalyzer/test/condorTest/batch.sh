#! /bin/bash

sbox=$PWD
cd /afs/cern.ch/work/t/taroni/private/newDeadCh102X/src/ZmmAnalyser/ZmmAnalyser/test/fromPrasanna/RecvsCentral/condorTest/
eval `scramv1 runtime -sh`
export XRD_NETWORKSTACK=IPv4
export X509_USER_PROXY=x509up_u29820
cd $sbox

iRun=$@
#echo 'arguments' $iRun
cp /afs/cern.ch/work/t/taroni/private/newDeadCh102X/src/ZmmAnalyser/ZmmAnalyser/test/fromPrasanna/RecvsCentral/condorTest/prasanna_File$iRun.txt.
edmCopyPickMerge inputFiles_load=zeeSkim$iRun.txt outputFile=pickeventsPrasanna$iRun.root maxSize=10000000
cp pickeventsPrasanna$iRun.root /eos/cms/store/user/taroni/DoubleMuRecovPrasanna/condorTest/.

