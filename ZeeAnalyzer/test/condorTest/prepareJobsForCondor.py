import os 

#first job
i=0
numberOfJobs=100
lastjob=numberOfJobs+i

while i<lastjob:
    cfg="""import math
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RERECO',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.MessageLogger.cerr.FwkReport.reportEvery = 50

readFiles=cms.untracked.vstring()
# Input source
inputFiles = readFiles
process.source = cms.Source("PoolSource",
    fileNames = inputFiles       
)
events= cms.untracked.VEventRange()
#for line in open('recoveredEvt2018Cv3.txt'):
#    events.append(line.rstrip("\\n"))
#process.source.eventsToProcess = cms.untracked.VEventRange(events)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_PromptLike_v10', '')

process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string('EcalChannelStatusRcd'),
           tag = cms.string('EcalChannelStatus_isolated_deadchannels_data'),
           connect = cms.string("frontier://FrontierPrep/CMS_CONDITIONS")
           ),
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

from Configuration.EventContent.EventContent_cff import RAWRECOEventContent
process.skimContent = process.RAWRECOEventContent.clone()
process.load("DPGAnalysis.Skims.filterRecHitsRecovery_cfi")
process.recoveryfilter = cms.Path(process.recHitRecoveryFilter)
from RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi import ecalRecHit
process.ecalRecHit.singleChannelRecoveryThreshold=0.7
process.ecalRecHit.sum8ChannelRecoveryThreshold=20.

# Output definition
outputFile = "testZeeC_recov_"""+str(i)+""".root"


process.RAWRECOoutput = cms.OutputModule("PoolOutputModule",
                             
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RAW-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(outputFile),
    outputCommands = process.RAWRECOEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0),
                                         )
process.RAWRECOoutput.outputCommands.append('drop *_*_*_RECO')


process.ntupler = cms.EDAnalyzer(
    'ElectronTree',
    beamSpot = cms.InputTag('offlineBeamSpot'),
    genEventInfoProduct = cms.InputTag('generator'),
    electrons    = cms.InputTag("gedGsfElectrons"),
    genParticles = cms.InputTag("genParticles"),
    vertices     = cms.InputTag("offlinePrimaryVertices"),
    conversions  = cms.InputTag('allConversions'),
    inputRecHitsEB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    isMC         = cms.bool(False)
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("eleTree_"+outputFile)
                                   )

process.load("DPGAnalysis/Skims/ZElectronSkim_cff") 
process.p = cms.Path(process.recoverySequence*process.zdiElectronSequence*process.ntupler)

# Path and EndPath definitions

process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWRECOoutput_step = cms.EndPath(process.RAWRECOoutput)

##process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step, process.RAWRECOoutput_step , process.endjob_step, process.p)


#number of input file per job: 44
mf = open(\"ZeeFiles.txt\")
inputfiles=[x.strip(\"\\n\") for x in mf]
filePerJob=int(math.ceil(float(len(inputfiles))/float("""+str(numberOfJobs)+""")))
start=filePerJob*"""+str(i)+"""
end=filePerJob+filePerJob*"""+str(i)+"""
count=0
print len(inputfiles), """+str(numberOfJobs)+""", start, end

for ii in inputfiles:
    if (count >= start and count < end):
        process.source.fileNames.append("root://xrootd-cms.infn.it/"+ii)
        
    count+=1

"""

    cfg_file = open("ZeeRecovery_"+str(i)+"_cfg.py","w")
    cfg_file.write(cfg)
    cfg_file.close()


    #crea lo script per lanciare su cluster
    sh = """#!/bin/tcsh -f
    
            set W_DIR = \"/afs/cern.ch/work/t/taroni/private/JulyDeadCh/src/ZeeAnalyzer/ZeeAnalyzer/test/condorTest\"
            set CFG = \"/afs/cern.ch/work/t/taroni/private/JulyDeadCh/src/ZeeAnalyzer/ZeeAnalyzer/test/condorTest/ZeeRecovery_"""+str(i)+"""_cfg.py\"
            set FILELIST = \"/afs/cern.ch/user/t/taroni/work/private/JulyDeadCh/src/ZeeAnalyzer/ZeeAnalyzer/test/condorTest/ZeeFiles.txt\"
            cp $FILELIST .
            cd $W_DIR
            eval `scram runtime -csh`
            cd -
            cmsRun $CFG
            cp testZeeC_recov_"""+str(i)+""".root /eos/cms/store/user/taroni/condorTest/.
            cp eleTree_testZeeC_recov_"""+str(i)+""".root /eos/cms/store/user/taroni/condorTest/.
            
    """
    sh_file = open("batchZee_"+str(i)+".sh","w")
    sh_file.write(sh)
    sh_file.close()
    
    #make it executable
    os.popen("chmod a+x batchZee_"+str(i)+".sh" )

    i+=1

cndsh="""
use_x509userproxy = true
x509userproxy = /afs/cern.ch/user/t/taroni/work/private/JulyDeadCh/src/ZeeAnalyzer/ZeeAnalyzer/test/condorTest/x509up_u29820
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
getenv = True
executable = /afs/cern.ch/user/t/taroni/work/private/JulyDeadCh/src/ZeeAnalyzer/ZeeAnalyzer/test/condorTest/batchZee_$(ProcId).sh
+MaxRuntime = 21600
#the queue to use: list of the que in slide 55 of the tutorial https://indico.cern.ch/event/783287/ 
+JobFlavour = "tomorrow"
Output = testZee_$(ProcId).out
Error = testZee_$(ProcId).err
Log = testZee_$(ProcId).log
#the number after 'queue' is the number of jobs
queue """+str(numberOfJobs)+"""
"""
cndsh_file=open("condor.jdl", "w")
cndsh_file.write(cndsh)
cndsh_file.close()

