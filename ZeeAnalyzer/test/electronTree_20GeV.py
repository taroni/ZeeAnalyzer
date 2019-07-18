

import FWCore.ParameterSet.Config as cms
process = cms.Process("TestElectrons")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string('EcalChannelStatusRcd'),
           tag = cms.string('EcalChannelStatus_isolated_deadchannels_data'),
           connect = cms.string("frontier://FrontierPrep/CMS_CONDITIONS")
           ),
)
readFiles=cms.untracked.vstring()
count=0
start=50*0
end=50+50*0

#lines=[]
#for line in open('pierreOutputFiles.txt'):
#    if (count >= start and count < end):
#        lines.append('file:%s' %(line.rstrip('\n')))
#    count+=1

readFiles.extend([
#'file:/eos/cms/store/user/taroni/EGamma/ZEEDeadCh2018_2_pierreTag_Cv3/181219_151536/step3_Zee_sum8gt0.root',
'file:/eos/cms/store/user/taroni/EGamma/ZeeDeadCh2018_sum8gt20_pierre/EGamma_sum8gt20_pierre_2.root',
'file:/eos/cms/store/user/taroni/EGamma/ZeeDeadCh2018_sum8gt20_pierre/EGamma_sum8gt20_pierre_1.root',
'file:/eos/cms/store/user/taroni/EGamma/ZeeDeadCh2018_sum8gt20_pierre/EGamma_sum8gt20_pierre_0.root',
])

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
inputFiles = readFiles
outputFile = "electronTreeZEE_PierreTag_sum8gt20_xtalInfo.root"
process.source = cms.Source ("PoolSource", fileNames = inputFiles )                             

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

import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt').getVLuminosityBlockRange()


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

process.load("DPGAnalysis/Skims/ZElectronSkim_cff") 
process.p = cms.Path(process.zdiElectronSequence*process.ntupler)


