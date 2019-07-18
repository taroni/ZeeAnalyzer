

import FWCore.ParameterSet.Config as cms
process = cms.Process("TestElectrons")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15', '')
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
#'file:/eos/cms/store/group/dpg_ecal/comm_ecal/taroni/DYToEE_M-50_NNPDF31_TuneCP5_13TeV-powheg-pythia8/recoveryMC_pierreTag_15GeV/181123_143703/0000/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_15GeV_81.root'
        #'file:test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_15GeV_onefile.root'
        #'file:/eos/cms/store/group/dpg_ecal/comm_ecal/taroni/DYToEE_M-50_NNPDF31_TuneCP5_13TeV-powheg-pythia8/recoveryMC2_pierreTag_15GeV/181128_114055/0000/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_15GeV_81.root'
        #'file:/eos/cms/store/group/dpg_ecal/comm_ecal/taroni/DYToEE_M-50_NNPDF31_TuneCP5_13TeV-powheg-pythia8/recoveryMC2_pierreTag_15GeV/181128_114055/0000/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_15GeV_81.root'
'file:/eos/cms/store/group/dpg_ecal/comm_ecal/taroni/DYToEE_M-50_NNPDF31_TuneCP5_13TeV-powheg-pythia8/recoveryMC100j_pierreTag_15GeV/181128_145459/0000/test_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_15GeV_81.root'
])

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
inputFiles = readFiles
outputFile = "electronTree_PierreTag_15GeV_81.root"
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

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

process.load("DPGAnalysis/Skims/ZElectronSkim_cff") 
process.p = cms.Path(process.zdiElectronSequence*process.ntupler)


