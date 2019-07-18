

import FWCore.ParameterSet.Config as cms
process = cms.Process("TestElectrons")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15', '')
process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string('EcalChannelStatusRcd'),
           tag = cms.string('EcalChannelStatus_isolated_deadchannels_mc'),
           connect = cms.string("frontier://FrontierPrep/CMS_CONDITIONS")
           ),
)
readFiles=cms.untracked.vstring()

readFiles.extend([
'file:/eos/cms/store/user/taroni/DYToEE_M-50_NNPDF31_TuneCP5_13TeV-powheg-pythia8/DYEE_pierreTag/newtest_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO_PU_noRecov.root'
])

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
inputFiles = readFiles
outputFile = "electronTree_PierreTag_MCPU_noRecov_reweight.root"
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
    isMC         = cms.bool(True)
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

process.load("DPGAnalysis/Skims/ZElectronSkim_cff") 
process.p = cms.Path(process.zdiElectronSequence*process.ntupler)


