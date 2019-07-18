from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('outname', 'test.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
)
options.parseArguments()

import FWCore.ParameterSet.Config as cms
process = cms.Process("TestElectrons")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Prompt_v9', '')
process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string('EcalChannelStatusRcd'),
           tag = cms.string('EcalChannelStatus_isolated_deadchannels_data'),
           connect = cms.string("frontier://FrontierPrep/CMS_CONDITIONS")
           ),
)
readFiles=cms.untracked.vstring()

readFiles.extend([
#'file:/eos/cms/store/user/taroni/EGamma/ZEEDeadCh2018_2_pierreTag_Cv3/181219_151536/0000/step3_2018_95.root',
#'file:/eos/cms/store/user/taroni/EGamma/ZEEDeadCh2018_2_pierreTag_Cv3/181219_151536/0000/step3_2018_96.root',
#'file:/eos/cms/store/user/taroni/EGamma/ZEEDeadCh2018_2_pierreTag_Cv3/181219_151536/0000/step3_2018_97.root',
#'file:/eos/cms/store/user/taroni/EGamma/ZEEDeadCh2018_2_pierreTag_Cv3/181219_151536/0000/step3_2018_98.root',
#'file:/eos/cms/store/user/taroni/EGamma/ZEEDeadCh2018_2_pierreTag_Cv3/181219_151536/0000/step3_2018_99.root',
'root://xrootd-cms.infn.it//store/data/Run2018C/EGamma/RAW-RECO/ZElectron-17Sep2018-v1/70000/2087231C-0B8D-CF44-972D-06ADF3129788.root',
'root://xrootd-cms.infn.it//store/data/Run2018C/EGamma/RAW-RECO/ZElectron-17Sep2018-v1/70000/D2AABEF7-A296-444C-95F1-00577878FC24.root',
'root://xrootd-cms.infn.it//store/data/Run2018C/EGamma/RAW-RECO/ZElectron-17Sep2018-v1/70000/A22B58D7-DDAA-4A4E-B4A6-867497DF9343.root',
'root://xrootd-cms.infn.it//store/data/Run2018C/EGamma/RAW-RECO/ZElectron-17Sep2018-v1/70000/8AD4DAF9-0D49-4146-879A-14CCA7CC4EEF.root',
'root://xrootd-cms.infn.it//store/data/Run2018C/EGamma/RAW-RECO/ZElectron-17Sep2018-v1/70000/5B5869B8-8CE0-B94B-B94B-E43F98B6936C.root',
])

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )
inputFiles = readFiles


outputFile = "testZeeD_official.root" if options.outname=='' else options.outname
process.source = cms.Source ("PoolSource", fileNames = inputFiles )  
print outputFile
                 
import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt').getVLuminosityBlockRange()

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


