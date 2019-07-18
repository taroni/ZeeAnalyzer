import FWCore.ParameterSet.Config as cms
process = cms.Process("TestElectrons")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Prompt_v9', '')
readFiles=cms.untracked.vstring()



readFiles.extend([
'file:/eos/cms/store/user/taroni/EGamma/noChangePierreTag102X/EGamma_Zee_102X_pierre_0.root',
'file:/eos/cms/store/user/taroni/EGamma/noChangePierreTag102X/EGamma_Zee_102X_pierre_1.root',
'file:/eos/cms/store/user/taroni/EGamma/noChangePierreTag102X/EGamma_Zee_102X_pierre_2.root',
])
#lines=[]
#for line in open('filelist_Official.txt'):
#    lines.append('%s' %(line.rstrip('\n')))
#readFiles.extend(lines)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
inputFiles = readFiles

outputFile = "electronTreeZEE_PierreTag_noChange.root"
process.source = cms.Source ("PoolSource", fileNames = inputFiles )  
events= cms.untracked.VEventRange()
#for line in open('eventList.txt'):
#    events.append(line.rstrip("\n"))
#process.source.eventsToProcess = cms.untracked.VEventRange(events)
    
                 
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


