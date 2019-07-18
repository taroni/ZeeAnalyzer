

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
'file:/eos/cms/store/user/taroni/EGamma/EGamma_Zee_official.root'
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/369/00000/E2C382D4-7192-E811-945E-FA163E6A4783.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/344/00000/AA6954A2-5092-E811-B750-FA163E9C3E47.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/340/00000/7C736028-4792-E811-BDC5-02163E017F3A.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/338/00000/74F20D4E-4892-E811-A61A-02163E00F6AC.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/328/00000/4616BBA2-2D92-E811-A63A-FA163E6888C0.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/321/00000/0A21A116-0F92-E811-8BFC-FA163E02F83C.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/285/00000/D0EC2C76-D391-E811-A088-FA163E754652.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/277/00000/0E082D17-CA91-E811-BC32-02163E01A132.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/195/00000/68345347-6B91-E811-A911-FA163E761FF5.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/192/00000/80DADA95-4A91-E811-948A-FA163E92F260.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/191/00000/8CD591B1-1B91-E811-9347-FA163E9C2CB6.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/190/00000/DE4D608F-1391-E811-B429-FA163ED4037D.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/189/00000/D83C665C-0C91-E811-A380-FA163ED5778A.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/186/00000/D66E03DA-0291-E811-BFDB-FA163EE6485C.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/090/00000/8ADFBE9C-0B90-E811-919B-FA163E439461.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/077/00000/64DB6FED-0190-E811-9559-FA163EA8ED32.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/F6B8CCAE-F18F-E811-8846-A4BF0112F7B8.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/F41870D2-F48F-E811-A278-FA163EBFD315.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/F255F355-0D90-E811-A313-02163E019FD7.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/F23A19C3-F48F-E811-9784-FA163E33047A.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/EE8922DC-EB8F-E811-A875-FA163EA82112.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/ECF95C73-E88F-E811-B9B1-FA163EF0B127.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/ECD69B33-EE8F-E811-A2C2-FA163E5234B3.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/ECBF9B08-0890-E811-BB0E-02163E015FB2.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/EC974828-E68F-E811-806D-FA163EDFFB85.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/EC8C61E6-F18F-E811-B9E2-02163E015FB2.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/EC782898-F18F-E811-A5BC-FA163E664457.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/EAF5E5CD-F18F-E811-BB61-FA163EF4E3E9.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/E87926B1-F18F-E811-9052-FA163E200265.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/E6E728DC-0790-E811-A08D-FA163EE76A9C.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/E631CC40-F28F-E811-B5F8-FA163EB29A29.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/E61BA664-E78F-E811-A4EC-FA163E104BF4.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/E4C0125C-F48F-E811-B724-FA163EC26C98.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/E48387B2-E88F-E811-8379-FA163E80361E.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/DED852EC-F38F-E811-BD09-FA163E4BA778.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/DEC95051-F58F-E811-A336-02163E00ADE7.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/DE468499-E88F-E811-B718-02163E01A048.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/DE1A7EE1-EB8F-E811-B9F3-FA163E05C2A9.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/DAB9F85B-F78F-E811-B252-FA163E91365C.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/DA6C81BC-F18F-E811-9F87-FA163E75056D.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/DA33D911-0590-E811-A6F0-FA163E4BB6D2.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/D80BFA1A-EC8F-E811-9910-FA163EE3AE81.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/D6BC68AC-F18F-E811-9A27-FA163E4A6585.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/D4A894E4-EE8F-E811-AD23-02163E01770B.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/D2FF0403-F28F-E811-8B6E-FA163E14986C.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/D2F8D0F8-EB8F-E811-95BD-FA163E046877.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/D07CC71B-EC8F-E811-AFD1-FA163EA27FD6.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/CE01711C-EC8F-E811-A0CD-FA163E199BE5.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/CCC83D52-F78F-E811-95F2-FA163ED2F018.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/C8F9583E-F58F-E811-B674-02163E016512.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/C846982C-0290-E811-9721-FA163E355512.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/C6A2840A-F78F-E811-BE1A-FA163E564EBD.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/C4E93806-E98F-E811-8252-FA163E49B7CC.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/C076856E-E28F-E811-8BC5-FA163EF033E8.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/BE8397B0-F18F-E811-87C0-FA163EDE4F34.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/BE3E27A6-E88F-E811-B64F-FA163E43D7B9.root',
##'/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/065/00000/BC890B92-E88F-E811-8DEB-FA163E9681D4.root'
])
#lines=[]
#for line in open('filelist_Official.txt'):
#    lines.append('%s' %(line.rstrip('\n')))
#readFiles.extend(lines)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
inputFiles = readFiles

outputFile = "electronTreeZEE_official_xtalInfo.root"
process.source = cms.Source ("PoolSource", fileNames = inputFiles )  
events= cms.untracked.VEventRange()
for line in open('eventList.txt'):
    events.append(line.rstrip("\n"))
process.source.eventsToProcess = cms.untracked.VEventRange(events)
    
                 
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


