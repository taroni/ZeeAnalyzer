

import FWCore.ParameterSet.Config as cms
process = cms.Process("TestElectrons")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Prompt_v9', '')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
inputFiles=cms.untracked.vstring()
lines=[]
for line in open('newReRecoFiles.txt'):
    lines.append('file:%s' %(line.rstrip('\n')))
    
inputFiles.extend(lines)


outputFile = "electronPlots_allZee_test.root"
process.source = cms.Source ("PoolSource", fileNames = inputFiles )                             

process.ntupler = cms.EDAnalyzer(
    'ElectronPlots',
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
process.p = cms.Path(process.ntupler)


