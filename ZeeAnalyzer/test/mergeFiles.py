import FWCore.ParameterSet.Config as cms
process = cms.Process("mergingFilesgt15")

process.load("FWCore.MessageService.MessageLogger_cfi")

lines=[]
readFiles=cms.untracked.vstring()

for line in open('filesZeesum8gt0.txt'):
    lines.append('file:%s' %(line.rstrip('\n')))

readFiles.extend(lines)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
inputFiles = readFiles
outputFile = "/tmp/taroni/step3_Zee_sum8gt0.root"
process.source = cms.Source ("PoolSource", fileNames = inputFiles )                             

process.Out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string (outputFile)
)
process.end = cms.EndPath(process.Out)
