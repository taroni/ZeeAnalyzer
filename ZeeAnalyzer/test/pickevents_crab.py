
## Edited By Raman Khurana
##
## CRAB documentation : https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCrab
##
## CRAB 3 parameters : https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile#CRAB_configuration_parameters
##
## Once you are happy with this file, please run
## crab submit

## In CRAB3 the configuration file is in Python language. It consists of creating a Configuration object imported from the WMCore library: 

from WMCore.Configuration import Configuration
config = Configuration()

##  Once the Configuration object is created, it is possible to add new sections into it with corresponding parameters
config.section_("General")
config.General.requestName = 'pickEvents'
config.General.workArea = 'crab_pickevents_20190320_165826'


config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw/CMSSW_10_2_3/src/PhysicsTools/Utilities/configuration/copyPickMerge_cfg.py'
config.JobType.pyCfgParams = ['eventsToProcess_load=pickevents_runEvents.txt', 'outputFile=pickevents.root']

config.section_("Data")
config.Data.inputDataset = '/EGamma/Run2018C-ZElectron-PromptReco-v3/RAW-RECO'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 5
config.Data.lumiMask = 'pickevents.json'
#config.Data.publication = True
#config.Data.publishDbsUrl = 'phys03'
#config.Data.publishDataName = 'CRAB3_CSA_DYJets'
#config.JobType.allowNonProductionCMSSW=True

config.section_("Site")
## Change site name accordingly
config.Site.storageSite = 'T3_US_NotreDame'
config.Data.ignoreLocality = True
config.Site.ignoreGlobalBlacklist = True
config.Site.whitelist = ['T3_US_FNALLPC', 'T3_US_Kansas', 'T3_RU_FIAN', 'T3_US_MIT',  'T3_US_UCD', 'T3_CO_Uniandes', 'T3_US_NotreDame', 'T1_US_FNAL', 'T2_IT_Rome', 'T3_IN_PUHEP', 'T2_CH_CERN_HLT', 'T2_AT_Vienna', 'T3_IN_TIFRCloud', 'T3_GR_IASA', 'T3_CN_PKU',  'T2_RU_ITEP', 'T3_US_JHU', 'T3_BY_NCPHEP', 'T3_US_FSU', 'T3_KR_UOS', 'T3_CH_PSI']

