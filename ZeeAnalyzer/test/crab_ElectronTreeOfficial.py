# from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'ZeeOfficial_noRecov'

config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'electronTree_Zee_official.py'

config.section_("Data")
#config.Data.inputDataset = '/DYToEE_M-50_NNPDF31_TuneCP5_13TeV-powheg-pythia8/taroni-DY_sum8gt0_pierreTag-240b766fa4a3ff5cac77f0390e9bd2e8/USER'
#config.Data.inputDBS = 'phys03'
config.Data.inputDataset ='/EGamma/Run2018C-ZElectron-PromptReco-v3/RAW-RECO'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
config.Data.runRange = '319840-319915'
config.Data.splitting = 'Automatic'
config.Data.totalUnits = -1
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.publication = True
config.Data.outputDatasetTag =  'ZeeOfficial_noRecov'

#config.Site.storageSite = 'T2_CH_CERN'
config.Site.storageSite = 'T3_US_NotreDame'
config.Data.ignoreLocality = True
config.Site.ignoreGlobalBlacklist = True

config.Site.whitelist = ['T3_US_FNALLPC', 'T3_US_Kansas', 'T3_RU_FIAN', 'T3_US_MIT',  'T3_US_UCD', 'T3_CO_Uniandes', 'T3_US_NotreDame', 'T1_US_FNAL', 'T2_IT_Rome', 'T3_IN_PUHEP', 'T2_CH_CERN_HLT', 'T2_AT_Vienna', 'T3_IN_TIFRCloud', 'T3_GR_IASA', 'T3_CN_PKU',  'T2_RU_ITEP', 'T3_US_JHU', 'T3_BY_NCPHEP', 'T3_US_FSU', 'T3_KR_UOS', 'T3_CH_PSI']
#all the configuration parameters https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile
#all crab commands https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3Commands
