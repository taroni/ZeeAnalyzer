# from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'ZeeMC_sum8gt0'

config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'electronTree.py'

config.section_("Data")
config.Data.inputDataset = '/DYToEE_M-50_NNPDF31_TuneCP5_13TeV-powheg-pythia8/taroni-DY_sum8gt0_pierreTag-240b766fa4a3ff5cac77f0390e9bd2e8/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'Automatic'
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.publication = True
config.Data.outputDatasetTag =  'ZeeMC_sum8gt0'

config.Site.storageSite = 'T2_CH_CERN'
config.Site.whitelist = ['T3_US_NotreDame']
#all the configuration parameters https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile
#all crab commands https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3Commands
