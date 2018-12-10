from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'Diamond_crab2_test'
#config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'test.py'
#config.JobType.disableAutomaticOutputCollection = False
config.JobType.outputFiles = ['output.root']

#config.Data.inputDataset = '/SingleMuon/Run2017E-PromptReco-v1/RECO'
#config.Data.inputDataset = '/SingleElectron/Run2017E-PromptReco-v1/RECO'
config.Data.inputDataset = '/ZeroBias/Run2018C-PromptReco-v1/AOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 5000
config.Data.runRange = '319337'

config.Data.outLFNDirBase = '/afs/cern.ch/work/e/ebossini/Analysis/CMSSW_10_2_0_pre6/src/CrabOut'
config.Site.storageSite = 'T2_CH_CERN'