from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = 'Charmonium2017C_mini'
# request name is the name of the folder where crab log is saved

config.General.workArea = 'crab3_projects'
config.General.transferOutputs = True

config.section_('JobType')
# set here the path of your working area
config.JobType.psetName = 'cfg_recoWmhis.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['his.root']
config.JobType.allowUndistributedCMSSW = True

config.section_('Data')

# dataset to process, run range and data quality JSON file
config.Data.inputDataset = '/Charmonium/Run2017C-PromptReco-v2/MINIAOD'
#config.Data.runRange = '299368-299368'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/PromptReco/Cert_294927-301141_13TeV_PromptReco_Collisions17_JSON.txt'

#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 1
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 15

# set here the path of a storage area you can write to
config.Data.outLFNDirBase = '/store/user/paolocms'
config.Data.publication = False

##############

#config.Data.ignoreLocality = True

### patch for /Charmonium/Run2017C-PromptReco-v2/MINIAOD stored at T1_FR_CCIN2P3_Disk
config.Data.ignoreLocality = True

config.section_("Site")

# set here a storage site you can write to
config.Site.storageSite = 'T2_IT_Legnaro'

### patch for /Charmonium/Run2017C-PromptReco-v2/MINIAOD stored at T1_FR_CCIN2P3_Disk
config.Site.whitelist = ['T2_FR_*']

