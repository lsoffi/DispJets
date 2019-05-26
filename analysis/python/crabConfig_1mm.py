from CRABClient.UserUtilities import config
config = config()

## Name of output directory ##
config.General.requestName = 'DispJets_m50_ct1mm'
config.General.workArea    = 'crab_projects'

config.JobType.allowUndistributedCMSSW = True
## Input analyzer pset ## 
config.JobType.pluginName  = 'analysis'
config.JobType.psetName    = 'crabConfigTemplate.py'
#config.JobType.pyCfgParams = ['globalTag=92X_dataRun2_Prompt_v7','useOOTPhotons=True']
#config.JobType.inputFiles  = ['HLTpaths.txt','HLTfilters.txt']

## Input Data ##
#config.Data.inputDataset   = ''
config.Data.userInputFiles = open('inputfiles_1mm.txt').readlines()
config.Data.unitsPerJob    = 1
config.Data.splitting      = 'FileBased' 

## Where to run ##
config.Site.whitelist     = ['T1_US_FNAL']

## Output Data ##
config.Data.outputPrimaryDataset = 'XXQQQQ_m50'
config.Data.publication   = False
config.Site.storageSite   = 'T2_CH_CERN'
config.Data.outLFNDirBase = '/store/group/phys_egamma/soffi/XXQQQQ-GEN-SIM/DisplJetAnalyzer'
