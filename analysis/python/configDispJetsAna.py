import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("dispjets")
process.load("FWCore.MessageService.MessageLogger_cfi")

# set any input variables
opts = VarParsing.VarParsing('analysis')
opts.register('debug',                                      # option name
              False,                                        # default value
              VarParsing.VarParsing.multiplicity.singleton, # singleton or list
              VarParsing.VarParsing.varType.bool,           # type: string, int, float, bool
              "Print out debug info")                       # description
opts.register('sampleID',
              100,                                          # default value
              VarParsing.VarParsing.multiplicity.singleton, # singleton or list
              VarParsing.VarParsing.varType.int,            # type: string, int, float, bool
              "SampleID [1-99] Bkg, [100+] Sig")            # description
opts.parseArguments()

# max events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

# get input file
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # local file needs format 'file:mylocalfile.root'
        # nonlocal file just needs full DAS path
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer15GS/XXTo4J_M-100_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/GEN-SIM/MCRUN2_71_V1-v1/70000/14117C30-7FC2-E511-A4D6-00266CF2679C.root',
	#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer15GS/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/GEN-SIM/MCRUN2_71_V1-v1/60003/FECD584C-1560-E511-BCA4-02163E015537.root',
	#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17wmLHEGS/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/GEN-SIM/92X_upgrade2017_realistic_v11-v2/00016/B817B949-BBAB-E711-AF6A-24BE05C6B701.root',
    )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple_dispjets.root")
)

# check event content
process.content  = cms.EDAnalyzer('EventContentAnalyzer')

# run displaced jet analyzer
process.dispjets = cms.EDAnalyzer('DisplacedJetsAnalyzer',
                                verbose         = cms.untracked.bool(opts.debug),
				sampleID	= cms.untracked.int32(opts.sampleID),
				generatorInfo	= cms.InputTag("generator"),
				genjets		= cms.untracked.InputTag("ak4GenJets", "", "SIM"),
                                genparticles    = cms.untracked.InputTag("genParticles", "", "SIM"),
				pileupInfo	= cms.untracked.InputTag("slimmedAddPileupInfo"),
				vertices	= cms.untracked.InputTag("g4SimHits", "", "SIM"),
)

if opts.debug: process.p = cms.Path( process.content*process.dispjets )
else:          process.p = cms.Path( process.dispjets )
