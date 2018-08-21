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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# get input file
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( opts.inputFiles
        # local file needs format 'file:mylocalfile.root'
        # nonlocal file just needs full DAS path
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
