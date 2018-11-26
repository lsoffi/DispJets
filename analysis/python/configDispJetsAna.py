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
	# xrootd files FNAL:            'root://cmsxrootd.fnal.gov//PATH/file.root'
        # xrootd files GLOBAL REDIRECT: 'root://cms-xrd-global.cern.ch//PATH/file.root'
	#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer15GS/QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6/GEN-SIM/MCRUN2_71_V1-v1/40000/448C5AAE-F16B-E511-AA32-3417EBE51CDF.root'
	'root://cmsxrootd.fnal.gov//store/user/mzientek/DisplacedJetSamples_GENSIM/XXQQQQ_m50_ctau100mm/ggF-H-S1S2-Sdecay-qqnunu_UFO_ctau-10p0_384049_10444854_GENSIM.root'
	#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer15GS/XXTo4J_M-50_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/GEN-SIM/MCRUN2_71_V1-v1/10000/F80A1959-95C2-E511-8FA6-001EC9ADE758.root'
	#'/store/mc/RunIISummer15GS/XXTo4J_M-50_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/GEN-SIM/MCRUN2_71_V1-v1/10000/F80A1959-95C2-E511-8FA6-001EC9ADE758.root'
        #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer15GS/XXTo4J_M-100_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/GEN-SIM/MCRUN2_71_V1-v1/70000/14117C30-7FC2-E511-A4D6-00266CF2679C.root',
        #'/store/mc/RunIISummer15GS/XXTo4J_M-100_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/GEN-SIM/MCRUN2_71_V1-v1/70000/14117C30-7FC2-E511-A4D6-00266CF2679C.root',
	#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer15GS/XXTo4J_M-50_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/GEN-SIM/MCRUN2_71_V1-v1/10000/F80A1959-95C2-E511-8FA6-001EC9ADE758.root'
	#'root://cmsxrootd.fnal.gov//store/user/mzientek/DisplacedJetSamples_GENSIM/XXQQQQ_m50_ctau10mm/ggF-H-S1S2-Sdecay-qqnunu_UFO_ctau-1p0_713309_17563488_GENSIM.root'
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
				tracks		= cms.untracked.InputTag("g4SimHits", "", "SIM"),
)

if opts.debug: process.p = cms.Path( process.content*process.dispjets )
else:          process.p = cms.Path( process.dispjets )
