import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('inputFile',
                '/store/mc/RunIISummer19UL18MiniAOD/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/280000/BF4D011D-4171-B749-8F28-EA1E912E7131.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "input file to process"
                 )
options.register('maxEvts',
                 -1,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 'Number of events to process'
                )
                  
options.register('globalTag',
                 "",
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 'Global Tag'
                )
                
options.register('wantSummary',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                'Want summary report'
                )

options.register('isMC',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 'Running on MC' 
                )


### Events to process: 'maxEvents' is already registered by the framework
options.setDefault('maxEvents', 10)

options.parseArguments()


process = cms.Process("Demo")
#message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# set input to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvts) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFile),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
    )
    
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = options.globalTag

#process.load('PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi')

process.demo = cms.EDAnalyzer('Analysis',
	isMC=cms.untracked.bool(options.isMC),
	Cand_Pt_THR=cms.double(1.0),
	Jet_Pt_THR=cms.double(10.0),
   )

process.TFileService = cms.Service("TFileService",
   fileName = cms.string("Output.root")
   )

process.p=cms.Path(process.demo)
