import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('inputFile',
                 "",
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
#options.setDefault('maxEvts', 10)

options.parseArguments()

isTTbar=False
#ttbar_MC_input="root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/TTJets_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/00000/000EEBFE-D1D0-E511-B584-0025905964BA.root"
ttbar_MC_input="root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/002253C9-DFB8-E511-8B0A-001A648F1C42.root"
qcd_MC_input="root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/0402F447-08BA-E511-B4FF-008CFA002184.root"

outfileName=cms.string("")

if options.isMC:
	options.globalTag="76X_mcRun2_asymptotic_RunIIFall15DR76_v1"	
	options.inputFile= ttbar_MC_input if isTTbar else qcd_MC_input
	outfileName=cms.string("OutFile_MC.root")

else:
	options.globalTag="76X_dataRun2_v15"
	options.inputFile="root://eospublic.cern.ch//eos/opendata/cms/Run2015D/JetHT/MINIAOD/16Dec2015-v1/00000/301A497D-70B0-E511-9630-002590D0AFA8.root"
	outfileName=cms.string("OutFile_data.root")



process = cms.Process("Demo")
#message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 100

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

process.TFileService = cms.Service("TFileService",
   fileName = outfileName
   )

process.demo = cms.EDAnalyzer('Analysis',
	isMC=cms.untracked.bool(options.isMC),
	Cand_Pt_THR=cms.double(1.0),
	Jet_Pt_THR=cms.double(10.0),
   )

process.p=cms.Path(process.demo)
