import FWCore.ParameterSet.Config as cms
from HiggsAnalysis.HeavyChHiggsToTauNu.HChOptions import getOptionsDataVersion

dataVersion = "44XmcS6"

options, dataVersion = getOptionsDataVersion(dataVersion, useDefaultSignalTrigger=False)

process = cms.Process("TauEmbeddingAnalysis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string(dataVersion.getGlobalTag())

process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        "/store/group/local/HiggsChToTauNuFullyHadronic/pattuples/CMSSW_4_4_X/TTJets_TuneZ2_Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/Fall11_PU_S6_START44_V9B_v1_AODSIM_tauembedding_gentauskim_v44_5/9ecb3a23e436fc2ffd8a803eac2a3a15/pattuple_1012_1_LSv.root",
    ),
)

process.load("HiggsAnalysis.HeavyChHiggsToTauNu.HChCommon_cfi")
# Fragment to run PAT on the fly if requested from command line
from HiggsAnalysis.HeavyChHiggsToTauNu.HChPatTuple import addPatOnTheFly
patArgs = {"doPatTrigger": False,
#           "doPatTaus": False,
#           "doHChTauDiscriminators": False,
           "doPatElectronID": True,
           "doTauHLTMatching": False,
           }
process.commonSequence, additionalCounters = addPatOnTheFly(process, options, dataVersion, patArgs=patArgs,
                                                            doHBHENoiseFilter=False,
                                                            )
import HiggsAnalysis.HeavyChHiggsToTauNu.AnalysisConfiguration as AnalysisConfiguration
dataEras = [
    "Run2011AB",
    "Run2011A",
    "Run2011B",
]
puWeights = AnalysisConfiguration.addPuWeightProducers(dataVersion, process, process.commonSequence, dataEras)

# Add GenTau skim counters
import HiggsAnalysis.HeavyChHiggsToTauNu.CustomGenTauSkim as tauSkim
additionalCounters = tauSkim.getCounters() + additionalCounters

# Add configuration information to histograms.root
import HiggsAnalysis.HeavyChHiggsToTauNu.HChTools as HChTools
process.infoPath = HChTools.addConfigInfo(process, options, dataVersion)

# Jet selection
import HiggsAnalysis.HeavyChHiggsToTauNu.tauEmbedding.muonSelectionPF as muonSelection
muonSelection.addMuonSelectionForEmbedding(process)
process.commonSequence += (
    process.goodJets +
    process.goodJetFilter +
    process.muonSelectionJets
)
additionalCounters.append("muonSelectionJets")

# Configuration
import HiggsAnalysis.HeavyChHiggsToTauNu.HChSignalAnalysisParameters_cff as param
analyzer = cms.EDAnalyzer("HPlusEmbeddingDebugTauAnalyzer",
    jetSrc = cms.untracked.InputTag("goodJets"),
    genSrc = cms.untracked.InputTag("genParticles"),

    tauPtCut = cms.untracked.double(40),
    tauEtaCut = cms.untracked.double(2.1),

    eventCounter = param.eventCounter.clone(),
    histogramAmbientLevel = cms.untracked.string("Informative"),
)

HChTools.addAnalysis(process, "debugAnalyzer", analyzer,
                     preSequence=process.commonSequence,
                     additionalCounters=additionalCounters)
process.debugAnalyzer.eventCounter.printMainCounter = True


f = open("configDumpEmbeddingDebugTauAnalysis.py", "w")
f.write(process.dumpPython())
f.close()
