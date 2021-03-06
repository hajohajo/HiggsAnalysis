import FWCore.ParameterSet.Config as cms
import HiggsAnalysis.MiniAOD2TTree.tools.git as git #HiggsAnalysis.HeavyChHiggsToTauNu.tools.git as git

process = cms.Process("TTreeDump")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:ttjets_for_testing.root' # tmp
#	'/store/mc/RunIISpring15DR74/QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/2A98D4CF-F9FE-E411-AB1A-047D7BD6DD44.root',
#        'file:miniAOD-prod_PAT_TT_RECO_721_11112014.root'
#	'file:miniAOD-prod_PAT_TT_RECO_740p1_02122014.root'
#	'file:miniAOD-prod_PAT_Hp200_RECO_740p1_09122014.root'
#	'file:miniAOD-prod_PAT_SingleMu_Run2012D_v1_RECO.root'
#	'file:PYTHIA6_Tauola_TTbar_H160_taunu_13TeV_cff_py_GEN.root'
    )
)

# Set up electron ID (VID framework)
# https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
# define which IDs we want to produce and add them to the VID producer
for idmod in ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_PHYS14_PU20bx25_nonTrig_V1_cff']:
    setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)

process.dump = cms.EDFilter('MiniAOD2TTreeFilter',
    OutputFileName = cms.string("miniaod2tree.root"),
    CodeVersion = cms.string(git.getCommitId()),
    DataVersion = cms.string("74Xmc"),
    CMEnergy = cms.int32(13),
    EventInfo = cms.PSet(
	PileupSummaryInfoSrc = cms.InputTag("addPileupInfo"),
#	LHESrc = cms.InputTag(""),
	OfflinePrimaryVertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    ),
    Trigger = cms.PSet(
	TriggerResults = cms.InputTag("TriggerResults::HLT"),
	TriggerBits = cms.vstring(
	    "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v1",
#	    "HLT_IsoMu24_IterTrk02_v1"
        ),
	L1Extra = cms.InputTag("l1extraParticles::MET"),
	TriggerObjects = cms.InputTag("selectedPatTrigger"),
	filter = cms.untracked.bool(False)
    ),
    METNoiseFilter = cms.PSet(
        triggerResults = cms.InputTag("TriggerResults","PAT"),
        printTriggerResultsList = cms.untracked.bool(False),
        CSCTightHaloFilter = cms.string(""),
        goodVerticesFilter = cms.string(""),
        EEBadScFilter = cms.string("")
    ),
    Taus = cms.VPSet(
        cms.PSet(
            branchname = cms.untracked.string("Taus"),
            src = cms.InputTag("slimmedTaus"),
            discriminators = cms.vstring(
                #"againstElectronLoose",
                "againstElectronLooseMVA5",
                "againstElectronMVA5category",
                #"againstElectronMVA5raw",
                #"againstElectronMedium",
                "againstElectronMediumMVA5",
                #"againstElectronTight",
                "againstElectronTightMVA5",
                "againstElectronVLooseMVA5",
                "againstElectronVTightMVA5",
                #"againstMuonLoose",
                #"againstMuonLoose2",
                "againstMuonLoose3",
                #"againstMuonLooseMVA",
                #"againstMuonMVAraw",
                #"againstMuonMedium",
                #"againstMuonMedium2",
                #"againstMuonMediumMVA",
                #"againstMuonTight",
                #"againstMuonTight2",
                "againstMuonTight3",
                #"againstMuonTightMVA",
                "byCombinedIsolationDeltaBetaCorrRaw3Hits",
                "byIsolationMVA3newDMwLTraw",
                "byIsolationMVA3newDMwoLTraw",
                "byIsolationMVA3oldDMwLTraw",
                "byIsolationMVA3oldDMwoLTraw",
                "byLooseCombinedIsolationDeltaBetaCorr3Hits",
                "byLooseIsolationMVA3newDMwLT",
                "byLooseIsolationMVA3newDMwoLT",
                "byLooseIsolationMVA3oldDMwLT",
                "byLooseIsolationMVA3oldDMwoLT",
                "byMediumCombinedIsolationDeltaBetaCorr3Hits",
                "byMediumIsolationMVA3newDMwLT",
                "byMediumIsolationMVA3newDMwoLT",
                "byMediumIsolationMVA3oldDMwLT",
                "byMediumIsolationMVA3oldDMwoLT",
                "byTightCombinedIsolationDeltaBetaCorr3Hits",
                "byTightIsolationMVA3newDMwLT",
                "byTightIsolationMVA3newDMwoLT",
                "byTightIsolationMVA3oldDMwLT",
                "byTightIsolationMVA3oldDMwoLT",
                "byVLooseIsolationMVA3newDMwLT",
                "byVLooseIsolationMVA3newDMwoLT",
                "byVLooseIsolationMVA3oldDMwLT",
                "byVLooseIsolationMVA3oldDMwoLT",
                "byVTightIsolationMVA3newDMwLT",
                "byVTightIsolationMVA3newDMwoLT",
                "byVTightIsolationMVA3oldDMwLT",
                "byVTightIsolationMVA3oldDMwoLT",
                "byVVTightIsolationMVA3newDMwLT",
                "byVVTightIsolationMVA3newDMwoLT",
                "byVVTightIsolationMVA3oldDMwLT",
                "byVVTightIsolationMVA3oldDMwoLT",
                "chargedIsoPtSum",
                "decayModeFinding",
                "decayModeFindingNewDMs", # only for analysis involving high pT taus
                "neutralIsoPtSum",
                "puCorrPtSum"
	    ),
            filter = cms.untracked.bool(False),
            TESvariation = cms.untracked.double(0.03),
            TESvariationExtreme = cms.untracked.double(0.10)
        )
    ),
    Electrons = cms.VPSet(
        cms.PSet(
            branchname = cms.untracked.string("Electrons"),
            src = cms.InputTag("slimmedElectrons"),
            rhoSource = cms.InputTag("fixedGridRhoFastjetAll"), # for PU mitigation in isolation
            IDprefix = cms.vstring("egmGsfElectronIDs"),
            discriminators = cms.vstring("mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp80",
                                         "mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp90")
        )
    ),
    Muons = cms.VPSet(   
        cms.PSet(
            branchname = cms.untracked.string("Muons"),   
            src = cms.InputTag("slimmedMuons"),    
            discriminators = cms.vstring() 
        )   
    ),
    Jets = cms.VPSet(      
        cms.PSet(
            branchname = cms.untracked.string("Jets"),       
            src = cms.InputTag("slimmedJets"),      
            discriminators = cms.vstring(
                "jetBProbabilityBJetTags",
                "jetProbabilityBJetTags",
                "trackCountingHighPurBJetTags", 
                "trackCountingHighEffBJetTags",
                "simpleSecondaryVertexHighEffBJetTags",
                "simpleSecondaryVertexHighPurBJetTags",
                "combinedSecondaryVertexBJetTags",
                "combinedInclusiveSecondaryVertexBJetTags",
                "combinedInclusiveSecondaryVertexV2BJetTags",
            ),
	    userFloats = cms.vstring(
		"pileupJetId:fullDiscriminant"
	    ),
        )
    ),
    METs = cms.VPSet(
        cms.PSet(
            branchname = cms.untracked.string("MET_Type1"),
            src = cms.InputTag("slimmedMETs")
        ),
    )
)

# module execution
process.runEDFilter = cms.Path(process.egmGsfElectronIDSequence*process.dump)

#process.output = cms.OutputModule("PoolOutputModule",
    #outputCommands = cms.untracked.vstring(
        #"keep *",
    #),
    #fileName = cms.untracked.string("tmp.root")
#)
#process.out_step = cms.EndPath(process.output)
