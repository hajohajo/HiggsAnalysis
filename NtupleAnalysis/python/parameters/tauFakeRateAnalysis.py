#!/usr/bin/env python

from HiggsAnalysis.NtupleAnalysis.main import PSet
import HiggsAnalysis.NtupleAnalysis.parameters.scaleFactors as scaleFactors
import HiggsAnalysis.NtupleAnalysis.parameters.jsonReader as jsonReader

#================================================================================================
# General parameters
#================================================================================================
verbose               = True
histogramAmbientLevel = "Informative"  # ("Systematics", "Vital", "Informative", "Debug")

#================================================================================================
# Trigger [scanned in range _v1--_v100 (=>remove the '_v' suffix)]
#================================================================================================
trigger = PSet(
    triggerOR = [
        # SingleMuon Primary Dataset (PD)
        "HLT_IsoMu24",
        "HLT_IsoTkMu24",
        # # SingleElectron Primary Dataset (PD)
        # "HLT_Ele27_WPTight_Gsf",
        # # DoubleMuon Primary Dataset (PD)
        # "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
        # "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
        # "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
        # "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
        # # DoubleEG Primary Dataset (PD)
        # "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        # "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        # "HLT_DoubleEle33_CaloIdL",
        # # MuonEG Primary Dataset (PD)
        # "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
        # "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        # "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
        # "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        # "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
        # "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        ],
    triggerOR2 = [],
    )

#================================================================================================
# MET filter
#================================================================================================
metFilter = PSet(
    discriminators = [
        "Flag_HBHENoiseFilter",
        "Flag_HBHENoiseIsoFilter",
        "Flag_EcalDeadCellTriggerPrimitiveFilter",
        "Flag_eeBadScFilter",
        "Flag_goodVertices",
        #"Flag_globalTightHalo2016Filter",
        "Flag_globalSuperTightHalo2016Filter",
        "badPFMuonFilter",
        "badChargedCandidateFilter"]
    )

#================================================================================================
# Electron Selection
#================================================================================================
electronSelection = PSet(
    applyTriggerMatching = False,
    triggerMatchingCone  = 0.1,   # DeltaR for matching offline tau with trigger tau
    electronPtCut        = 10.0,    # [default: 10.0]
    electronEtaCut       = 2.4,     # [default: 2.4]
    electronIDType       = "MVA",   # [default: "MVA] ("default", "MVA")
    electronID           = "cutBasedElectronID_Spring15_25ns_V1_standalone_veto",
    electronMVA          = "ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values",
    electronMVACut       = "Loose", # [default: "Loose"]
    electronIsolation    = "veto",  # [default: "veto"] ("veto", "tight")
    electronIsolType     = "mini",  # [default: "mini"] ("mini", "default")
    )

#================================================================================================
# Muon Selection
#================================================================================================
muonSelection = PSet(
    applyTriggerMatching = True ,       # HLTMuon_Pt Missing from tree!
    triggerMatchingCone  = 0.1,         # DeltaR for matching offline tau with trigger tau
    muonPtCut            = 26.0,        # [default: 10.0]
    muonEtaCut           = 2.4,         # [default: 2.4]
    muonID               = "muIDTight", # [default: "muIDTight"] ("muIDLoose", "muIDMedium", "muIDTight")
    muonIsolation        = "tight",     # [default: "tight"]
    muonIsolType         = "default",   # [default: "default"]
)

#================================================================================================
# Tau selection (sync with HToTauNu analysis)
#================================================================================================
tauSelection = PSet(
    applyTriggerMatching = False, # [default: False]
    triggerMatchingCone  =   0.1, # [default: False]
    tauPtCut             =  20.0, # [default: 20.0]
    tauEtaCut            =   2.3, # [default: 2.3]
    tauLdgTrkPtCut       =   0.0, # [default: 0.0]
    prongs               =  -1,   # [default: -1] (1, 2, 3, 12, 13, 23, 123 or -1 (all))
    rtau                 =   0.0, # [default: 0.0] (to disable set to 0.0)
    againstElectronDiscr = "againstElectronTightMVA6",
    againstMuonDiscr     = "againstMuonLoose3",
    isolationDiscr       = "byVLooseIsolationMVArun2v1DBoldDMwLT",
    #isolationDiscr       = "byLooseIsolationMVArun2v1DBoldDMwLT",
    #isolationDiscr       = "byMediumIsolationMVArun2v1DBoldDMwLT", # Default
    #isolationDiscr       = "byTightIsolationMVArun2v1DBoldDMwLT",
    # isolationDiscr       = "byLooseCombinedIsolationDeltaBetaCorr3Hits",
    # isolationDiscr       = "byMediumCombinedIsolationDeltaBetaCorr3Hits",
    )


# tau mis-id fake rate measurement
looseTauSelection = PSet(
  applyTriggerMatching = False, # no effect now
   triggerMatchingCone = 0.1,   # DeltaR for matching offline tau with trigger tau
              tauPtCut = 20.0,
             tauEtaCut = 2.3,
        tauLdgTrkPtCut = 1.0,
                prongs = -1,    # options: 1, 2, 3, 12, 13, 23, 123 or -1 (all)
                  rtau = 0.0,   # to disable set to 0.0
  againstElectronDiscr = "againstElectronVLooseMVA6",
      againstMuonDiscr = "againstMuonLoose3", #"againstMuonTight3", #"againstMuonLoose3",
        isolationDiscr = "byVLooseIsolationMVArun2v1DBoldDMwLT", #"byMediumIsolationMVArun2v1DBnewDMwLT",
)

# tau identification scale factors
scaleFactors.assignTauIdentificationSF(tauSelection)
scaleFactors.assignTauIdentificationSF(looseTauSelection)

# tau misidentification scale factor
scaleFactors.assignTauMisidentificationSF(tauSelection, "eToTau", "nominal")
scaleFactors.assignTauMisidentificationSF(tauSelection, "muToTau", "nominal")
scaleFactors.assignTauMisidentificationSF(tauSelection, "jetToTau",  "nominal")

scaleFactors.assignTauMisidentificationSF(looseTauSelection, "eToTau", "nominal")
scaleFactors.assignTauMisidentificationSF(looseTauSelection, "muToTau", "nominal")
scaleFactors.assignTauMisidentificationSF(looseTauSelection, "jetToTau",  "nominal")
'''
newer version still unsuported in this branch
scaleFactors.assignTauMisidentificationSF(tauSelection, "eToTau", "full", "nominal")
scaleFactors.assignTauMisidentificationSF(tauSelection, "muToTau", "full", "nominal")
scaleFactors.assignTauMisidentificationSF(tauSelection, "jetToTau", "full", "nominal")

scaleFactors.assignTauMisidentificationSF(looseTauSelection, "eToTau", "full", "nominal")
scaleFactors.assignTauMisidentificationSF(looseTauSelection, "muToTau", "full", "nominal")
scaleFactors.assignTauMisidentificationSF(looseTauSelection, "jetToTau", "full", "nominal")
'''

#================================================================================================
# Jet selection
#================================================================================================
jetSelection = PSet(
    jetType                  = "Jets",    # [default: "jets"] ("Jets", "JetsPuppi")
    jetPtCuts                = [30.0],    # [default: [30.0]]
    jetEtaCuts               = [4.7],     # [default: [4.7]]
    numberOfJetsCutValue     = 0,         # [default: 3]
    numberOfJetsCutDirection = ">=",      # [default: ">="] (==, !=, <, <=, >, >=)
    jetIDDiscr               = "IDloose", # [default: "IDloose"] ("IDloose", "IDtight", "IDtightLeptonVeto")
    jetPUIDDiscr             = "",        # [default: ""]
    tauMatchingDeltaR        = 0.4,       # [default: 0.4]
    HTCutValue               = 0.0,       # [default: 500.0]
    HTCutDirection           = ">=",      # [default: ">="]
    JTCutValue               = 0.0,       # [default: 0.0]
    JTCutDirection           = ">=",      # [default: ">="]
    MHTCutValue              = 0.0,       # [default: 0.0]
    MHTCutDirection          = ">=",      # [default: ">="]
)

#================================================================================================
# B-jet selection
#================================================================================================
bjetSelection = PSet(
    triggerMatchingApply      = False,    # [default: False]
    triggerMatchingCone       = 0.1,      # [default: 0.1 ]
    jetPtCuts                 = [30.0],   # [default: [20.]]
    jetEtaCuts                = [2.4],    # [default: [2.4]]
    bjetDiscr                 = "pfCombinedInclusiveSecondaryVertexV2BJetTags", # default
    bjetDiscrWorkingPoint     = "Medium", # [default: "Medium"] ("Medium", "Tight")
    numberOfBJetsCutValue     = 0,        # [default: 1]
    numberOfBJetsCutDirection = ">=",     # [default: ">="] (==, !=, <, <=, >, >=)
)

#=================================================================================================
# Fat jet selection
#=================================================================================================
fatjetVeto = PSet(
    fatjetType                  = "FatJets", # [default: "FatJets"]  
    fatjetPtCuts                = [450.0],   # [default: [450.0] ]
    fatjetEtaCuts               = [2.4],     # [default: [2.4] ]
    fatjetIDDiscr               = "IDloose", # [default: "IDLoose"] ("IDloose", "IDtight", "IDtightLeptonVeto")
    fatjetPUIDDiscr             = "",        # [default: ""]
    topMatchDeltaR              = 0.8,       # [default: 0.8]
    topMatchTypes               = [1],       # [default: 1]   (kJJB=1, kJJ=2, kJB=3, kJJBorJJ=4, kJJBorJB=5, kJJorJB=6, kAll=7, any = -1)
    numberOfFatJetsCutValue     = 0,         # [default: 0]
    numberOfFatJetsCutDirection = ">=",      # [default: "=="] (TO DISABLE: >=0)
)

#================================================================================================
# MET selection
#================================================================================================
metSelection = PSet(
    METCutValue                 = 0.0,         # [default: 30]
    METCutDirection             = ">=",        # [default: ">"] (==, !=, <, <=, >, >=)
    METSignificanceCutValue     = -1000.0,     # [default: -1000.0]
    METSignificanceCutDirection = ">",         # [default: ">"] (==, !=, <, <=, >, >=)
    METType                     = "MET_Type1", # [default: "MET_Type1"] ("MET_Type1_NoHF", "MET_Puppi", "GenMET", "L1MET", "HLTMET", "CaloMET"))
    applyPhiCorrections         = False        #  [default: False]
    )

#================================================================================================
# Scale Factors (SFs)
#================================================================================================
if bjetSelection.bjetDiscr == "pfCombinedInclusiveSecondaryVertexV2BJetTags":
    scaleFactors.setupBtagSFInformation(btagPset               = bjetSelection, 
                                        btagPayloadFilename    = "CSVv2.csv",
                                        #btagEfficiencyFilename = "btageff_hybrid_HToTB.json",
                                        btagEfficiencyFilename = "btageff_HToTB.json",
                                        direction              = "nominal")
elif bjetSelection.bjetDiscr == "pfCombinedMVAV2BJetTags":
    scaleFactors.setupBtagSFInformation(btagPset               = bjetSelection, 
                                        btagPayloadFilename    = "cMVAv2_Moriond17_B_H.csv", # use this for MVA b-tagging
                                        btagEfficiencyFilename = "btageff_Hybrid_TT+WJetsHT.json", # use with taunu analysis and WJetsHT samples
                                        direction              = "nominal")
else:
    raise Exception("This should never be reached!")


#================================================================================================
# Common plots options
#================================================================================================
commonPlotsOptions = PSet(
    histogramSplitting         = [],    # Splitting of histograms as function of one or more parameters
    enableGenuineBHistograms   = False,
    enablePUDependencyPlots    = True,  # Enable/Disable some debug-level plots
    # Bin settings (final bin setting done in datacardGenerator, there also variable bin width is supported)
    nVerticesBins     = PSet(nBins = 100, axisMin =  0.0, axisMax =  200.0),
    ptBins            = PSet(nBins =  50, axisMin =  0.0, axisMax =  500.0),
    etaBins           = PSet(nBins =  50, axisMin = -5.0, axisMax =    5.0),
    phiBins           = PSet(nBins =  64, axisMin = -3.2, axisMax =    3.2),
    deltaEtaBins      = PSet(nBins = 100, axisMin =  0.0, axisMax =   10.0),
    #deltaPhiBins      = PSet(nBins =  32, axisMin =  0.0, axisMax =    3.2),
    deltaPhiBins      = PSet(nBins =  36, axisMin =  0.0, axisMax =  180.0), # 5 degress
    deltaRBins        = PSet(nBins = 100, axisMin =  0.0, axisMax =   10.0),
    rtauBins          = PSet(nBins =  55, axisMin =  0.0, axisMax =    1.1), # HToTauNu
    njetsBins         = PSet(nBins =  18, axisMin =  0.0, axisMax =   18.0),
    metBins           = PSet(nBins = 160, axisMin =  0.0, axisMax =  800.0), # 5 GeV bin width
    htBins            = PSet(nBins = 400, axisMin =  0.0, axisMax = 2000.0), # 5 GeV bin width 
    bjetDiscrBins     = PSet(nBins = 120, axisMin =  0.0, axisMax =    1.2),
    angularCuts1DBins = PSet(nBins =  52, axisMin =  0.0, axisMax =  260.0), 
    topMassBins       = PSet(nBins = 200, axisMin =  0.0, axisMax = 1000.0), #  5 GeV bin width 
    wMassBins         = PSet(nBins = 200, axisMin =  0.0, axisMax = 1000.0), #  5 GeV bin width 
    mtBins            = PSet(nBins = 200, axisMin =  0.0, axisMax = 2000.0), # 10 GeV bin width
    invMassBins       = PSet(nBins = 600, axisMin =  0.0, axisMax = 3000.0), #  5 GeV bin width
)

#================================================================================================
# Build all selections group
#================================================================================================
allSelections = PSet(
    Verbose               = verbose,
    Trigger               = trigger,
    METFilter             = metFilter,
    ElectronSelection     = electronSelection,
    MuonSelection         = muonSelection,
    TauSelection          = tauSelection,
    JetSelection          = jetSelection,
    BJetSelection         = bjetSelection,
    METSelection          = metSelection,
    # FatJetSelection       = fatjetVeto,
    CommonPlots           = commonPlotsOptions,
    HistogramAmbientLevel = histogramAmbientLevel,
)
