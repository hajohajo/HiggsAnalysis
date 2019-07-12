#!/usr/bin/env python

from HiggsAnalysis.NtupleAnalysis.main import PSet
import HiggsAnalysis.NtupleAnalysis.parameters.scaleFactors as scaleFactors


#====== General parameters
histoLevel = "Debug"  # Options: Systematics, Vital, Informative, Debug

#====== Trigger
trg = PSet(
  # No need to specify version numbers, they are automatically scanned in range 1--100 (remove the '_v' suffix)
  TautriggerEfficiencyJsonName = "tauLegTriggerEfficiency_2016_fit.json",
  METtriggerEfficiencyJsonName = "metLegTriggerEfficiency_2016_MET90_fit.json",
  L1ETM = 80,
  triggerOR = ["HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90"
               ],
  triggerOR2 = [
	        ],
)

#====== MET filter
metFilter = PSet(
  discriminators = [#"hbheNoiseTokenRun2Loose", # Loose is recommended
#                    "hbheIsoNoiseToken", # under scrutiny
                    "Flag_HBHENoiseFilter",
                    "Flag_HBHENoiseIsoFilter",
                    "Flag_EcalDeadCellTriggerPrimitiveFilter",
#                    "Flag_CSCTightHaloFilter",
                    "Flag_eeBadScFilter",
                    "Flag_goodVertices",
                    "Flag_globalTightHalo2016Filter",
                    "badPFMuonFilter",
                    "badChargedCandidateFilter"],
    runOnlyData = ["Flag_eeBadScFilter"]
)

#====== Tau selection
tauSelection = PSet(
  applyTriggerMatching = True,
   triggerMatchingCone = 0.1,   # DeltaR for matching offline tau with trigger tau
              tauPtCut = 50.0,
             tauEtaCut = 2.1,
        tauLdgTrkPtCut = 30.0,
#                prongs = 13,    # options: 1, 2, 3, 12, 13, 23, 123 or -1 (all)
                prongs = 1,    # options: 1, 2, 3, 12, 13, 23, 123 or -1 (all)
                  rtau = 0.75,   # to disable set to 0.0
  againstElectronDiscr = "againstElectronTightMVA6",
#  againstElectronDiscr = "",
      againstMuonDiscr = "againstMuonLoose3",
        isolationDiscr = "byLooseIsolationMVArun2v1DBoldDMwLT", # MVA (default)
#        isolationDiscr = "byLooseCombinedIsolationDeltaBetaCorr3Hits",  # cut-based
)
# tau identification scale factors
scaleFactors.assignTauIdentificationSF(tauSelection)
# tau misidentification scale factors
scaleFactors.assignTauMisidentificationSF(tauSelection, "eToTau", "nominal")
scaleFactors.assignTauMisidentificationSF(tauSelection, "muToTau", "nominal")
scaleFactors.assignTauMisidentificationSF(tauSelection, "jetToTau", "nominal")
# tau trigger SF

scaleFactors.assignTauTriggerSF(tauSelection, "nominal", trg.TautriggerEfficiencyJsonName)

#====== Electron veto
eVeto = PSet(
    electronPtCut = 15.0,
    electronEtaCut = 2.5,
#            electronID = "mvaEleID_PHYS14_PU20bx25_nonTrig_V1_wp90", # highest (wp90) for vetoing (2012: wp95)
    electronID = "cutBasedElectronID_Spring15_25ns_V1_standalone_veto",
    electronIDType    = "MVA",  # options: "default", "MVA"
    electronMVA       = "ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values",
    electronMVACut    = "Loose",
    electronIsolation = "veto", # loosest possible for vetoing ("veto"), "tight" for selecting
    electronIsolType  = "mini", # options: "mini", "default"
)

#====== Muon veto
muVeto = PSet(
    muonPtCut = 10.0,
    muonEtaCut = 2.5,
    muonID = "muIDLoose", # loosest option for vetoing (options: muIDLoose, muIDMedium, muIDTight)
    muonIsolation = "veto", # loosest possible for vetoing ("veto"), "tight" for selecting
    muonIsolType      = "mini",      # options: "mini", "default" 
)

#====== Muon selection (for embedding)
muForEmbedding = PSet(
             muonPtCut = 40.0,
            muonEtaCut = 2.5,
                muonID = "muIDTight", # options: muIDLoose, muIDMedium, muIDTight
         muonIsolation = "tight", # for selecting, not vetoing
)

#====== Jet selection
jetSelection = PSet(
               jetType  = "Jets", # options: Jets (AK4PFCHS), JetsPuppi (AK4Puppi)
              jetPtCuts = [30.0],
             jetEtaCuts = [4.7],
     tauMatchingDeltaR  = 0.4,
  numberOfJetsCutValue  = 3,
  numberOfJetsCutDirection = ">=", # options: ==, !=, <, <=, >, >=
            jetIDDiscr = "IDloose", # options: IDloose, IDtight, IDtightLeptonVeto
          jetPUIDDiscr = "", # does not work at the moment 
            HTCutValue = 0.0,
    HTCutDirection     = ">=",
            JTCutValue = 0.0,
    JTCutDirection     = ">=",
           MHTCutValue = 0.0,
    MHTCutDirection    = ">=",
)
 
#====== Angular cuts / collinear
angularCutsCollinear = PSet(
       nConsideredJets = 3,    # Number of highest-pt jets to consider (excluding jet corresponding to tau)
enableOptimizationPlots = True, # 2D histograms for optimizing angular cuts
        cutValueJet1 = 0.0,   # Cut value in degrees (circular cut)
        cutValueJet2 = 0.0,   # Cut value in degrees (circular cut)
        cutValueJet3 = 0.0,   # Cut value in degrees (circular cut)
        cutValueJet4 = 0.0,   # Cut value in degrees (circular cut)
)
 
#====== B-jet selection
bjetSelection = PSet(
    triggerMatchingApply= False,
    triggerMatchingCone = 0.1,  # DeltaR for matching offline bjet with trigger::TriggerBjet 
              jetPtCuts = [30.0],
             jetEtaCuts = [2.4],
             bjetDiscr  = "pfCombinedInclusiveSecondaryVertexV2BJetTags", # default
#             bjetDiscr  = "pfCombinedMVAV2BJetTags", # use this for MVA b-tagging
 bjetDiscrWorkingPoint  = "Medium", #optimal for CSVv2
# bjetDiscrWorkingPoint  = "Tight", #optimal for CMVAv2
 numberOfBJetsCutValue  = 1,
 numberOfBJetsCutDirection = ">=", # options: ==, !=, <, <=, >, >=
)

scaleFactors.setupBtagSFInformation(btagPset=bjetSelection, 
                                    btagPayloadFilename="CSVv2.csv",
                                    #btagPayloadFilename="cMVAv2_Moriond17_B_H.csv", # use this for MVA b-tagging
                                    #btagEfficiencyFilename="btageff_TTJets.json",
                                    #btagEfficiencyFilename="btageff_WJetsHT.json",
                                    #btagEfficiencyFilename="btageff_hybrid.json",
                                    #btagEfficiencyFilename="btageff_hybrid_HToTB.json",
                                    btagEfficiencyFilename="btageff_Hybrid_TT+WJetsHT.json", # use with taunu analysis and WJetsHT samples
                                    direction="nominal")

#====== MET selection
metSelection = PSet(
           METCutValue = 90.0,
       METCutDirection = ">", # options: ==, !=, <, <=, >, >=
  METSignificanceCutValue = -1000.0,
  METSignificanceCutDirection = ">", # options: ==, !=, <, <=, >, >=
               METType = "MET_Type1", # options: MET_Type1, MET_Type1_NoHF, MET_Puppi, GenMET, L1MET, HLTMET, CaloMET
   applyPhiCorrections = False  # FIXME: no effect yet
)
# MET trigger SF
scaleFactors.assignMETTriggerSF(metSelection, bjetSelection.bjetDiscrWorkingPoint, "nominal", trg.METtriggerEfficiencyJsonName)
#====== Angular cuts / back-to-back
angularCutsBackToBack = PSet(
       nConsideredJets = 3,    # Number of highest-pt jets to consider (excluding jet corresponding to tau)
enableOptimizationPlots = True, # 2D histograms for optimizing angular cuts
        cutValueJet1 = 40.0,   # Cut value in degrees (circular cut)
        cutValueJet2 = 40.0,   # Cut value in degrees (circular cut)
        cutValueJet3 = 40.0,   # Cut value in degrees (circular cut)
        cutValueJet4 = 40.0,   # Cut value in degrees (circular cut)
)
#====== Experimental


#================================================================================================
# Top selection BDT
#================================================================================================
topSelectionBDT = PSet(
        AnyTopMVACutValue      = -0.95,   # [default: -1.0] NOTE: Defines StandardSelections
        AnyTopMVACutDirection  =  ">",    # [default: ">"]
        TopMVACutValue         =  0.40,   # [default: 0.40] NOTE: Only use numbers with 2 decimals (e.g 0.40, 0.30, 0.00)
        TopMVACutDirection     =  ">=",   # [default: ">="]
        TopMassLowCutValue     =   0.0,   # [default: 0.0]
        TopMassLowCutDirection =  ">=",   # [default: ">="]
        TopMassUppCutValue     =  400.0,  # [default: 400.0]  # Do not evaluate top candidate if top mass greater than this cut (600 takes TOO long!)
        TopMassUppCutDirection =  "<=",   # [default: "<"]
        CSV_bDiscCutValue      = 0.8484,  # [default: 0.8484] # Do not evaluate top candidate if b-jet assigned as b from top fails this cut
        CSV_bDiscCutDirection  = ">=",    # [default: ">="]
        WeightFile             = "BDTG_DeltaR0p3_DeltaPtOverPt0p32_BJetPt40_noTopPtRew_24Oct2018.weights.xml",  # Disabled top-pt reweighting
)

# top-tagging (json files available for: defaut, fatJet, ldgJet)
MVAstring = "%.2f" % topSelectionBDT.TopMVACutValue
# Determine which top JSON files to use depending on the BDT trainigh weightfile used
if "noDeltaRqq_noTopPtRew" in topSelectionBDT.WeightFile:
    # dR(q,q') > 0.8 removed from training (q,q': partons from top decay)
    topMisID     = "topMisID_BDT0p40_TopMassCut400_BDTGnoDRqq_noTopPtRew.json"
    topTagEff    = "toptagEff_BDT0p40_GenuineTT_TopMassCut400_BDTGnoDRqq_noTopPtRew.json"
    topTagEffUnc = "toptagEffUncert_BDT0p40_GenuineTT_TopMassCut400_BDTGnoDRqq_noTopPtRew.json"
elif "noDeltaRqq" in topSelectionBDT.WeightFile:
    # dR(q,q') > 0.8 removed from training (q,q': partons from top decay)
    topMisID     = "topMisID_BDT0p40_TopMassCut400_BDTGnoDRqq.json"
    topTagEff    = "toptagEff_BDT0p40_GenuineTT_TopMassCut400_BDTGnoDRqq.json"
    topTagEffUnc = "toptagEffUncert_BDT0p40_GenuineTT_TopMassCut400_BDTGnoDRqq.json"
elif "noTopPtRew" in topSelectionBDT.WeightFile:
    # Disabled top-pt reweighting
    topMisID     = "topMisID_BDT0p40_TopMassCut400_noTopPtRew.json"
    topTagEff    = "toptagEff_BDT0p40_GenuineTT_TopMassCut400_noTopPtRew.json"
    topTagEffUnc = "toptagEffUncert_BDT0p40_GenuineTT_TopMassCut400_noTopPtRew.json"
else:
    # Defaut
    topMisID     = "topMisID_BDT%s_TopMassCut400.json" % MVAstring.replace(".", "p").replace("-", "m")
    topTagEff    = "toptagEff_BDT%s_GenuineTT_TopMassCut400.json" % MVAstring.replace(".", "p").replace("-", "m")
    topTagEffUnc = "toptagEffUncert_BDT%s_GenuineTT_TopMassCut400.json" % MVAstring.replace(".", "p").replace("-", "m")
scaleFactors.setupToptagSFInformation(topTagPset                     = topSelectionBDT,
                                      topTagMisidFilename            = topMisID,
                                      topTagEfficiencyFilename       = topTagEff,
                                      topTagEffUncertaintiesFilename = topTagEffUnc,
                                      direction                      = "nominal",
                                      variationInfo                  = None)

jetCorrelations = PSet (

)


def setAngularCutsWorkingPoint(pset, workingPoint):
    if workingPoint == "NoCut":
        pset.cutValueJet1 = 0.0
        pset.cutValueJet2 = 0.0
        pset.cutValueJet3 = 0.0
        pset.cutValueJet4 = 0.0
    elif workingPoint == "Loose":
        pset.cutValueJet1 = 40.0
        pset.cutValueJet2 = 40.0
        pset.cutValueJet3 = 40.0
        pset.cutValueJet4 = 40.0
    elif workingPoint == "Medium":
        pset.cutValueJet1 = 60.0
        pset.cutValueJet2 = 60.0
        pset.cutValueJet3 = 60.0
        pset.cutValueJet4 = 60.0
    elif workingPoint == "Tight":
        pset.cutValueJet1 = 80.0
        pset.cutValueJet2 = 80.0
        pset.cutValueJet3 = 80.0
        pset.cutValueJet4 = 80.0
    elif workingPoint == "VTight":
        pset.cutValueJet1 = 100.0
        pset.cutValueJet2 = 100.0
        pset.cutValueJet3 = 100.0
        pset.cutValueJet4 = 100.0
    else:
        raise Exception("Error: Unknown working point '%s' requested!"%workingPoint)

#====== Common plots options
commonPlotsOptions = PSet(
  # Splitting of histograms as function of one or more parameters
  # Example: histogramSplitting = [PSet(label="tauPt", binLowEdges=[60, 70, 80, 100, 120], useAbsoluteValues=False)],
  histogramSplitting = [],
  # By default, inclusive (i.e. fake tau+genuine tau) and fake tau histograms are produced. Set to true to also produce genuine tau histograms (Note: will slow down running and enlarge resulting files).
  enableGenuineTauHistograms = False, 
  # Bin settings (final bin setting done in datacardGenerator, there also variable bin width is supported)
       nVerticesBins = PSet(nBins=60, axisMin=0., axisMax=60.),
              ptBins = PSet(nBins=500, axisMin=0., axisMax=5000.),
             etaBins = PSet(nBins=100, axisMin=-5.0, axisMax=5.0),
             phiBins = PSet(nBins=72, axisMin=-3.1415926, axisMax=3.1415926),
        deltaEtaBins = PSet(nBins=50, axisMin=0., axisMax=10.0),
        deltaPhiBins = PSet(nBins=18, axisMin=0., axisMax=180.), # used in 2D plots, i.e. putting high number of bins here will cause troubles
        deltaRBins   = PSet(nBins=50, axisMin=0., axisMax=10.),
            rtauBins = PSet(nBins=55, axisMin=0., axisMax=1.1),
           njetsBins = PSet(nBins=20, axisMin=0., axisMax=20.),
             metBins = PSet(nBins=80, axisMin=0., axisMax=800.), # please use 10 GeV bin width because of QCD measurement
             htBins = PSet(nBins=240, axisMin=0., axisMax=2400.), 
       bjetDiscrBins = PSet(nBins=20, axisMin=-1.0, axisMax=1.0),
   angularCuts1DBins = PSet(nBins=52, axisMin=0., axisMax=260.),
         topMassBins = PSet(nBins=60, axisMin=0., axisMax=600.),
           wMassBins = PSet(nBins=60, axisMin=0., axisMax=300.),
              mtBins = PSet(nBins=2000, axisMin=0., axisMax=10000.), # 5 GeV bin width for tail fitter
         invMassBins = PSet(nBins=500, axisMin=0., axisMax=5000.),
  # Enable/Disable some debug-level plots
       enablePUDependencyPlots = True,
)

#====== Build all selections group
allSelections = PSet(
 histogramAmbientLevel = histoLevel,
               Trigger = trg,
             METFilter = metFilter,
          TauSelection = tauSelection,
     ElectronSelection = eVeto,
         MuonSelection = muVeto,
      MuonForEmbedding = muForEmbedding,
          JetSelection = jetSelection,
  AngularCutsCollinear = angularCutsCollinear,
         BJetSelection = bjetSelection,
          METSelection = metSelection,
       TopSelectionBDT = topSelectionBDT,
 AngularCutsBackToBack = angularCutsBackToBack,
       JetCorrelations = jetCorrelations,
           CommonPlots = commonPlotsOptions,
)


## Parses command line parameters and returns suffix for analysis
def obtainAnalysisSuffix(argv):
    suffix = "" 
    if "1prong" in argv or "1pr" in argv:
        suffix = "1pr"
        print "Running on 1-prong taus"
    elif "2prong" in argv or "2pr" in argv:
        suffix = "2pr"
        print "Running on 2-prong taus"
    elif "3prong" in argv or "3pr" in argv:
        suffix = "3pr"
        print "Running on 3-prong taus"
    return suffix

## Parses command line parameters and adjusts the parameters accordingly
def applyAnalysisCommandLineOptions(argv, config):
    if len(argv) < 3:
        return
    print "Applying command line options"
    if "1prong" in argv or "1pr" in argv:
        config.TauSelection.prongs = 1
    elif "2prong" in argv or "2pr" in argv:
        config.TauSelection.prongs = 2
    elif "3prong" in argv or "3pr" in argv:
        config.TauSelection.prongs = 3
    scaleFactors.assignTauTriggerSF(config.TauSelection, "nominal",config.Trigger.TautriggerEfficiencyJsonName)
