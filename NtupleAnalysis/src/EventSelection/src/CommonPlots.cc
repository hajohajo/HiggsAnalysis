#include "EventSelection/interface/CommonPlots.h"
#include "EventSelection/interface/TransverseMass.h"
#include "EventSelection/interface/PUDependencyPlots.h"
#include "DataFormat/interface/Event.h"

CommonPlots::CommonPlots(const ParameterSet& config, const AnalysisType type, HistoWrapper& histoWrapper)
  : fEnableGenuineTauHistograms(true), // Needed always for limits
    //fEnableGenuineTauHistograms(config.getParameter<bool>("enableGenuineTauHistograms")),
    // Analysis type
    fAnalysisType(type),
    // HistoWrapper
    fHistoWrapper(histoWrapper),
    // Histogram splitter
    fHistoSplitter(config, histoWrapper),
    // Settings for histogram binning
    fNVerticesBinSettings(config.getParameter<ParameterSet>("nVerticesBins")),
    fPtBinSettings(config.getParameter<ParameterSet>("ptBins")),
    fEtaBinSettings(config.getParameter<ParameterSet>("etaBins")),
    fPhiBinSettings(config.getParameter<ParameterSet>("phiBins")),
    fDeltaEtaBinSettings(config.getParameter<ParameterSet>("deltaEtaBins")),
    fDeltaPhiBinSettings(config.getParameter<ParameterSet>("deltaPhiBins")),
    fDeltaRBinSettings(config.getParameter<ParameterSet>("deltaRBins")),
    fRtauBinSettings(config.getParameter<ParameterSet>("rtauBins")),
    fNjetsBinSettings(config.getParameter<ParameterSet>("njetsBins")),
    fMetBinSettings(config.getParameter<ParameterSet>("metBins")),
    fHtBinSettings(config.getParameter<ParameterSet>("htBins")),
    fBJetDiscriminatorBinSettings(config.getParameter<ParameterSet>("bjetDiscrBins")),
    fAngularCuts1DSettings(config.getParameter<ParameterSet>("angularCuts1DBins")),
    fDnnSelectionBinSettings(config.getParameter<ParameterSet>("dnnSelectionBins")),
    fWMassBinSettings(config.getParameter<ParameterSet>("wMassBins")),
    fTopMassBinSettings(config.getParameter<ParameterSet>("topMassBins")),
    fInvmassBinSettings(config.getParameter<ParameterSet>("invMassBins")),
    fMtBinSettings(config.getParameter<ParameterSet>("mtBins")),
    hNSelectedVsRunNumber(nullptr)
{ 
  // Create CommonPlotsBase objects
  bool enableStatus = config.getParameter<bool>("enablePUDependencyPlots");
  if (fAnalysisType == kQCDNormalizationSystematicsSignalRegion || fAnalysisType == kQCDNormalizationSystematicsControlRegion) {
    enableStatus = false;
  }
  fPUDependencyPlots = new PUDependencyPlots(histoWrapper, enableStatus, fNVerticesBinSettings);
  fBaseObjects.push_back(fPUDependencyPlots);
}


CommonPlots::CommonPlots(const ParameterSet& config, const AnalysisType type, HistoWrapper& histoWrapper, bool test)
  : fEnableGenuineTauHistograms(config.getParameter<bool>("enableGenuineBHistograms")), // means GenuineB for Htb
    fAnalysisType(type),
    fHistoWrapper(histoWrapper),
    fHistoSplitter(config, histoWrapper),
    fNVerticesBinSettings(config.getParameter<ParameterSet>("nVerticesBins")),
    fPtBinSettings(config.getParameter<ParameterSet>("ptBins")),
    fEtaBinSettings(config.getParameter<ParameterSet>("etaBins")),
    fPhiBinSettings(config.getParameter<ParameterSet>("phiBins")),
    fDeltaEtaBinSettings(config.getParameter<ParameterSet>("deltaEtaBins")),
    fDeltaPhiBinSettings(config.getParameter<ParameterSet>("deltaPhiBins")),
    fDeltaRBinSettings(config.getParameter<ParameterSet>("deltaRBins")),
    fRtauBinSettings(config.getParameter<ParameterSet>("rtauBins")),
    fNjetsBinSettings(config.getParameter<ParameterSet>("njetsBins")),
    fMetBinSettings(config.getParameter<ParameterSet>("metBins")),
    fHtBinSettings(config.getParameter<ParameterSet>("htBins")),
    fBJetDiscriminatorBinSettings(config.getParameter<ParameterSet>("bjetDiscrBins")),
    fAngularCuts1DSettings(config.getParameter<ParameterSet>("angularCuts1DBins")),
    fDnnSelectionBinSettings(config.getParameter<ParameterSet>("dnnSelectionBins")),
    fWMassBinSettings(config.getParameter<ParameterSet>("wMassBins")),
    fTopMassBinSettings(config.getParameter<ParameterSet>("topMassBins")),
    fInvmassBinSettings(config.getParameter<ParameterSet>("invMassBins")),
    fMtBinSettings(config.getParameter<ParameterSet>("mtBins")),
    hNSelectedVsRunNumber(nullptr)
{ 
  // Create CommonPlotsBase objects
  bool enableStatus = config.getParameter<bool>("enablePUDependencyPlots");
  if (fAnalysisType == kQCDNormalizationSystematicsSignalRegion || fAnalysisType == kQCDNormalizationSystematicsControlRegion) {
    enableStatus = false;
  }
  fPUDependencyPlots = new PUDependencyPlots(histoWrapper, enableStatus, fNVerticesBinSettings);
  fBaseObjects.push_back(fPUDependencyPlots);
}


CommonPlots::~CommonPlots() {
  fHistoSplitter.deleteHistograms(hCtrlNjets);
  fHistoSplitter.deleteHistograms(hCtrlNjetsAfterJetSelectionAndMETSF);
  fHistoSplitter.deleteHistograms(hCtrlCollinearAngularCutsMinimum);
  fHistoSplitter.deleteHistograms(hCtrlCollinearAngularCutsJet1);
  fHistoSplitter.deleteHistograms(hCtrlCollinearAngularCutsJet2);
  fHistoSplitter.deleteHistograms(hCtrlCollinearAngularCutsJet3);
  fHistoSplitter.deleteHistograms(hCtrlCollinearAngularCutsJet4);
  fHistoSplitter.deleteHistograms(hCtrlAngularCutsDeltaPhiTaus);
  fHistoSplitter.deleteHistograms(hCtrlAngularCutsDeltaPhiTauMET);
  fHistoSplitter.deleteHistograms(hCtrlAngularCutsDeltaPhiLdgTauMET);
  fHistoSplitter.deleteHistograms(hCtrlAngularCutsDeltaPhiSubldgTauMET);
  fHistoSplitter.deleteHistograms(hCtrlAngularCutsDeltaPhiLdgJetMET);
  fHistoSplitter.deleteHistograms(hCtrlAngularCutsDeltaPhiSubldgJetMET);
  fHistoSplitter.deleteHistograms(hCtrlAngularCutsDeltaPhiMuonMET);
  fHistoSplitter.deleteHistograms(hCtrlAngularCutsDeltaPhiLdgTauMuon);
  fHistoSplitter.deleteHistograms(hCtrlAngularCutsDeltaPhiSubldgTauMuon);
  fHistoSplitter.deleteHistograms(hCtrlAngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiSubldgTauMET);
  fHistoSplitter.deleteHistograms(hCtrlAngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiMuonMET);
  fHistoSplitter.deleteHistograms(hCtrlAngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiLdgBJetMET);
  fHistoSplitter.deleteHistograms(hCtrlAngularCutsDeltaPhiLdgTauMuon_Vs_DeltaPhiSubldgTauMuon);
  fHistoSplitter.deleteHistograms(hCtrlNVerticesAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauPtAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauEtaAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauPhiAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauEtaPhiAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauLdgTrkPtAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauDecayModeAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauNProngsAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauRtauAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauSourceAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSubldgTauPtAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSubldgTauEtaAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSubldgTauPhiAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSubldgTauEtaPhiAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSubldgTauLdgTrkPtAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSubldgTauDecayModeAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSubldgTauNProngsAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedMuonPtAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedMuonEtaAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedMuonPhiAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedMuonEtaPhiAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlNJetsAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlJetPtAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlJetEtaAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlJetEtaPhiAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlMETAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlMETPhiAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlDeltaPhiTauMetAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlDeltaPhiMuonMetAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlNBJetsAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlBJetPtAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlBJetEtaAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlBDiscriminatorAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlHTAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlMHTAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlNTopsAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopPtAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopMVAAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopDijetPtAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopDijetMassAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopMassAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopMassWMassRatioAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopPt_Vs_TopDijetPtAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopBJetPtAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopBJetEtaAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopBJetBdiscAfterStdSelections);
  fHistoSplitter.deleteHistograms(hCtrlMET);
  fHistoSplitter.deleteHistograms(hCtrlMETPhi);
  fHistoSplitter.deleteHistograms(hCtrlNBJets);
  fHistoSplitter.deleteHistograms(hCtrlBJetPt);
  fHistoSplitter.deleteHistograms(hCtrlBJetEta);
  fHistoSplitter.deleteHistograms(hCtrlBDiscriminator);
  fHistoSplitter.deleteHistograms(hCtrlBackToBackAngularCutsMinimum);
  fHistoSplitter.deleteHistograms(hCtrlBackToBackAngularCutsJet1);
  fHistoSplitter.deleteHistograms(hCtrlBackToBackAngularCutsJet2);
  fHistoSplitter.deleteHistograms(hCtrlBackToBackAngularCutsJet3);
  fHistoSplitter.deleteHistograms(hCtrlBackToBackAngularCutsJet4);
  fHistoSplitter.deleteHistograms(hCtrlDnnSelection);
  fHistoSplitter.deleteHistograms(hCtrlDnnSelectionLoose);
  fHistoSplitter.deleteHistograms(hCtrlDnnSelectionMedium);
  fHistoSplitter.deleteHistograms(hCtrlDnnSelectionTight);
  fHistoSplitter.deleteHistograms(hCtrlNVerticesAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTausPtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTausEtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTausPhiAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauPtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauEtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauPhiAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauEtaPhiAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauLdgTrkPtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauDecayModeAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauNProngsAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauRtauAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauSourceAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedTauIPxyAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSubldgTauPtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSubldgTauEtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSubldgTauPhiAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSubldgTauEtaPhiAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSubldgTauLdgTrkPtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSubldgTauDecayModeAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSubldgTauNProngsAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSubldgTauIPxyAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedMuonPtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedMuonEtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedMuonPhiAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedElectronPtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedElectronEtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlSelectedElectronPhiAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlNJetsAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlJetPtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlJetEtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlJetEtaPhiAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlHTAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlMHTAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlMinDeltaPhiJetMHTAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlMaxDeltaPhiJetMHTAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlMinDeltaRJetMHTAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlMinDeltaRReversedJetMHTAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlCollinearAngularCutsMinimumAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlMETAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlMETPhiAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlNBJetsAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlBJetPtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlBJetEtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlJet1PtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlJet2PtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlJet3PtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlJet4PtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlJet5PtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlJet6PtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlJet7PtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlJet1EtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlJet2EtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlJet3EtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlJet4EtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlJet5EtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlJet6EtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlJet7EtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlBJet1PtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlBJet2PtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlBJet3PtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlBJet4PtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlBJet1EtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlBJet2EtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlBJet3EtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlBJet4EtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlBDiscriminatorAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlBackToBackAngularCutsMinimumAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlDnnSelectionAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlDeltaPhiTauMetAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlDeltaPhiMuonMetAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlNTopsAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopPtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopMVAAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopDijetPtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopDijetMassAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopMassAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopMassWMassRatioAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopPt_Vs_TopDijetPtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopBJetPtAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopBJetEtaAfterAllSelections);
  fHistoSplitter.deleteHistograms(hCtrlTopBJetBdiscAfterAllSelections);
  fHistoSplitter.deleteHistograms(hShapeTransverseMass);
  fHistoSplitter.deleteHistograms(hShapeProbabilisticBtagTransverseMass);
  if (hNSelectedVsRunNumber != nullptr) delete hNSelectedVsRunNumber;
  for (auto p: fBaseObjects) delete p;
}

void CommonPlots::book(TDirectory *dir, bool isData) { 
  fHistoSplitter.bookHistograms(dir);
  // Create directories for data driven control plots
  std::string tausOrB = "Taus";
  if ( (fAnalysisType == kFakeBMeasurement) || (fAnalysisType == kHplus2tbAnalysis) ) tausOrB = "B";

  std::string myLabel = "ForDataDrivenCtrlPlots";
  std::string myFakeLabel = "ForDataDrivenCtrlPlotsEWKFake" + tausOrB;
  std::string myGenuineLabel = "ForDataDrivenCtrlPlotsEWKGenuine" + tausOrB;
  
  if (fAnalysisType == kQCDNormalizationSystematicsSignalRegion) {
    myLabel += "QCDNormalizationSignal";
    myFakeLabel += "QCDNormalizationSignalEWKFake" + tausOrB;
    myGenuineLabel += "QCDNormalizationSignalEWKGenuine" + tausOrB;
  }
  if (fAnalysisType == kQCDNormalizationSystematicsControlRegion) {
    myLabel += "QCDNormalizationControl";
    myFakeLabel += "QCDNormalizationControlEWKFake" + tausOrB;
    myGenuineLabel += "QCDNormalizationControlEWKGenuine" + tausOrB;
  }

  // Create the directories
  TDirectory* myCtrlDir            = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myLabel);
  TDirectory* myCtrlEWKFakeTausDir = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myFakeLabel);
  TDirectory* myCtrlGenuineTausDir = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myGenuineLabel);
  std::vector<TDirectory*> myDirs2 = {myCtrlDir, myCtrlEWKFakeTausDir};
  std::vector<TDirectory*> myDirs3 = {myCtrlDir, myCtrlEWKFakeTausDir, myCtrlGenuineTausDir};
  std::vector<TDirectory*> myDirs;
  
  /// Needed for TauIDSyst
  auto dirName = dir->GetName();
  std::string str(dirName);

  tauIDup=false;
  tauIDdown=false;

  if(str.std::string::find("TauIDSystPlus") != std::string::npos){
    tauIDup=true;
  }
  if(str.std::string::find("TauIDSystMinus")!= std::string::npos){
    tauIDdown=true;
  }

  if (fEnableGenuineTauHistograms) {
    for (auto& p: myDirs3)
      myDirs.push_back(p);
  } else {
    for (auto& p: myDirs2)
      myDirs.push_back(p);
  }
    
  // Create histograms
  const bool hplus2tb     = ( (fAnalysisType == kFakeBMeasurement) || (fAnalysisType == kHplus2tbAnalysis) );
  const bool hplus2hw     = ( (fAnalysisType == kHplus2hwAnalysis) || (fAnalysisType == kHplus2hwAnalysisWithTop) || (fAnalysisType == kQCDMeasurement_mmt) || (fAnalysisType == kHplus2hwAnalysis_mmt));
  const bool hplus2hw_top = ( (fAnalysisType == kHplus2hwAnalysisWithTop) );
  const bool hplus2hw_ele = ( (fAnalysisType == kHplus2hw_ele_Analysis) || (fAnalysisType == kQCDMeasurement_eet) || (fAnalysisType == kHplus2hwAnalysis_eet));

  // vertex

  // tau selection

  // tau trigger SF

  // veto tau selection
  
  // electron veto
  
  // muon veto

  // tau veto


  //==========================================
  // jet selection
  //==========================================
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlNjets, 
						   "Njets", ";Number of selected jets;N_{events}",
						   fNjetsBinSettings.bins(), fNjetsBinSettings.min(), fNjetsBinSettings.max());
  
  //==========================================
  // MET trigger SF
  //==========================================
  if (!hplus2tb)
    {
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlNjetsAfterJetSelectionAndMETSF,
						       "NjetsAfterJetSelectionAndMETSF", ";Number of selected jets;N_{events}",
						       fNjetsBinSettings.bins(), fNjetsBinSettings.min(), fNjetsBinSettings.max());
    }// if (!hplus2tb)

  
  //==========================================     
  // collinear angular cuts 
  //==========================================
  if (!hplus2tb)
    {
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlCollinearAngularCutsMinimum, 
						       "CollinearAngularCutsMinimum", ";min(#sqrt{#Delta#phi(#tau,MET)^{2}+(180^{#circ}-#Delta#phi(jet_{1..n},MET))^{2}}), ^{#circ};N_{events}", 
						       fAngularCuts1DSettings.bins(), fAngularCuts1DSettings.min(), fAngularCuts1DSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlCollinearAngularCutsJet1, 
						       "CollinearAngularCutsJet1", ";min(#sqrt{#Delta#phi(#tau,MET)^{2}+(180^{#circ}-#Delta#phi(jet_{1},MET))^{2}}), ^{#circ};N_{events}", 
						       fAngularCuts1DSettings.bins(), fAngularCuts1DSettings.min(), fAngularCuts1DSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlCollinearAngularCutsJet2, 
						       "CollinearAngularCutsJet2", ";min(#sqrt{#Delta#phi(#tau,MET)^{2}+(180^{#circ}-#Delta#phi(jet_{2},MET))^{2}}), ^{#circ};N_{events}", 
						       fAngularCuts1DSettings.bins(), fAngularCuts1DSettings.min(), fAngularCuts1DSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlCollinearAngularCutsJet3, 
						       "CollinearAngularCutsJet3", ";min(#sqrt{#Delta#phi(#tau,MET)^{2}+(180^{#circ}-#Delta#phi(jet_{3},MET))^{2}}), ^{#circ};N_{events}", 
						       fAngularCuts1DSettings.bins(), fAngularCuts1DSettings.min(), fAngularCuts1DSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlCollinearAngularCutsJet4, 
						       "CollinearAngularCutsJet4", ";min(#sqrt{#Delta#phi(#tau,MET)^{2}+(180^{#circ}-#Delta#phi(jet_{4},MET))^{2}}), ^{#circ};N_{events}", 
						       fAngularCuts1DSettings.bins(), fAngularCuts1DSettings.min(), fAngularCuts1DSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlAngularCutsDeltaPhiTaus,
						       "AngularCutsDeltaPhiTaus", ";#Delta#phi(#tau_{1},#tau_{2}) (^{#circ});N_{events}", 
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlAngularCutsDeltaPhiTauMET,
						       "AngularCutsDeltaPhiTauMET", ";#Delta#phi(#tau_{1}+#tau_{2}, MET) (^{#circ});N_{events}", 
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlAngularCutsDeltaPhiLdgTauMET,
						       "AngularCutsDeltaPhiLdgTauMET", ";#Delta#phi(#tau_{1}, MET) (^{#circ});N_{events}", 
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlAngularCutsDeltaPhiSubldgTauMET,
						       "AngularCutsDeltaPhiSubldgTauMET", ";#Delta#phi(#tau_{2}, MET) (^{#circ});N_{events}", 
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlAngularCutsDeltaPhiLdgJetMET,
						       "AngularCutsDeltaPhiLdgJetMET", ";#Delta#phi(jet_{1}, MET) (^{#circ});N_{events}", 
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlAngularCutsDeltaPhiSubldgJetMET,
						       "AngularCutsDeltaPhiSubldgJetMET", ";#Delta#phi(jet_{1}, MET) (^{#circ});N_{events}", 
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlAngularCutsDeltaPhiMuonMET,
						       "AngularCutsDeltaPhiMuonMET", ";#Delta#phi(jet_{1}, MET) (^{#circ});N_{events}", 
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlAngularCutsDeltaPhiLdgTauMuon,
						       "AngularCutsDeltaPhiLdgTauMuon", ";#Delta#phi(#tau_{1}, #mu) (^{#circ});N_{events}", 
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlAngularCutsDeltaPhiSubldgTauMuon,
						       "AngularCutsDeltaPhiSubldgTauMuon", ";#Delta#phi(#tau_{2}, #mu) (^{#circ});N_{events}", 
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH2F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlAngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiSubldgTauMET,
						       "AngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiSubldgTauMET" , ";#Delta#phi(#tau_{1}, MET) (^{#circ});#Delta#phi(#tau_{2}, MET) (^{#circ})",
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max(),
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH2F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlAngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiMuonMET,
						       "AngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiMuonMET" , ";#Delta#phi(#tau_{1}, MET) (^{#circ});#Delta#phi(#mu, MET) (^{#circ})",
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max(),
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH2F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlAngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiLdgBJetMET,
						       "AngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiLdgBJetMET" , ";#Delta#phi(#tau_{1}, MET) (^{#circ});#Delta#phi(b jet_{1}, MET) (^{#circ})",
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max(),
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH2F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlAngularCutsDeltaPhiLdgTauMuon_Vs_DeltaPhiSubldgTauMuon,
						       "AngularCutsDeltaPhiLdgTauMuon_Vs_DeltaPhiSubldgTauMuon" , ";#Delta#phi(#tau_{1}, #mu) (^{#circ});#Delta#phi(#tau_{2}, #mu) (^{#circ})",
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max(),
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max());
    }// if (!hplus2tb)


  //==========================================     
  // standard selections
  //==========================================     
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlNVerticesAfterStdSelections, 
						   "NVertices_AfterStandardSelections", ";N_{vertices};N_{events}",
						   fNVerticesBinSettings.bins(), fNVerticesBinSettings.min(), fNVerticesBinSettings.max());

  //==========================================     
  // standard selections: tau 
  //==========================================     
  if (!hplus2tb)
    {
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTauPtAfterStdSelections, 
						       "SelectedTau_pT_AfterStandardSelections", ";#tau p_{T}, GeV/c;N_{events}",
						       fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTauEtaAfterStdSelections, 
						       "SelectedTau_eta_AfterStandardSelections", ";#tau #eta;N_{events}",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTauPhiAfterStdSelections, 
						       "SelectedTau_phi_AfterStandardSelections", ";#tau #phi;N_{events}",
						       fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH2F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlSelectedTauEtaPhiAfterStdSelections, 
						       "SelectedTau_etaphi_AfterStandardSelections", ";#tau #eta;#tau #phi;",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max(),
						       fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTauLdgTrkPtAfterStdSelections, 
						       "SelectedTau_ldgTrkPt_AfterStandardSelections", ";#tau ldg. trk p_{T}, GeV/c;N_{events}",
						       fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTauDecayModeAfterStdSelections, 
						       "SelectedTau_DecayMode_AfterStandardSelections", ";#tau decay mode;N_{events}",
						       20, 0, 20);
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTauNProngsAfterStdSelections, 
						       "SelectedTau_Nprongs_AfterStandardSelections", ";N_{prongs};N_{events}",
						       10, 0, 10);
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTauRtauAfterStdSelections, 
						       "SelectedTau_Rtau_AfterStandardSelections", ";R_{#tau};N_{events}",
						       fRtauBinSettings.bins(), fRtauBinSettings.min(), fRtauBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTauSourceAfterStdSelections, 
						       "SelectedTau_source_AfterStandardSelections", ";;N_{events}",
						       fHelper.getTauSourceBinCount(), 0, fHelper.getTauSourceBinCount());
      
      for (int i = 0; i < fHelper.getTauSourceBinCount(); ++i) 
	{
	  fHistoSplitter.SetBinLabel(hCtrlSelectedTauSourceAfterStdSelections, i+1, fHelper.getTauSourceBinLabel(i));
	}


      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSubldgTauPtAfterStdSelections, 
						       "SubldgTau_pT_AfterStandardSelections", ";#tau p_{T}, GeV/c;N_{events}",
						       fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSubldgTauEtaAfterStdSelections, 
						       "SubldgTau_eta_AfterStandardSelections", ";#tau #eta;N_{events}",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSubldgTauPhiAfterStdSelections, 
						       "SubldgTau_phi_AfterStandardSelections", ";#tau #phi;N_{events}",
						       fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH2F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlSubldgTauEtaPhiAfterStdSelections, 
						       "SubldgTau_etaphi_AfterStandardSelections", ";#tau #eta;#tau #phi;",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max(),
						       fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSubldgTauLdgTrkPtAfterStdSelections, 
						       "SubldgTau_ldgTrkPt_AfterStandardSelections", ";#tau ldg. trk p_{T}, GeV/c;N_{events}",
						       fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSubldgTauDecayModeAfterStdSelections, 
						       "SubldgTau_DecayMode_AfterStandardSelections", ";#tau decay mode;N_{events}",
						       20, 0, 20);
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSubldgTauNProngsAfterStdSelections, 
						       "SubldgTau_Nprongs_AfterStandardSelections", ";N_{prongs};N_{events}",
						       10, 0, 10);
    } // if (!hplus2tb)

  //========================================== 
  // standard selections: muon
  //==========================================
  if (!hplus2tb)
    {
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedMuonPtAfterStdSelections, 
						       "SelectedMu_pT_AfterStandardSelections", ";#mu p_{T}, GeV/c;N_{events}",
						       fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedMuonEtaAfterStdSelections, 
						       "SelectedMu_eta_AfterStandardSelections", ";#mu #eta;N_{events}",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedMuonPhiAfterStdSelections, 
						       "SelectedMu_phi_AfterStandardSelections", ";#mu #phi;N_{events}",
						       fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH2F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlSelectedMuonEtaPhiAfterStdSelections, 
						       "SelectedMu_etaphi_AfterStandardSelections", ";#mu #eta;#mu #phi;",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max(),
						       fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());
    }// if (!hplus2tb)

  //==========================================
  // standard selections: jets  
  //==========================================
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlNJetsAfterStdSelections, 
						   "Njets_AfterStandardSelections", ";Number of selected jets;N_{events}",
						   fNjetsBinSettings.bins(), fNjetsBinSettings.min(), fNjetsBinSettings.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJetPtAfterStdSelections, 
						   "JetPt_AfterStandardSelections", ";Selected jets p_{T}, GeV/c;N_{events}",
						   fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJetEtaAfterStdSelections, 
						   "JetEta_AfterStandardSelections", ";Selected jets #eta;N_{events}",
						   fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

  fHistoSplitter.createShapeHistogramTriplet<TH2F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlJetEtaPhiAfterStdSelections, 
						   "JetEtaPhi_AfterStandardSelections", ";Selected jets #eta;Selected jets #phi",
						   fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max(),
						   fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());

  if (hplus2tb || hplus2hw_top)
    {
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlNBJetsAfterStdSelections, 
						       "NBjets_AfterStandardSelections", ";Number of selected b jets;N_{events}",
						       fNjetsBinSettings.bins(), fNjetsBinSettings.min(), fNjetsBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBJetPtAfterStdSelections, 
						       "BjetPt_AfterStandardSelections", ";Selected b jets p_{T}, GeV/c;N_{events}",
						       fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBJetEtaAfterStdSelections, 
						       "BjetEta_AfterStandardSelections", ";Selected b jets #eta;N_{events}",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBDiscriminatorAfterStdSelections, 
						       "BtagDiscriminator_AfterStandardSelections", ";b tag discriminator;N_{events}",
						       fBJetDiscriminatorBinSettings.bins(), fBJetDiscriminatorBinSettings.min(), fBJetDiscriminatorBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlHTAfterStdSelections, 
						       "HT_AfterStandardSelections", ";H_{T}, GeV;N_{events}",
						       fHtBinSettings.bins(), fHtBinSettings.min(), fHtBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlMHTAfterStdSelections, 
						       "MHT_AfterStandardSelections", ";MHT, GeV;N_{events}",
						       fMetBinSettings.bins(), fMetBinSettings.min(), fMetBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlNTopsAfterStdSelections,
						       "NTops_AfterStandardSelections", ";top multiplicity;N_{events}",
						       fNjetsBinSettings.bins()*10, fNjetsBinSettings.min(), fNjetsBinSettings.max()*10);

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopPtAfterStdSelections,
						       "TopPt_AfterStandardSelections", ";p_{T} (GeV/c);N_{events}", 						       
						       2*fPtBinSettings.bins(), fPtBinSettings.min(), 2*fPtBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopMVAAfterStdSelections,
						       "TopMVA_AfterStandardSelections", ";MVA score;N_{events}", 40, -1.0, +1.0);

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopDijetPtAfterStdSelections,
						       "TopDijetPt_AfterStandardSelections", ";p_{T} (GeV/c);N_{events}",
						       2*fPtBinSettings.bins(), fPtBinSettings.min(), 2*fPtBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopDijetMassAfterStdSelections,
						       "TopDijetMass_AfterStandardSelections", ";m_{jjb} (GeV/c^{2});N_{events}",
						       fTopMassBinSettings.bins(), fTopMassBinSettings.min(), fTopMassBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopMassAfterStdSelections,
						       "TopMass_AfterStandardSelections", ";m_{jjb} (GeV/c^{2});N_{events}",
						       fTopMassBinSettings.bins(), fTopMassBinSettings.min(), fTopMassBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopMassWMassRatioAfterStdSelections,
						       "TopMassWMassRatio_AfterStandardSelections", ";R_{32}", 100 , 0.0, 10.0);

      fHistoSplitter.createShapeHistogramTriplet<TH2F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopPt_Vs_TopDijetPtAfterStdSelections, 
						       "TopPt_Vs_TopDijetPtAfterStandardSelections" ,";p_{T} (GeV/c);p_{T} (GeV/c)",
						       2*fPtBinSettings.bins(), fPtBinSettings.min(), 2*fPtBinSettings.max(), 
						       2*fPtBinSettings.bins(), fPtBinSettings.min(), 2*fPtBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopBJetPtAfterStdSelections,
						       "TopBjetPt_AfterStandardSelections", ";p_{T} (GeV/c);N_{events}",
						       2*fPtBinSettings.bins(), fPtBinSettings.min(), 2*fPtBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopBJetEtaAfterStdSelections,
						       "TopBjetEta_AfterStandardSelections", ";#eta;N_{events}",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopBJetBdiscAfterStdSelections,
						       "TopBjetBdisc_AfterStandardSelections", ";b-tag discriminator;N_{events}",
						       fBJetDiscriminatorBinSettings.bins(), fBJetDiscriminatorBinSettings.min(), fBJetDiscriminatorBinSettings.max());

      // All Selections
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlNTopsAfterAllSelections,
						       "NTops_AfterAllSelections", ";top multiplicity;N_{events}",
						       fNjetsBinSettings.bins()*10, fNjetsBinSettings.min(), fNjetsBinSettings.max()*10);

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopPtAfterAllSelections,
						       "TopPt_AfterAllSelections", ";p_{T} (GeV/c);N_{events}", 						       
						       2*fPtBinSettings.bins(), fPtBinSettings.min(), 2*fPtBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopMVAAfterAllSelections,
						       "TopMVA_AfterAllSelections", ";MVA score;N_{events}", 40, -1.0, +1.0);

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopDijetPtAfterAllSelections,
						       "TopDijetPt_AfterAllSelections", ";p_{T} (GeV/c);N_{events}", 						       
						       2*fPtBinSettings.bins(), fPtBinSettings.min(), 2*fPtBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopDijetMassAfterAllSelections,
						       "TopDijetMass_AfterAllSelections", ";m_{jjb} (GeV/c^{2});N_{events}",
						       fTopMassBinSettings.bins(), fTopMassBinSettings.min(), fTopMassBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopMassAfterAllSelections,
						       "TopMass_AfterAllSelections", ";m_{jjb} (GeV/c^{2});N_{events}",
						       fTopMassBinSettings.bins(), fTopMassBinSettings.min(), fTopMassBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopMassWMassRatioAfterAllSelections,
						       "TopMassWMassRatioAfterAllSelections", ";R_{32}", 100 , 0.0, 10.0);

      fHistoSplitter.createShapeHistogramTriplet<TH2F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopPt_Vs_TopDijetPtAfterAllSelections, 
						       "TopPt_Vs_TopDijetPtAfterAllSelections" ,";p_{T} (GeV/c);p_{T} (GeV/c)",
						       2*fPtBinSettings.bins(), fPtBinSettings.min(), 2*fPtBinSettings.max(), 
						       2*fPtBinSettings.bins(), fPtBinSettings.min(), 2*fPtBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopBJetPtAfterAllSelections,
						       "TopBjetPt_AfterAllSelections", ";p_{T} (GeV/c);N_{events}",
						       2*fPtBinSettings.bins(), fPtBinSettings.min(), 2*fPtBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopBJetEtaAfterAllSelections,
						       "TopBjetEta_AfterAllSelections", ";#eta;N_{events}",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlTopBJetBdiscAfterAllSelections,
						       "TopBjetBdisc_AfterAllSelections", ";b-tag discriminator;N_{events}",
						       fBJetDiscriminatorBinSettings.bins(), fBJetDiscriminatorBinSettings.min(), fBJetDiscriminatorBinSettings.max());
    }// if (hplus2tb || hplus2hw_top)
  
  //==========================================
  // standard selections: MET
  //==========================================
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlMETAfterStdSelections, 
						   "MET_AfterStandardSelections", ";MET, GeV;N_{events}", 
						   fMetBinSettings.bins(), fMetBinSettings.min(), fMetBinSettings.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlMETPhiAfterStdSelections, 
						   "METPhi_AfterStandardSelections", ";MET #phi;N_{events}", 
						   fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());

  if (!hplus2tb)
    {
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlDeltaPhiTauMetAfterStdSelections, 
						       "DeltaPhiTauMet_AfterStandardSelections", ";#Delta#phi(#tau,MET), {}^{#circ};N_{events}", 
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlDeltaPhiMuonMetAfterStdSelections, 
						       "DeltaPhiMuonMet_AfterStandardSelections", ";#Delta#phi(#mu,MET), {}^{#circ};N_{events}", 
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max());
    }// if (!hplus2tb)
  
  //==========================================
  // MET selection
  //==========================================
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlMET, 
						   "MET", ";MET, GeV;N_{events}", 
						   fHtBinSettings.bins(), fHtBinSettings.min(), fHtBinSettings.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlMETPhi, 
						   "METPhi", ";MET #phi;N_{events}", 
						   fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());
  
  //==========================================
  // b tagging
  //==========================================
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlNBJets, 
						   "NBjets", ";Number of selected b jets;N_{events}",
						   fNjetsBinSettings.bins(), fNjetsBinSettings.min(), fNjetsBinSettings.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBJetPt, 
						   "BJetPt", ";Selected b jets p_{T}, GeV/c;N_{events}",
						   fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBJetEta, 
						   "BJetEta", ";Selected b jets #eta;N_{events}",
						   fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBDiscriminator, 
						   "BtagDiscriminator", ";b tag discriminator;N_{events}",
						   fBJetDiscriminatorBinSettings.bins(), fBJetDiscriminatorBinSettings.min(), fBJetDiscriminatorBinSettings.max());
  
  //==========================================
  // back-to-back angular cuts
  //==========================================
  if (!hplus2tb)
    {      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBackToBackAngularCutsMinimum, 
						       "BackToBackAngularCutsMinimum", ";min(#sqrt{(180^{#circ}-#Delta#phi(#tau,MET))^{2}+#Delta#phi(jet_{1..n},MET)^{2}}), ^{#circ};N_{events}", 
						       fAngularCuts1DSettings.bins(), fAngularCuts1DSettings.min(), fAngularCuts1DSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlBackToBackAngularCutsJet1, 
						       "BackToBackAngularCutsJet1", ";min(#sqrt{#Delta#phi(#tau,MET)^{2}+(180^{#circ}-#Delta#phi(jet_{1},MET))^{2}}), ^{#circ};N_{events}", 
						       fAngularCuts1DSettings.bins(), fAngularCuts1DSettings.min(), fAngularCuts1DSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlBackToBackAngularCutsJet2, 
						       "BackToBackAngularCutsJet2", ";min(#sqrt{#Delta#phi(#tau,MET)^{2}+(180^{#circ}-#Delta#phi(jet_{2},MET))^{2}}), ^{#circ};N_{events}", 
						       fAngularCuts1DSettings.bins(), fAngularCuts1DSettings.min(), fAngularCuts1DSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlBackToBackAngularCutsJet3, 
						       "BackToBackAngularCutsJet3", ";min(#sqrt{#Delta#phi(#tau,MET)^{2}+(180^{#circ}-#Delta#phi(jet_{3},MET))^{2}}), ^{#circ};N_{events}", 
						       fAngularCuts1DSettings.bins(), fAngularCuts1DSettings.min(), fAngularCuts1DSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlBackToBackAngularCutsJet4, 
						       "BackToBackAngularCutsJet4", ";min(#sqrt{#Delta#phi(#tau,MET)^{2}+(180^{#circ}-#Delta#phi(jet_{4},MET))^{2}}), ^{#circ};N_{events}", 
						       fAngularCuts1DSettings.bins(), fAngularCuts1DSettings.min(), fAngularCuts1DSettings.max());
    }// if (!hplus2tb)

  //==========================================
  // Dnn selection cuts
  //==========================================

  if (!hplus2tb)
    {
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlDnnSelection,
                               "DnnSelection", ";Dnn output value", fDnnSelectionBinSettings.bins(), fDnnSelectionBinSettings.min(), fDnnSelectionBinSettings.max());
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlDnnSelectionLoose,
                               "DnnSelection_MtAfterLoose", "m_{T};Events", fMtBinSettings.bins(), fMtBinSettings.min(), fMtBinSettings.max());
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlDnnSelectionMedium,
                               "DnnSelection_MtAfterMedium", "m_{T};Events", fMtBinSettings.bins(), fMtBinSettings.min(), fMtBinSettings.max());
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlDnnSelectionTight,
                               "DnnSelection_MtAfterTight", "m_{T};Events", fMtBinSettings.bins(), fMtBinSettings.min(), fMtBinSettings.max());

    }


  //==========================================
  // all selections: control plots
  //==========================================
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlNVerticesAfterAllSelections, 
						   "NVertices_AfterAllSelections", ";N_{vertices};N_{events}",
						   fNVerticesBinSettings.bins(), fNVerticesBinSettings.min(), fNVerticesBinSettings.max());

  if (!hplus2tb)
    {  
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTausPtAfterAllSelections, 
                                                       "SelectedTaus_pT_AfterAllSelections", ";#taus p_{T}, GeV/c;N_{events}",
                                                       fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTausEtaAfterAllSelections, 
                                                       "SelectedTaus_eta_AfterAllSelections", ";#taus #eta;N_{events}",
                                                       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTausPhiAfterAllSelections, 
                                                       "SelectedTaus_phi_AfterAllSelections", ";#taus #phi;N_{events}",
                                                       fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTauPtAfterAllSelections, 
						       "SelectedTau_pT_AfterAllSelections", ";#tau p_{T}, GeV/c;N_{events}",
						       fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTauEtaAfterAllSelections, 
						       "SelectedTau_eta_AfterAllSelections", ";#tau #eta;N_{events}",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTauPhiAfterAllSelections, 
						       "SelectedTau_phi_AfterAllSelections", ";#tau #phi;N_{events}",
						       fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH2F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlSelectedTauEtaPhiAfterAllSelections, 
						       "SelectedTau_etaphi_AfterAllSelections", ";#tau #eta;#tau #phi;",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max(),
						       fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTauLdgTrkPtAfterAllSelections, 
						       "SelectedTau_ldgTrkPt_AfterAllSelections", ";#tau ldg. trk p_{T}, GeV/c;N_{events}",
						       fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTauDecayModeAfterAllSelections, 
						       "SelectedTau_DecayMode_AfterAllSelections", ";#tau decay mode;N_{events}",
						       20, 0, 20);
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTauNProngsAfterAllSelections, 
						       "SelectedTau_Nprongs_AfterAllSelections", ";N_{prongs};N_{events}",
						       10, 0, 10);
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTauRtauAfterAllSelections, 
						       "SelectedTau_Rtau_AfterAllSelections", ";R_{#tau};N_{events}",
						       fRtauBinSettings.bins(), fRtauBinSettings.min(), fRtauBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTauSourceAfterAllSelections, 
						       "SelectedTau_source_AfterAllSelections", ";;N_{events}",
						       fHelper.getTauSourceBinCount(), 0, fHelper.getTauSourceBinCount());
      
      for (int i = 0; i < fHelper.getTauSourceBinCount(); ++i) 
	{
	  fHistoSplitter.SetBinLabel(hCtrlSelectedTauSourceAfterAllSelections, i+1, fHelper.getTauSourceBinLabel(i));
	}
      

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedTauIPxyAfterAllSelections, 
						       "SelectedTau_IPxy_AfterAllSelections", ";IP_{T} (cm);N_{events}",
						       100, 0, 0.2);
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSubldgTauPtAfterAllSelections, 
						       "SubldgTau_pT_AfterAllSelections", ";#tau p_{T}, GeV/c;N_{events}",
						       fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSubldgTauEtaAfterAllSelections, 
						       "SubldgTau_eta_AfterAllSelections", ";#tau #eta;N_{events}",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSubldgTauPhiAfterAllSelections, 
						       "SubldgTau_phi_AfterAllSelections", ";#tau #phi;N_{events}",
						       fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH2F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlSubldgTauEtaPhiAfterAllSelections, 
						       "SubldgTau_etaphi_AfterAllSelections", ";#tau #eta;#tau #phi;",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max(),
						       fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSubldgTauLdgTrkPtAfterAllSelections, 
						       "SubldgTau_ldgTrkPt_AfterAllSelections", ";#tau ldg. trk p_{T}, GeV/c;N_{events}",
						       fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSubldgTauDecayModeAfterAllSelections, 
						       "SubldgTau_DecayMode_AfterAllSelections", ";#tau decay mode;N_{events}",
						       20, 0, 20);
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSubldgTauNProngsAfterAllSelections, 
						       "SubldgTau_Nprongs_AfterAllSelections", ";N_{prongs};N_{events}",
						       10, 0, 10);
            
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSubldgTauIPxyAfterAllSelections, 
						       "SubldgTau_IPxy_AfterAllSelections", ";IP_{T} (cm);N_{events}",
						       100, 0, 0.2);
    }// if (!hplus2tb)   
  
    if (hplus2hw || fAnalysisType == kQCDMeasurement_muon)
    {
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedMuonPtAfterAllSelections, 
                                                       "SelectedMu_pT_AfterAllSelections", ";#mu p_{T}, GeV/c;N_{events}",
                                                       fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedMuonEtaAfterAllSelections, 
                                                       "SelectedMu_eta_AfterAllSelections", ";#mu #eta;N_{events}",
                                                       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedMuonPhiAfterAllSelections, 
                                                       "SelectedMu_phi_AfterAllSelections", ";#mu #phi;N_{events}",
                                                       fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());

    }// if (hplushw_muon)


    if (hplus2hw_ele || fAnalysisType == kQCDMeasurement_ele)
    {
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedElectronPtAfterAllSelections, 
                                                       "SelectedEle_pT_AfterAllSelections", ";e p_{T}, GeV/c;N_{events}",
                                                       fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedElectronEtaAfterAllSelections, 
                                                       "SelectedEle_eta_AfterAllSelections", ";e #eta;N_{events}",
                                                       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlSelectedElectronPhiAfterAllSelections, 
                                                       "SelectedEle_phi_AfterAllSelections", ";e #phi;N_{events}",
                                                       fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());

    }// if (hplushw_ele)

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlNJetsAfterAllSelections, 
						   "Njets_AfterAllSelections", ";Number of selected jets;N_{events}",
						   fNjetsBinSettings.bins(), fNjetsBinSettings.min(), fNjetsBinSettings.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJetPtAfterAllSelections, 
						   "JetPt_AfterAllSelections", ";Selected jets p_{T}, GeV/c;N_{events}",
						   2*fPtBinSettings.bins(), fPtBinSettings.min(), 2*fPtBinSettings.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJetEtaAfterAllSelections, 
						   "JetEta_AfterAllSelections", ";Selected jets #eta;N_{events}",
						   fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

  fHistoSplitter.createShapeHistogramTriplet<TH2F>(fEnableGenuineTauHistograms, HistoLevel::kInformative, myDirs, hCtrlJetEtaPhiAfterAllSelections, 
						   "JetEtaPhi_AfterAllSelections", ";Selected jets #eta;Selected jets #phi",
						   fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max(),
						   fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());
  
  //==========================================
  /// all selections: Experimental
  //==========================================
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlHTAfterAllSelections, 
						   "HT_AfterAllSelections", ";H_{T}, GeV;N_{events}",
						   fHtBinSettings.bins(), fHtBinSettings.min(), fHtBinSettings.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlMHTAfterAllSelections, 
						   "MHT_AfterAllSelections", ";MHT, GeV;N_{events}",
						   fMetBinSettings.bins(), fMetBinSettings.min(), fMetBinSettings.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlMinDeltaPhiJetMHTAfterAllSelections, 
						   "MinDeltaPhiJetMHT_AfterAllSelections", ";min(#Delta#phi(jet_{i}, MHT-jet_{i}));N_{events}",
						   fPhiBinSettings.bins(), 0.0, fPhiBinSettings.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlMaxDeltaPhiJetMHTAfterAllSelections, 
						   "MaxDeltaPhiJetMHT_AfterAllSelections", ";max(#Delta#phi(jet_{i}, MHT-jet_{i}));N_{events}",
						   fPhiBinSettings.bins(), 0.0, fPhiBinSettings.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlMinDeltaRJetMHTAfterAllSelections, 
						   "MinDeltaRJetMHT_AfterAllSelections", ";min(#DeltaR(jet_{i}, MHT-jet_{i}));N_{events}",
						   fDeltaRBinSettings.bins(), 0.0, fDeltaRBinSettings.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlMinDeltaRReversedJetMHTAfterAllSelections, 
						   "MinDeltaRJetMHTReversed_AfterAllSelections", ";min(#DeltaR(-jet_{i}, MHT-jet_{i}));N_{events}",
						   fDeltaRBinSettings.bins(), 0.0, fDeltaRBinSettings.max());

  //==========================================
  /// after btag SF
  //==========================================

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlNJetsAfterBtagSF, 
						   "Njets_AfterBtagSF", ";Number of selected jets;N_{events}",
						   fNjetsBinSettings.bins(), fNjetsBinSettings.min(), fNjetsBinSettings.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJetPtAfterBtagSF, 
						   "JetPt_AfterBtagSF", ";Selected jets p_{T}, GeV/c;N_{events}",
						   fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBJetPtAfterBtagSF,
						   "BJetPt_AfterBtagSF", ";Selected b jets p_{T}, GeV/c;N_{events}",
						   fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());


  //==========================================
  // all selections
  //========================================== 
  if (!hplus2tb)
    {
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlCollinearAngularCutsMinimumAfterAllSelections, 
						       "CollinearAngularCutsMinimum_AfterAllSelections", ";min(#sqrt{#Delta#phi(#tau,MET)^{2}+(180^{#circ}-#Delta#phi(jet_{1..n},MET))^{2}}), ^{#circ};N_{events}", 
						       fAngularCuts1DSettings.bins(), fAngularCuts1DSettings.min(), fAngularCuts1DSettings.max());
    }// if (!hplus2tb)
      
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlMETAfterAllSelections, 
						   "MET_AfterAllSelections", ";MET, GeV;N_{events}", 
						   fMetBinSettings.bins(), fMetBinSettings.min(), fMetBinSettings.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlMETPhiAfterAllSelections,
						   "METPhi_AfterAllSelections", ";MET #phi;N_{events}", 
						   fPhiBinSettings.bins(), fPhiBinSettings.min(), fPhiBinSettings.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlNBJetsAfterAllSelections,
						   "NBjets_AfterAllSelections", ";Number of selected b jets;N_{events}",
						   fNjetsBinSettings.bins(), fNjetsBinSettings.min(), fNjetsBinSettings.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBJetPtAfterAllSelections,
						   "BJetPt_AfterAllSelections", ";Selected b jets p_{T}, GeV/c;N_{events}",
						   2*fPtBinSettings.bins(), fPtBinSettings.min(), 2*fPtBinSettings.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBJetEtaAfterAllSelections,
						   "BJetEta_AfterAllSelections", ";Selected b jets #eta;N_{events}",
						   fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBDiscriminatorAfterAllSelections, 
						   "BtagDiscriminator_AfterAllSelections", ";b tag discriminator;N_{events}",
						   fBJetDiscriminatorBinSettings.bins(), fBJetDiscriminatorBinSettings.min(), fBJetDiscriminatorBinSettings.max());

  if (hplus2tb)
    {

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJet1PtAfterAllSelections,
						       "Jet1Pt_AfterAllSelections", ";p_{T} (GeV/c);Events",
						       3*fPtBinSettings.bins(), 3*fPtBinSettings.min(), 3*fPtBinSettings.max());
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJet2PtAfterAllSelections,
						       "Jet2Pt_AfterAllSelections", ";p_{T} (GeV/c);Events",
						       3*fPtBinSettings.bins(), 3*fPtBinSettings.min(), 3*fPtBinSettings.max());
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJet3PtAfterAllSelections,
						       "Jet3Pt_AfterAllSelections", ";p_{T} (GeV/c);Events",
						       3*fPtBinSettings.bins(), 3*fPtBinSettings.min(), 3*fPtBinSettings.max());
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJet4PtAfterAllSelections,
						       "Jet4Pt_AfterAllSelections", ";p_{T} (GeV/c);Events",
						       2*fPtBinSettings.bins(), 2*fPtBinSettings.min(), 2*fPtBinSettings.max());
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJet5PtAfterAllSelections,
						       "Jet5Pt_AfterAllSelections", ";p_{T} (GeV/c);Events",
						       2*fPtBinSettings.bins(), 2*fPtBinSettings.min(), 2*fPtBinSettings.max());
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJet6PtAfterAllSelections,
						       "Jet6Pt_AfterAllSelections", ";p_{T} (GeV/c);Events",
						       fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJet7PtAfterAllSelections,
						       "Jet7Pt_AfterAllSelections", ";p_{T} (GeV/c);Events",
						       fPtBinSettings.bins(), fPtBinSettings.min(), fPtBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBJet1PtAfterAllSelections,
						       "BJet1Pt_AfterAllSelections", ";p_{T} (GeV/c);Events",
						       3*fPtBinSettings.bins(), 3*fPtBinSettings.min(), 3*fPtBinSettings.max());
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBJet2PtAfterAllSelections,
						       "BJet2Pt_AfterAllSelections", ";p_{T} (GeV/c);Events",
						       3*fPtBinSettings.bins(), 3*fPtBinSettings.min(), 3*fPtBinSettings.max());
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBJet3PtAfterAllSelections,
						       "BJet3Pt_AfterAllSelections", ";p_{T} (GeV/c);Events",
						       3*fPtBinSettings.bins(), 3*fPtBinSettings.min(), 3*fPtBinSettings.max());
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBJet4PtAfterAllSelections,
						       "BJet4Pt_AfterAllSelections", ";p_{T} (GeV/c);Events",
						       2*fPtBinSettings.bins(), 2*fPtBinSettings.min(), 2*fPtBinSettings.max());
    
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJet1EtaAfterAllSelections,
						       "Jet1Eta_AfterAllSelections", ";#eta;Events",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJet2EtaAfterAllSelections,
						       "Jet2Eta_AfterAllSelections", ";#eta;Events",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJet3EtaAfterAllSelections,
						       "Jet3Eta_AfterAllSelections", ";#eta;Events",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJet4EtaAfterAllSelections,
						       "Jet4Eta_AfterAllSelections", ";#eta;Events",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJet5EtaAfterAllSelections,
						       "Jet5Eta_AfterAllSelections", ";#eta;Events",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJet6EtaAfterAllSelections,
						       "Jet6Eta_AfterAllSelections", ";#eta;Events",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlJet7EtaAfterAllSelections,
						       "Jet7Eta_AfterAllSelections", ";#eta;Events",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBJet1EtaAfterAllSelections,
						       "BJet1Eta_AfterAllSelections", ";#eta;Events",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBJet2EtaAfterAllSelections,
						       "BJet2Eta_AfterAllSelections", ";#eta;Events",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBJet3EtaAfterAllSelections,
						       "BJet3Eta_AfterAllSelections", ";#eta;Events",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBJet4EtaAfterAllSelections,
						       "BJet4Eta_AfterAllSelections", ";#eta;Events",
						       fEtaBinSettings.bins(), fEtaBinSettings.min(), fEtaBinSettings.max());
    }  

  if (!hplus2tb)
    {

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlBackToBackAngularCutsMinimumAfterAllSelections, 
						       "BackToBackAngularCutsMinimum_AfterAllSelections", ";min(#sqrt{(180^{#circ}-#Delta#phi(#tau,MET))^{2}+#Delta#phi(jet_{1..n},MET)^{2}}), ^{#circ};N_{events}", 
						       fAngularCuts1DSettings.bins(), fAngularCuts1DSettings.min(), fAngularCuts1DSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlDeltaPhiTauMetAfterAllSelections, 
						       "DeltaPhiTauMet_AfterAllSelections", ";#Delta#phi(#tau,MET), {}^{#circ};N_{events}", 
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max());

      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlDeltaPhiMuonMetAfterAllSelections, 
						       "DeltaPhiMuonMet_AfterAllSelections", ";#Delta#phi(#mu,MET), {}^{#circ};N_{events}", 
						       fDeltaPhiBinSettings.bins(), fDeltaPhiBinSettings.min(), fDeltaPhiBinSettings.max());

    }// if (!hplus2tb)
    
  //========================================== 
  // Dnn selection cuts
  //==========================================
  if (!hplus2tb)
    {
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, hCtrlDnnSelectionAfterAllSelections,
                               "DnnSelection_AfterAllSelections", ";Dnn output value", fDnnSelectionBinSettings.bins(), fDnnSelectionBinSettings.min(), fDnnSelectionBinSettings.max());
    }

  //==========================================  
  // all selections: shape plots
  //==========================================
  if (!hplus2tb)
    {
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(true, HistoLevel::kSystematics, myDirs3, hShapeTransverseMass, 
						       "shapeTransverseMass", ";m_{T}(tau,MET), GeV/c^{2};N_{events}",
						       fMtBinSettings.bins(), fMtBinSettings.min(), fMtBinSettings.max());
      
      fHistoSplitter.createShapeHistogramTriplet<TH1F>(true, HistoLevel::kSystematics, myDirs3, hShapeProbabilisticBtagTransverseMass, 
						       "shapeTransverseMassProbabilisticBTag", ";m_{T}(tau,MET), GeV/c^{2};N_{events}",
						       fMtBinSettings.bins(), fMtBinSettings.min(), fMtBinSettings.max());
    }
  

  if (isData) 
    {
      hNSelectedVsRunNumber = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, 
							 "NSelectedVsRunNumber", "NSelectedVsRunNumber;Run number;N_{events}", 14000, 246000, 260000);
    }
 
  for (auto& p: fBaseObjects) p->book(dir, isData);
  return;
}

void CommonPlots::initialize() {
  iVertices = -1;
  fTauData  = TauSelection::Data();
  //FakeTauIdentifier::Data fFakeTauData;
  bIsGenuineTau = false; // generic variables. Could be renamed as "bIsGenuineBkg"
  fElectronData = ElectronSelection::Data();
  fMuonData     = MuonSelection::Data();
  fJetData      = JetSelection::Data();
  fCollinearAngularCutsData = AngularCutsBackToBack::Data();
  fBJetData = BJetSelection::Data();
  fMETData  = METSelection::Data();
  // fQGLRData = QuarkGluonLikelihoodRatio::Data();
  fTopData  = TopSelectionMVA::Data();
  fBackToBackAngularCutsData = AngularCutsCollinear::Data();
  fDnnSelectionData = DnnSelection::Data();
  // fFatJetData = FatJetSelection::Data();
  // fFatJetSoftDropData = FatJetSoftDropSelection::Data();
  fHistoSplitter.initialize();
  
  for (auto& p: fBaseObjects) p->reset();
  return;
}

//===== unique filling methods (to be called inside the event selection routine only)
void CommonPlots::fillControlPlotsAtVertexSelection(const Event& event) {
  for (auto& p: fBaseObjects) {
    p->fillControlPlotsAtVertexSelection(event);
  } 
}

void CommonPlots::fillControlPlotsAtElectronSelection(const Event& event, const ElectronSelection::Data& data) {
  fElectronData = data;
  for (auto& p: fBaseObjects) {
    p->fillControlPlotsAtElectronSelection(event, data);
  }
}

void CommonPlots::fillControlPlotsAtMuonSelection(const Event& event, const MuonSelection::Data& data) {
  fMuonData = data;
  for (auto& p: fBaseObjects) {
    p->fillControlPlotsAtMuonSelection(event, data);
  }
}

void CommonPlots::fillControlPlotsAtTauSelection(const Event& event, const TauSelection::Data& data) {
  fTauData = data;
  for (auto& p: fBaseObjects) {
    p->fillControlPlotsAtTauSelection(event, data);
  }
}

void CommonPlots::fillControlPlotsAtJetSelection(const Event& event, const JetSelection::Data& data) {
  fJetData = data;
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlNjets, bIsGenuineTau, fJetData.getNumberOfSelectedJets());
  for (auto& p: fBaseObjects) {
    p->fillControlPlotsAtJetSelection(event, data);
  }
}

void CommonPlots::fillControlPlotsAtAngularCutsCollinear(const Event& event, const AngularCutsCollinear::Data& data) {
  fCollinearAngularCutsData = data;
  // if std::cout << "bIsGenuineTau= " << bIsGenuineTau << ", fCollinearAngularCutsData.getMinimumCutValue() = " << fCollinearAngularCutsData.getMinimumCutValue() 
  // 	       << " fCollinearAngularCutsData.get1DCutVariable(0) = " << fCollinearAngularCutsData.get1DCutVariable(0) 
  // 	       << " fCollinearAngularCutsData.get1DCutVariable(1) = " << fCollinearAngularCutsData.get1DCutVariable(1) 
  // 	       << " fCollinearAngularCutsData.get1DCutVariable(2) = " << fCollinearAngularCutsData.get1DCutVariable(2) 
  // 	       << " fCollinearAngularCutsData.get1DCutVariable(3) = " << fCollinearAngularCutsData.get1DCutVariable(3) 
  // 	       << std::endl;
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlCollinearAngularCutsMinimum, bIsGenuineTau, fCollinearAngularCutsData.getMinimumCutValue()); 
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlCollinearAngularCutsJet1, bIsGenuineTau, fCollinearAngularCutsData.get1DCutVariable(0));
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlCollinearAngularCutsJet2, bIsGenuineTau, fCollinearAngularCutsData.get1DCutVariable(1));
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlCollinearAngularCutsJet3, bIsGenuineTau, fCollinearAngularCutsData.get1DCutVariable(2));
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlCollinearAngularCutsJet4, bIsGenuineTau, fCollinearAngularCutsData.get1DCutVariable(3));

  if (fAnalysisType == kHplus2hwAnalysisWithTop)
    {
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiTaus         , bIsGenuineTau, fCollinearAngularCutsData.getDeltaPhiTaus() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiTauMET       , bIsGenuineTau, fCollinearAngularCutsData.getDeltaPhiTauMET() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiLdgTauMET    , bIsGenuineTau, fCollinearAngularCutsData.getDeltaPhiLdgTauMET() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiSubldgTauMET , bIsGenuineTau, fCollinearAngularCutsData.getDeltaPhiSubldgTauMET() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiLdgJetMET    , bIsGenuineTau, fCollinearAngularCutsData.getDeltaPhiLdgJetMET() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiSubldgJetMET , bIsGenuineTau, fCollinearAngularCutsData.getDeltaPhiSubldgJetMET() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiMuonMET      , bIsGenuineTau, fCollinearAngularCutsData.getDeltaPhiMuonMET() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiLdgTauMuon   , bIsGenuineTau, fCollinearAngularCutsData.getDeltaPhiLdgTauMuon() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiSubldgTauMuon, bIsGenuineTau, fCollinearAngularCutsData.getDeltaPhiSubldgTauMuon() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiSubldgTauMET  , bIsGenuineTau, 
					       fCollinearAngularCutsData.getDeltaPhiLdgTauMET(), fCollinearAngularCutsData.getDeltaPhiSubldgTauMET() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiMuonMET       , bIsGenuineTau, 
					       fCollinearAngularCutsData.getDeltaPhiLdgTauMET(), fCollinearAngularCutsData.getDeltaPhiMuonMET() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiLdgBJetMET    , bIsGenuineTau, 
					       fCollinearAngularCutsData.getDeltaPhiLdgTauMET(), fCollinearAngularCutsData.getDeltaPhiLdgBJetMET() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiLdgTauMuon_Vs_DeltaPhiSubldgTauMuon, bIsGenuineTau, 
					       fCollinearAngularCutsData.getDeltaPhiLdgTauMuon(), fCollinearAngularCutsData.getDeltaPhiSubldgTauMuon() );
    }

  for (auto& p: fBaseObjects) 
    {
      p->fillControlPlotsAtAngularCutsCollinear(event, data);
    }
}

void CommonPlots::fillControlPlotsAtBtagging(const Event& event, const BJetSelection::Data& data) {
  fBJetData = data;
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlNBJets, bIsGenuineTau, fBJetData.getNumberOfSelectedBJets());
  for (auto& p: fJetData.getSelectedJets()) {
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlBDiscriminator, bIsGenuineTau, p.bjetDiscriminator());
  }
  for (auto& p: fBJetData.getSelectedBJets()) {
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlBJetPt, bIsGenuineTau, p.pt());
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlBJetEta, bIsGenuineTau, p.eta());
  }
  for (auto& p: fBaseObjects) {
    p->fillControlPlotsAtBtagging(event, data);
  }
}

void CommonPlots::fillControlPlotsAtMETSelection(const Event& event, const METSelection::Data& data) {
  fMETData = data;
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlMET, bIsGenuineTau, fMETData.getMET().R());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlMETPhi, bIsGenuineTau, fMETData.getMET().phi());
  for (auto& p: fBaseObjects) {
    p->fillControlPlotsAtMETSelection(event, data);
  }
}

void CommonPlots::fillControlPlotsAtAngularCutsBackToBack(const Event& event, const AngularCutsBackToBack::Data& data) {
  fBackToBackAngularCutsData = data;
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlBackToBackAngularCutsMinimum, bIsGenuineTau, fBackToBackAngularCutsData.getMinimumCutValue());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlBackToBackAngularCutsJet1, bIsGenuineTau, fBackToBackAngularCutsData.get1DCutVariable(0));
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlBackToBackAngularCutsJet2, bIsGenuineTau, fBackToBackAngularCutsData.get1DCutVariable(1));
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlBackToBackAngularCutsJet3, bIsGenuineTau, fBackToBackAngularCutsData.get1DCutVariable(2));
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlBackToBackAngularCutsJet4, bIsGenuineTau, fBackToBackAngularCutsData.get1DCutVariable(3));

  // All histograms below are already filled in CommonPlots::fillControlPlotsAtAngularCutsCollinear() - fixme (decide on correct implementation later)
  // ==> Comment out to avoid double filling
  // if (fAnalysisType == kHplus2hwAnalysisWithTop)
  //   {
  //     fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiTaus         , bIsGenuineTau, fBackToBackAngularCutsData.getDeltaPhiTaus() );
  //     fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiTauMET       , bIsGenuineTau, fBackToBackAngularCutsData.getDeltaPhiTauMET() );
  //     fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiLdgTauMET    , bIsGenuineTau, fBackToBackAngularCutsData.getDeltaPhiLdgTauMET() );
  //     fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiSubldgTauMET , bIsGenuineTau, fBackToBackAngularCutsData.getDeltaPhiSubldgTauMET() );
  //     fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiLdgJetMET    , bIsGenuineTau, fBackToBackAngularCutsData.getDeltaPhiLdgJetMET() );
  //     fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiSubldgJetMET , bIsGenuineTau, fBackToBackAngularCutsData.getDeltaPhiSubldgJetMET() );
  //     fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiMuonMET      , bIsGenuineTau, fBackToBackAngularCutsData.getDeltaPhiMuonMET() );
  //     fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiLdgTauMuon   , bIsGenuineTau, fBackToBackAngularCutsData.getDeltaPhiLdgTauMuon() );
  //     fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiSubldgTauMuon, bIsGenuineTau, fBackToBackAngularCutsData.getDeltaPhiSubldgTauMuon() );
  //     fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiSubldgTauMET  , bIsGenuineTau,
  // 					       fBackToBackAngularCutsData.getDeltaPhiLdgTauMET(), fBackToBackAngularCutsData.getDeltaPhiSubldgTauMET() );
  //     fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiMuonMET       , bIsGenuineTau, 
  // 					       fBackToBackAngularCutsData.getDeltaPhiLdgTauMET(), fBackToBackAngularCutsData.getDeltaPhiMuonMET() );
  //     fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiLdgBJetMET    , bIsGenuineTau, 
  // 					       fBackToBackAngularCutsData.getDeltaPhiLdgTauMET(), fBackToBackAngularCutsData.getDeltaPhiLdgBJetMET() );
  //     fHistoSplitter.fillShapeHistogramTriplet(hCtrlAngularCutsDeltaPhiLdgTauMuon_Vs_DeltaPhiSubldgTauMuon, bIsGenuineTau, 
  // 					       fBackToBackAngularCutsData.getDeltaPhiLdgTauMuon(), fBackToBackAngularCutsData.getDeltaPhiSubldgTauMuon() );
  //   }


  for (auto& p: fBaseObjects) 
    {
    p->fillControlPlotsAtAngularCutsBackToBack(event, data);
  }
  
  return;
}

void CommonPlots::fillControlPlotsAtDnnSelection(const Event& event, const DnnSelection::Data& data){
  fDnnSelectionData = data;
  float myTransverseMass = data.getMt();

  fHistoSplitter.fillShapeHistogramTriplet(hCtrlDnnSelection, bIsGenuineTau, data.getDnnOutput());
  if(data.passedLooseSelection()){
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlDnnSelectionLoose, bIsGenuineTau, myTransverseMass);
  }
  if(data.passedMediumSelection()){
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlDnnSelectionMedium, bIsGenuineTau, myTransverseMass);
  }
  if(data.passedTightSelection()){
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlDnnSelectionTight, bIsGenuineTau, myTransverseMass);
  } 
}

//===== unique filling methods (to be called AFTER return statement from analysis routine)
void CommonPlots::fillControlPlotsAfterTrigger(const Event& event) {
  for (auto& p: fBaseObjects) {
    p->fillControlPlotsAfterTrigger(event);
  } 
}

void CommonPlots::fillControlPlotsAfterMETFilter(const Event& event) {
  for (auto& p: fBaseObjects) {
    p->fillControlPlotsAfterMETFilter(event);
  } 
}

void CommonPlots::fillControlPlotsAfterTauSelection(const Event& event, const TauSelection::Data& data) {
  // Code logic: if there is no identified tau (or anti-isolated tau for QCD), the code will for sure crash later
  // This piece of code is called from TauSelection, so there one cannot judge if things go right or not, 
  // that kind of check needs to be done in the analysis code (i.e. cut away event if tau selection is not passed)
  for (auto& p: fBaseObjects) {
    p->fillControlPlotsAfterTauSelection(event, data);
  }
  fTauData = data;
  if (event.isData()) {
    bIsGenuineTau = true;
    return;
  }
  if (usesAntiIsolatedTaus()) {
    if (data.hasAntiIsolatedTaus()) {
      if (data.getAntiIsolatedTaus().size()==1) {
//        bIsGenuineTau = data.getAntiIsolatedTauIsGenuineTau();
	bIsGenuineTau = data.getAntiIsolatedTaus()[0].isGenuineTau();
      } else if (data.getAntiIsolatedTaus().size() >= 2 && (data.getAntiIsolatedTaus()[0].isGenuineTau() || data.getAntiIsolatedTaus()[1].isGenuineTau())) {
	bIsGenuineTau = true;
      } else {
	bIsGenuineTau = false;
      }
    }
  } else {
    if (data.hasIdentifiedTaus()) {
//      bIsGenuineTau = data.isGenuineTau();
      if (data.getSelectedTaus().size()==1) {
        bIsGenuineTau = data.getSelectedTaus()[0].isGenuineTau();
      } else if (data.getSelectedTaus()[0].isGenuineTau() || data.getSelectedTaus()[1].isGenuineTau()) {
        bIsGenuineTau = true;
      } else {
	bIsGenuineTau = false;
      }
    }
  }
}

void CommonPlots::fillControlPlotsAfterBjetSelection(const Event& event, const BJetSelection::Data& data) {
  fBJetData = data;
  for (auto& p: fBaseObjects) p->fillControlPlotsAfterBjetSelection(event, fBJetData);
  // bIsGenuineB = fBJetData.isGenuineB();
  return;
}

void CommonPlots::fillControlPlotsAfterBtagSF(const Event& event, const JetSelection::Data& jetData, const BJetSelection::Data& bjetData) {
  fJetData = jetData;
  fBJetData = bjetData;
  // pT of all selected b jets
  for (auto& p: fBJetData.getSelectedBJets()) 
    {
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlBJetPtAfterBtagSF, bIsGenuineTau, p.pt() );
    }
  // pT of all selected jets
  for (auto& p: fJetData.getSelectedJets()) 
    {
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlJetPtAfterBtagSF, bIsGenuineTau, p.pt() );
    }
  // Jet multiplicity 
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlNJetsAfterBtagSF, bIsGenuineTau, fJetData.getSelectedJets().size());

  for (auto& p: fBaseObjects) p->fillControlPlotsAfterBtagSF(event, fJetData, fBJetData);
  return;
}

void CommonPlots::fillControlPlotsAfterAntiIsolatedTauSelection(const Event& event, const TauSelection::Data& data) {

  for (auto& p: fBaseObjects) {
//    std::cout << p->isNull() << "\n";
    p->fillControlPlotsAfterAntiIsolatedTauSelection(event, data);
  }
}

void CommonPlots::fillControlPlotsAfterMETTriggerScaleFactor(const Event& event) {
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlNjetsAfterJetSelectionAndMETSF, bIsGenuineTau, fJetData.getNumberOfSelectedJets());
}

void CommonPlots::fillControlPlotsAfterStandardSelections(const Event& event, 
							  const JetSelection::Data& jetData,
							  const BJetSelection::Data& bjetData,
							  const METSelection::Data& METData,
							  const TopSelectionMVA::Data& topData){
  fJetData      = jetData;
  fBJetData     = bjetData;
  fTopData      = topData;
  fMETData      = METData;
  // WARNING! Remember to set value to "bIsGenuineTau" in your analyzer through:
  // fCommonPlots.setGenuineBkgStatus(isGenuineTau);

  // Fill Histograms
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlNVerticesAfterStdSelections, bIsGenuineTau, iVertices);

  // Taus
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauPtAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTau().pt());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauEtaAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTau().eta());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauPhiAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTau().phi());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauEtaPhiAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTau().eta(), fTauData.getSelectedTau().phi());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauLdgTrkPtAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTau().lChTrkPt());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauDecayModeAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTau().decayMode());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauNProngsAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTau().nProngs());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauRtauAfterStdSelections, bIsGenuineTau, fTauData.getRtauOfSelectedTau());
  for (auto& p: fHelper.getTauSourceData(!event.isMC(), fTauData.getSelectedTau()))
    {
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauSourceAfterStdSelections, bIsGenuineTau, p);
    }

  if (fTauData.getSelectedTaus().size() > 1)
    {
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSubldgTauPtAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].pt());
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSubldgTauEtaAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].eta());
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSubldgTauPhiAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].phi());
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSubldgTauEtaPhiAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].eta(), fTauData.getSelectedTaus()[1].phi());
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSubldgTauLdgTrkPtAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].lChTrkPt());
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSubldgTauDecayModeAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].decayMode());
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSubldgTauNProngsAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].nProngs());
    }
  
  // Only available if collinearData are filled
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlDeltaPhiTauMetAfterStdSelections, bIsGenuineTau, fBackToBackAngularCutsData.getDeltaPhiTauMET());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlDeltaPhiMuonMetAfterStdSelections, bIsGenuineTau, fBackToBackAngularCutsData.getDeltaPhiMuonMET());
  // fHistoSplitter.fillShapeHistogramTriplet(hCtrlDeltaPhiTauMetAfterStdSelections, bIsGenuineTau, fCollinearAngularCutsData.getDeltaPhiTauMET());
  // fHistoSplitter.fillShapeHistogramTriplet(hCtrlDeltaPhiTauMetAfterStdSelections, bIsGenuineTau, fCollinearAngularCutsData.getDeltaPhiMuonMET());
  
  // Muons  
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedMuonPtAfterStdSelections, bIsGenuineTau, fMuonData.getHighestSelectedMuonPt());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedMuonEtaAfterStdSelections, bIsGenuineTau, fMuonData.getHighestSelectedMuonEta());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedMuonPhiAfterStdSelections, bIsGenuineTau, fMuonData.getHighestSelectedMuonPhi());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedMuonEtaPhiAfterStdSelections, bIsGenuineTau, fMuonData.getHighestSelectedMuonEta(), fMuonData.getHighestSelectedMuonPhi());

  // Hadronic jets  
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlNJetsAfterStdSelections, bIsGenuineTau, fJetData.getNumberOfSelectedJets());
  for (auto& p: fJetData.getSelectedJets()) {
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlJetPtAfterStdSelections, bIsGenuineTau, p.pt());
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlJetEtaAfterStdSelections, bIsGenuineTau, p.eta());
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlJetEtaPhiAfterStdSelections, bIsGenuineTau, p.eta(), p.phi());
  }
  // For-loop: All selected jets
  for (auto& p: fJetData.getSelectedJets()) 
    {
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlJetPtAfterStdSelections         , bIsGenuineTau, p.pt() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlJetEtaAfterStdSelections        , bIsGenuineTau, p.eta() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlJetEtaPhiAfterStdSelections     , bIsGenuineTau, p.eta(), p.phi() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlBDiscriminatorAfterStdSelections, bIsGenuineTau, p.bjetDiscriminator() );
    }
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlHTAfterStdSelections , bIsGenuineTau, fJetData.HT());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlMHTAfterStdSelections, bIsGenuineTau, std::sqrt(fJetData.MHT().perp2()));

  // For-loop: All selected bjets
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlNBJetsAfterStdSelections, bIsGenuineTau, fBJetData.getSelectedBJets().size());
  for (auto& p: fBJetData.getSelectedBJets()) 
    {
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlBJetPtAfterStdSelections        , bIsGenuineTau, p.pt() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlBJetEtaAfterStdSelections       , bIsGenuineTau, p.eta() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlBDiscriminatorAfterStdSelections, bIsGenuineTau, p.bjetDiscriminator() );
    }

  // MET
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlMETAfterStdSelections   , bIsGenuineTau, fMETData.getMET().R() );
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlMETPhiAfterStdSelections, bIsGenuineTau, fMETData.getMET().Phi() );

  // QGLR
  // fHistoSplitter.fillShapeHistogramTriplet(hCtrlQGLRAfterStdSelections, bIsGenuineTau, fQGLRData.getQGLR());
  // fHistoSplitter.fillShapeHistogramTriplet(hCtrlQGLRNJetsAfterStdSelections,bIsGenuineTau, fQGLRData.getNumberOfJetsForQGLR());
  // fHistoSplitter.fillShapeHistogramTriplet(hCtrlQGLRNLightJetsAfterStdSelections,bIsGenuineTau, fQGLRData.getNumberOfLightJets());
  // fHistoSplitter.fillShapeHistogramTriplet(hCtrlQGLRNGluonJetsAfterStdSelections,bIsGenuineTau, fQGLRData.getNumberOfGluonJets());					   

  // TopSelection histograms
  if (fAnalysisType == kHplus2hwAnalysisWithTop)
    {
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlNTopsAfterStdSelections, bIsGenuineTau, fTopData.getSelectedCleanedTopsMVA().size() );
      
      if (fTopData.getAllCleanedTopsSize() > 0)
	{
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopPtAfterStdSelections        , bIsGenuineTau, fTopData.getTop().pt() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopMVAAfterStdSelections       , bIsGenuineTau, fTopData.getTopMVA() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopDijetPtAfterStdSelections   , bIsGenuineTau, fTopData.getTopDijet().pt() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopDijetMassAfterStdSelections , bIsGenuineTau, fTopData.getTopDijet().mass() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopMassAfterStdSelections      , bIsGenuineTau, fTopData.getTop().mass() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopMassWMassRatioAfterStdSelections, bIsGenuineTau, fTopData.getTopMassWMassRatio() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopPt_Vs_TopDijetPtAfterStdSelections, bIsGenuineTau, fTopData.getTop().pt(), fTopData.getTopDijet().pt());
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopBJetPtAfterStdSelections    , bIsGenuineTau, fTopData.getTopBJet().p4().pt() ); 
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopBJetEtaAfterStdSelections   , bIsGenuineTau, fTopData.getTopBJet().p4().eta() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopBJetBdiscAfterStdSelections , bIsGenuineTau, fTopData.getTopBJet().bjetDiscriminator() );
	}
    }

  return;
}

void CommonPlots::fillControlPlotsAfterAllSelections(const Event& event, int isGenuineB) {
  // NB: Call only after fillControlPlotsAfterStandardSelections() has been called
  // Variables fJetData, fBJetData, fQGLRData, fTopData, fMETData, bIsGenuineB already set!
  
  // Fill Histogram Triplets
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlNVerticesAfterAllSelections, isGenuineB, iVertices);
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlNVerticesAfterAllSelections, isGenuineB, iVertices);
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlNJetsAfterAllSelections    , isGenuineB, fJetData.getSelectedJets().size());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlNBJetsAfterAllSelections   , isGenuineB, fBJetData.getSelectedBJets().size());

  // For-loop: All selected jets
  for (auto& p: fJetData.getSelectedJets()) 
    {
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlJetPtAfterAllSelections         , isGenuineB, p.pt() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlJetEtaAfterAllSelections        , isGenuineB, p.eta() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlJetEtaPhiAfterAllSelections     , isGenuineB, p.eta(), p.phi() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlBDiscriminatorAfterAllSelections, isGenuineB, p.bjetDiscriminator() );
    }
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlHTAfterAllSelections , isGenuineB, fJetData.HT());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlMHTAfterAllSelections, isGenuineB, std::sqrt(fJetData.MHT().perp2()));
    

  // For-loop: All selected bjets
  for (auto& p: fBJetData.getSelectedBJets()) 
    {
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlBJetPtAfterAllSelections        , isGenuineB, p.pt() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlBJetEtaAfterAllSelections       , isGenuineB, p.eta() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlBDiscriminatorAfterAllSelections, isGenuineB, p.bjetDiscriminator() );
    }

  // MET
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlMETAfterAllSelections   , isGenuineB, fMETData.getMET().R() );
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlMETPhiAfterAllSelections, isGenuineB, fMETData.getMET().Phi() );

  // Leading jets
  unsigned int index = -1;
  for (auto jet: fJetData.getSelectedJets())
    {
      index++;
      if (index > 6) break;
      if (index == 0)
	{
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlJet1PtAfterAllSelections, isGenuineB, jet.pt() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlJet1EtaAfterAllSelections, isGenuineB, jet.eta() );
	}

      if (index == 1)
	{
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlJet2PtAfterAllSelections, isGenuineB, jet.pt() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlJet2EtaAfterAllSelections, isGenuineB, jet.eta() );
	}

      if (index == 2)
	{
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlJet3PtAfterAllSelections, isGenuineB, jet.pt() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlJet3EtaAfterAllSelections, isGenuineB, jet.eta() );
	}

      if (index == 3)
	{
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlJet4PtAfterAllSelections, isGenuineB, jet.pt() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlJet4EtaAfterAllSelections, isGenuineB, jet.eta() );
	}

      if (index == 4)
	{
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlJet5PtAfterAllSelections, isGenuineB, jet.pt() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlJet5EtaAfterAllSelections, isGenuineB, jet.eta() );
	}

      if (index == 5)
	{
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlJet6PtAfterAllSelections, isGenuineB, jet.pt() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlJet6EtaAfterAllSelections, isGenuineB, jet.eta() );
	}

      if (index == 6)
	{
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlJet7PtAfterAllSelections, isGenuineB, jet.pt() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlJet7EtaAfterAllSelections, isGenuineB, jet.eta() );
	}
    }

  // Leading b-jets
  index = -1;
  for (auto jet: fBJetData.getSelectedBJets())
    {
      index++;
      if (index > 3) break;
      if (index == 0)
	{
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlBJet1PtAfterAllSelections, isGenuineB, jet.pt() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlBJet1EtaAfterAllSelections, isGenuineB, jet.eta() );
	}

      if (index == 1)
	{
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlBJet2PtAfterAllSelections, isGenuineB, jet.pt() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlBJet2EtaAfterAllSelections, isGenuineB, jet.eta() );
	}

      if (index == 2)
	{
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlBJet3PtAfterAllSelections, isGenuineB, jet.pt() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlBJet3EtaAfterAllSelections, isGenuineB, jet.eta() );
	}

      if (index == 3)
	{
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlBJet4PtAfterAllSelections, isGenuineB, jet.pt() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlBJet4EtaAfterAllSelections, isGenuineB, jet.eta() );
	}
    }

  // QGLR
  // fHistoSplitter.fillShapeHistogramTriplet(hCtrlQGLRAfterAllSelections, isGenuineB, fQGLRData.getQGLR());
  // fHistoSplitter.fillShapeHistogramTriplet(hCtrlQGLRNJetsAfterAllSelections,isGenuineB, fQGLRData.getNumberOfJetsForQGLR());
  // fHistoSplitter.fillShapeHistogramTriplet(hCtrlQGLRNLightJetsAfterAllSelections,isGenuineB, fQGLRData.getNumberOfLightJets());
  // fHistoSplitter.fillShapeHistogramTriplet(hCtrlQGLRNGluonJetsAfterAllSelections,isGenuineB, fQGLRData.getNumberOfGluonJets());

  // TopSelection histograms
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlNTopsAfterAllSelections, isGenuineB, fTopData.getSelectedCleanedTopsMVA().size() );

  if (fTopData.getAllCleanedTopsSize() > 0)
    {
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopPtAfterAllSelections        , isGenuineB, fTopData.getTop().pt() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopMVAAfterAllSelections       , isGenuineB, fTopData.getTopMVA() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopDijetPtAfterAllSelections   , isGenuineB, fTopData.getTopDijet().pt() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopDijetMassAfterAllSelections , isGenuineB, fTopData.getTopDijet().mass() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopMassAfterAllSelections      , isGenuineB, fTopData.getTop().mass() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopMassWMassRatioAfterAllSelections, isGenuineB, fTopData.getTopMassWMassRatio() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopPt_Vs_TopDijetPtAfterAllSelections, isGenuineB, fTopData.getTop().pt(), fTopData.getTopDijet().pt());
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopBJetPtAfterAllSelections    , isGenuineB, fTopData.getTopBJet().p4().pt() ); 
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopBJetEtaAfterAllSelections   , isGenuineB, fTopData.getTopBJet().p4().eta() );
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopBJetBdiscAfterAllSelections , isGenuineB, fTopData.getTopBJet().bjetDiscriminator());
    }
      return;
}


void CommonPlots::fillControlPlotsAfterTopologicalSelections(const Event& event, bool withoutTau, bool withMu) {
  // I.e. plots after standard selections
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlNVerticesAfterStdSelections, bIsGenuineTau, iVertices);

  if (withoutTau == false)
    {
      if (usesAntiIsolatedTaus()) {
	fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauPtAfterStdSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().pt());
	fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauEtaAfterStdSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().eta());
	fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauPhiAfterStdSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().phi());
	fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauEtaPhiAfterStdSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().eta(), fTauData.getAntiIsolatedTau().phi());
	fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauLdgTrkPtAfterStdSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().lChTrkPt());
	fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauDecayModeAfterStdSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().decayMode());
	fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauNProngsAfterStdSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().nProngs());
	fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauRtauAfterStdSelections, bIsGenuineTau, fTauData.getRtauOfAntiIsolatedTau());
	for (auto& p: fHelper.getTauSourceData(!event.isMC(), fTauData.getAntiIsolatedTau())) {
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauSourceAfterStdSelections, bIsGenuineTau, p);
	}
      } else {
	fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauPtAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTau().pt());
	fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauEtaAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTau().eta());
	fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauPhiAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTau().phi());
	fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauEtaPhiAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTau().eta(), fTauData.getSelectedTau().phi());
	fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauLdgTrkPtAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTau().lChTrkPt());
	fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauDecayModeAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTau().decayMode());
	fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauNProngsAfterStdSelections, bIsGenuineTau, fTauData.getSelectedTau().nProngs());
	fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauRtauAfterStdSelections, bIsGenuineTau, fTauData.getRtauOfSelectedTau());
	for (auto& p: fHelper.getTauSourceData(!event.isMC(), fTauData.getSelectedTau())) {
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauSourceAfterStdSelections, bIsGenuineTau, p);
	}
      }
    }// if (withoutTau == false)
  
  if (withMu) // muons for embedding studies
    {
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedMuonPtAfterStdSelections, bIsGenuineTau, fMuonData.getHighestSelectedMuonPt());
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedMuonEtaAfterStdSelections, bIsGenuineTau, fMuonData.getHighestSelectedMuonEta());
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedMuonPhiAfterStdSelections, bIsGenuineTau, fMuonData.getHighestSelectedMuonPhi());
      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedMuonEtaPhiAfterStdSelections, bIsGenuineTau, fMuonData.getHighestSelectedMuonEta(), fMuonData.getHighestSelectedMuonPhi());
    }

  fHistoSplitter.fillShapeHistogramTriplet(hCtrlNJetsAfterStdSelections, bIsGenuineTau, fJetData.getNumberOfSelectedJets());
  for (auto& p: fJetData.getSelectedJets()) {
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlJetPtAfterStdSelections, bIsGenuineTau, p.pt());
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlJetEtaAfterStdSelections, bIsGenuineTau, p.eta());
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlJetEtaPhiAfterStdSelections, bIsGenuineTau, p.eta(), p.phi());
  }

  return;
}

void CommonPlots::fillControlPlotsAfterAllSelections(const Event& event, bool withoutTau, int first, int second) {
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlNVerticesAfterAllSelections, bIsGenuineTau, iVertices);

  if (withoutTau == false)
    {
      if (usesAntiIsolatedTaus()) {
        if (fTauData.getAntiIsolatedTaus().size() >= 2 && fTauData.getSelectedTaus().size() == 0) {
          fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausPtAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTaus()[first].pt());
          fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausEtaAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTaus()[first].eta());
          fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausPhiAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTaus()[first].phi());
          fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausPtAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTaus()[second].pt());
          fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausEtaAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTaus()[second].eta());
          fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausPhiAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTaus()[second].phi());

          fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauPtAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().pt());
          fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauEtaAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().eta());
          fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauPhiAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().phi());
          fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauEtaPhiAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().eta(), fTauData.getAntiIsolatedTau().phi());
          fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauLdgTrkPtAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().lChTrkPt());
          fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauDecayModeAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().decayMode());
          fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauNProngsAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().nProngs());
          fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauRtauAfterAllSelections, bIsGenuineTau, fTauData.getRtauOfAntiIsolatedTau());
          for (auto& p: fHelper.getTauSourceData(!event.isMC(), fTauData.getAntiIsolatedTau())) {
            fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauSourceAfterAllSelections, bIsGenuineTau, p);
          }
          fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauIPxyAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().IPxy());

        }
	else
	  {
	    if (fTauData.getAntiIsolatedTaus().size() >= 1 && fTauData.getSelectedTaus().size() == 1)
	      {
		fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausPtAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTaus()[first].pt());
		fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausEtaAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTaus()[first].eta());
		fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausPhiAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTaus()[first].phi());
		fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausPtAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[0].pt());
		fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausEtaAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[0].eta());
		fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausPhiAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[0].phi());
	      }
	    else if (fTauData.getAntiIsolatedTaus().size() == 0 && fTauData.getSelectedTaus().size() == 1)
              {
                fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausPtAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[0].pt());
                fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausEtaAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[0].eta());
                fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausPhiAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[0].phi());
              }
            else if (fTauData.getAntiIsolatedTaus().size() >= 1 && fTauData.getSelectedTaus().size() == 0)
              {
	        fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauPtAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().pt());
	        fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauEtaAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().eta());
	        fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauPhiAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().phi());
	        fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauEtaPhiAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().eta(), fTauData.getAntiIsolatedTau().phi());
	        fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauLdgTrkPtAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().lChTrkPt());
	        fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauDecayModeAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().decayMode());
	        fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauNProngsAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().nProngs());
	        fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauRtauAfterAllSelections, bIsGenuineTau, fTauData.getRtauOfAntiIsolatedTau());
              }
//	    for (auto& p: fHelper.getTauSourceData(!event.isMC(), fTauData.getAntiIsolatedTau())) 
//	      {
//		fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauSourceAfterAllSelections, bIsGenuineTau, p);
//	      }
//	    fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauIPxyAfterAllSelections, bIsGenuineTau, fTauData.getAntiIsolatedTau().IPxy());
	  }
      } // antiIsolatedTaus
      else
	{
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausPtAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[0].pt());
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausEtaAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[0].eta());
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausPhiAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[0].phi());
	  
	  if(fTauData.getSelectedTaus().size() >= 2)
	    {
	      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausPtAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].pt());
	      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausEtaAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].eta());
	      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTausPhiAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].phi());

	      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSubldgTauPtAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].pt());
	      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSubldgTauEtaAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].eta());
	      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSubldgTauPhiAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].phi());
	      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSubldgTauEtaPhiAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].eta(), fTauData.getSelectedTaus()[1].phi());
	      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSubldgTauLdgTrkPtAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].lChTrkPt());
	      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSubldgTauDecayModeAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].decayMode());
	      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSubldgTauNProngsAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].nProngs());
	      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSubldgTauIPxyAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTaus()[1].IPxy());
	    }
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauPtAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTau().pt());
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauEtaAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTau().eta());
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauPhiAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTau().phi());
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauEtaPhiAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTau().eta(), fTauData.getSelectedTau().phi());
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauLdgTrkPtAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTau().lChTrkPt());
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauDecayModeAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTau().decayMode());
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauNProngsAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTau().nProngs());
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauRtauAfterAllSelections, bIsGenuineTau, fTauData.getRtauOfSelectedTau());
	  
	  for (auto& p: fHelper.getTauSourceData(!event.isMC(), fTauData.getSelectedTau())) 
	    {
	      fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauSourceAfterAllSelections, bIsGenuineTau, p);
	    }
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedTauIPxyAfterAllSelections, bIsGenuineTau, fTauData.getSelectedTau().IPxy());
	}
      
    } // if (withoutTau == false)
  
  // Muons
  if(fMuonData.hasIdentifiedMuons()){
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedMuonPtAfterAllSelections, bIsGenuineTau, fMuonData.getSelectedMuons()[0].pt());
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedMuonEtaAfterAllSelections, bIsGenuineTau, fMuonData.getSelectedMuons()[0].eta());
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedMuonPhiAfterAllSelections, bIsGenuineTau, fMuonData.getSelectedMuons()[0].phi());
  }

  // Electrons
  if(fElectronData.hasIdentifiedElectrons()){
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedElectronPtAfterAllSelections, bIsGenuineTau, fElectronData.getSelectedElectrons()[0].pt());
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedElectronEtaAfterAllSelections, bIsGenuineTau, fElectronData.getSelectedElectrons()[0].eta());
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlSelectedElectronPhiAfterAllSelections, bIsGenuineTau, fElectronData.getSelectedElectrons()[0].phi());
  }

  fHistoSplitter.fillShapeHistogramTriplet(hCtrlNJetsAfterAllSelections, bIsGenuineTau, fJetData.getNumberOfSelectedJets());
  for (auto& p: fJetData.getSelectedJets()) {
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlJetPtAfterAllSelections, bIsGenuineTau, p.pt());
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlJetEtaAfterAllSelections, bIsGenuineTau, p.eta());
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlJetEtaPhiAfterAllSelections, bIsGenuineTau, p.eta(), p.phi());
  }
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlHTAfterAllSelections, bIsGenuineTau, fJetData.HT());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlMHTAfterAllSelections, bIsGenuineTau, std::sqrt(fJetData.MHT().perp2()));
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlMinDeltaPhiJetMHTAfterAllSelections, bIsGenuineTau, fJetData.minDeltaPhiJetMHT());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlMaxDeltaPhiJetMHTAfterAllSelections, bIsGenuineTau, fJetData.maxDeltaPhiJetMHT());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlMinDeltaRJetMHTAfterAllSelections, bIsGenuineTau, fJetData.minDeltaRJetMHT());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlMinDeltaRReversedJetMHTAfterAllSelections, bIsGenuineTau, fJetData.minDeltaRReversedJetMHT());
  
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlCollinearAngularCutsMinimumAfterAllSelections, bIsGenuineTau, fCollinearAngularCutsData.getMinimumCutValue());

  fHistoSplitter.fillShapeHistogramTriplet(hCtrlMETAfterAllSelections, bIsGenuineTau, fMETData.getMET().R());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlMETPhiAfterAllSelections, bIsGenuineTau, fMETData.getMET().phi());
  
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlNBJetsAfterAllSelections, bIsGenuineTau, fBJetData.getNumberOfSelectedBJets());
  for (auto& p: fBJetData.getSelectedBJets()) {
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlBJetPtAfterAllSelections, bIsGenuineTau, p.pt());
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlBJetEtaAfterAllSelections, bIsGenuineTau, p.eta());
    fHistoSplitter.fillShapeHistogramTriplet(hCtrlBDiscriminatorAfterAllSelections, bIsGenuineTau, p.bjetDiscriminator());
  }
  
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlBackToBackAngularCutsMinimumAfterAllSelections, bIsGenuineTau, fBackToBackAngularCutsData.getMinimumCutValue());

  fHistoSplitter.fillShapeHistogramTriplet(hCtrlDnnSelectionAfterAllSelections, bIsGenuineTau, fDnnSelectionData.getDnnOutput());
  
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlDeltaPhiTauMetAfterAllSelections, bIsGenuineTau, fBackToBackAngularCutsData.getDeltaPhiTauMET());
  fHistoSplitter.fillShapeHistogramTriplet(hCtrlDeltaPhiMuonMetAfterAllSelections , bIsGenuineTau, fBackToBackAngularCutsData.getDeltaPhiMuonMET());

  double myTransverseMass = -1.0;


  if (withoutTau == false)
    {
      if (usesAntiIsolatedTaus()) {
        if (fTauData.getAntiIsolatedTaus().size()>0) {
	  if (fAnalysisType == kQCDMeasurement ) myTransverseMass = TransverseMass::reconstruct(fTauData.getAntiIsolatedTau(), fMETData.getMET());
        }
        if (fTauData.getAntiIsolatedTaus().size()>=2 && fTauData.getSelectedTaus().size() == 0) {
          if (fAnalysisType == kHplus2hwAnalysis || fAnalysisType == kQCDMeasurement_muon) myTransverseMass = TransverseMass::reconstruct(fTauData.getAntiIsolatedTaus()[first],fTauData.getAntiIsolatedTaus()[second],fMuonData.getSelectedMuons()[0], fMETData.getMET());
          if (fAnalysisType == kHplus2hw_ele_Analysis || fAnalysisType == kQCDMeasurement_ele) myTransverseMass = TransverseMass::reconstruct(fTauData.getAntiIsolatedTaus()[first],fTauData.getAntiIsolatedTaus()[second],fElectronData.getSelectedElectrons()[0], fMETData.getMET());
        }
        if (fTauData.getAntiIsolatedTaus().size()>=1 && fTauData.getSelectedTaus().size() == 1) {
          if (fAnalysisType == kHplus2hwAnalysis || fAnalysisType == kQCDMeasurement_muon) myTransverseMass = TransverseMass::reconstruct(fTauData.getSelectedTaus()[0],fTauData.getAntiIsolatedTaus()[first],fMuonData.getSelectedMuons()[0], fMETData.getMET());
          if (fAnalysisType == kHplus2hw_ele_Analysis || fAnalysisType == kQCDMeasurement_ele) myTransverseMass = TransverseMass::reconstruct(fTauData.getSelectedTaus()[0],fTauData.getAntiIsolatedTaus()[first],fElectronData.getSelectedElectrons()[0], fMETData.getMET());
        }
        if (fTauData.getSelectedTaus().size() == 0 && fTauData.getAntiIsolatedTaus().size() >= 1 && fMuonData.getSelectedMuons().size() == 1 && fMuonData.getAntiIsolatedMuons().size() >= 1) {
          if (fAnalysisType == kQCDMeasurement_mmt) myTransverseMass = TransverseMass::reconstruct(fTauData.getAntiIsolatedTaus()[first],fMuonData.getSelectedMuons()[0],fMuonData.getAntiIsolatedMuons()[0], fMETData.getMET());
        }
        if (fTauData.getAntiIsolatedTaus().size() >= 1 && fMuonData.getSelectedMuons().size() == 2) {
          if (fAnalysisType == kQCDMeasurement_mmt) myTransverseMass = TransverseMass::reconstruct(fTauData.getAntiIsolatedTaus()[first],fMuonData.getSelectedMuons()[0],fMuonData.getSelectedMuons()[1], fMETData.getMET());
        }
        if (fTauData.getAntiIsolatedTaus().size() >= 1) {
          if (fAnalysisType == kQCDMeasurement_eet) myTransverseMass = TransverseMass::reconstruct(fTauData.getAntiIsolatedTaus()[first],fElectronData.getSelectedElectrons()[0],fElectronData.getSelectedElectrons()[1], fMETData.getMET());
        }

//        myTransverseMass = TransverseMass::reconstruct(fLooseTauData.getSelectedTaus()[0],fLooseTauData.getSelectedTaus()[0],fMuonData.getSelectedMuons()[0], fMETData.getMET());

      } else {
       if ((fAnalysisType == kSignalAnalysis) || fAnalysisType == kQCDNormalizationSystematicsSignalRegion) myTransverseMass = TransverseMass::reconstruct(fTauData.getSelectedTau(), fMETData.getMET());
       if ( (fAnalysisType == kHplus2hwAnalysis) || (fAnalysisType == kHplus2hwAnalysisWithTop) )
	 {
	   if (fTauData.getSelectedTaus().size() > 1) myTransverseMass = TransverseMass::reconstruct(fTauData.getSelectedTaus()[0],fTauData.getSelectedTaus()[1],fMuonData.getSelectedMuons()[0], fMETData.getMET());
	   else myTransverseMass = -1.0;
	 }
       if (fAnalysisType == kHplus2hw_ele_Analysis) myTransverseMass = TransverseMass::reconstruct(fTauData.getSelectedTaus()[0],fTauData.getSelectedTaus()[1],fElectronData.getSelectedElectrons()[0], fMETData.getMET());
       if (fAnalysisType == kHplus2hwAnalysis_mmt) myTransverseMass = TransverseMass::reconstruct(fTauData.getSelectedTaus()[0],fMuonData.getSelectedMuons()[0],fMuonData.getSelectedMuons()[1], fMETData.getMET());
       if (fAnalysisType == kHplus2hwAnalysis_eet) myTransverseMass = TransverseMass::reconstruct(fTauData.getSelectedTaus()[0],fElectronData.getSelectedElectrons()[0],fElectronData.getSelectedElectrons()[1],fMETData.getMET());
      }

      // Create the up and down variation for tau ID shape
      // Could probably be done in a nicer way elsewhere...
      if(tauIDup && fTauData.hasIdentifiedTaus() &&fTauData.getSelectedTau().pt()>=200){
	fHistoSplitter.fillShapeHistogramTriplet(hShapeTransverseMass, bIsGenuineTau, myTransverseMass, (hShapeTransverseMass[0]->UnprotectedGetWeight()*(1.0+0.05*fTauData.getSelectedTau().pt()/1000.0)));
      }else if(tauIDdown && fTauData.hasIdentifiedTaus() && fTauData.getSelectedTau().pt()>=200){
	fHistoSplitter.fillShapeHistogramTriplet(hShapeTransverseMass, bIsGenuineTau, myTransverseMass, (hShapeTransverseMass[0]->UnprotectedGetWeight()*(1.0-0.35*fTauData.getSelectedTau().pt()/1000.0)));
      }else{
//	fHistoSplitter.fillShapeHistogramTriplet(hShapeTransverseMass, bIsGenuineTau, myTransverseMass);
	fHistoSplitter.fillShapeHistogramTriplet(hShapeTransverseMass, bIsGenuineTau, myTransverseMass);
      }
    }// if (withoutTau == false)

  if (event.isData()) {
    hNSelectedVsRunNumber->Fill(event.eventID().run());
  }
  for (auto& p: fBaseObjects) {
    p->fillControlPlotsAfterAllSelections(event,withoutTau,first,second);
  }

  // TopSelection histograms
  if (fAnalysisType == kHplus2hwAnalysisWithTop)
    {

      fHistoSplitter.fillShapeHistogramTriplet(hCtrlNTopsAfterAllSelections, bIsGenuineTau, fTopData.getSelectedCleanedTopsMVA().size() );
      
      if (fTopData.getAllCleanedTopsSize() > 0)
	{
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopPtAfterAllSelections        , bIsGenuineTau, fTopData.getTop().pt() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopMVAAfterAllSelections       , bIsGenuineTau, fTopData.getTopMVA() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopDijetPtAfterAllSelections   , bIsGenuineTau, fTopData.getTopDijet().pt() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopDijetMassAfterAllSelections , bIsGenuineTau, fTopData.getTopDijet().mass() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopMassAfterAllSelections      , bIsGenuineTau, fTopData.getTop().mass() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopMassWMassRatioAfterAllSelections, bIsGenuineTau, fTopData.getTopMassWMassRatio() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopPt_Vs_TopDijetPtAfterAllSelections, bIsGenuineTau, fTopData.getTop().pt(), fTopData.getTopDijet().pt());
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopBJetPtAfterAllSelections    , bIsGenuineTau, fTopData.getTopBJet().p4().pt() ); 
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopBJetEtaAfterAllSelections   , bIsGenuineTau, fTopData.getTopBJet().p4().eta() );
	  fHistoSplitter.fillShapeHistogramTriplet(hCtrlTopBJetBdiscAfterAllSelections , bIsGenuineTau, fTopData.getTopBJet().bjetDiscriminator() );
	}
    }

  return;
}

void CommonPlots::fillControlPlotsAfterAllSelectionsWithProbabilisticBtag(const Event& event, const METSelection::Data& metData, double btagWeight) {
  double myTransverseMass = -1.0;
  if (usesAntiIsolatedTaus()) {
    myTransverseMass = TransverseMass::reconstruct(fTauData.getAntiIsolatedTau(), metData.getMET());
  } else {
    myTransverseMass = TransverseMass::reconstruct(fTauData.getSelectedTau(), metData.getMET());
  }
  fHistoSplitter.fillShapeHistogramTriplet(hShapeProbabilisticBtagTransverseMass, bIsGenuineTau, myTransverseMass, btagWeight);
  for (auto& p: fBaseObjects) {
    p->fillControlPlotsAfterAllSelectionsWithProbabilisticBtag(event, metData, btagWeight);
  }
}

//===== Filling of control plots for determining QCD shape uncertainty
void CommonPlots::fillControlPlotsForQCDShapeUncertainty(const Event& event,
                                                         const AngularCutsBackToBack::Data& collinearAngularCutsData,
                                                         const BJetSelection::Data& bJetData,
                                                         const METSelection::Data& metData,
                                                         const AngularCutsCollinear::Data& backToBackAngularCutsData) {
  fillControlPlotsAfterTopologicalSelections(event);
  // Note that the following methods store the data object as members
  fillControlPlotsAtAngularCutsCollinear(event, collinearAngularCutsData);
  fillControlPlotsAtMETSelection(event, metData);
  fillControlPlotsAtBtagging(event, bJetData);
  fillControlPlotsAtAngularCutsBackToBack(event, backToBackAngularCutsData);
  // Fill plots after final selection
  fillControlPlotsAfterAllSelections(event);
}
