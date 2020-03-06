#include "EventSelection/interface/CommonPlots_lt.h"
#include "EventSelection/interface/TransverseMass.h"
#include "EventSelection/interface/PUDependencyPlots.h"
#include "DataFormat/interface/Event.h"

CommonPlots_lt::CommonPlots_lt(const ParameterSet& config, const AnalysisType type, HistoWrapper& histoWrapper)
  : fEnableGenuineTauHistograms(true), // Needed always for limits
    fAnalysisType(type),
    fHistoWrapper(histoWrapper),
    fHistoSplitter(config, histoWrapper),
    fVerticesBins(config.getParameter<ParameterSet>("nVerticesBins")),
    fPtBins(config.getParameter<ParameterSet>("ptBins")),
    fEtaBins(config.getParameter<ParameterSet>("etaBins")),
    fPhiBins(config.getParameter<ParameterSet>("phiBins")),
    fDeltaEtaBins(config.getParameter<ParameterSet>("deltaEtaBins")),
    fDeltaPhiBins(config.getParameter<ParameterSet>("deltaPhiBins")),
    fDeltaRBins(config.getParameter<ParameterSet>("deltaRBins")),
    fNjetsBins(config.getParameter<ParameterSet>("njetsBins")),
    fMetBins(config.getParameter<ParameterSet>("metBins")),
    fHtBins(config.getParameter<ParameterSet>("htBins")),
    fBjetDiscrBins(config.getParameter<ParameterSet>("bjetDiscrBins")),
    fMtBins(config.getParameter<ParameterSet>("mtBins")),
    fInvmassBins(config.getParameter<ParameterSet>("invMassBins")),
    hNSelectedVsRunNumber(nullptr)
{ 
  // Create CommonPlotsBase objects
  bool enableStatus = config.getParameter<bool>("enablePUDependencyPlots");
  if (fAnalysisType == kHToHW_wSystematics) enableStatus = false;
  fPUDependencyPlots = new PUDependencyPlots(histoWrapper, enableStatus, fVerticesBins);
  fBaseObjects.push_back(fPUDependencyPlots);
}

CommonPlots_lt::~CommonPlots_lt() 
{
  // AtJetSelection
  fHistoSplitter.deleteHistograms(hJetN_AtJetSelection);
  fHistoSplitter.deleteHistograms(hJetPt_AtJetSelection);
  fHistoSplitter.deleteHistograms(hJetEta_AtJetSelection);
  // AtBtagging
  fHistoSplitter.deleteHistograms(hBjetN_AtBtagging);
  fHistoSplitter.deleteHistograms(hBjetPt_AtBtagging);
  fHistoSplitter.deleteHistograms(hBjetEta_AtBtagging);
  fHistoSplitter.deleteHistograms(hBjetDiscr_AtBtagging);
  // AtMETSelection
  fHistoSplitter.deleteHistograms(hMET_AtMETSelection);
  fHistoSplitter.deleteHistograms(hMETPhi_AtMETSelection);
  // AfterBtagSF
  fHistoSplitter.deleteHistograms(hJetN_AfterBtagSF);
  fHistoSplitter.deleteHistograms(hJetPt_AfterBtagSF);
  fHistoSplitter.deleteHistograms(hBjetN_AfterBtagSF);
  fHistoSplitter.deleteHistograms(hBjetPt_AfterBtagSF);

  // AfterPreselections
  fHistoSplitter.deleteHistograms(hVerticesN_AfterPreselections);
  fHistoSplitter.deleteHistograms(hTauN_AfterPreselections);
  fHistoSplitter.deleteHistograms(hTauPt_AfterPreselections);
  fHistoSplitter.deleteHistograms(hTauEta_AfterPreselections);
  fHistoSplitter.deleteHistograms(hTauPhi_AfterPreselections);
  fHistoSplitter.deleteHistograms(hTauLdgTrkPt_AfterPreselections);
  fHistoSplitter.deleteHistograms(hTauDM_AfterPreselections);
  fHistoSplitter.deleteHistograms(hTauNProngs_AfterPreselections);
  fHistoSplitter.deleteHistograms(hTauSource_AfterPreselections);
  fHistoSplitter.deleteHistograms(hTauIPxy_AfterPreselections);
  fHistoSplitter.deleteHistograms(hMuonN_AfterPreselections);
  fHistoSplitter.deleteHistograms(hMuonPt_AfterPreselections);
  fHistoSplitter.deleteHistograms(hMuonPhi_AfterPreselections);
  fHistoSplitter.deleteHistograms(hJetN_AfterPreselections);
  fHistoSplitter.deleteHistograms(hJetPt_AfterPreselections);
  fHistoSplitter.deleteHistograms(hJetEta_AfterPreselections);
  fHistoSplitter.deleteHistograms(hBjetN_AfterPreselections);
  fHistoSplitter.deleteHistograms(hBjetPt_AfterPreselections);
  fHistoSplitter.deleteHistograms(hBjetEta_AfterPreselections);
  fHistoSplitter.deleteHistograms(hBjetDiscr_AfterPreselections);
  fHistoSplitter.deleteHistograms(hMET_AfterPreselections);
  fHistoSplitter.deleteHistograms(hMETPhi_AfterPreselections);
  fHistoSplitter.deleteHistograms(hDeltaPhiTauMet_AfterPreselections);
  fHistoSplitter.deleteHistograms(hDeltaPhiMuonMet_AfterPreselections);
  fHistoSplitter.deleteHistograms(hHT_AfterPreselections);
  fHistoSplitter.deleteHistograms(hMHT_AfterPreselections);
  fHistoSplitter.deleteHistograms(hTransverseMass_AfterPreselections);

  // AfterSelections
  fHistoSplitter.deleteHistograms(hVerticesN_AfterSelections);
  fHistoSplitter.deleteHistograms(hTauN_AfterSelections);
  fHistoSplitter.deleteHistograms(hTauPt_AfterSelections);
  fHistoSplitter.deleteHistograms(hTauEta_AfterSelections);
  fHistoSplitter.deleteHistograms(hTauPhi_AfterSelections);
  fHistoSplitter.deleteHistograms(hTauLdgTrkPt_AfterSelections);
  fHistoSplitter.deleteHistograms(hTauDM_AfterSelections);
  fHistoSplitter.deleteHistograms(hTauNProngs_AfterSelections);
  fHistoSplitter.deleteHistograms(hTauSource_AfterSelections);
  fHistoSplitter.deleteHistograms(hTauIPxy_AfterSelections);
  fHistoSplitter.deleteHistograms(hMuonN_AfterSelections);
  fHistoSplitter.deleteHistograms(hMuonPt_AfterSelections);
  fHistoSplitter.deleteHistograms(hMuonPhi_AfterSelections);
  fHistoSplitter.deleteHistograms(hJetN_AfterSelections);
  fHistoSplitter.deleteHistograms(hJetPt_AfterSelections);
  fHistoSplitter.deleteHistograms(hJetEta_AfterSelections);
  fHistoSplitter.deleteHistograms(hBjetN_AfterSelections);
  fHistoSplitter.deleteHistograms(hBjetPt_AfterSelections);
  fHistoSplitter.deleteHistograms(hBjetEta_AfterSelections);
  fHistoSplitter.deleteHistograms(hBjetDiscr_AfterSelections);
  fHistoSplitter.deleteHistograms(hMET_AfterSelections);
  fHistoSplitter.deleteHistograms(hMETPhi_AfterSelections);
  fHistoSplitter.deleteHistograms(hDeltaPhiTauMet_AfterSelections);
  fHistoSplitter.deleteHistograms(hDeltaPhiMuonMet_AfterSelections);
  fHistoSplitter.deleteHistograms(hHT_AfterSelections);
  fHistoSplitter.deleteHistograms(hMHT_AfterSelections);
  fHistoSplitter.deleteHistograms(hTransverseMass_AfterSelections);

  if (hNSelectedVsRunNumber != nullptr) delete hNSelectedVsRunNumber;
  for (auto p: fBaseObjects) delete p;

}

void CommonPlots_lt::book(TDirectory *dir, bool isData) 
{ 
  fHistoSplitter.bookHistograms(dir);
  
  // Create directories for data driven control plots
  std::string myLabel        = "ForDataDrivenCtrlPlots";
  std::string myFakeLabel    = myLabel + "FakeTau";    // "ForDataDrivenCtrlPlotsFakeTau";
  std::string myGenuineLabel = myLabel + "GenuineTau"; // "ForDataDrivenCtrlPlotsGenuineTau";
  
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

  tauIDup   = false;
  tauIDdown = false;

  if(str.std::string::find("TauIDSystPlus") != std::string::npos)
    {
    tauIDup=true;
    }
  if(str.std::string::find("TauIDSystMinus")!= std::string::npos)
    {
      tauIDdown=true;
    }

  if (fEnableGenuineTauHistograms) 
    {
      for (auto& p: myDirs3) myDirs.push_back(p);
    } 
  else
    {
      for (auto& p: myDirs2) myDirs.push_back(p);
    }
    

  //================================================================================================
  // AtJetSelection
  //================================================================================================
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs,
						   hJetN_AtJetSelection, "JetN_AtJetSelection", 
						   ";jet multiplicity;Events", fNjetsBins.bins(), fNjetsBins.min(), fNjetsBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hJetPt_AtJetSelection, "JetPt_AtJetSelection", 
						   ";jet p_{T} (GeV);Events", fPtBins.bins(), fPtBins.min(), fPtBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hJetEta_AtJetSelection, "JetEta_AtJetSelection", 
						   ";jet #eta;Events", fEtaBins.bins(), fEtaBins.min(), fEtaBins.max());

  //================================================================================================
  // AtBtagging
  //================================================================================================
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs,
						   hBjetN_AtBtagging, "BjetN_AtBtagging", 
						   ";b-jet multiplicity;Events", fNjetsBins.bins(), fNjetsBins.min(), fNjetsBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hBjetPt_AtBtagging, "BjetPt_AtBtagging", 
						   ";b-jet p_{T} (GeV);Events", fPtBins.bins(), fPtBins.min(), fPtBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hBjetEta_AtBtagging, "BjetEta_AtBtagging", 
						   ";jet #eta;Events", fEtaBins.bins(), fEtaBins.min(), fEtaBins.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hBjetDiscr_AtBtagging, "BjetDiscr_AtBtagging", 
						   ";b-jet discriminator;Events", fBjetDiscrBins.bins(), fBjetDiscrBins.min(), fBjetDiscrBins.max());

  //================================================================================================
  // AfterBtagSF
  //================================================================================================
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs,
						   hJetN_AfterBtagSF, "JetN_AfterBtagSF", 
						   ";jet multiplicity;Events", fNjetsBins.bins(), fNjetsBins.min(), fNjetsBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hJetPt_AfterBtagSF, "JetPt_AfterBtagSF", 
						   ";jet p_{T} (GeV);Events", fPtBins.bins(), fPtBins.min(), fPtBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs,
						   hBjetN_AfterBtagSF, "BjetN_AfterBtagSF", 
						   ";b-jet multiplicity;Events", fNjetsBins.bins(), fNjetsBins.min(), fNjetsBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hBjetPt_AfterBtagSF, "BjetPt_AfterBtagSF", 
						   ";b-jet p_{T} (GeV);Events", fPtBins.bins(), fPtBins.min(), fPtBins.max());
  
  //================================================================================================
  // AtMETSelection
  //================================================================================================
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hMET_AtMETSelection, "MET_AtMETSelection", 
						   ";E_{T}^{miss} (GeV);Events", fMetBins.bins(), fMetBins.min(), fMetBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hMETPhi_AtMETSelection, "METPhi_AtMETSelection", 
						   ";E_{T}^{miss} #phi (rads);Events", fPhiBins.bins(), fPhiBins.min(), fPhiBins.max());

  //================================================================================================
  // Preselections
  //================================================================================================
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hVerticesN_AfterPreselections, "VerticesN_AfterPreselections", 
						   ";vertex multiplicity;Events", fVerticesBins.bins(), fVerticesBins.min(), fVerticesBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs,
						   hTauN_AfterPreselections, "TauN_AfterPreselections",
						   ";#tau_{h} multiplicity;Events", fNjetsBins.bins(), fNjetsBins.min(), fNjetsBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs,
						   hTauPt_AfterPreselections, "TauPt_AfterPreselections",
						   ";#tau_{h} p_{T} (GeV);Events", fPtBins.bins(), fPtBins.min(), fPtBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hTauEta_AfterPreselections, "TauEta_AfterPreselections", 
						   ";#tau_{h} #eta;Events", fEtaBins.bins(), fEtaBins.min(), fEtaBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hTauPhi_AfterPreselections, "TauPhi_AfterPreselections", 
						   ";#tau_{h} #phi (rads);Events", fPhiBins.bins(), fPhiBins.min(), fPhiBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hTauLdgTrkPt_AfterPreselections, "TauLdgTrkPt_AfterPreselections", 
						   ";#tau_{h} ldg. trk p_{T} (GeV);Events", fPtBins.bins(), fPtBins.min(), fPtBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hTauDM_AfterPreselections, "TauDM_AfterPreselections",
						   ";#tau_{h} decay mode;Events", 20, 0, 20);

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hTauNProngs_AfterPreselections, "TauNProngs_AfterPreselections", 
						   ";#tau_{h} charged particle multiplicity;Events", 11, 0, 11);
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hTauSource_AfterPreselections, "TauSource_AfterPreselections", // fixme: needed? works?
						   ";#tau_{h} source;Events", fHelper.getTauSourceBinCount(), 0, fHelper.getTauSourceBinCount());
  // For-loop: All tau source bin labels (std::string)
  for (int i = 0; i < fHelper.getTauSourceBinCount(); ++i) 
    {
      fHistoSplitter.SetBinLabel(hTauSource_AfterPreselections, i+1, fHelper.getTauSourceBinLabel(i));
    }
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hTauIPxy_AfterPreselections, "TauIPxy_AfterPreselections", 
						   ";IP_{T} (cm);Events", 100, 0, 0.2);
    
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hMuonN_AfterPreselections, "MuonN_AfterPreselections", 
						   ";#mu multiplicity;Events", fNjetsBins.bins(), fNjetsBins.min(), fNjetsBins.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hMuonPt_AfterPreselections, "MuonPt_AfterPreselections", 
						   ";#mu p_{T} (GeV);Events", fPtBins.bins(), fPtBins.min(), fPtBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hMuonEta_AfterPreselections, "MuonEta_AfterPreselections", 
						   ";#mu #eta;Events", fEtaBins.bins(), fEtaBins.min(), fEtaBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hMuonPhi_AfterPreselections, "MuonPhi_AfterPreselections", 
						   ";#mu #phi (rads);Events", fPhiBins.bins(), fPhiBins.min(), fPhiBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs,
						   hJetN_AfterPreselections, "JetN_AfterPreselections", 
						   ";jet multiplicity;Events", fNjetsBins.bins(), fNjetsBins.min(), fNjetsBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hJetPt_AfterPreselections, "JetPt_AfterPreselections", 
						   ";jet p_{T} (GeV);Events", fPtBins.bins(), fPtBins.min(), fPtBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hJetEta_AfterPreselections, "JetEta_AfterPreselections", 
						   ";jet #eta;Events", fEtaBins.bins(), fEtaBins.min(), fEtaBins.max());
        
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hBjetN_AfterPreselections, "BjetN_AfterPreselections", 
						   ";b-jet multiplicity;Events", fNjetsBins.bins(), fNjetsBins.min(), fNjetsBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hBjetPt_AfterPreselections, "BjetPt_AfterPreselections", 
						   ";b-jet p_{T} (GeV);Events", fPtBins.bins(), fPtBins.min(), fPtBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hBjetEta_AfterPreselections, "BjetEta_AfterPreselections", 
						   ";b-jet #eta;Events", fEtaBins.bins(), fEtaBins.min(), fEtaBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hBjetDiscr_AfterPreselections, "BjetDiscr_AfterPreselections", 
						   ";b-jet discriminator;Events", fBjetDiscrBins.bins(), fBjetDiscrBins.min(), fBjetDiscrBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hMET_AfterPreselections, "MET_AfterPreselections", 
						   ";E_{T}^{miss} (GeV);Events", fMetBins.bins(), fMetBins.min(), fMetBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hMETPhi_AfterPreselections, "METPhi_AfterPreselections", 
						   ";E_{T}^{miss} #phi (rads);Events", fPhiBins.bins(), fPhiBins.min(), fPhiBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hDeltaPhiTauMet_AfterPreselections, "DeltaPhiTauMet_AfterPreselections", 
						   ";#Delta#phi(#tau_{h},E_{T}^{miss}) ({}^{#circ});Events", fDeltaPhiBins.bins(), fDeltaPhiBins.min(), fDeltaPhiBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hDeltaPhiMuonMet_AfterPreselections, "DeltaPhiMuonMet_AfterPreselections", 
						   ";#Delta#phi(#mu,E_{T}^{miss}) ({}^{#circ});Events", fDeltaPhiBins.bins(), fDeltaPhiBins.min(), fDeltaPhiBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hHT_AfterPreselections, "HT_AfterPreselections", 
						   ";H_{T} (GeV);Events", fHtBins.bins(), fHtBins.min(), fHtBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hMHT_AfterPreselections, "MHT_AfterPreselections", 
						   ";MHT (GeV);Events", fMetBins.bins(), fMetBins.min(), fMetBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(true, HistoLevel::kSystematics, myDirs3, 
						   hTransverseMass_AfterPreselections, "TransverseMass_AfterPreselections", 
						   ";m_{T}(#tau_{h}, #ell, j_{1}, j_{2}) (GeV);Events", fMtBins.bins(), fMtBins.min(), fMtBins.max());

  //================================================================================================
  // Selections
  //================================================================================================
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hVerticesN_AfterSelections, "VerticesN_AfterSelections", 
						   ";vertex multiplicity;Events", fVerticesBins.bins(), fVerticesBins.min(), fVerticesBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs,
						   hTauN_AfterSelections, "TauN_AfterSelections",
						   ";#tau_{h} multiplicity;Events", fNjetsBins.bins(), fNjetsBins.min(), fNjetsBins.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs,
						   hTauPt_AfterSelections, "TauPt_AfterSelections",
						   ";#tau_{h} p_{T} (GeV);Events", fPtBins.bins(), fPtBins.min(), fPtBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hTauEta_AfterSelections, "TauEta_AfterSelections", 
						   ";#tau_{h} #eta;Events", fEtaBins.bins(), fEtaBins.min(), fEtaBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hTauPhi_AfterSelections, "TauPhi_AfterSelections", 
						   ";#tau_{h} #phi (rads);Events", fPhiBins.bins(), fPhiBins.min(), fPhiBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hTauLdgTrkPt_AfterSelections, "TauLdgTrkPt_AfterSelections", 
						   ";#tau_{h} ldg. trk p_{T} (GeV);Events", fPtBins.bins(), fPtBins.min(), fPtBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hTauDM_AfterSelections, "TauDM_AfterSelections",
						   ";#tau_{h} decay mode;Events", 20, 0, 20);

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hTauNProngs_AfterSelections, "TauNProngs_AfterSelections", 
						   ";#tau_{h} charged particle multiplicity;Events", 11, 0, 11);
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hTauSource_AfterSelections, "TauSource_AfterSelections", // fixme: needed? works?
						   ";#tau_{h} source;Events", fHelper.getTauSourceBinCount(), 0, fHelper.getTauSourceBinCount());
  // For-loop: All tau source bin labels (std::string)
  for (int i = 0; i < fHelper.getTauSourceBinCount(); ++i) 
    {
      fHistoSplitter.SetBinLabel(hTauSource_AfterSelections, i+1, fHelper.getTauSourceBinLabel(i));
    }
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hTauIPxy_AfterSelections, "TauIPxy_AfterSelections", 
						   ";IP_{T} (cm);Events", 100, 0, 0.2);
    
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hMuonN_AfterSelections, "MuonN_AfterSelections", 
						   ";#mu multiplicity;Events", fNjetsBins.bins(), fNjetsBins.min(), fNjetsBins.max());

  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hMuonPt_AfterSelections, "MuonPt_AfterSelections", 
						   ";#mu p_{T} (GeV);Events", fPtBins.bins(), fPtBins.min(), fPtBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hMuonEta_AfterSelections, "MuonEta_AfterSelections", 
						   ";#mu #eta;Events", fEtaBins.bins(), fEtaBins.min(), fEtaBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hMuonPhi_AfterSelections, "MuonPhi_AfterSelections", 
						   ";#mu #phi (rads);Events", fPhiBins.bins(), fPhiBins.min(), fPhiBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs,
						   hJetN_AfterSelections, "JetN_AfterSelections", 
						   ";jet multiplicity;Events", fNjetsBins.bins(), fNjetsBins.min(), fNjetsBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hJetPt_AfterSelections, "JetPt_AfterSelections", 
						   ";jet p_{T} (GeV);Events", fPtBins.bins(), fPtBins.min(), fPtBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hJetEta_AfterSelections, "JetEta_AfterSelections", 
						   ";jet #eta;Events", fEtaBins.bins(), fEtaBins.min(), fEtaBins.max());
        
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hBjetN_AfterSelections, "BjetN_AfterSelections", 
						   ";b-jet multiplicity;Events", fNjetsBins.bins(), fNjetsBins.min(), fNjetsBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hBjetPt_AfterSelections, "BjetPt_AfterSelections", 
						   ";b-jet p_{T} (GeV);Events", fPtBins.bins(), fPtBins.min(), fPtBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hBjetEta_AfterSelections, "BjetEta_AfterSelections", 
						   ";b-jet #eta;Events", fEtaBins.bins(), fEtaBins.min(), fEtaBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hBjetDiscr_AfterSelections, "BjetDiscr_AfterSelections", 
						   ";b-jet discriminator;Events", fBjetDiscrBins.bins(), fBjetDiscrBins.min(), fBjetDiscrBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hMET_AfterSelections, "MET_AfterSelections", 
						   ";E_{T}^{miss} (GeV);Events", fMetBins.bins(), fMetBins.min(), fMetBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hMETPhi_AfterSelections, "METPhi_AfterSelections", 
						   ";E_{T}^{miss} #phi (rads);Events", fPhiBins.bins(), fPhiBins.min(), fPhiBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hDeltaPhiTauMet_AfterSelections, "DeltaPhiTauMet_AfterSelections", 
						   ";#Delta#phi(#tau_{h},E_{T}^{miss}) ({}^{#circ});Events", fDeltaPhiBins.bins(), fDeltaPhiBins.min(), fDeltaPhiBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hDeltaPhiMuonMet_AfterSelections, "DeltaPhiMuonMet_AfterSelections", 
						   ";#Delta#phi(#mu,E_{T}^{miss}) ({}^{#circ});Events", fDeltaPhiBins.bins(), fDeltaPhiBins.min(), fDeltaPhiBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hHT_AfterSelections, "HT_AfterSelections", 
						   ";H_{T} (GeV);Events", fHtBins.bins(), fHtBins.min(), fHtBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(fEnableGenuineTauHistograms, HistoLevel::kSystematics, myDirs, 
						   hMHT_AfterSelections, "MHT_AfterSelections", 
						   ";MHT (GeV);Events", fMetBins.bins(), fMetBins.min(), fMetBins.max());
  
  fHistoSplitter.createShapeHistogramTriplet<TH1F>(true, HistoLevel::kSystematics, myDirs3, 
						   hTransverseMass_AfterSelections, "TransverseMass_AfterSelections", 
						   ";m_{T}(#tau_{h}, #ell, j_{1}, j_{2});Events", fMtBins.bins(), fMtBins.min(), fMtBins.max());

  if (isData) 
    {
      hNSelectedVsRunNumber = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "NSelectedVsRunNumber",
							 "NSelectedVsRunNumber;Run number;Events", 
							 14000, 246000, 260000);
    }
 
  for (auto& p: fBaseObjects) p->book(dir, isData);
  return;
}

void CommonPlots_lt::initialize() 
{
  bIsGenuineTau = false;
  nVertices     = -1;
  nElectrons    = -1;
  nMuons        = -1;
  nTaus         = -1;
  nJets         = -1;
  nBJets        = -1;
  fElectronData = ElectronSelection::Data();
  fMuonData     = MuonSelection::Data();
  fTauData      = TauSelection::Data();
  fJetData      = JetSelection::Data();
  fBJetData     = BJetSelection::Data();
  fMETData      = METSelection::Data();
  // fTopData      = TopSelectionMVA::Data();
  // fBackToBackAngularCutsData = AngularCutsCollinear::Data();
  // fCollinearAngularCutsData  = AngularCutsBackToBack::Data();
  // fFatJetData                = FatJetSelection::Data();
  // fFatJetSoftDropData        = FatJetSoftDropSelection::Data();
  fHistoSplitter.initialize();
  
  for (auto& p: fBaseObjects) p->reset();

  return;
}

//===== unique filling methods (to be called inside the event selection routine only)
void CommonPlots_lt::fillControlPlotsAtVertexSelection(const Event& event) 
{
  for (auto& p: fBaseObjects) p->fillControlPlotsAtVertexSelection(event);

  return;
}

void CommonPlots_lt::fillControlPlotsAtElectronSelection(const Event& event, const ElectronSelection::Data& data) 
{
  fElectronData = data;

  for (auto& p: fBaseObjects) p->fillControlPlotsAtElectronSelection(event, data);

  return;
}

void CommonPlots_lt::fillControlPlotsAtMuonSelection(const Event& event, const MuonSelection::Data& data)
{
  fMuonData = data;

  for (auto& p: fBaseObjects) p->fillControlPlotsAtMuonSelection(event, data);

  return;
}

void CommonPlots_lt::fillControlPlotsAtTauSelection(const Event& event, const TauSelection::Data& data) 
{
  
  fTauData = data;

  for (auto& p: fBaseObjects) p->fillControlPlotsAtTauSelection(event, data);

  return;
}

void CommonPlots_lt::fillControlPlotsAtJetSelection(const Event& event, const JetSelection::Data& data)
{
  fJetData = data;

  fHistoSplitter.fillShapeHistogramTriplet(hJetN_AtJetSelection, bIsGenuineTau, fJetData.getNumberOfSelectedJets());
  // For-loop: All jets
  for (auto& p: fJetData.getSelectedJets()) 
    {
      fHistoSplitter.fillShapeHistogramTriplet(hJetPt_AtJetSelection , bIsGenuineTau, p.pt());
      fHistoSplitter.fillShapeHistogramTriplet(hJetEta_AtJetSelection, bIsGenuineTau, p.eta());
    }



  for (auto& p: fBaseObjects) p->fillControlPlotsAtJetSelection(event, data);

  return;
}

void CommonPlots_lt::fillControlPlotsAtBtagging(const Event& event, const BJetSelection::Data& data)
{
  fBJetData = data;

  fHistoSplitter.fillShapeHistogramTriplet(hBjetN_AtBtagging, bIsGenuineTau, fBJetData.getNumberOfSelectedBJets());
  // For-loop: All b-jets
  for (auto& p: fBJetData.getSelectedBJets()) 
    {
      fHistoSplitter.fillShapeHistogramTriplet(hBjetPt_AtBtagging   , bIsGenuineTau, p.pt());
      fHistoSplitter.fillShapeHistogramTriplet(hBjetEta_AtBtagging  , bIsGenuineTau, p.eta());
      fHistoSplitter.fillShapeHistogramTriplet(hBjetDiscr_AtBtagging, bIsGenuineTau, p.bjetDiscriminator());
    }

  for (auto& p: fBaseObjects) p->fillControlPlotsAtBtagging(event, data);

  return;
}

void CommonPlots_lt::fillControlPlotsAtMETSelection(const Event& event, const METSelection::Data& data) 
{
  fMETData = data;

  fHistoSplitter.fillShapeHistogramTriplet(hMET_AtMETSelection   , bIsGenuineTau, fMETData.getMET().R());
  fHistoSplitter.fillShapeHistogramTriplet(hMETPhi_AtMETSelection, bIsGenuineTau, fMETData.getMET().phi());
  
  for (auto& p: fBaseObjects) p->fillControlPlotsAtMETSelection(event, data);
  
  return;
}

void CommonPlots_lt::fillControlPlotsAfterTrigger(const Event& event) 
{
  for (auto& p: fBaseObjects) p->fillControlPlotsAfterTrigger(event);

  return;
}

void CommonPlots_lt::fillControlPlotsAfterMETFilter(const Event& event)
{
  for (auto& p: fBaseObjects) p->fillControlPlotsAfterMETFilter(event);

  return;
}

void CommonPlots_lt::fillControlPlotsAfterTauSelection(const Event& event, const TauSelection::Data& data) 
{
  // Code logic: if there is no identified tau (or anti-isolated tau for QCD), the code will for sure crash later
  // This piece of code is called from TauSelection, so there one cannot judge if things go right or not, 
  // that kind of check needs to be done in the analysis code (i.e. cut away event if tau selection is not passed)
  for (auto& p: fBaseObjects) p->fillControlPlotsAfterTauSelection(event, data);

  fTauData = data;
  if (event.isData()) 
    {
      bIsGenuineTau = true;
      return;
    }

  if (usesAntiIsolatedTaus()) 
    {

      if (data.hasAntiIsolatedTaus())
	{
	  if (data.getAntiIsolatedTaus().size()==1)
	    {
	      bIsGenuineTau = data.getAntiIsolatedTaus()[0].isGenuineTau();
	    } 
	  else if (data.getAntiIsolatedTaus().size() >= 2 && (data.getAntiIsolatedTaus()[0].isGenuineTau() || data.getAntiIsolatedTaus()[1].isGenuineTau())) 
	    {
	      bIsGenuineTau = true;
	    } 
	  else 
	    {
	      bIsGenuineTau = false;
	    }
	}
    } 
  else {
    if (data.hasIdentifiedTaus()) 
      {
	if (data.getSelectedTaus().size()==1) 
	  {
	    bIsGenuineTau = data.getSelectedTaus()[0].isGenuineTau();
	  } 
	else if (data.getSelectedTaus()[0].isGenuineTau() || data.getSelectedTaus()[1].isGenuineTau()) 
	  {
	    bIsGenuineTau = true;
	  } 
	else {
	  bIsGenuineTau = false;
	}
      }
  }
  
  return;
}

void CommonPlots_lt::fillControlPlotsAfterBjetSelection(const Event& event, const BJetSelection::Data& data) 
{
  fBJetData = data;

  for (auto& p: fBaseObjects) p->fillControlPlotsAfterBjetSelection(event, fBJetData);

  return;
}

void CommonPlots_lt::fillControlPlotsAfterBtagSF(const Event& event, const JetSelection::Data& jetData, const BJetSelection::Data& bjetData) 
{
  fJetData  = jetData;
  fBJetData = bjetData;

  // For-loop: Selected jets
  for (auto& p: fJetData.getSelectedJets()) 
    {
      fHistoSplitter.fillShapeHistogramTriplet(hJetPt_AfterBtagSF, bIsGenuineTau, p.pt() );
    }

  // For-loop: Selected b-jets
  for (auto& p: fBJetData.getSelectedBJets()) 
    {
      fHistoSplitter.fillShapeHistogramTriplet(hBjetPt_AfterBtagSF, bIsGenuineTau, p.pt() );
    }

  // Jet multiplicities
  fHistoSplitter.fillShapeHistogramTriplet(hJetN_AfterBtagSF, bIsGenuineTau, fJetData.getSelectedJets().size());
  fHistoSplitter.fillShapeHistogramTriplet(hBjetN_AfterBtagSF, bIsGenuineTau, fBJetData.getSelectedBJets().size());

  for (auto& p: fBaseObjects) p->fillControlPlotsAfterBtagSF(event, fJetData, fBJetData);

  return;
}

void CommonPlots_lt::fillControlPlotsAfterAntiIsolatedTauSelection(const Event& event, const TauSelection::Data& data) 
{
  for (auto& p: fBaseObjects) p->fillControlPlotsAfterAntiIsolatedTauSelection(event, data);

  return;
}

void CommonPlots_lt::fillControlPlotsAfterPreselections(const Event& event, 
							const JetSelection::Data& jetData,
							const BJetSelection::Data& bjetData,
							const METSelection::Data& METData)
{

  // NOTE: This assumes that fillControlPlotsfillControlPlotsAtMuonSelection() and fillControlPlotsAtMuonSelection()
  //       have already been called, and hence all relevant "Data" structures (fElectronData, fMuonData) have already been assigned

  // NOTE: Remember to set value to "bIsGenuineTau" in your analyzer through:
  // fCommonPlots.setGenuineBkgStatus(isGenuineTau);


  // Assign data structures
  fJetData      = jetData;
  fBJetData     = bjetData;
  fMETData      = METData;

  // Assign object multiplicities
  nElectrons = fElectronData.getSelectedElectrons().size();
  nMuons     = fMuonData.getSelectedMuons().size();
  nTaus      = fTauData.getSelectedTaus().size();
  nJets      = fJetData.getNumberOfSelectedJets();
  nBJets     = fBJetData.getSelectedBJets().size();

  // Fill Histograms
  fHistoSplitter.fillShapeHistogramTriplet(hVerticesN_AfterPreselections, bIsGenuineTau, nVertices);

  // Muons    
  fHistoSplitter.fillShapeHistogramTriplet(hMuonN_AfterPreselections, bIsGenuineTau, nMuons);
  if (nMuons > 0)
    {
      fHistoSplitter.fillShapeHistogramTriplet(hMuonPt_AfterPreselections , bIsGenuineTau, fMuonData.getHighestSelectedMuonPt());
      fHistoSplitter.fillShapeHistogramTriplet(hMuonEta_AfterPreselections, bIsGenuineTau, fMuonData.getHighestSelectedMuonEta());
      fHistoSplitter.fillShapeHistogramTriplet(hMuonPhi_AfterPreselections, bIsGenuineTau, fMuonData.getHighestSelectedMuonPhi());

      double deltaPhiMuonMet = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(fMuonData.getSelectedMuon().p4(), fMETData.getMET())*57.29578);
      fHistoSplitter.fillShapeHistogramTriplet(hDeltaPhiMuonMet_AfterPreselections, bIsGenuineTau, deltaPhiMuonMet);
    }


  // Taus
  fHistoSplitter.fillShapeHistogramTriplet(hTauN_AfterPreselections, bIsGenuineTau, nTaus);
  if (nTaus > 0)
    {
      fHistoSplitter.fillShapeHistogramTriplet(hTauPt_AfterPreselections      , bIsGenuineTau, fTauData.getSelectedTau().pt());
      fHistoSplitter.fillShapeHistogramTriplet(hTauEta_AfterPreselections     , bIsGenuineTau, fTauData.getSelectedTau().eta());
      fHistoSplitter.fillShapeHistogramTriplet(hTauPhi_AfterPreselections     , bIsGenuineTau, fTauData.getSelectedTau().phi());
      fHistoSplitter.fillShapeHistogramTriplet(hTauLdgTrkPt_AfterPreselections, bIsGenuineTau, fTauData.getSelectedTau().lChTrkPt());
      fHistoSplitter.fillShapeHistogramTriplet(hTauDM_AfterPreselections      , bIsGenuineTau, fTauData.getSelectedTau().decayMode());
      fHistoSplitter.fillShapeHistogramTriplet(hTauNProngs_AfterPreselections , bIsGenuineTau, fTauData.getSelectedTau().nProngs());
      for (auto& p: fHelper.getTauSourceData(!event.isMC(), fTauData.getSelectedTau()))
	{
	  fHistoSplitter.fillShapeHistogramTriplet(hTauSource_AfterPreselections, bIsGenuineTau, p);
	}
      
      double deltaPhiTauMet = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(fTauData.getSelectedTau().p4(), fMETData.getMET())*57.29578);
      fHistoSplitter.fillShapeHistogramTriplet(hDeltaPhiTauMet_AfterPreselections, bIsGenuineTau, deltaPhiTauMet);
    }


  // Hadronic jets  
  fHistoSplitter.fillShapeHistogramTriplet(hJetN_AfterPreselections, bIsGenuineTau, nJets);
  for (auto& p: fJetData.getSelectedJets()) 
    {
      fHistoSplitter.fillShapeHistogramTriplet(hJetPt_AfterPreselections, bIsGenuineTau, p.pt());
      fHistoSplitter.fillShapeHistogramTriplet(hJetEta_AfterPreselections, bIsGenuineTau, p.eta());
    }
  fHistoSplitter.fillShapeHistogramTriplet(hHT_AfterPreselections , bIsGenuineTau, fJetData.HT());
  fHistoSplitter.fillShapeHistogramTriplet(hMHT_AfterPreselections, bIsGenuineTau, std::sqrt(fJetData.MHT().perp2()));


  // b-jets
  fHistoSplitter.fillShapeHistogramTriplet(hBjetN_AfterPreselections, bIsGenuineTau, nBJets);
  for (auto& p: fBJetData.getSelectedBJets()) 
    {
      fHistoSplitter.fillShapeHistogramTriplet(hBjetPt_AfterPreselections  , bIsGenuineTau, p.pt() );
      fHistoSplitter.fillShapeHistogramTriplet(hBjetEta_AfterPreselections , bIsGenuineTau, p.eta() );
      fHistoSplitter.fillShapeHistogramTriplet(hBjetDiscr_AfterPreselections, bIsGenuineTau, p.bjetDiscriminator() );
    }


  // MET
  fHistoSplitter.fillShapeHistogramTriplet(hMET_AfterPreselections   , bIsGenuineTau, fMETData.getMET().R() );
  fHistoSplitter.fillShapeHistogramTriplet(hMETPhi_AfterPreselections, bIsGenuineTau, fMETData.getMET().Phi() );


  // Transverse Mass
  double mT = -1.0;
  if (nTaus >= 1 && nMuons >= 1 && nJets >= 2)
    { 
      mT = TransverseMass::reconstruct(fTauData.getSelectedTau().p2(), fMuonData.getSelectedMuon().p2(), fJetData.getSelectedJets()[0].p2(), fJetData.getSelectedJets()[1].p2(), fMETData.getMET());
    }
  fHistoSplitter.fillShapeHistogramTriplet(hTransverseMass_AfterPreselections, bIsGenuineTau, mT);

  return;
}


void CommonPlots_lt::fillControlPlotsAfterSelections(const Event& event)
{

  // NOTE: This assumes that fillControlPlots_AfterPreselections() has already been called,
  //       and hence all "Data" structures have been assigned

  // Fill Histograms
  fHistoSplitter.fillShapeHistogramTriplet(hVerticesN_AfterSelections, bIsGenuineTau, nVertices);

  // Muons    
  fHistoSplitter.fillShapeHistogramTriplet(hMuonN_AfterSelections, bIsGenuineTau, nMuons);
  if (nMuons > 0)
    {
      fHistoSplitter.fillShapeHistogramTriplet(hMuonPt_AfterSelections , bIsGenuineTau, fMuonData.getHighestSelectedMuonPt());
      fHistoSplitter.fillShapeHistogramTriplet(hMuonEta_AfterSelections, bIsGenuineTau, fMuonData.getHighestSelectedMuonEta());
      fHistoSplitter.fillShapeHistogramTriplet(hMuonPhi_AfterSelections, bIsGenuineTau, fMuonData.getHighestSelectedMuonPhi());

      double deltaPhiMuonMet = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(fMuonData.getSelectedMuon().p4(), fMETData.getMET())*57.29578);
      fHistoSplitter.fillShapeHistogramTriplet(hDeltaPhiMuonMet_AfterSelections, bIsGenuineTau, deltaPhiMuonMet);
    }


  // Taus
  fHistoSplitter.fillShapeHistogramTriplet(hTauN_AfterSelections, bIsGenuineTau, nTaus);
  if (nTaus > 0)
    {
      fHistoSplitter.fillShapeHistogramTriplet(hTauPt_AfterSelections      , bIsGenuineTau, fTauData.getSelectedTau().pt());
      fHistoSplitter.fillShapeHistogramTriplet(hTauEta_AfterSelections     , bIsGenuineTau, fTauData.getSelectedTau().eta());
      fHistoSplitter.fillShapeHistogramTriplet(hTauPhi_AfterSelections     , bIsGenuineTau, fTauData.getSelectedTau().phi());
      fHistoSplitter.fillShapeHistogramTriplet(hTauLdgTrkPt_AfterSelections, bIsGenuineTau, fTauData.getSelectedTau().lChTrkPt());
      fHistoSplitter.fillShapeHistogramTriplet(hTauDM_AfterSelections      , bIsGenuineTau, fTauData.getSelectedTau().decayMode());
      fHistoSplitter.fillShapeHistogramTriplet(hTauNProngs_AfterSelections , bIsGenuineTau, fTauData.getSelectedTau().nProngs());
      for (auto& p: fHelper.getTauSourceData(!event.isMC(), fTauData.getSelectedTau()))
	{
	  fHistoSplitter.fillShapeHistogramTriplet(hTauSource_AfterSelections, bIsGenuineTau, p);
	}
      
      double deltaPhiTauMet = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(fTauData.getSelectedTau().p4(), fMETData.getMET())*57.29578);
      fHistoSplitter.fillShapeHistogramTriplet(hDeltaPhiTauMet_AfterSelections, bIsGenuineTau, deltaPhiTauMet);
    }


  // Hadronic jets  
  fHistoSplitter.fillShapeHistogramTriplet(hJetN_AfterSelections, bIsGenuineTau, nJets);
  for (auto& p: fJetData.getSelectedJets()) 
    {
      fHistoSplitter.fillShapeHistogramTriplet(hJetPt_AfterSelections, bIsGenuineTau, p.pt());
      fHistoSplitter.fillShapeHistogramTriplet(hJetEta_AfterSelections, bIsGenuineTau, p.eta());
    }
  fHistoSplitter.fillShapeHistogramTriplet(hHT_AfterSelections , bIsGenuineTau, fJetData.HT());
  fHistoSplitter.fillShapeHistogramTriplet(hMHT_AfterSelections, bIsGenuineTau, std::sqrt(fJetData.MHT().perp2()));


  // b-jets
  fHistoSplitter.fillShapeHistogramTriplet(hBjetN_AfterSelections, bIsGenuineTau, nBJets);
  for (auto& p: fBJetData.getSelectedBJets()) 
    {
      fHistoSplitter.fillShapeHistogramTriplet(hBjetPt_AfterSelections  , bIsGenuineTau, p.pt() );
      fHistoSplitter.fillShapeHistogramTriplet(hBjetEta_AfterSelections , bIsGenuineTau, p.eta() );
      fHistoSplitter.fillShapeHistogramTriplet(hBjetDiscr_AfterSelections, bIsGenuineTau, p.bjetDiscriminator() );
    }


  // MET
  fHistoSplitter.fillShapeHistogramTriplet(hMET_AfterSelections   , bIsGenuineTau, fMETData.getMET().R() );
  fHistoSplitter.fillShapeHistogramTriplet(hMETPhi_AfterSelections, bIsGenuineTau, fMETData.getMET().Phi() );


  // Transverse Mass
  double mT = -1.0;
  if (nTaus >= 1 && nMuons >= 1 && nJets >= 2)
    { 
      mT = TransverseMass::reconstruct(fTauData.getSelectedTau().p2(), fMuonData.getSelectedMuon().p2(), fJetData.getSelectedJets()[0].p2(), fJetData.getSelectedJets()[1].p2(), fMETData.getMET());
    }
  fHistoSplitter.fillShapeHistogramTriplet(hTransverseMass_AfterSelections, bIsGenuineTau, mT );

  return;
}
