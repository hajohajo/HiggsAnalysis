// -*- c++ -*-
#ifndef EventSelection_CommonPlots_lt_h
#define EventSelection_CommonPlots_lt_h

#include "EventSelection/interface/EventSelections.h"
#include "EventSelection/interface/CommonPlotsHelper.h"
#include "EventSelection/interface/CommonPlotsBase.h"
#include "EventSelection/interface/PUDependencyPlots.h"
#include "Framework/interface/ParameterSet.h"
#include "Framework/interface/HistogramSettings.h"
#include "Framework/interface/HistoSplitter.h"
#include "Framework/interface/HistoWrapper.h"

#include "TDirectory.h"

#include <vector>

class CommonPlots_lt {
public:
  enum AnalysisType {
    kHToHW,
    kHToHW_wSystematics,
    kTauFakeRateMeasurement,
  };

  CommonPlots_lt(const ParameterSet& config, const CommonPlots_lt::AnalysisType type, HistoWrapper& histoWrapper);
  CommonPlots_lt(const ParameterSet& config, const CommonPlots_lt::AnalysisType type, HistoWrapper& histoWrapper, bool test);
  ~CommonPlots_lt();

  //Tau ID syst switches
  bool tauIDup;
  bool tauIDdown;
  
  void book(TDirectory* dir, bool isData);
  
  /// Initialize (call this at the beginning of each event; prevents double-counting of events)
  void initialize();
  /// Sets factorisation bin (call this for each event before filling the first histogram!)
  void setFactorisationBinForEvent(const std::vector<float>& values=std::vector<float>{}) { fHistoSplitter.setFactorisationBinForEvent(values); }
  
  /// Returns the histogram splitter object (usecase: QCD measurement)
  HistoSplitter& getHistoSplitter() { return fHistoSplitter; }
  const HistogramSettings& getPtBinSettings() const { return fPtBins; }
  const HistogramSettings& getEtaBinSettings() const { return fEtaBins; }
  const HistogramSettings& getMetBinSettings() const { return fMetBins; }
  const HistogramSettings& getHtBinSettings() const { return fHtBins; }
  const HistogramSettings& getBJetDiscBinSettings() const { return fBjetDiscrBins;}
  const HistogramSettings& getNjetsBinSettings() const { return fNjetsBins;}
  const HistogramSettings& getNVtxBinSettings() const { return fVerticesBins; }
  const HistogramSettings& getPhiBinSettings() const { return fPhiBins; }
  const HistogramSettings& getDeltaEtaBinSettings() const { return fDeltaEtaBins; }
  const HistogramSettings& getDeltaPhiBinSettings() const { return fDeltaPhiBins; }
  const HistogramSettings& getDeltaRBinSettings() const { return fDeltaRBins; }
  const HistogramSettings& getInvMassBinSettings() const { return fInvmassBins; }
  const HistogramSettings& getMtBinSettings() const { return fMtBins; }

  void setGenuineTauStatus(const bool isGenuineTau) { bIsGenuineTau = isGenuineTau; };
  
  // unique filling methods (to be called inside the event selection routine only (BEFORE passing decision is done)
  void fillControlPlotsAtVertexSelection(const Event& event);
  void fillControlPlotsAtElectronSelection(const Event& event, const ElectronSelection::Data& data);
  void fillControlPlotsAtMuonSelection(const Event& event, const MuonSelection::Data& data);
  void fillControlPlotsAtTauSelection(const Event& event, const TauSelection::Data& data);
  void fillControlPlotsAtJetSelection(const Event& event, const JetSelection::Data& data);
  void fillControlPlotsAtBtagging(const Event& event, const BJetSelection::Data& data);
  void fillControlPlotsAtMETSelection(const Event& event, const METSelection::Data& data);
  
  // unique filling methods (AFTER return statement from analysis routine)
  void setNvertices(int vtx) { nVertices = vtx; fPUDependencyPlots->setNvtx(vtx); }
  void fillControlPlotsAfterTrigger(const Event& event);
  void fillControlPlotsAfterMETFilter(const Event& event);
  void fillControlPlotsAfterTauSelection(const Event& event, const TauSelection::Data& data);
  void fillControlPlotsAfterAntiIsolatedTauSelection(const Event& event, const TauSelection::Data& data);
  void fillControlPlotsAfterBjetSelection(const Event& event, const BJetSelection::Data& data);
  void fillControlPlotsAfterBtagSF(const Event& event,const JetSelection::Data& jetData ,const BJetSelection::Data& bjetData);
  void fillControlPlotsAfterPreselections(const Event& event, const JetSelection::Data& jetData, const BJetSelection::Data& bjetData, const METSelection::Data& METData);
  void fillControlPlotsAfterSelections(const Event& event);

  // Filling of control plots for determining QCD shape uncertainty
  void fillControlPlotsForQCDShapeUncertainty(const Event& event,
                                              const AngularCutsBackToBack::Data& collinearAngularCutsData,
                                              const BJetSelection::Data& bJetData,
                                              const METSelection::Data& metData,
                                              const AngularCutsCollinear::Data& backToBackAngularCutsData);
  
  /// Getter for all vertices
  int getVertices() const { return nVertices; }

private:
  /// Returns true if anti-isolated taus need to be used
  const bool usesAntiIsolatedTaus() const { return fAnalysisType == kTauFakeRateMeasurement; }
  
private:
  /// Config params
  const bool fEnableGenuineTauHistograms;
  
  /// Analysis type
  const AnalysisType fAnalysisType;

  /// HistoWrapper;
  HistoWrapper fHistoWrapper;
  
  /// Histogram splitter
  HistoSplitter fHistoSplitter;

  /// Settings for histogram binning
  const HistogramSettings fVerticesBins;
  const HistogramSettings fPtBins;
  const HistogramSettings fEtaBins;
  const HistogramSettings fPhiBins;
  const HistogramSettings fDeltaEtaBins;
  const HistogramSettings fDeltaPhiBins;
  const HistogramSettings fDeltaRBins;
  const HistogramSettings fNjetsBins;
  const HistogramSettings fMetBins;
  const HistogramSettings fHtBins;
  const HistogramSettings fBjetDiscrBins;
  const HistogramSettings fMtBins;
  const HistogramSettings fInvmassBins;

  /// Histograms
  // NOTE: think before adding a histogram - they do slow down the analysis a lot
  // NOTE: the histograms with the prefix hCtrl are used as data driven control plots
  // NOTE: the histograms with the prefix hShape are used as shape histograms
  // NOTE: histogram triplets contain the inclusive and events with fake tau histograms


  // AtJetSelection  
  HistoSplitter::SplittedTripletTH1s hJetN_AtJetSelection;
  HistoSplitter::SplittedTripletTH1s hJetPt_AtJetSelection;
  HistoSplitter::SplittedTripletTH1s hJetEta_AtJetSelection;
  // AtBtagging
  HistoSplitter::SplittedTripletTH1s hBjetN_AtBtagging;
  HistoSplitter::SplittedTripletTH1s hBjetPt_AtBtagging;
  HistoSplitter::SplittedTripletTH1s hBjetEta_AtBtagging;
  HistoSplitter::SplittedTripletTH1s hBjetDiscr_AtBtagging;
  // AtMETSelection
  HistoSplitter::SplittedTripletTH1s hMET_AtMETSelection;
  HistoSplitter::SplittedTripletTH1s hMETPhi_AtMETSelection;
  // AfterBtagSF
  HistoSplitter::SplittedTripletTH1s hJetN_AfterBtagSF;
  HistoSplitter::SplittedTripletTH1s hJetPt_AfterBtagSF;
  HistoSplitter::SplittedTripletTH1s hBjetN_AfterBtagSF;
  HistoSplitter::SplittedTripletTH1s hBjetPt_AfterBtagSF;

  // AfterPreselections
  HistoSplitter::SplittedTripletTH1s hVerticesN_AfterPreselections;  
  HistoSplitter::SplittedTripletTH1s hTauN_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hTauPt_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hTauEta_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hTauPhi_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hTauLdgTrkPt_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hTauDM_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hTauNProngs_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hTauSource_AfterPreselections; 
  HistoSplitter::SplittedTripletTH1s hTauIPxy_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hMuonN_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hMuonPt_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hMuonEta_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hMuonPhi_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hJetN_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hJetPt_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hJetEta_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hBjetN_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hBjetPt_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hBjetEta_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hBjetDiscr_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hMET_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hMETPhi_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hDeltaPhiTauMet_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hDeltaPhiMuonMet_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hHT_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hMHT_AfterPreselections;
  HistoSplitter::SplittedTripletTH1s hTransverseMass_AfterPreselections;

  // AfterSelections
  HistoSplitter::SplittedTripletTH1s hVerticesN_AfterSelections;  
  HistoSplitter::SplittedTripletTH1s hTauN_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hTauPt_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hTauEta_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hTauPhi_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hTauLdgTrkPt_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hTauDM_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hTauNProngs_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hTauSource_AfterSelections; 
  HistoSplitter::SplittedTripletTH1s hTauIPxy_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hMuonN_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hMuonPt_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hMuonEta_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hMuonPhi_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hJetN_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hJetPt_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hJetEta_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hBjetN_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hBjetPt_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hBjetEta_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hBjetDiscr_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hMET_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hMETPhi_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hDeltaPhiTauMet_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hDeltaPhiMuonMet_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hHT_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hMHT_AfterSelections;
  HistoSplitter::SplittedTripletTH1s hTransverseMass_AfterSelections;
  
  // Other plots
  WrappedTH1* hNSelectedVsRunNumber; // For data only
  
  // Plots from base class
  std::vector<CommonPlotsBase*> fBaseObjects;
  PUDependencyPlots* fPUDependencyPlots;
  
  // Global variables
  bool bIsGenuineTau;
  int nVertices;
  int nElectrons;
  int nMuons;
  int nTaus;
  int nJets;
  int nBJets;

  // Data cache (Cached data objects from silent analyze)
  ElectronSelection::Data fElectronData;
  MuonSelection::Data fMuonData;
  TauSelection::Data fTauData;
  JetSelection::Data fJetData;
  BJetSelection::Data fBJetData;
  METSelection::Data fMETData;
  // TopSelectionMVA::Data fTopData;
  // AngularCutsBackToBack::Data fCollinearAngularCutsData;
  // AngularCutsCollinear::Data fBackToBackAngularCutsData;
  // FatJetSelection::Data fFatJetData;
  // FatJetSoftDropSelection::Data fFatJetSoftDropData;

  /// Helper
  CommonPlotsHelper fHelper;
};

#endif
