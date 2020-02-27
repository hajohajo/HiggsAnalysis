// -*- c++ -*-
#ifndef EventSelection_CommonPlots_ttm_h
#define EventSelection_CommonPlots_ttm_h

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

class CommonPlots_ttm {
public:
  enum AnalysisType {
    kSignalAnalysis,
    kHplus2tbAnalysis,
    kHplus2hwAnalysis,
    kHplus2hwAnalysis_mmt,
    kHplus2hwAnalysis_eet,
    kHplus2hwAnalysisWithTop,
    kHplus2hw_ele_Analysis,
    kBTagEfficiencyAnalysis,
    kTauAnalysis,
    kTopReco,
    kMuAnalysis,
    kQCDMeasurement,
    kQCDMeasurement_ele,
    kQCDMeasurement_muon,
    kQCDMeasurement_mmt,
    kQCDMeasurement_eet,
    kFakeBMeasurement,
    kQCDNormalizationSystematicsSignalRegion, // Needed for obtaining normalization systematics to data-driven control plots
    kQCDNormalizationSystematicsControlRegion // Needed for obtaining normalization systematics to data-driven control plots
  };

  CommonPlots_ttm(const ParameterSet& config, const CommonPlots_ttm::AnalysisType type, HistoWrapper& histoWrapper);
  CommonPlots_ttm(const ParameterSet& config, const CommonPlots_ttm::AnalysisType type, HistoWrapper& histoWrapper, bool test);
//  CommonPlots();
  ~CommonPlots_ttm();

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
  /// Returns the histogram settings for pt histograms (usecase: QCD measurement)
  const HistogramSettings& getPtBinSettings() const { return fPtBinSettings; }
  /// Returns the histogram settings for eta histograms (usecase: QCD measurement)
  const HistogramSettings& getEtaBinSettings() const { return fEtaBinSettings; }
  /// Returns the histogram settings for MET bins (usecase: QCD measurement)
  const HistogramSettings& getMetBinSettings() const { return fMetBinSettings; }
  /// Returns the histogram settings for HT bins (usecase: Htb analysis)
  const HistogramSettings& getHtBinSettings() const { return fHtBinSettings; }
  /// Returns the histogram settings for Mt bins (usecase: QCD measurement)
  const HistogramSettings& getMtBinSettings() const { return fMtBinSettings; }
  /// Returns the histogram settings for Mt bins (usecase: QCD measurement)
  const HistogramSettings& getBJetDiscBinSettings() const { return fBJetDiscriminatorBinSettings;}
  /// Returns the histogram settings for Njets bins (usecase: FakeB measurement)
  const HistogramSettings& getNjetsBinSettings() const { return fNjetsBinSettings;}
  /// Returns the histogram settings for # Vertices historams
  const HistogramSettings& getNVtxBinSettings() const { return fNVerticesBinSettings; }
  /// Returns the histogram settings for Phi histograms
  const HistogramSettings& getPhiBinSettings() const { return fPhiBinSettings; }
  /// Returns the histogram settings for DeltaEta histograms
  const HistogramSettings& getDeltaEtaBinSettings() const { return fDeltaEtaBinSettings; }
  /// Returns the histogram settings for DeltaPhi histograms
  const HistogramSettings& getDeltaPhiBinSettings() const { return fDeltaPhiBinSettings; }
  /// Returns the histogram settings for DeltaR histograms
  const HistogramSettings& getDeltaRBinSettings() const { return fDeltaRBinSettings; }
  /// Returns the histogram settings for Rtau histograms
  const HistogramSettings& getRtauBinSettings() const { return fRtauBinSettings; }
  /// Returns the histogram settings for Angular Cuts (1D) histograms
  const HistogramSettings& getAngularCuts1DBinSettings() const { return fAngularCuts1DSettings; }
  /// Returns the histogram settings for W-mass histograms
  const HistogramSettings& getWMassBinSettings() const { return fWMassBinSettings; }
  /// Returns the histogram settings for Topmass histograms
  const HistogramSettings& getTopMassBinSettings() const { return fTopMassBinSettings; }
  /// Returns the histogram settings for InvMass histograms
  const HistogramSettings& getInvMassBinSettings() const { return fInvmassBinSettings; }

  /** Special method for setting genuine tau status 
   * (it is usually set through TauSelection via CommonPlots::fillControlPlotsAfterTauSelection)
   */
  void setGenuineTauStatus(const bool isGenuineTau) { bIsGenuineTau = isGenuineTau; };
  
  //===== unique filling methods (to be called inside the event selection routine only, i.e. (before a passing decision is done))
  void fillControlPlotsAtVertexSelection(const Event& event);
  //void fillControlPlotsAtVetoTauSelection(const Event& event, const VetoTauSelection::Data& tauVetoData);
  void fillControlPlotsAtElectronSelection(const Event& event, const ElectronSelection::Data& data);
  void fillControlPlotsAtMuonSelection(const Event& event, const MuonSelection::Data& data);
  void fillControlPlotsAtTauSelection(const Event& event, const TauSelection::Data& data);
  void fillControlPlotsAtJetSelection(const Event& event, const JetSelection::Data& data);
  void fillControlPlotsAtAngularCutsCollinear(const Event& event, const AngularCutsCollinear::Data& data);
  void fillControlPlotsAtMETSelection(const Event& event, const METSelection::Data& data);
  void fillControlPlotsAtBtagging(const Event& event, const BJetSelection::Data& data);
  void fillControlPlotsAtAngularCutsBackToBack(const Event& event, const AngularCutsBackToBack::Data& data);
  //void fillControlPlotsAtTopSelection(const Event& event, const TopSelectionManager::Data& data);
  //void fillControlPlotsAtEvtTopology(const Event& event, const EvtTopology::Data& data);
  
  //===== unique filling methods (to be called AFTER return statement from analysis routine)
  void setNvertices(int vtx) { iVertices = vtx; fPUDependencyPlots->setNvtx(vtx); }
  void fillControlPlotsAfterTrigger(const Event& event);
  void fillControlPlotsAfterMETFilter(const Event& event);
  void fillControlPlotsAfterTauSelection(const Event& event, const TauSelection::Data& data);
  void fillControlPlotsAfterAntiIsolatedTauSelection(const Event& event, const TauSelection::Data& data);
  //void fillControlPlotsAfterMuonSelection(const Event& event);
  void fillControlPlotsAfterMETTriggerScaleFactor(const Event& event);
  void fillControlPlotsAfterBjetSelection(const Event& event, const BJetSelection::Data& data);
  void fillControlPlotsAfterBtagSF(const Event& event,const JetSelection::Data& jetData ,const BJetSelection::Data& bjetData);
  void fillControlPlotsAfterStandardSelections(const Event& event, 
					       const JetSelection::Data& jetData, 
					       const BJetSelection::Data& bjetData, 
					       const METSelection::Data& METData, 
					       const TopSelectionMVA::Data& topData);
  void fillControlPlotsAfterTopologicalSelections(const Event& event, bool withoutTau=false, bool withMu=false);
  void fillControlPlotsAfterAllSelections(const Event& event, bool withoutTau=false, int i=0, int j=1);
  void fillControlPlotsAfterAllSelections(const Event& event, int isGenuineB);  //HToTB-specific
  void fillControlPlotsAfterAllSelectionsWithProbabilisticBtag(const Event& event, const METSelection::Data& metData, double btagWeight);
  //void fillControlPlotsAfterAllSelectionsWithFullMass(const Event& event, FullHiggsMassCalculator::Data& data);

  //===== Filling of control plots for determining QCD shape uncertainty
  void fillControlPlotsForQCDShapeUncertainty(const Event& event,
                                              const AngularCutsBackToBack::Data& collinearAngularCutsData,
                                              const BJetSelection::Data& bJetData,
                                              const METSelection::Data& metData,
                                              const AngularCutsCollinear::Data& backToBackAngularCutsData);
  
  /// Getter for all vertices
  int nVertices() const { return iVertices; }

private:
  /// Returns true if anti-isolated taus need to be used (QCD measurement)
  const bool usesAntiIsolatedTaus() const { return fAnalysisType == kQCDMeasurement || fAnalysisType == kQCDNormalizationSystematicsControlRegion || fAnalysisType == kQCDMeasurement_ele || fAnalysisType == kQCDMeasurement_muon || fAnalysisType == kQCDMeasurement_mmt; }
  
private:
  ///===== Config params
  const bool fEnableGenuineTauHistograms;
  
  ///===== Analysis type
  const AnalysisType fAnalysisType;

  ///===== HistoWrapper;
  HistoWrapper fHistoWrapper;
  
  ///===== Histogram splitter
  HistoSplitter fHistoSplitter;

  ///===== Settings for histogram binning
  const HistogramSettings fNVerticesBinSettings;
  const HistogramSettings fPtBinSettings;
  const HistogramSettings fEtaBinSettings;
  const HistogramSettings fPhiBinSettings;
  const HistogramSettings fDeltaEtaBinSettings;
  const HistogramSettings fDeltaPhiBinSettings;
  const HistogramSettings fDeltaRBinSettings;
  const HistogramSettings fRtauBinSettings;
  const HistogramSettings fNjetsBinSettings;
  const HistogramSettings fMetBinSettings;
  const HistogramSettings fHtBinSettings;
  const HistogramSettings fBJetDiscriminatorBinSettings;
  const HistogramSettings fAngularCuts1DSettings;
  const HistogramSettings fWMassBinSettings;
  const HistogramSettings fTopMassBinSettings;
  const HistogramSettings fInvmassBinSettings;
  const HistogramSettings fMtBinSettings;

  ///===== Histograms
  // NOTE: think before adding a histogram - they do slow down the analysis a lot
  // NOTE: the histograms with the prefix hCtrl are used as data driven control plots
  // NOTE: the histograms with the prefix hShape are used as shape histograms
  // NOTE: histogram triplets contain the inclusive and events with fake tau histograms
  
  // jet selection
  HistoSplitter::SplittedTripletTH1s hCtrlNjets;
  
  // MET trigger SF
  HistoSplitter::SplittedTripletTH1s hCtrlNjetsAfterJetSelectionAndMETSF;
  
  // collinear angular cuts
  HistoSplitter::SplittedTripletTH1s hCtrlCollinearAngularCutsMinimum;
  HistoSplitter::SplittedTripletTH1s hCtrlCollinearAngularCutsJet1;
  HistoSplitter::SplittedTripletTH1s hCtrlCollinearAngularCutsJet2;
  HistoSplitter::SplittedTripletTH1s hCtrlCollinearAngularCutsJet3;
  HistoSplitter::SplittedTripletTH1s hCtrlCollinearAngularCutsJet4;
  // Collinear or B2B? Decide later
  HistoSplitter::SplittedTripletTH1s hCtrlAngularCutsDeltaPhiTaus;
  HistoSplitter::SplittedTripletTH1s hCtrlAngularCutsDeltaPhiTauMET;
  HistoSplitter::SplittedTripletTH1s hCtrlAngularCutsDeltaPhiLdgTauMET;
  HistoSplitter::SplittedTripletTH1s hCtrlAngularCutsDeltaPhiSubldgTauMET;
  HistoSplitter::SplittedTripletTH1s hCtrlAngularCutsDeltaPhiLdgJetMET;
  HistoSplitter::SplittedTripletTH1s hCtrlAngularCutsDeltaPhiSubldgJetMET;
  HistoSplitter::SplittedTripletTH1s hCtrlAngularCutsDeltaPhiMuonMET;
  HistoSplitter::SplittedTripletTH1s hCtrlAngularCutsDeltaPhiLdgTauMuon;
  HistoSplitter::SplittedTripletTH1s hCtrlAngularCutsDeltaPhiSubldgTauMuon;
  HistoSplitter::SplittedTripletTH2s hCtrlAngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiSubldgTauMET;
  HistoSplitter::SplittedTripletTH2s hCtrlAngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiMuonMET;
  HistoSplitter::SplittedTripletTH2s hCtrlAngularCutsDeltaPhiLdgTauMET_Vs_DeltaPhiLdgBJetMET;
  HistoSplitter::SplittedTripletTH2s hCtrlAngularCutsDeltaPhiLdgTauMuon_Vs_DeltaPhiSubldgTauMuon;
  
  // this is the point of "standard selections"
  HistoSplitter::SplittedTripletTH1s hCtrlNVerticesAfterStdSelections;  
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTauPtAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTauEtaAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTauPhiAfterStdSelections;
  HistoSplitter::SplittedTripletTH2s hCtrlSelectedTauEtaPhiAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTauLdgTrkPtAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTauDecayModeAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTauNProngsAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTauRtauAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTauSourceAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSubldgTauPtAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSubldgTauEtaAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSubldgTauPhiAfterStdSelections;
  HistoSplitter::SplittedTripletTH2s hCtrlSubldgTauEtaPhiAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSubldgTauLdgTrkPtAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSubldgTauDecayModeAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSubldgTauNProngsAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedMuonPtAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedMuonEtaAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedMuonPhiAfterStdSelections;
  HistoSplitter::SplittedTripletTH2s hCtrlSelectedMuonEtaPhiAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlNJetsAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJetPtAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJetEtaAfterStdSelections;
  HistoSplitter::SplittedTripletTH2s hCtrlJetEtaPhiAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlNBJetsAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlBJetPtAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlBJetEtaAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlBDiscriminatorAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlMETAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlMETPhiAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlDeltaPhiTauMetAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlDeltaPhiMuonMetAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlHTAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlMHTAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlNTopsAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopPtAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopBDTAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopDijetPtAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopDijetMassAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopMassAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopMassWMassRatioAfterStdSelections;
  HistoSplitter::SplittedTripletTH2s hCtrlTopPt_Vs_TopDijetPtAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopBJetPtAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopBJetEtaAfterStdSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopBJetBdiscAfterStdSelections;

  // MET
  HistoSplitter::SplittedTripletTH1s hCtrlMET;
  HistoSplitter::SplittedTripletTH1s hCtrlMETPhi;
  
  // b tagging
  HistoSplitter::SplittedTripletTH1s hCtrlNBJets;
  HistoSplitter::SplittedTripletTH1s hCtrlBJetPt;
  HistoSplitter::SplittedTripletTH1s hCtrlBJetEta;
  HistoSplitter::SplittedTripletTH1s hCtrlBDiscriminator;

  // after b-tag SF
  HistoSplitter::SplittedTripletTH1s hCtrlNJetsAfterBtagSF;
  HistoSplitter::SplittedTripletTH1s hCtrlJetPtAfterBtagSF;
  HistoSplitter::SplittedTripletTH1s hCtrlBJetPtAfterBtagSF;
  
  // back-to-back angular cuts
  HistoSplitter::SplittedTripletTH1s hCtrlBackToBackAngularCutsMinimum;
  HistoSplitter::SplittedTripletTH1s hCtrlBackToBackAngularCutsJet1;
  HistoSplitter::SplittedTripletTH1s hCtrlBackToBackAngularCutsJet2;
  HistoSplitter::SplittedTripletTH1s hCtrlBackToBackAngularCutsJet3;
  HistoSplitter::SplittedTripletTH1s hCtrlBackToBackAngularCutsJet4;

  // control plots after all selections
  HistoSplitter::SplittedTripletTH1s hCtrlNVerticesAfterAllSelections;  
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTausPtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTausEtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTausPhiAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTauPtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTauEtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTauPhiAfterAllSelections;
  HistoSplitter::SplittedTripletTH2s hCtrlSelectedTauEtaPhiAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTauLdgTrkPtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTauDecayModeAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTauNProngsAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTauRtauAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTauSourceAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedTauIPxyAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSubldgTauPtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSubldgTauEtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSubldgTauPhiAfterAllSelections;
  HistoSplitter::SplittedTripletTH2s hCtrlSubldgTauEtaPhiAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSubldgTauLdgTrkPtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSubldgTauDecayModeAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSubldgTauNProngsAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSubldgTauIPxyAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedMuonPtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedMuonEtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedMuonPhiAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedElectronPtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedElectronEtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlSelectedElectronPhiAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlNJetsAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJetPtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJetEtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH2s hCtrlJetEtaPhiAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlHTAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlMHTAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlMinDeltaPhiJetMHTAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlMaxDeltaPhiJetMHTAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlMinDeltaRJetMHTAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlMinDeltaRReversedJetMHTAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlCollinearAngularCutsMinimumAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlMETAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlMETPhiAfterAllSelections;  
  HistoSplitter::SplittedTripletTH1s hCtrlNBJetsAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlBJetPtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlBJetEtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJet1PtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJet2PtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJet3PtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJet4PtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJet5PtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJet6PtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJet7PtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJet1EtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJet2EtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJet3EtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJet4EtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJet5EtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJet6EtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlJet7EtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlBJet1PtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlBJet2PtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlBJet3PtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlBJet4PtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlBJet1EtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlBJet2EtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlBJet3EtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlBJet4EtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlBDiscriminatorAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlBackToBackAngularCutsMinimumAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlDeltaPhiTauMetAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlDeltaPhiMuonMetAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlNTopsAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopPtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopBDTAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopDijetPtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopDijetMassAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopMassAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopMassWMassRatioAfterAllSelections;
  HistoSplitter::SplittedTripletTH2s hCtrlTopPt_Vs_TopDijetPtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopBJetPtAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopBJetEtaAfterAllSelections;
  HistoSplitter::SplittedTripletTH1s hCtrlTopBJetBdiscAfterAllSelections;
  
  // shape plots after all selections
  HistoSplitter::SplittedTripletTH1s hShapeTransverseMass;
  HistoSplitter::SplittedTripletTH1s hShapeProbabilisticBtagTransverseMass;

  // Other plots
  WrappedTH1* hNSelectedVsRunNumber; // For data only
  
  //====== Plots from base class
  std::vector<CommonPlotsBase*> fBaseObjects;
  PUDependencyPlots* fPUDependencyPlots;
  
  //====== Data cache
  /// Cached data objects from silent analyze
  //VertexSelection::Data fVertexData;
  int iVertices;
  TauSelection::Data fTauData;
  //FakeTauIdentifier::Data fFakeTauData;
  bool bIsGenuineTau;
  ElectronSelection::Data fElectronData;
  MuonSelection::Data fMuonData;
  JetSelection::Data fJetData;
  AngularCutsBackToBack::Data fCollinearAngularCutsData;
  BJetSelection::Data fBJetData;
  METSelection::Data fMETData;
  TopologySelection::Data fTopologyData;
  // QuarkGluonLikelihoodRatio::Data fQGLRData;
  TopSelectionMVA::Data fTopData;
  AngularCutsCollinear::Data fBackToBackAngularCutsData;
  // FatJetSelection::Data fFatJetData;
  // FatJetSoftDropSelection::Data fFatJetSoftDropData;

  /// Helper
  CommonPlotsHelper fHelper;
};

#endif