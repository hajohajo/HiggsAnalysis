// -*- c++ -*-
#include "Framework/interface/BaseSelector.h"
#include "Framework/interface/makeTH.h"

#include "EventSelection/interface/CommonPlots.h"
#include "EventSelection/interface/EventSelections.h"

#include "TDirectory.h"

class QGLAnalysis: public BaseSelector {
public:
  explicit QGLAnalysis(const ParameterSet& config, const TH1* skimCounters);
  virtual ~QGLAnalysis() {}

  /// Books histograms
  virtual void book(TDirectory *dir) override;
  /// Sets up branches for reading the TTree
  virtual void setupBranches(BranchManager& branchManager) override;
  /// Called for each event
  virtual void process(Long64_t entry) override;

private:
  
  // Common plots
  CommonPlots fCommonPlots;
  // Event selection classes and event counters (in same order like they are applied)
  Count cAllEvents;
  Count cTrigger;
  METFilterSelection fMETFilterSelection;
  Count cVertexSelection;
  ElectronSelection fElectronSelection;
  MuonSelection fMuonSelection;
  TauSelection fTauSelection;
  JetSelection fJetSelection;
  BJetSelection fBJetSelection;
  Count cBTaggingSFCounter;
  METSelection fMETSelection;
  TopSelectionMVA fTopSelection;
  Count cSelected;
  
  // Non-common histograms
  WrappedTH1 *hCtrlMETAfterJetSelection;
  WrappedTH1 *hCtrlNVerticesAfterJetSelection;
  WrappedTH1 *hCtrlNJetsAfterJetSelection;
  WrappedTH1 *hCtrlJetsPtAfterJetSelection;
  WrappedTH1 *hCtrlJetsEtaAfterJetSelection;
  WrappedTH1 *hCtrlJetsPhiAfterJetSelection;
  WrappedTH1 *hCtrlJetsQGLAfterJetSelection;
  WrappedTH1 *hCtrlJetsBDiscAfterJetSelection;
  WrappedTH1 *hCtrlJetsHadronFlavourAfterJetSelection;
  WrappedTH1 *hCtrlJetsPartonFlavourAfterJetSelection;
  WrappedTH1 *hCtrlHTAfterJetSelection;
  
  WrappedTH1 *hCtrlNBJets;
  WrappedTH1 *hCtrlBJetsPt;
  WrappedTH1 *hCtrlBJetsEta;
  WrappedTH1 *hCtrlBJetsPhi;
  WrappedTH1 *hCtrlBJetsBDisc;
  WrappedTH1 *hCtrlBJetsQGL;
  
  WrappedTH1 *hCtrlNCJets;
  WrappedTH1 *hCtrlCJetsPt;
  WrappedTH1 *hCtrlCJetsEta;
  WrappedTH1 *hCtrlCJetsPhi;
  WrappedTH1 *hCtrlCJetsBDisc;
  WrappedTH1 *hCtrlCJetsQGL;
  
  WrappedTH1 *hCtrlNUDSGJets;
  WrappedTH1 *hCtrlUDSGJetsPt;
  WrappedTH1 *hCtrlUDSGJetsEta;
  WrappedTH1 *hCtrlUDSGJetsPhi;
  WrappedTH1 *hCtrlUDSGJetsBDisc;
  WrappedTH1 *hCtrlUDSGJetsQGL;
  WrappedTH1 *hCtrlNGJets;
  WrappedTH1 *hCtrlGJetsPt;
  WrappedTH1 *hCtrlGJetsEta;
  WrappedTH1 *hCtrlGJetsPhi;
  WrappedTH1 *hCtrlGJetsBDisc;
  WrappedTH1 *hCtrlGJetsQGL;
  WrappedTH1 *hCtrlNUDSJets;
  WrappedTH1 *hCtrlUDSJetsPt;
  WrappedTH1 *hCtrlUDSJetsEta;
  WrappedTH1 *hCtrlUDSJetsPhi;
  WrappedTH1 *hCtrlUDSJetsBDisc;
  WrappedTH1 *hCtrlUDSJetsQGL;
  
  WrappedTH2 *hBJets_QGLvsPt;
  WrappedTH2 *hCJets_QGLvsPt;  
  WrappedTH2 *hUDSGJets_QGLvsPt;
  WrappedTH2 *hUDSJets_QGLvsPt;
  WrappedTH2 *hGJets_QGLvsPt;
  
  WrappedTH1 *hGJetsQGL_30pt35;
  WrappedTH1 *hGJetsQGL_35pt40;
  WrappedTH1 *hGJetsQGL_40pt45;
  WrappedTH1 *hGJetsQGL_45pt50;
  WrappedTH1 *hGJetsQGL_50pt55;
  WrappedTH1 *hGJetsQGL_55pt60;
  WrappedTH1 *hGJetsQGL_60pt65;
  WrappedTH1 *hGJetsQGL_65pt70;
  WrappedTH1 *hGJetsQGL_70pt75;
  WrappedTH1 *hGJetsQGL_75pt80;
  WrappedTH1 *hGJetsQGL_80pt90;
  WrappedTH1 *hGJetsQGL_90pt100;
  WrappedTH1 *hGJetsQGL_100pt110;
  WrappedTH1 *hGJetsQGL_110pt120;
  WrappedTH1 *hGJetsQGL_120pt140;
  WrappedTH1 *hGJetsQGL_140pt160;
  WrappedTH1 *hGJetsQGL_160pt180;
  WrappedTH1 *hGJetsQGL_180pt200;
  WrappedTH1 *hGJetsQGL_200pt250;
  WrappedTH1 *hGJetsQGL_250pt300;
  WrappedTH1 *hGJetsQGL_300pt350;
  WrappedTH1 *hGJetsQGL_350pt400;
  WrappedTH1 *hGJetsQGL_400pt450;
  WrappedTH1 *hGJetsQGL_450pt500;
  WrappedTH1 *hGJetsQGL_500pt550;
  WrappedTH1 *hGJetsQGL_550pt600;
  WrappedTH1 *hGJetsQGL_600pt700;
  WrappedTH1 *hGJetsQGL_700pt800;
  WrappedTH1 *hGJetsQGL_800pt1000;
  WrappedTH1 *hGJetsQGL_1000ptInf;
  
  // UDS Jets
  WrappedTH1 *hUDSJetsQGL_30pt35;
  WrappedTH1 *hUDSJetsQGL_35pt40;
  WrappedTH1 *hUDSJetsQGL_40pt45;
  WrappedTH1 *hUDSJetsQGL_45pt50;
  WrappedTH1 *hUDSJetsQGL_50pt55;
  WrappedTH1 *hUDSJetsQGL_55pt60;
  WrappedTH1 *hUDSJetsQGL_60pt65;
  WrappedTH1 *hUDSJetsQGL_65pt70;
  WrappedTH1 *hUDSJetsQGL_70pt75;
  WrappedTH1 *hUDSJetsQGL_75pt80;
  WrappedTH1 *hUDSJetsQGL_80pt90;
  WrappedTH1 *hUDSJetsQGL_90pt100;
  WrappedTH1 *hUDSJetsQGL_100pt110;
  WrappedTH1 *hUDSJetsQGL_110pt120;
  WrappedTH1 *hUDSJetsQGL_120pt140;
  WrappedTH1 *hUDSJetsQGL_140pt160;
  WrappedTH1 *hUDSJetsQGL_160pt180;
  WrappedTH1 *hUDSJetsQGL_180pt200;
  WrappedTH1 *hUDSJetsQGL_200pt250;
  WrappedTH1 *hUDSJetsQGL_250pt300;
  WrappedTH1 *hUDSJetsQGL_300pt350;
  WrappedTH1 *hUDSJetsQGL_350pt400;
  WrappedTH1 *hUDSJetsQGL_400pt450;
  WrappedTH1 *hUDSJetsQGL_450pt500;
  WrappedTH1 *hUDSJetsQGL_500pt550;
  WrappedTH1 *hUDSJetsQGL_550pt600;
  WrappedTH1 *hUDSJetsQGL_600pt700;
  WrappedTH1 *hUDSJetsQGL_700pt800;
  WrappedTH1 *hUDSJetsQGL_800pt1000;
  WrappedTH1 *hUDSJetsQGL_1000ptInf;
  
};

#include "Framework/interface/SelectorFactory.h"
REGISTER_SELECTOR(QGLAnalysis);

QGLAnalysis::QGLAnalysis(const ParameterSet& config, const TH1* skimCounters)
  : BaseSelector(config, skimCounters),
    fCommonPlots(config.getParameter<ParameterSet>("CommonPlots"), CommonPlots::kHplus2tbAnalysis, fHistoWrapper),
    cAllEvents(fEventCounter.addCounter("all events")),
    cTrigger(fEventCounter.addCounter("passed trigger")),
    fMETFilterSelection(config.getParameter<ParameterSet>("METFilter"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cVertexSelection(fEventCounter.addCounter("passed PV")),
    fElectronSelection(config.getParameter<ParameterSet>("ElectronSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    fMuonSelection(config.getParameter<ParameterSet>("MuonSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    fTauSelection(config.getParameter<ParameterSet>("TauSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    fJetSelection(config.getParameter<ParameterSet>("JetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fBJetSelection(config.getParameter<ParameterSet>("BJetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cBTaggingSFCounter(fEventCounter.addCounter("b tag SF")),
    // fMETSelection(config.getParameter<ParameterSet>("METSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fMETSelection(config.getParameter<ParameterSet>("METSelection")), // no subcounter in main counter
    fTopSelection(config.getParameter<ParameterSet>("TopSelectionMVA"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cSelected(fEventCounter.addCounter("Selected Events"))
{ }


void QGLAnalysis::book(TDirectory *dir) {

  
  // Book common plots histograms
  fCommonPlots.book(dir, isData());

  // Book histograms in event selection classes
  fMETFilterSelection.bookHistograms(dir);
  fElectronSelection.bookHistograms(dir);
  fMuonSelection.bookHistograms(dir);
  fTauSelection.bookHistograms(dir);
  fJetSelection.bookHistograms(dir);
  fBJetSelection.bookHistograms(dir);
  fMETSelection.bookHistograms(dir);
  fTopSelection.bookHistograms(dir);
  
  // Fixed Binning
  const int nBinsQGL = 100;
  const float minQGL = 0.0;
  const float maxQGL = 1.0;
  
  const int nPhiBins = 50;
  const float fPhiMin = -5.0;
  const float fPhiMax = +5.0;
  
  // Obtain binning
  const int  nPtBins = 2*fCommonPlots.getPtBinSettings().bins();
  const float fPtMin = fCommonPlots.getPtBinSettings().min();
  const float fPtMax = 2*fCommonPlots.getPtBinSettings().max();
  
  const int nMetBins = fCommonPlots.getMetBinSettings().bins();
  const float fMetMin = fCommonPlots.getMetBinSettings().min();
  const float fMetMax = fCommonPlots.getMetBinSettings().max();
  
  const int  nHtBins = fCommonPlots.getHtBinSettings().bins();
  const float fHtMin = fCommonPlots.getHtBinSettings().min();
  const float fHtMax = fCommonPlots.getHtBinSettings().max();
  
  const int nEtaBins = fCommonPlots.getEtaBinSettings().bins();
  const float fEtaMin = fCommonPlots.getEtaBinSettings().min();
  const float fEtaMax = fCommonPlots.getEtaBinSettings().max();

  const int nBinsBDisc = fCommonPlots.getBJetDiscBinSettings().bins();
  const float minBDisc = fCommonPlots.getBJetDiscBinSettings().min();
  const float maxBDisc = fCommonPlots.getBJetDiscBinSettings().max();
  
  // Book non-common histograms
  
  // Control Plots
  hCtrlMETAfterJetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlMETAfterJetSelection", ";E_{T}^{miss};Occur / %.1f", nMetBins, fMetMin, fMetMax);
  hCtrlNVerticesAfterJetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlNVerticesAfterJetSelection", ";Number of Vertices;Occur / %.2f", 150, 0, 150.0);
  hCtrlNJetsAfterJetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlNJetsAfterJetSelection", "Jets Multiplicity", 20, 0, 20.0);
  hCtrlJetsPtAfterJetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlJetsPtAfterJetSelection", "Jets p_{T} (GeV)", nPtBins, fPtMin, fPtMax);
  hCtrlJetsEtaAfterJetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlJetsEtaAfterJetSelection", "Jets eta;#eta", nEtaBins, fEtaMin, fEtaMax);
  hCtrlJetsPhiAfterJetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlJetsPhiAfterJetSelection", "Jets phi;#phi", nPhiBins, fPhiMin, fPhiMax);
  hCtrlJetsBDiscAfterJetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlJetsBDiscrAfterJetSelection", "Jets b-discriminator", nBinsBDisc, minBDisc, maxBDisc);
  hCtrlJetsQGLAfterJetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlJetsQGLAfterJetSelection", "Jets QGL", nBinsBDisc, minBDisc, maxBDisc);
  hCtrlJetsHadronFlavourAfterJetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlJetsHadronFlavourAfterJetSelection", "Jets Hadron Flavour", 22, 0, 22.0);
  hCtrlJetsPartonFlavourAfterJetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlJetsPartonFlavourAfterJetSelection", "Jets Parton Flavour", 22, 0, 22.0);
  hCtrlHTAfterJetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlHTAfterJetSelection", "H_{T} (GeV)", nHtBins, fHtMin, fHtMax);
  hCtrlNBJets = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlNBJets", "Jets Multiplicity", 20, 0, 20.0);
  hCtrlBJetsPt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlBJetsPt", "Jets p_{T} (GeV)", nPtBins, fPtMin, fPtMax);
  hCtrlBJetsEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlBJetsEta", "Jets eta;#eta", nEtaBins, fEtaMin, fEtaMax);
  hCtrlBJetsPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlBJetsPhi", "Jets phi;#phi", nPhiBins, fPhiMin, fPhiMax);
  hCtrlBJetsBDisc = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlBJetsBDisc", "Jets b-discriminator", nBinsBDisc, minBDisc, maxBDisc);
  hCtrlBJetsQGL = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlBJetsQGL", "Jets QGL", nBinsBDisc, minBDisc, maxBDisc);
  hCtrlNCJets = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlNCJets", "Jets Multiplicity", 20, 0, 20.0);
  hCtrlCJetsPt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlCJetsPt", "Jets p_{T} (GeV)", nPtBins, fPtMin, fPtMax);
  hCtrlCJetsEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlCJetsEta", "Jets eta;#eta", nEtaBins, fEtaMin, fEtaMax);
  hCtrlCJetsPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlCJetsPhi", "Jets phi;#phi", nPhiBins, fPhiMin, fPhiMax);
  hCtrlCJetsBDisc = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlCJetsBDisc", "Jets b-discriminator", nBinsBDisc, minBDisc, maxBDisc);
  hCtrlCJetsQGL = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlCJetsQGL", "Jets QGL", nBinsBDisc, minBDisc, maxBDisc);
  hCtrlNUDSGJets = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlNUDSGJets", "Jets Multiplicity", 20, 0, 20.0);
  hCtrlUDSGJetsPt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlUDSGJetsPt", "Jets p_{T} (GeV)", nPtBins, fPtMin, fPtMax);
  hCtrlUDSGJetsEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlUDSGJetsEt", "Jets eta;#eta", nEtaBins, fEtaMin, fEtaMax);
  hCtrlUDSGJetsPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlUDSGJetsPhi", "Jets phi;#phi", nPhiBins, fPhiMin, fPhiMax);
  hCtrlUDSGJetsBDisc = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlUDSGJetsBDisc", "Jets b-discriminator", nBinsBDisc, minBDisc, maxBDisc);
  hCtrlUDSGJetsQGL = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlUDSGJetsQGL", "Jets QGL", nBinsBDisc, minBDisc, maxBDisc);
  hCtrlNGJets = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlNGJets", "Jets Multiplicity", 20, 0, 20.0);
  hCtrlGJetsPt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlGJetsPt", "Jets p_{T} (GeV)", nPtBins, fPtMin, fPtMax);
  hCtrlGJetsEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlGJetsEta", "Jets eta;#eta", nEtaBins, fEtaMin, fEtaMax);
  hCtrlGJetsPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlGJetsPhi", "Jets phi;#phi", nPhiBins, fPhiMin, fPhiMax);
  hCtrlGJetsBDisc = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlGJetsBDisc", "Jets b-discriminator", nBinsBDisc, minBDisc, maxBDisc);
  hCtrlGJetsQGL = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlGJetsQGL", "Jets QGL", nBinsBDisc, minBDisc, maxBDisc);
  hCtrlNUDSJets = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlNUDSJets", "Jets Multiplicity", 20, 0, 20.0);
  hCtrlUDSJetsPt = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlUDSJetsPt", "Jets p_{T} (GeV)", nPtBins, fPtMin, fPtMax);
  hCtrlUDSJetsEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlUDSJetsEta", "Jets eta;#eta", nEtaBins, fEtaMin, fEtaMax);
  hCtrlUDSJetsPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlUDSJetsPhi", "Jets phi;#phi", nPhiBins, fPhiMin, fPhiMax);
  hCtrlUDSJetsBDisc = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlUDSJetsBDisc", "Jets b-discriminator", nBinsBDisc, minBDisc, maxBDisc);
  hCtrlUDSJetsQGL = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "CtrlUDSJetsQGL", "Jets QGL", nBinsBDisc, minBDisc, maxBDisc);
  
  hBJets_QGLvsPt = fHistoWrapper.makeTH<TH2F>(HistoLevel::kInformative, dir, "BJets_QGLvsPt", "BJets QGL vs p_{T}", nBinsQGL, minQGL, maxQGL, nPtBins, fPtMin, fPtMax);
  hCJets_QGLvsPt = fHistoWrapper.makeTH<TH2F>(HistoLevel::kInformative, dir, "CJets_QGLvsPt", "CJets QGL vs p_{T}", nBinsQGL, minQGL, maxQGL, nPtBins, fPtMin, fPtMax);
  hUDSGJets_QGLvsPt = fHistoWrapper.makeTH<TH2F>(HistoLevel::kInformative, dir, "UDSGJets_QGLvsPt", "UDSG Jets QGL vs p_{T}", nBinsQGL, minQGL, maxQGL, nPtBins, fPtMin, fPtMax);
  hUDSJets_QGLvsPt = fHistoWrapper.makeTH<TH2F>(HistoLevel::kInformative, dir, "UDSJets_QGLvsPt", "UDS Jets QGL vs p_{T}", nBinsQGL, minQGL, maxQGL, nPtBins, fPtMin, fPtMax);
  hGJets_QGLvsPt = fHistoWrapper.makeTH<TH2F>(HistoLevel::kInformative, dir, "GJets_QGLvsPt", "G Jets QGL vs p_{T}", nBinsQGL, minQGL, maxQGL, nPtBins, fPtMin, fPtMax);
  
  // Gluon Jets pT bins
  hGJetsQGL_30pt35 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_30pt35", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_35pt40 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_35pt40", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_40pt45 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_40pt45", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_45pt50 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_45pt50", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_50pt55 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_50pt55", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_55pt60 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_55pt60", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_60pt65 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_60pt65", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_65pt70 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_65pt70", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_70pt75 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_70pt75", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_75pt80 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_75pt80", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_80pt90 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_80pt90", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_90pt100 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_90pt100", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_100pt110 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_100pt110", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_110pt120 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_110pt120", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_120pt140 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_120pt140", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_140pt160 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_140pt160", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_160pt180 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_160pt180", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_180pt200 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_180pt200", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_200pt250 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_200pt250", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_250pt300 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_250pt300", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_300pt350 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_300pt350", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_350pt400 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_350pt400", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_400pt450 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_400pt450", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_450pt500 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_450pt500", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_500pt550 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_500pt550", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_550pt600 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_550pt600", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_600pt700 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_600pt700", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_700pt800 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_700pt800", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_800pt1000 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_800pt1000", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hGJetsQGL_1000ptInf = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "GJetsQGL_1000ptInf", "QGL disc", nBinsQGL, minQGL, maxQGL);
  
  // UDS Jets pT bins
  hUDSJetsQGL_30pt35 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_30pt35", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_35pt40 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_35pt40", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_40pt45 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_40pt45", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_45pt50 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_45pt50", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_50pt55 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_50pt55", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_55pt60 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_55pt60", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_60pt65 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_60pt65", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_65pt70 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_65pt70", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_70pt75 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_70pt75", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_75pt80 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_75pt80", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_80pt90 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_80pt90", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_90pt100 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_90pt100", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_100pt110 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_100pt110", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_110pt120 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_110pt120", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_120pt140 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_120pt140", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_140pt160 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_140pt160", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_160pt180 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_160pt180", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_180pt200 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_180pt200", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_200pt250 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_200pt250", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_250pt300 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_250pt300", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_300pt350 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_300pt350", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_350pt400 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_350pt400", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_400pt450 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_400pt450", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_450pt500 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_450pt500", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_500pt550 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_500pt550", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_550pt600 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_550pt600", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_600pt700 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_600pt700", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_700pt800 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_700pt800", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_800pt1000 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_800pt1000", "QGL disc", nBinsQGL, minQGL, maxQGL);
  hUDSJetsQGL_1000ptInf = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "UDSJetsQGL_1000ptInf", "QGL disc", nBinsQGL, minQGL, maxQGL);
  
  return;
}


void QGLAnalysis::setupBranches(BranchManager& branchManager) {
  fEvent.setupBranches(branchManager);
  return;
}


void QGLAnalysis::process(Long64_t entry) {
  //====== Initialize
  fCommonPlots.initialize();
  fCommonPlots.setFactorisationBinForEvent(std::vector<float> {});
  cAllEvents.increment();

  //================================================================================================   
  // 1) Apply trigger 
  //================================================================================================   
  if (0) std::cout << "=== Trigger" << std::endl;
  if ( !(fEvent.passTriggerDecision()) ) return;
  
  cTrigger.increment();
  int nVertices = fEvent.vertexInfo().value();
  fCommonPlots.setNvertices(nVertices);
  fCommonPlots.fillControlPlotsAfterTrigger(fEvent);

  //================================================================================================   
  // 2) MET filters (to remove events with spurious sources of fake MET)
  //================================================================================================   
  if (0) std::cout << "=== MET Filter" << std::endl;
  const METFilterSelection::Data metFilterData = fMETFilterSelection.analyze(fEvent);
  if (!metFilterData.passedSelection()) return;
  fCommonPlots.fillControlPlotsAfterMETFilter(fEvent);  

  //================================================================================================   
  // 3) Primarty Vertex (Check that a PV exists)
  //================================================================================================   
  if (0) std::cout << "=== Vertices" << std::endl;
  if (nVertices < 1) return;
  cVertexSelection.increment();
  fCommonPlots.fillControlPlotsAtVertexSelection(fEvent);
  
  //================================================================================================   
  // 4) Electron veto (Fully hadronic + orthogonality)
  //================================================================================================   
  if (0) std::cout << "=== Electron veto" << std::endl;
  const ElectronSelection::Data eData = fElectronSelection.analyze(fEvent);
  if (eData.hasIdentifiedElectrons()) return;

  //================================================================================================
  // 5) Muon veto (Fully hadronic + orthogonality)
  //================================================================================================
  if (0) std::cout << "=== Muon veto" << std::endl;
  const MuonSelection::Data muData = fMuonSelection.analyze(fEvent);
  if (muData.hasIdentifiedMuons()) return;

  //================================================================================================   
  // 6) Tau Veto (HToTauNu Orthogonality)
  //================================================================================================   
  if (0) std::cout << "=== Tau Veto" << std::endl;
  const TauSelection::Data tauData = fTauSelection.analyze(fEvent);
  if (tauData.hasIdentifiedTaus() ) return;
  
  //================================================================================================
  // 7) MET selection
  //================================================================================================
  if (0) std::cout << "=== MET selection" << std::endl;
  const METSelection::Data METData = fMETSelection.silentAnalyze(fEvent, nVertices);
  // if (!METData.passedSelection()) return;
  double MET = METData.getMET().R();
  
  //================================================================================================
  // 8) Jet selection
  //================================================================================================
  if (0) std::cout << "=== Jet selection" << std::endl;
  const JetSelection::Data jetData = fJetSelection.analyzeWithoutTau(fEvent);
  if (!jetData.passedSelection()) return;
  fCommonPlots.fillControlPlotsAfterTopologicalSelections(fEvent, true);
  
  //================================================================================================
  // Fill Control Plots AfterJetSelection
  //================================================================================================
  
  hCtrlMETAfterJetSelection -> Fill(MET);
  hCtrlNVerticesAfterJetSelection -> Fill(nVertices);
  hCtrlNJetsAfterJetSelection -> Fill(jetData.getNumberOfSelectedJets());
  hCtrlHTAfterJetSelection -> Fill(jetData.HT());
  for(const Jet& jet: jetData.getSelectedJets()) {
    hCtrlJetsPtAfterJetSelection -> Fill(jet.pt());
    hCtrlJetsEtaAfterJetSelection -> Fill(jet.eta());
    hCtrlJetsPhiAfterJetSelection -> Fill(jet.phi());
    hCtrlJetsQGLAfterJetSelection -> Fill(jet.QGTaggerAK4PFCHSqgLikelihood());
    hCtrlJetsBDiscAfterJetSelection -> Fill(jet.bjetDiscriminator());
    hCtrlJetsHadronFlavourAfterJetSelection -> Fill(jet.hadronFlavour());
    hCtrlJetsPartonFlavourAfterJetSelection -> Fill(std::abs(jet.partonFlavour()));
  }
  
  //================================================================================================
  // 9) B-Jet selection
  //================================================================================================
  //if (0) std::cout << "B-Jet selection" << std::endl;
  //const BJetSelection::Data bjetData = fBJetSelection.silentAnalyze(fEvent, jetData);
  // if (!bjetData.passedSelection()) return;
  

  //================================================================================================
  // 9) QGL selection
  //================================================================================================
  
  int nCJets = 0;
  int nBJets = 0;
  int nUDSGJets = 0;
  int nGJets = 0;
  int nUDSJets = 0;
  
  // Loop over selected jets
  for(const Jet& jet: jetData.getSelectedJets()) {
    
    double QGL = jet.QGTaggerAK4PFCHSqgLikelihood();
    double pt  = jet.pt();
    
    if (std::abs(jet.hadronFlavour()) == 4 && std::abs(jet.partonFlavour()) == 4)
      {
	nCJets++;
	hCtrlCJetsPt -> Fill(jet.pt());
	hCtrlCJetsEta -> Fill(jet.eta());
	hCtrlCJetsPhi -> Fill(jet.phi());
	hCtrlCJetsBDisc -> Fill(jet.bjetDiscriminator());
	hCtrlCJetsQGL -> Fill(QGL);
	hCJets_QGLvsPt -> Fill(QGL, pt);
      }
    else if (std::abs(jet.hadronFlavour()) == 5 && std::abs(jet.partonFlavour()) == 5)
      {
	nBJets++;
	hCtrlBJetsPt -> Fill(jet.pt());
	hCtrlBJetsEta -> Fill(jet.eta());
	hCtrlBJetsPhi -> Fill(jet.phi());
	hCtrlBJetsBDisc -> Fill(jet.bjetDiscriminator());
	hCtrlBJetsQGL -> Fill(QGL);
	hBJets_QGLvsPt -> Fill(QGL, pt);
      }
    
    //=== Reject jets consistent with b or c
    if (jet.hadronFlavour() != 0) continue;
    
    const short jetPartonFlavour = std::abs(jet.partonFlavour());
    
    //=== Keep only jets considered as light
    if (jetPartonFlavour != 21 && jetPartonFlavour != 1 && jetPartonFlavour != 2 && jetPartonFlavour != 3) continue;
    
    nUDSGJets++;
    hCtrlUDSGJetsPt    -> Fill(jet.pt());
    hCtrlUDSGJetsEta   -> Fill(jet.eta());
    hCtrlUDSGJetsPhi   -> Fill(jet.phi());
    hCtrlUDSGJetsBDisc -> Fill(jet.bjetDiscriminator());
    hCtrlUDSGJetsQGL   -> Fill(QGL);
    hUDSGJets_QGLvsPt -> Fill(QGL, pt);
    
    // Gluon Jets
    if (jetPartonFlavour == 21)
      {
	nGJets++;
	hCtrlGJetsPt -> Fill(jet.pt());
	hCtrlGJetsEta -> Fill(jet.eta());
	hCtrlGJetsPhi -> Fill(jet.phi());
	hCtrlGJetsBDisc -> Fill(jet.bjetDiscriminator());
	hCtrlGJetsQGL -> Fill(QGL);
	hGJets_QGLvsPt -> Fill(QGL, pt);
	
	if (pt >= 30 && pt < 35) hGJetsQGL_30pt35->Fill(QGL);
	else if (pt >= 35 && pt < 40) hGJetsQGL_35pt40 -> Fill(QGL);
	else if (pt >= 40 && pt < 45) hGJetsQGL_40pt45 -> Fill(QGL);
	else if (pt >= 45 && pt < 50) hGJetsQGL_45pt50 -> Fill(QGL);
	else if (pt >= 50 && pt < 55) hGJetsQGL_50pt55 -> Fill(QGL);
	else if (pt >= 55 && pt < 60) hGJetsQGL_55pt60 -> Fill(QGL);
	else if (pt >= 60 && pt < 65) hGJetsQGL_60pt65 -> Fill(QGL);
	else if (pt >= 65 && pt < 70) hGJetsQGL_65pt70 -> Fill(QGL);
	else if (pt >= 70 && pt < 75) hGJetsQGL_70pt75 -> Fill(QGL);
	else if (pt >= 75 && pt < 80) hGJetsQGL_75pt80 -> Fill(QGL);
	else if (pt >= 80 && pt < 90) hGJetsQGL_80pt90 -> Fill(QGL);
	else if (pt >= 90 && pt < 100) hGJetsQGL_90pt100 -> Fill(QGL);
	else if (pt >= 100 && pt < 110) hGJetsQGL_100pt110 -> Fill(QGL);
	else if (pt >= 110 && pt < 120) hGJetsQGL_110pt120 -> Fill(QGL);
	else if (pt >= 120 && pt < 140) hGJetsQGL_120pt140 -> Fill(QGL);
	else if (pt >= 140 && pt < 160) hGJetsQGL_140pt160 -> Fill(QGL);
	else if (pt >= 160 && pt < 180) hGJetsQGL_160pt180 -> Fill(QGL);
	else if (pt >= 180 && pt < 200) hGJetsQGL_180pt200 -> Fill(QGL);
	else if (pt >= 200 && pt < 250) hGJetsQGL_200pt250 -> Fill(QGL);
	else if (pt >= 250 && pt < 300) hGJetsQGL_250pt300 -> Fill(QGL);
	else if (pt >= 300 && pt < 350) hGJetsQGL_300pt350 -> Fill(QGL);
	else if (pt >= 350 && pt < 400) hGJetsQGL_350pt400 -> Fill(QGL);
	else if (pt >= 400 && pt < 450) hGJetsQGL_400pt450 -> Fill(QGL);
	else if (pt >= 450 && pt < 500) hGJetsQGL_450pt500 -> Fill(QGL);
	else if (pt >= 500 && pt < 550) hGJetsQGL_500pt550 -> Fill(QGL);
	else if (pt >= 550 && pt < 600) hGJetsQGL_550pt600 -> Fill(QGL);
	else if (pt >= 600 && pt < 700) hGJetsQGL_600pt700 -> Fill(QGL);
	else if (pt >= 700 && pt < 800) hGJetsQGL_700pt800 -> Fill(QGL);
	else if (pt >= 800 && pt < 1000) hGJetsQGL_800pt1000 -> Fill(QGL);
	else if (pt >= 1000) hGJetsQGL_1000ptInf -> Fill(QGL);
      }
    
    
    // Light Jets
    if (jetPartonFlavour == 1 || jetPartonFlavour == 2 || jetPartonFlavour == 3)
      {
	nUDSJets++;
	hUDSJets_QGLvsPt -> Fill(QGL, pt);
	hCtrlUDSJetsPt -> Fill(jet.pt());
	hCtrlUDSJetsEta -> Fill(jet.eta());
	hCtrlUDSJetsPhi -> Fill(jet.phi());
	hCtrlUDSJetsQGL -> Fill(QGL);
	hCtrlUDSJetsBDisc -> Fill(jet.bjetDiscriminator());
	
	if (pt >= 30 && pt < 35) hUDSJetsQGL_30pt35->Fill(QGL);
	else if (pt >= 35 && pt < 40) hUDSJetsQGL_35pt40 -> Fill(QGL);
	else if (pt >= 40 && pt < 45) hUDSJetsQGL_40pt45 -> Fill(QGL);
	else if (pt >= 45 && pt < 50) hUDSJetsQGL_45pt50 -> Fill(QGL);
	else if (pt >= 50 && pt < 55) hUDSJetsQGL_50pt55 -> Fill(QGL);
	else if (pt >= 55 && pt < 60) hUDSJetsQGL_55pt60 -> Fill(QGL);
	else if (pt >= 60 && pt < 65) hUDSJetsQGL_60pt65 -> Fill(QGL);
	else if (pt >= 65 && pt < 70) hUDSJetsQGL_65pt70 -> Fill(QGL);
	else if (pt >= 70 && pt < 75) hUDSJetsQGL_70pt75 -> Fill(QGL);
	else if (pt >= 75 && pt < 80) hUDSJetsQGL_75pt80 -> Fill(QGL);
	else if (pt >= 80 && pt < 90) hUDSJetsQGL_80pt90 -> Fill(QGL);
	else if (pt >= 90 && pt < 100) hUDSJetsQGL_90pt100 -> Fill(QGL);
	else if (pt >= 100 && pt < 110) hUDSJetsQGL_100pt110 -> Fill(QGL);
	else if (pt >= 110 && pt < 120) hUDSJetsQGL_110pt120 -> Fill(QGL);
	else if (pt >= 120 && pt < 140) hUDSJetsQGL_120pt140 -> Fill(QGL);
	else if (pt >= 140 && pt < 160) hUDSJetsQGL_140pt160 -> Fill(QGL);
	else if (pt >= 160 && pt < 180) hUDSJetsQGL_160pt180 -> Fill(QGL);
	else if (pt >= 180 && pt < 200) hUDSJetsQGL_180pt200 -> Fill(QGL);
	else if (pt >= 200 && pt < 250) hUDSJetsQGL_200pt250 -> Fill(QGL);
	else if (pt >= 250 && pt < 300) hUDSJetsQGL_250pt300 -> Fill(QGL);
	else if (pt >= 300 && pt < 350) hUDSJetsQGL_300pt350 -> Fill(QGL);
	else if (pt >= 350 && pt < 400) hUDSJetsQGL_350pt400 -> Fill(QGL);
	else if (pt >= 400 && pt < 450) hUDSJetsQGL_400pt450 -> Fill(QGL);
	else if (pt >= 450 && pt < 500) hUDSJetsQGL_450pt500 -> Fill(QGL);
	else if (pt >= 500 && pt < 550) hUDSJetsQGL_500pt550 -> Fill(QGL);
	else if (pt >= 550 && pt < 600) hUDSJetsQGL_550pt600 -> Fill(QGL);
	else if (pt >= 600 && pt < 700) hUDSJetsQGL_600pt700 -> Fill(QGL);
	else if (pt >= 700 && pt < 800) hUDSJetsQGL_700pt800 -> Fill(QGL);
	else if (pt >= 800 && pt < 1000) hUDSJetsQGL_800pt1000 -> Fill(QGL);
	else if (pt >= 1000) hUDSJetsQGL_1000ptInf -> Fill(QGL);
      }
  }
  hCtrlNCJets -> Fill(nCJets);
  hCtrlNBJets -> Fill(nBJets);
  hCtrlNUDSGJets -> Fill(nUDSGJets);
  hCtrlNGJets -> Fill(nGJets);
  hCtrlNUDSJets -> Fill(nUDSJets);
  
  //================================================================================================
  // Finalize
  //================================================================================================
  fEventSaver.save();
  
  return;
}
