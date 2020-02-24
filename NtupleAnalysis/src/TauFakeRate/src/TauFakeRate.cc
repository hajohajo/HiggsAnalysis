// -*- c++ -*-
#include "Framework/interface/BaseSelector.h"
#include "Framework/interface/makeTH.h"

#include "EventSelection/interface/CommonPlots.h"
#include "EventSelection/interface/EventSelections.h"
#include "EventSelection/interface/TransverseMass.h"


#include "TDirectory.h"

class TauFakeRate: public BaseSelector {
public:
  explicit TauFakeRate(const ParameterSet& config, const TH1* skimCounters);
  virtual ~TauFakeRate() {}

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
  Count cLeptonOSCounter;
  Count cLeptonMassCounter;
  TauSelection fTauSelection;
  TauSelection fLooseTauSelection;
  Count cTauNCounter;
  Count cTauSFCounter;
  Count cFakeTauSFCounter;
  JetSelection fJetSelection;
  BJetSelection fBJetSelection;
  Count cBTaggingSFCounter;
  METSelection fMETSelection;
  // FatJetSelection fFatJetSelection;
  Count cSelected;
    
  void doLooseTaus(const Event& event, const TauSelection::Data& tauData, const int nVertices);
  void doTightTaus(const Event& event, const TauSelection::Data& looseTauData, const int nVertices);

  // Non-common histograms
  WrappedTH1 *hTauPt_num_dm0;
  WrappedTH1 *hTauPt_num_dm1;
  WrappedTH1 *hTauPt_num_dm10;
  WrappedTH1 *hTauPt_den_dm0;
  WrappedTH1 *hTauPt_den_dm1;
  WrappedTH1 *hTauPt_den_dm10;
  WrappedTH1 *hTauPt_num_g_dm0;
  WrappedTH1 *hTauPt_num_g_dm1;
  WrappedTH1 *hTauPt_num_g_dm10;
  WrappedTH1 *hTauPt_den_g_dm0;
  WrappedTH1 *hTauPt_den_g_dm1;
  WrappedTH1 *hTauPt_den_g_dm10;
  WrappedTH1 *hTauEta_num_dm0;
  WrappedTH1 *hTauEta_num_dm1;
  WrappedTH1 *hTauEta_num_dm10;
  WrappedTH1 *hTauEta_den_dm0;
  WrappedTH1 *hTauEta_den_dm1;
  WrappedTH1 *hTauEta_den_dm10;

  // Barrel ( |eta| < 1.5)
  WrappedTH1 *hTauPt_num_dm0_barrel;
  WrappedTH1 *hTauPt_num_dm1_barrel;
  WrappedTH1 *hTauPt_num_dm10_barrel;
  WrappedTH1 *hTauPt_den_dm0_barrel;
  WrappedTH1 *hTauPt_den_dm1_barrel;
  WrappedTH1 *hTauPt_den_dm10_barrel;
  WrappedTH1 *hTauPt_num_g_dm0_barrel;
  WrappedTH1 *hTauPt_num_g_dm1_barrel;
  WrappedTH1 *hTauPt_num_g_dm10_barrel;
  WrappedTH1 *hTauPt_den_g_dm0_barrel;
  WrappedTH1 *hTauPt_den_g_dm1_barrel;
  WrappedTH1 *hTauPt_den_g_dm10_barrel;
  WrappedTH1 *hTauEta_num_dm0_barrel;
  WrappedTH1 *hTauEta_num_dm1_barrel;
  WrappedTH1 *hTauEta_num_dm10_barrel;
  WrappedTH1 *hTauEta_den_dm0_barrel;
  WrappedTH1 *hTauEta_den_dm1_barrel;
  WrappedTH1 *hTauEta_den_dm10_barrel;

  // Endcap ( |eta| >= 1.5)
  WrappedTH1 *hTauPt_num_dm0_endcap;
  WrappedTH1 *hTauPt_num_dm1_endcap;
  WrappedTH1 *hTauPt_num_dm10_endcap;
  WrappedTH1 *hTauPt_den_dm0_endcap;
  WrappedTH1 *hTauPt_den_dm1_endcap;
  WrappedTH1 *hTauPt_den_dm10_endcap;
  WrappedTH1 *hTauPt_num_g_dm0_endcap;
  WrappedTH1 *hTauPt_num_g_dm1_endcap;
  WrappedTH1 *hTauPt_num_g_dm10_endcap;
  WrappedTH1 *hTauPt_den_g_dm0_endcap;
  WrappedTH1 *hTauPt_den_g_dm1_endcap;
  WrappedTH1 *hTauPt_den_g_dm10_endcap;
  WrappedTH1 *hTauEta_num_dm0_endcap;
  WrappedTH1 *hTauEta_num_dm1_endcap;
  WrappedTH1 *hTauEta_num_dm10_endcap;
  WrappedTH1 *hTauEta_den_dm0_endcap;
  WrappedTH1 *hTauEta_den_dm1_endcap;
  WrappedTH1 *hTauEta_den_dm10_endcap;

  WrappedTH1 *hDileptonMass_BeforeOnZSelection;

  WrappedTH1 *hLeptonN_AfterLeptonSelection;
  WrappedTH1 *hLeptonPt_AfterLeptonSelection;
  WrappedTH1 *hLeptonEta_AfterLeptonSelection;
  WrappedTH1 *hDileptonPt_AfterLeptonSelection;
  WrappedTH1 *hDileptonEta_AfterLeptonSelection;
  WrappedTH1 *hDileptonMass_AfterLeptonSelection;
  
  WrappedTH1 *hTauSrc_AfterTauSelection;
  WrappedTH1 *hNTau_AfterTauSelection;
  WrappedTH1 *hNTau_AfterJetSelection;
  WrappedTH1 *hNTau_AfterBJetSelection;
  WrappedTH1 *hNTau_AfterMetSelection;

  WrappedTH1 *hNTau;
  WrappedTH1 *hTauPt;
  WrappedTH1 *hTauEta;

  WrappedTH1 *hDileptonMass_AfterTauSelection;
  WrappedTH1 *hDileptonMass_AfterJetSelection;
  WrappedTH1 *hDileptonMass_AfterBJetSelection;
  WrappedTH1 *hDileptonMass_AfterMetSelection;

  WrappedTH1 *hNTau_AfterAllSelections;
  WrappedTH1 *hNJet_AfterAllSelections;
  WrappedTH1 *hNBjet_AfterAllSelections;
  WrappedTH1 *hMet_AfterAllSelections;
  WrappedTH1 *hDileptonPt_AfterAllSelections;
  WrappedTH1 *hDileptonEta_AfterAllSelections;
  WrappedTH1 *hDileptonMass_AfterAllSelections;

};

#include "Framework/interface/SelectorFactory.h"
REGISTER_SELECTOR(TauFakeRate);

TauFakeRate::TauFakeRate(const ParameterSet& config, const TH1* skimCounters)
  : BaseSelector(config, skimCounters),
    fCommonPlots(config.getParameter<ParameterSet>("CommonPlots"), CommonPlots::kHplus2hwAnalysisWithTop, fHistoWrapper),
    cAllEvents(fEventCounter.addCounter("all events")),
    cTrigger(fEventCounter.addCounter("passed trigger")),
    fMETFilterSelection(config.getParameter<ParameterSet>("METFilter"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cVertexSelection(fEventCounter.addCounter("passed PV")),
    fElectronSelection(config.getParameter<ParameterSet>("ElectronSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    fMuonSelection(config.getParameter<ParameterSet>("MuonSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cLeptonOSCounter(fEventCounter.addCounter("#mu OS")),
    cLeptonMassCounter(fEventCounter.addCounter("m_{#mu#mu} window")),
    // fTauSelection(config.getParameter<ParameterSet>("TauSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fTauSelection(config.getParameter<ParameterSet>("TauSelection")),
    fLooseTauSelection(config.getParameter<ParameterSet>("LooseTauSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cTauNCounter(fEventCounter.addCounter("#geq 1 loose #tau")),
    cTauSFCounter(fEventCounter.addCounter("#tau SF")),
    cFakeTauSFCounter(fEventCounter.addCounter("Fake #tau SF")),
    fJetSelection(config.getParameter<ParameterSet>("JetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fBJetSelection(config.getParameter<ParameterSet>("BJetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cBTaggingSFCounter(fEventCounter.addCounter("b-tag SF")),
    // fMETSelection(config.getParameter<ParameterSet>("METSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fMETSelection(config.getParameter<ParameterSet>("METSelection")),
    // fFatJetSelection(config.getParameter<ParameterSet>("FatJetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    cSelected(fEventCounter.addCounter("Selected Events"))
{ }


void TauFakeRate::book(TDirectory *dir) {

  if (0) std::cout << "=== TauFakeRate::book()" << std::endl;
  // Book common plots histograms
  fCommonPlots.book(dir, isData());

  // Book histograms in event selection classes
  fMETFilterSelection.bookHistograms(dir);
  fElectronSelection.bookHistograms(dir);
  fMuonSelection.bookHistograms(dir);
  fTauSelection.bookHistograms(dir);
  fLooseTauSelection.bookHistograms(dir);
  fJetSelection.bookHistograms(dir);
  fBJetSelection.bookHistograms(dir);
  fMETSelection.bookHistograms(dir);
  // fFatJetSelection.bookHistograms(dir);

  // Get binning from cfg file
  const int   ptN    = fCommonPlots.getPtBinSettings().bins();
  const float ptMin  = fCommonPlots.getPtBinSettings().min();
  const float ptMax  = fCommonPlots.getPtBinSettings().max();
  const int   etaN   = fCommonPlots.getEtaBinSettings().bins();
  const float etaMin = fCommonPlots.getEtaBinSettings().min();
  const float etaMax = fCommonPlots.getEtaBinSettings().max();
  const int   mN     = fCommonPlots.getInvMassBinSettings().bins();
  const float mMin   = fCommonPlots.getInvMassBinSettings().min();
  const float mMax   = fCommonPlots.getInvMassBinSettings().max();
  const int   nN     = fCommonPlots.getNjetsBinSettings().bins();
  const float nMin   = fCommonPlots.getNjetsBinSettings().min();
  const float nMax   = fCommonPlots.getNjetsBinSettings().max();
  const int   metN   = fCommonPlots.getMetBinSettings().bins();
  const float metMin = fCommonPlots.getMetBinSettings().min();
  const float metMax = fCommonPlots.getMetBinSettings().max();
  
  // Book non-common histograms 
  double bin[8] = {20,25,30,35,40,50,60,120}; // iro-fixme

  hTauPt_num_dm0 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm0", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_dm1 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm1", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_dm10 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm10", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm0 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm0", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm1 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm1", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm10 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm10", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm0 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm0", ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm1 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm1", ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm10 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm10", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm0 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm0", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm1 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm1", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm10 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm10", ";p_{T} (GeV)", 7, bin);
  hTauEta_num_dm0 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm0", ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_dm1 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm1", ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_dm10 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm10", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm0 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm0", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm1 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm1", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm10 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm10", ";#eta", etaN, etaMin, etaMax);

  // Barrel ( |eta| < 1.5)
  hTauPt_num_dm0_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm0_barrel", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_dm1_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm1_barrel", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_dm10_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm10_barrel", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm0_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm0_barrel", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm1_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm1_barrel", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm10_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm10_barrel", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm0_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm0_barrel", ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm1_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm1_barrel", ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm10_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm10_barrel", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm0_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm0_barrel", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm1_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm1_barrel", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm10_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm10_barrel", ";p_{T} (GeV)", 7, bin);
  hTauEta_num_dm0_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm0_barrel", ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_dm1_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm1_barrel", ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_dm10_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm10_barrel", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm0_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm0_barrel", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm1_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm1_barrel", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm10_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm10_barrel", ";#eta", etaN, etaMin, etaMax);

  // Endcap ( |eta| >= 1.5)
  hTauPt_num_dm0_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm0_endcap", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_dm1_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm1_endcap", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_dm10_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm10_endcap", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm0_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm0_endcap", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm1_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm1_endcap", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_dm10_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm10_endcap", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm0_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm0_endcap", ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm1_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm1_endcap", ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm10_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm10_endcap", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm0_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm0_endcap", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm1_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm1_endcap", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm10_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm10_endcap", ";p_{T} (GeV)", 7, bin);
  hTauEta_num_dm0_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm0_endcap", ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_dm1_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm1_endcap", ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_dm10_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm10_endcap", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm0_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm0_endcap", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm1_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm1_endcap", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm10_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm10_endcap", ";#eta", etaN, etaMin, etaMax);


  hDileptonMass_BeforeOnZSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_BeforeOnZSelection", ";m_{ll} (GeV)", mN, mMin, mMax);

  hLeptonN_AfterLeptonSelection      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LeptonN_AfterLeptonSelection"     , ";lepton multiplicity", nN  , nMin  , nMax );
  hLeptonPt_AfterLeptonSelection     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LeptonPt_AfterLeptonSelection"    , ";p_{T} (GeV)", ptN , ptMin , ptMax );
  hLeptonEta_AfterLeptonSelection    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LeptonEta_AfterLeptonSelection"   , ";#eta"       , etaN, etaMin, etaMax);
  hDileptonPt_AfterLeptonSelection   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonPt_AfterLeptonSelection"  , ";p_{T} (GeV)", ptN , ptMin , ptMax );
  hDileptonEta_AfterLeptonSelection  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonEta_AfterLeptonSelection" , ";#eta"       , etaN, etaMin, etaMax);
  hDileptonMass_AfterLeptonSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_AfterLeptonSelection", ";m_{ll} (GeV)", mN , mMin  , mMax  );

  hDileptonMass_AfterTauSelection    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_AfterTauSelection"  , ";m_{ll} (GeV)", mN, mMin, mMax);
  hDileptonMass_AfterJetSelection    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_AfterJetSelection"  , ";m_{ll} (GeV)", mN, mMin, mMax);
  hDileptonMass_AfterBJetSelection   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_AfterBJetSelection" , ";m_{ll} (GeV)", mN, mMin, mMax);
  hDileptonMass_AfterMetSelection    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_AfterMetSelection"  , ";m_{ll} (GeV)", mN, mMin, mMax);

  hTauSrc_AfterTauSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrc_AfterTauSelection", ";#tau source", nN, nMin, nMax  );
  hNTau_AfterTauSelection   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_AfterTauSelection"  , "; loose #tau-jet multiplicity" , nN  , nMin  , nMax  );
  hNTau_AfterJetSelection   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_AfterJetSelection"  , "; loose #tau-jet multiplicity" , nN  , nMin  , nMax  );
  hNTau_AfterBJetSelection  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_AfterBJetSelection" , "; loose #tau-jet multiplicity" , nN  , nMin  , nMax  );
  hNTau_AfterMetSelection   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_AfterMetSelection"  , "; loose #tau-jet multiplicity" , nN  , nMin  , nMax  );

  hNTau_AfterAllSelections         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_AfterAllSelections" , "; loose #tau-jet multiplicity" , nN  , nMin  , nMax  );
  hNJet_AfterAllSelections         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NJet_AfterAllSelections" , "; jet multiplicity"  , nN  , nMin  , nMax  );
  hNBjet_AfterAllSelections        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NBjet_AfterAllSelections", "; b-jet multiplicity", nN  , nMin  , nMax  );
  hMet_AfterAllSelections          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "Met_AfterAllSelections"  , "; E_{T}^{miss} (GeV)", metN, metMin, metMax);
  hDileptonPt_AfterAllSelections   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonPt_AfterAllSelections"  , ";p_{T} (GeV)" , ptN , ptMin , ptMax );
  hDileptonEta_AfterAllSelections  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonEta_AfterAllSelections" , ";#eta"        , etaN, etaMin, etaMax);
  hDileptonMass_AfterAllSelections = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_AfterAllSelections", ";m_{ll} (GeV)", mN  , mMin  , mMax  );

  return;
}


void TauFakeRate::setupBranches(BranchManager& branchManager) {
  fEvent.setupBranches(branchManager);
  return;
}


void TauFakeRate::process(Long64_t entry) {

  //====== Initialize
  fCommonPlots.initialize();
  fCommonPlots.setFactorisationBinForEvent(std::vector<float> {});
  cAllEvents.increment();

  //================================================================================================   
  // Apply trigger 
  //================================================================================================   
  if (0) std::cout << "=== Trigger" << std::endl;
  if ( !(fEvent.passTriggerDecision()) ) return;  
  cTrigger.increment();
  int nVertices = fEvent.vertexInfo().value();
  fCommonPlots.setNvertices(nVertices);

  // Fill histos
  fCommonPlots.fillControlPlotsAfterTrigger(fEvent);


  //================================================================================================   
  // MET filters (to remove events with spurious sources of fake MET)
  //================================================================================================   
  if (0) std::cout << "=== MET Filter" << std::endl;
  const METFilterSelection::Data metFilterData = fMETFilterSelection.analyze(fEvent);
  if (!metFilterData.passedSelection()) return;

  // Fill histos
  fCommonPlots.fillControlPlotsAfterMETFilter(fEvent);  


  //================================================================================================   
  // Primarty Vertex (Check that a PV exists)
  //================================================================================================   
  if (0) std::cout << "=== Vertices" << std::endl;
  if (nVertices < 1) return;
  cVertexSelection.increment();

  // Fill histos
  fCommonPlots.fillControlPlotsAtVertexSelection(fEvent);

  
  //================================================================================================   
  // Electron Selection
  //================================================================================================   
  if (0) std::cout << "=== Electron veto" << std::endl;
  const ElectronSelection::Data eData = fElectronSelection.analyze(fEvent);
  if (eData.hasIdentifiedElectrons()) return;

  // Fill histos
  fCommonPlots.fillControlPlotsAtElectronSelection(fEvent, eData);


  //================================================================================================
  // Muon Selection
  //================================================================================================
  if (0) std::cout << "=== Muon Selection" << std::endl;
  const MuonSelection::Data muData = fMuonSelection.analyze(fEvent);
  if(!(muData.hasIdentifiedMuons())) return;

  // Require exactly 2 muons
  if (muData.getSelectedMuons().size() != 2) return; // note: remember to disable trigger-matching option if using a single muon trigger

  // Calculate variables for dilepton system
  double dilepton_invMass = ( muData.getSelectedMuons()[0].p4() + muData.getSelectedMuons()[1].p4() ).M();
  double dilepton_pt      = ( muData.getSelectedMuons()[0].p4() + muData.getSelectedMuons()[1].p4() ).pt();
  double dilepton_eta     = ( muData.getSelectedMuons()[0].p4() + muData.getSelectedMuons()[1].p4() ).eta();
  const double cfg_massWindow = 20.0; // 15.0;

  // Apply muon trigger scale factors
  // https://gitlab.cern.ch/HPlus/HiggsAnalysis/blob/master/NtupleAnalysis/python/parameters/scaleFactors.py
  // https://gitlab.cern.ch/HPlus/HiggsAnalysis/blob/master/NtupleAnalysis/data/TriggerEfficiency/muonPAGEff.json
  if (fEvent.isMC()) 
    {
      fEventWeight.multiplyWeight(muData.getMuonTriggerSF());
    }

  // Apply muon ID scale factors
  if (fEvent.isMC()) 
    {
      fEventWeight.multiplyWeight(muData.getMuonIDSF());
    }
  

  // Require the muons to have opposite sign (OS)
  if( muData.getSelectedMuons()[0].charge() == muData.getSelectedMuons()[1].charge()) return;
  cLeptonOSCounter.increment(); 
  if (0) std::cout << "=== Muon Selection: OS requirement" << std::endl;

  hDileptonMass_BeforeOnZSelection->Fill(dilepton_invMass);

  // Apply on-Z mass requirement for dilepton pair
  if (dilepton_invMass < ( 91.1876 - cfg_massWindow ) ) return; // fixme. define in run.py 
  if (dilepton_invMass > ( 91.1876 + cfg_massWindow ) ) return; // fixme. define in run.py
  cLeptonMassCounter.increment(); 

  // Apply dR cut for dilepton system?
  //  dilepton_dR = ROOT::Math::VectorUtil::DeltaR(muData.getSelectedMuons()[0].p4(),muData.getSelectedMuons()[1].p4());
  //  if(dilepton_dR < 0.3) return;


  // Fill histos
  fCommonPlots.fillControlPlotsAtMuonSelection(fEvent, muData);
  hDileptonPt_AfterLeptonSelection ->Fill(dilepton_pt);
  hDileptonEta_AfterLeptonSelection->Fill(dilepton_eta);
  hDileptonMass_AfterLeptonSelection->Fill(dilepton_invMass);
  hLeptonN_AfterLeptonSelection ->Fill(muData.getSelectedMuons().size());
  // For-loop: All Selected leptons
  for(unsigned int i=0; i< muData.getSelectedMuons().size(); i++)
    {
      hLeptonPt_AfterLeptonSelection ->Fill(muData.getSelectedMuons()[i].pt());
      hLeptonEta_AfterLeptonSelection->Fill(muData.getSelectedMuons()[i].eta());
    }

  

  //================================================================================================   
  // Tau Selection
  //================================================================================================   
  if (0) std::cout << "=== Tau Selection" << std::endl;
  const TauSelection::Data tauData = fTauSelection.analyze(fEvent);
  const TauSelection::Data looseTauData = fLooseTauSelection.analyzeLoose(fEvent);

  if (!looseTauData.hasIdentifiedTaus() ) return;
  // if(looseTauData.getSelectedTaus().size() != 1) return; // same as line above (double-check)

  if (0) std::cout << "=== Tau Selection: Has at least 1 loose tau(s)" << std::endl;
  cTauNCounter.increment(); 


  // Determine if genuine or fake tau. If more than one tau present use the first one in the list (pt-sorted)
  bool bGenuineTau      = false; // 
  bool bFakeTauElectron = false; // electron -> tau fakes only
  bool bFakeTauMuon     = false; // muon     -> tau fakes only
  bool bFakeTauJets     = false; // jet      -> tau fakes only
  if ( fEvent.isMC() )
    {
      bGenuineTau      = looseTauData.isGenuineTau();
      bFakeTauElectron = looseTauData.getSelectedTau().isElectronToTau();
      bFakeTauMuon     = looseTauData.getSelectedTau().isMuonToTau();
    }
  bFakeTauJets = (!bGenuineTau && !(bFakeTauElectron || bFakeTauMuon) );
  int isGenuineTau = int(!(bFakeTauJets));
  int tauSrcBit    = pow(2,0)*bFakeTauElectron + pow(2,1)*bFakeTauMuon + pow(2,2)*bFakeTauJets + pow(2,3)*(bGenuineTau);

  // Overwrite default boolean bIsGenuineTau = (data.getSelectedTaus()[0].isGenuineTau() || data.getSelectedTaus()[1].isGenuineTau()) 
  // [NOTE: Filled when  "fillControlPlotsAfterTauSelection" is called()]
  if (0) fCommonPlots.setGenuineTauStatus(isGenuineTau);
  if (0) std::cout << "=== Tau Selection: isGenuineTau = " << isGenuineTau << " (src = " << tauSrcBit << ")" << std::endl;

  // Fill histos ( also sets value for boolean bIsGenuineTau
  hDileptonMass_AfterTauSelection->Fill(dilepton_invMass);
  hNTau_AfterTauSelection->Fill( looseTauData.getSelectedTaus().size() );
  hTauSrc_AfterTauSelection->Fill(tauSrcBit);
  fCommonPlots.fillControlPlotsAfterTauSelection(fEvent, looseTauData);

  // xenios

  //================================================================================================
  // Jet selection
  //================================================================================================
  if (0) std::cout << "=== Jet selection" << std::endl;
  const JetSelection::Data jetData = fJetSelection.analyze(fEvent, looseTauData.getSelectedTau());
  if (!jetData.passedSelection()) return;  
    
  // Fill histos
  fCommonPlots.fillControlPlotsAtJetSelection(fEvent, jetData);
  hDileptonMass_AfterJetSelection->Fill(dilepton_invMass);
  hNTau_AfterJetSelection->Fill( looseTauData.getSelectedTaus().size() );

 
  //================================================================================================  
  // BJet selection
  //================================================================================================
  if (0) std::cout << "=== BJet selection" << std::endl;
  const BJetSelection::Data bjetData = fBJetSelection.analyze(fEvent, jetData);
  if (!bjetData.passedSelection()) return;
  
  // Fill histos
  fCommonPlots.fillControlPlotsAtBtagging(fEvent, bjetData);

  if (fEvent.isMC()) 
    {
      fEventWeight.multiplyWeight(bjetData.getBTaggingScaleFactorEventWeight());
    }
  cBTaggingSFCounter.increment();

  // Fill histos
  fCommonPlots.fillControlPlotsAfterBtagSF(fEvent, jetData, bjetData);
  hDileptonMass_AfterBJetSelection->Fill(dilepton_invMass);
  hNTau_AfterBJetSelection->Fill( looseTauData.getSelectedTaus().size() );


  //================================================================================================
  // MET selection
  //================================================================================================
  if (0) std::cout << "=== MET selection" << std::endl;
  // const METSelection::Data metData = fMETSelection.silentAnalyze(fEvent, nVertices);
  const METSelection::Data metData = fMETSelection.analyze(fEvent, nVertices);
  if (!metData.passedSelection()) return;

  // Fill Histos
  fCommonPlots.fillControlPlotsAtMETSelection(fEvent, metData);
  hDileptonMass_AfterMetSelection->Fill(dilepton_invMass);
  hNTau_AfterMetSelection->Fill( looseTauData.getSelectedTaus().size() );

  //================================================================================================
  // All Selections
  //================================================================================================
  if (0) std::cout << "=== All Selections" << std::endl;
  cSelected.increment();

  // Fill histos
  fCommonPlots.fillControlPlotsAfterAllSelections(fEvent);

  hNTau_AfterAllSelections->Fill( looseTauData.getSelectedTaus().size() );
  hNJet_AfterAllSelections ->Fill( jetData.getSelectedJets().size() );
  hNBjet_AfterAllSelections->Fill( bjetData.getSelectedBJets().size() );
  hDileptonPt_AfterAllSelections->Fill(dilepton_pt);
  hDileptonEta_AfterAllSelections->Fill(dilepton_eta);
  hDileptonMass_AfterAllSelections->Fill(dilepton_invMass);
  hMet_AfterAllSelections->Fill(metData.getMET().R());

  // Tau stuff here
  if(looseTauData.getSelectedTaus().size() != 1) return;


  // "Tight" Tau
  if (tauData.hasIdentifiedTaus()) 
    {

      if (fEvent.isMC()) 
	{
	  // Apply "tight" tau ID scale factor (SF)
	  fEventWeight.multiplyWeight(tauData.getTauIDSF());
	  cTauSFCounter.increment(); 
	  
	  // Apply "tight" fake tau SF
	  fEventWeight.multiplyWeight(tauData.getTauMisIDSF());
	  cFakeTauSFCounter.increment();
	}

      // Do rest of event selection
      doTightTaus(fEvent, tauData, nVertices);

      if (fEvent.isMC())
	{ 
	  // Undo tau ID SF!
	  fEventWeight.multiplyWeight(1.0/tauData.getTauIDSF());

	  // Undo tau mis-ID SF!
	  fEventWeight.multiplyWeight(1.0/tauData.getTauMisIDSF());
	}
    }

  // "Loose" Tau
  if (looseTauData.hasIdentifiedTaus()) 
    {
      
      if (fEvent.isMC())
	{
	  // Apply "loose" tau ID SF
	  fEventWeight.multiplyWeight(looseTauData.getTauIDSF());
	  // Apply "loose" tau mis-ID SF
	  fEventWeight.multiplyWeight(looseTauData.getTauMisIDSF());
	}

      // Do rest of event selection
      doLooseTaus(fEvent, looseTauData, nVertices);
    }

  fEventSaver.save();

  return;
}


void TauFakeRate::doLooseTaus(const Event& event, const TauSelection::Data& tauData, const int nVertices) {


  // For-loop: All Selected taus
  for(unsigned int i=0; i<tauData.getSelectedTaus().size(); i++)
    {      
      
      // https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
      // DM = 5*(Nc-1)+Np  [Nc = # charged hadrons, Np = # of pi0s]. DN = 0-11 
      unsigned int tau_dm = tauData.getSelectedTaus()[i].decayMode();
      double tau_pt       = tauData.getSelectedTaus()[i].pt();
      double tau_eta      = tauData.getSelectedTaus()[i].eta();
      bool tau_inBarrel   = ( fabs(tau_eta) < 1.5);

      // Only if MC  and selected tau is genuine (not fake)
      if (event.isMC() && tauData.getSelectedTaus()[i].isGenuineTau()) 
	{
	  if (tau_dm == 0)  
	    {
	      hTauPt_den_g_dm0->Fill( tau_pt );
	      if (tau_inBarrel) hTauPt_den_g_dm0_barrel->Fill( tau_pt );
	      else hTauPt_den_g_dm0_endcap->Fill( tau_pt );
	    }
	  
	  if (tau_dm==1)
	    {
	      hTauPt_den_g_dm1->Fill( tau_pt );
	      if (tau_inBarrel) hTauPt_den_g_dm1_barrel->Fill( tau_pt );
	      else hTauPt_den_g_dm1_endcap->Fill( tau_pt );
	    }
	  
	  if (tau_dm==10) 
	    {
	      hTauPt_den_g_dm10->Fill( tau_pt );
	      if (tau_inBarrel) hTauPt_den_g_dm10_barrel->Fill( tau_pt );
	      else hTauPt_den_g_dm10_endcap->Fill( tau_pt );
	    }
	}
    
      // 1-prong decays; decay mode (DM) = 0
      if (tau_dm==0) 
	{
	  hTauPt_den_dm0 ->Fill( tau_pt  );
	  hTauEta_den_dm0->Fill( tau_eta );

	  if (tau_inBarrel)
	    {
	      hTauPt_den_dm0_barrel ->Fill( tau_pt  );
	      hTauEta_den_dm0_barrel->Fill( tau_eta );
	    }
	  else
	    {
	      hTauPt_den_dm0_endcap ->Fill( tau_pt  );
	      hTauEta_den_dm0_endcap->Fill( tau_eta );
	    }
	}
      
      // 2-prong decays; decay mode (DM) = 1
      if (tau_dm==1) 
	{
	  hTauPt_den_dm1 ->Fill( tau_pt  );
	  hTauEta_den_dm1->Fill( tau_eta );
	  
	  if (tau_inBarrel)
	    {
	      hTauPt_den_dm1_barrel ->Fill( tau_pt  );
	      hTauEta_den_dm1_barrel->Fill( tau_eta );
	    }
	  else
	    {
	      hTauPt_den_dm1_endcap ->Fill( tau_pt  );
	      hTauEta_den_dm1_endcap->Fill( tau_eta );
	    }
	}
      
      // 3-prong decays; decay mode (DM) = 10
      if (tau_dm==10) 
	{
	  hTauPt_den_dm10 ->Fill( tau_pt  );
	  hTauEta_den_dm10->Fill( tau_eta );
	  
	  if (tau_inBarrel)
	    {
	      hTauPt_den_dm10_barrel ->Fill( tau_pt  );
	      hTauEta_den_dm10_barrel->Fill( tau_eta );
	    }
	  else
	    {
	      hTauPt_den_dm10_endcap ->Fill( tau_pt  );
	      hTauEta_den_dm10_endcap->Fill( tau_eta );
	    }
	}
    }
  
  return;
}


void TauFakeRate::doTightTaus(const Event& event, const TauSelection::Data& tightTauData, const int nVertices) {

  // For-loop: All selected taus
  for(unsigned int i=0; i<tightTauData.getSelectedTaus().size(); i++)
    {
      
      // https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
      // DM = 5*(Nc-1)+Np  [Nc = # charged hadrons, Np = # of pi0s]. DN = 0-11 
      unsigned int tau_dm = tightTauData.getSelectedTaus()[i].decayMode();
      double tau_pt       = tightTauData.getSelectedTaus()[i].pt();
      double tau_eta      = tightTauData.getSelectedTaus()[i].eta();
      bool tau_inBarrel   = ( fabs(tau_eta) < 1.5);

      // Only if MC  and selected tau is genuine (not fake)
      if (event.isMC() && tightTauData.getSelectedTaus()[i].isGenuineTau()) 
	{
	  if (tau_dm==0) 
	    {
	      hTauPt_num_g_dm0->Fill( tau_pt );
	      if (tau_inBarrel) hTauPt_num_g_dm0_barrel->Fill( tau_pt );
	      else hTauPt_num_g_dm0_endcap->Fill( tau_pt );
	    }

	  if (tau_dm==1)
	    {
	      hTauPt_num_g_dm1->Fill( tau_pt );
	      if (tau_inBarrel) hTauPt_num_g_dm1_barrel->Fill( tau_pt );
	      else hTauPt_num_g_dm1_endcap->Fill( tau_pt );	   
	    }

	  if (tau_dm==10) 
	    {
	      hTauPt_num_g_dm10->Fill( tau_pt );
	      if (tau_inBarrel) hTauPt_num_g_dm10_barrel->Fill( tau_pt );
	      else hTauPt_num_g_dm10_endcap->Fill( tau_pt );
	    }
	}
      
      // 1-prong decays; decay mode (DM) = 0
      if (tau_dm==0) 
	{
	  hTauPt_num_dm0 ->Fill( tau_pt  );
	  hTauEta_num_dm0->Fill( tau_eta );

	  if (tau_inBarrel) 
	    {
	      hTauPt_num_dm0_barrel ->Fill( tau_pt  );
	      hTauEta_num_dm0_barrel->Fill( tau_eta );
	    }
	  else
	    {
	      hTauPt_num_dm0_endcap ->Fill( tau_pt  );
	      hTauEta_num_dm0_endcap->Fill( tau_eta );
	    }
	}
      
      // 2-prong decays; decay mode (DM) = 1
      if (tau_dm==1) 
	{
	  hTauPt_num_dm1 ->Fill( tau_pt  );
	  hTauEta_num_dm1->Fill( tau_eta );
	  if (tau_inBarrel) 
	    {
	      hTauPt_num_dm1_barrel ->Fill( tau_pt  );
	      hTauEta_num_dm1_barrel->Fill( tau_eta );
	    }
	  else
	    {
	      hTauPt_num_dm1_endcap ->Fill( tau_pt  );
	      hTauEta_num_dm1_endcap->Fill( tau_eta );
	    }
	}
      
      // 3-prong decays; decay mode (DM) = 10
      if (tau_dm==10) 
	{
	  hTauPt_num_dm10 ->Fill( tau_pt  );
	  hTauEta_num_dm10->Fill( tau_eta );
	  if (tau_inBarrel) 
	    {
	      hTauPt_num_dm10_barrel ->Fill( tau_pt  );
	      hTauEta_num_dm10_barrel->Fill( tau_eta );
	    }
	  else
	    {
	      hTauPt_num_dm10_endcap ->Fill( tau_pt  );
	      hTauEta_num_dm10_endcap->Fill( tau_eta );
	    }
	}
    }
  
  return;
}
