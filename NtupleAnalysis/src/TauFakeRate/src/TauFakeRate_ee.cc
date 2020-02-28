// -*- c++ -*-
#include "Framework/interface/BaseSelector.h"
#include "Framework/interface/makeTH.h"

#include "EventSelection/interface/CommonPlots.h"
#include "EventSelection/interface/EventSelections.h"
#include "EventSelection/interface/TransverseMass.h"


#include "TDirectory.h"

class TauFakeRate_ee: public BaseSelector {
public:
  explicit TauFakeRate_ee(const ParameterSet& config, const TH1* skimCounters);
  virtual ~TauFakeRate_ee() {}

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
    
  void doLooseTaus(const Event& event, const TauSelection::Data& tauData);
  void doTightTaus(const Event& event, const TauSelection::Data& looseTauData);

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

  WrappedTH1 *hTauSrc_AfterAllSelections;
  WrappedTH1 *hNTau_AfterAllSelections;
  WrappedTH1 *hNJet_AfterAllSelections;
  WrappedTH1 *hNBjet_AfterAllSelections;
  WrappedTH1 *hMet_AfterAllSelections;
  WrappedTH1 *hDileptonPt_AfterAllSelections;
  WrappedTH1 *hDileptonEta_AfterAllSelections;
  WrappedTH1 *hDileptonMass_AfterAllSelections;

};

#include "Framework/interface/SelectorFactory.h"
REGISTER_SELECTOR(TauFakeRate_ee);

TauFakeRate_ee::TauFakeRate_ee(const ParameterSet& config, const TH1* skimCounters)
  : BaseSelector(config, skimCounters),
    fCommonPlots(config.getParameter<ParameterSet>("CommonPlots"), CommonPlots::kHplus2hwAnalysisWithTop, fHistoWrapper),
    cAllEvents(fEventCounter.addCounter("all events")),
    cTrigger(fEventCounter.addCounter("passed trigger")),
    fMETFilterSelection(config.getParameter<ParameterSet>("METFilter"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cVertexSelection(fEventCounter.addCounter("passed PV")),
    fElectronSelection(config.getParameter<ParameterSet>("ElectronSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fMuonSelection(config.getParameter<ParameterSet>("MuonSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    cLeptonOSCounter(fEventCounter.addCounter("e OS")),
    cLeptonMassCounter(fEventCounter.addCounter("m_{ee} window")),
    fTauSelection(config.getParameter<ParameterSet>("TauSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    // fTauSelection(config.getParameter<ParameterSet>("TauSelection")), // Fixes "An object with name tauSelection_ exists already"
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


void TauFakeRate_ee::book(TDirectory *dir) {

  if (0) std::cout << "=== TauFakeRate_ee::book()" << std::endl;
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
  hNTau_AfterTauSelection   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_AfterTauSelection"  , ";loose #tau-jet multiplicity" , nN  , nMin  , nMax  );
  hNTau_AfterJetSelection   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_AfterJetSelection"  , ";loose #tau-jet multiplicity" , nN  , nMin  , nMax  );
  hNTau_AfterBJetSelection  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_AfterBJetSelection" , ";loose #tau-jet multiplicity" , nN  , nMin  , nMax  );
  hNTau_AfterMetSelection   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_AfterMetSelection"  , ";loose #tau-jet multiplicity" , nN  , nMin  , nMax  );

  hTauSrc_AfterAllSelections       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrc_AfterAllSelections"      , ";#tau source", nN, nMin, nMax  );
  hNTau_AfterAllSelections         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_AfterAllSelections"        , ";loose #tau-jet multiplicity", nN  , nMin  , nMax  );
  hNJet_AfterAllSelections         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NJet_AfterAllSelections"        , ";jet multiplicity"  , nN  , nMin  , nMax  );
  hNBjet_AfterAllSelections        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NBjet_AfterAllSelections"       , ";b-jet multiplicity", nN  , nMin  , nMax  );
  hMet_AfterAllSelections          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "Met_AfterAllSelections"         , ";E_{T}^{miss} (GeV)", metN, metMin, metMax);
  hDileptonPt_AfterAllSelections   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonPt_AfterAllSelections"  , ";p_{T} (GeV)" , ptN , ptMin , ptMax );
  hDileptonEta_AfterAllSelections  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonEta_AfterAllSelections" , ";#eta"        , etaN, etaMin, etaMax);
  hDileptonMass_AfterAllSelections = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_AfterAllSelections", ";m_{ll} (GeV)", mN  , mMin  , mMax  );

  return;
}


void TauFakeRate_ee::setupBranches(BranchManager& branchManager) {
  fEvent.setupBranches(branchManager);
  return;
}


void TauFakeRate_ee::process(Long64_t entry) {

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

  if (0) cout << "passed trigger" << endl;
  // Fill histos
  fCommonPlots.fillControlPlotsAfterTrigger(fEvent);


  //================================================================================================   
  // MET filters (to remove events with spurious sources of fake MET)
  //================================================================================================   
  if (0) std::cout << "=== MET Filter" << std::endl;
  const METFilterSelection::Data metFilterData = fMETFilterSelection.analyze(fEvent);
  if (!metFilterData.passedSelection()) return;

  if (0) cout << "passed MET" << endl;
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
  // Muon Selection
  //================================================================================================   
  if (0) std::cout << "=== Muon veto" << std::endl;
  const MuonSelection::Data muData = fMuonSelection.analyze(fEvent);
  if (muData.hasIdentifiedMuons()) return;

  // Fill histos
  fCommonPlots.fillControlPlotsAtMuonSelection(fEvent, muData);


  //================================================================================================
  // Electron Selection
  //================================================================================================
  if (0) std::cout << "=== Electron Selection" << std::endl;
  const ElectronSelection::Data eData = fElectronSelection.analyze(fEvent);
  if(!(eData.hasIdentifiedElectrons())) return;

  // Require exactly 2 electron
  if (eData.getSelectedElectrons().size() != 2) return; // note: remember to disable trigger-matching option if using a single muon trigger

  if (0) cout << "passed two electrons" << endl;
  // Calculate variables for dilepton system
  double dilepton_invMass = ( eData.getSelectedElectrons()[0].p4() + eData.getSelectedElectrons()[1].p4() ).M();
  double dilepton_pt      = ( eData.getSelectedElectrons()[0].p4() + eData.getSelectedElectrons()[1].p4() ).pt();
  double dilepton_eta     = ( eData.getSelectedElectrons()[0].p4() + eData.getSelectedElectrons()[1].p4() ).eta();
  const double cfg_massWindow = 20.0; // 15.0;

  if (0) cout << "dilepton_invMass: " << dilepton_invMass << endl;
  // Apply muon trigger scale factors
  // https://gitlab.cern.ch/HPlus/HiggsAnalysis/blob/master/NtupleAnalysis/python/parameters/scaleFactors.py
  // https://gitlab.cern.ch/HPlus/HiggsAnalysis/blob/master/NtupleAnalysis/data/TriggerEfficiency/muonPAGEff.json
  if (fEvent.isMC()) 
    {
      //fEventWeight.multiplyWeight(muData.getMuonTriggerSF());
      fEventWeight.multiplyWeight(eData.getElectronTriggerSF());
    }

  // Apply muon ID scale factors
  if (fEvent.isMC()) 
    {
      //fEventWeight.multiplyWeight(muData.getMuonIDSF());
      fEventWeight.multiplyWeight(eData.getElectronIDSF());
    }
  
  if (0) cout << "electron id SF applied" << endl;

  //cout << "eData.getSelectedElectrons()[0].getHLTElectronCharge(): " <<  eData.getSelectedElectrons()[0].getHLTElectronCharge() << endl; //fix me the charge collection

  // Require the muons to have opposite sign (OS)
  //if( eData.getSelectedElectrons()[0].charge() == eData.getSelectedElectrons()[1].charge()) return;
  //cLeptonOSCounter.increment(); 
  if (0) std::cout << "=== Electron Selection: OS requirement" << std::endl;

  hDileptonMass_BeforeOnZSelection->Fill(dilepton_invMass);

  // Apply on-Z mass requirement for dilepton pair
  if (dilepton_invMass < ( 91.1876 - cfg_massWindow ) ) return; // fixme. define in run.py 
  if (dilepton_invMass > ( 91.1876 + cfg_massWindow ) ) return; // fixme. define in run.py
  cLeptonMassCounter.increment(); 

  // Apply dR cut for dilepton system?
  //  dilepton_dR = ROOT::Math::VectorUtil::DeltaR(eData.getSelectedElectrons()[0].p4(),eData.getSelectedElectrons()[1].p4());
  //  if(dilepton_dR < 0.3) return;

  // Fill histos
  fCommonPlots.fillControlPlotsAtElectronSelection(fEvent, eData);
  hDileptonPt_AfterLeptonSelection->Fill(dilepton_pt);
  hDileptonEta_AfterLeptonSelection->Fill(dilepton_eta);
  hDileptonMass_AfterLeptonSelection->Fill(dilepton_invMass);
  hLeptonN_AfterLeptonSelection ->Fill(eData.getSelectedElectrons().size());
  // For-loop: All Selected leptons
  for(unsigned int i=0; i< eData.getSelectedElectrons().size(); i++)
    {
      hLeptonPt_AfterLeptonSelection ->Fill(eData.getSelectedElectrons()[i].pt());
      hLeptonEta_AfterLeptonSelection->Fill(eData.getSelectedElectrons()[i].eta());
    }

  

  //================================================================================================   
  // Tau Selection
  //================================================================================================   
  if (0) std::cout << "=== Tau Selection" << std::endl;
  const TauSelection::Data tauData = fTauSelection.analyze(fEvent);
  const TauSelection::Data looseTauData = fLooseTauSelection.analyzeLoose(fEvent);

  if (!looseTauData.hasIdentifiedTaus() ) return;
  //if(looseTauData.getSelectedTaus().size() != 1) return; // not the same as line above!

  if (0) std::cout << "=== Tau Selection: Has at least 1 loose tau(s)" << std::endl;
  cTauNCounter.increment(); 

  // Determine if genuine or fake tau. If more than one tau present use the first one in the list (pt-sorted)
  bool bEleToTau     = false; // ele  -> tau fakes only
  bool bMuonToTau    = false; // muon -> tau fakes only
  bool bJetToTau     = false; // jet  -> tau fakes only
  bool bLightQToTau  = false; // light quark-> tau fakes only (u,d,s)
  bool bHeavyQToTau  = false; // heavy quark-> tau fakes only (c,b)
  bool bGluonToTau   = false; // jet  -> tau fakes only
  bool bUnknownToTau = false; // unknown -> tau fakes only
  bool bGenuineTau   = false;
  if ( fEvent.isMC() )
    {
      bEleToTau     = looseTauData.getSelectedTau().isElectronToTau();
      bMuonToTau    = looseTauData.getSelectedTau().isMuonToTau();
      bJetToTau     = looseTauData.getSelectedTau().isJetToTau();
      bLightQToTau  = looseTauData.getSelectedTau().isQuarkToTau(1) || looseTauData.getSelectedTau().isQuarkToTau(2) || looseTauData.getSelectedTau().isQuarkToTau(3);
      bHeavyQToTau  = looseTauData.getSelectedTau().isQuarkToTau(4) || looseTauData.getSelectedTau().isQuarkToTau(5);
      bGluonToTau   = looseTauData.getSelectedTau().isGluonToTau();
      bUnknownToTau = looseTauData.getSelectedTau().isUnknownTauDecay();
      bGenuineTau   = looseTauData.getSelectedTau().isGenuineTau();
    }
  int tauSrcBit =  1*bEleToTau + 2*bMuonToTau + 3*bJetToTau + 4*bLightQToTau + 5*bHeavyQToTau + 6*bGluonToTau + 7*bUnknownToTau;// => genuineTau = 0
  // int tauSrcBit =  pow(2,0)*bEleToTau + pow(2,1)*bMuonToTau + pow(2,2)*bJetToTau + pow(2,3)*bQuarkToTau + pow(2,4)*bGluonToTau + pow(2,5)*bUnknownToTau;// => genuineTau = 0

  // Fake tau in our analysis is "not genuine tau and not e->tau and not mu->tau"
  bool bFakeTau = ( fEvent.isMC() && !bGenuineTau && !(bEleToTau || bMuonToTau) );
  
  // Fill histos ( also sets value for boolean bIsGenuineTau
  hDileptonMass_AfterTauSelection->Fill(dilepton_invMass);
  hNTau_AfterTauSelection->Fill( looseTauData.getSelectedTaus().size() );
  if (fEvent.isMC()) hTauSrc_AfterTauSelection->Fill(tauSrcBit);
  fCommonPlots.fillControlPlotsAfterTauSelection(fEvent, looseTauData); // uses: bIsGenuineTau = data.getSelectedTaus()[0].isGenuineTau();

  // Redefine what the "bGenuineTau" boolean is. N.B: Affects the genuineTau plots filled by control plots!
  if (0) std::cout << "=== Tau Selection: bFakeTau = " << bFakeTau << " (src = " << tauSrcBit << "). NTaus = " << looseTauData.getSelectedTaus().size() << std::endl;
  fCommonPlots.setGenuineTauStatus(!bFakeTau); // CommonPlots.cc (L1292) bIsGenuineTau = data.getSelectedTaus()[0].isGenuineTau();


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
  if (0) std::cout << "=== 1" << std::endl;
  //  fCommonPlots.fillControlPlotsAfterAllSelections(fEvent);
  if (0) std::cout << "=== 2" << std::endl;

  if (fEvent.isMC()) hTauSrc_AfterAllSelections->Fill(tauSrcBit);
  if (0) std::cout << "=== 3" << std::endl;
  hNTau_AfterAllSelections->Fill( looseTauData.getSelectedTaus().size() );
  if (0) std::cout << "=== 4" << std::endl;
  hNJet_AfterAllSelections ->Fill( jetData.getSelectedJets().size() );
  if (0) std::cout << "=== 5" << std::endl;
  hNBjet_AfterAllSelections->Fill( bjetData.getSelectedBJets().size() );
  if (0) std::cout << "=== 6" << std::endl;
  hDileptonPt_AfterAllSelections->Fill(dilepton_pt);
  if (0) std::cout << "=== 7" << std::endl;
  hDileptonEta_AfterAllSelections->Fill(dilepton_eta);
  if (0) std::cout << "=== 8" << std::endl;
  hDileptonMass_AfterAllSelections->Fill(dilepton_invMass);
  if (0) std::cout << "=== 9" << std::endl;
  hMet_AfterAllSelections->Fill(metData.getMET().R());
  if (0) std::cout << "=== 10" << std::endl;

  // Tau stuff here
  if(looseTauData.getSelectedTaus().size() != 1) return; // why is this necessary? (remove it and get crash)

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
      doTightTaus(fEvent, tauData);
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
      doLooseTaus(fEvent, looseTauData);
    }

  fEventSaver.save();

  return;
}


void TauFakeRate_ee::doLooseTaus(const Event& event, const TauSelection::Data& tauData) {


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
      bool bIsMC          = fEvent.isMC();
      bool bGenuineTau    = false;
      bool bEleToTau      = false;
      bool bMuonToTau     = false;
      bool bFakeTau       = false;
      if (bIsMC)
	{
	  bGenuineTau = tauData.getSelectedTaus()[i].isGenuineTau(); // tauData.isGenuineTau();
	  bEleToTau   = tauData.getSelectedTaus()[i].isElectronToTau();
	  bMuonToTau  = tauData.getSelectedTaus()[i].isMuonToTau();
	  bFakeTau    = (bIsMC && !bGenuineTau && !(bEleToTau || bMuonToTau) );
	}

      // Only if MC  and selected tau is genuine (not fake)
      // if (event.isMC() && tauData.getSelectedTaus()[i].isGenuineTau()) 
      if (bIsMC && !bFakeTau)
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
	}// (bIsMC && !bFakeTau)
    
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


void TauFakeRate_ee::doTightTaus(const Event& event, const TauSelection::Data& tightTauData) {

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
      bool bIsMC          = fEvent.isMC();
      bool bGenuineTau    = false;
      bool bEleToTau      = false;
      bool bMuonToTau     = false;
      bool bFakeTau       = false;
      if (bIsMC)
	{
	  bGenuineTau = tightTauData.getSelectedTaus()[i].isGenuineTau();
	  bEleToTau   = tightTauData.getSelectedTaus()[i].isElectronToTau();
	  bMuonToTau  = tightTauData.getSelectedTaus()[i].isMuonToTau();
	  bFakeTau    = (bIsMC && !bGenuineTau && !(bEleToTau || bMuonToTau) );
	}

      // Only if MC  and selected tau is genuine. (if not genuine tau, reject the events)
      //if (event.isMC() && tightTauData.getSelectedTaus()[i].isGenuineTau()) 
      if (bIsMC && !bFakeTau)
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
