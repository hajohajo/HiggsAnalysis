// -*- c++ -*-
#include "Framework/interface/BaseSelector.h"
#include "Framework/interface/makeTH.h"

#include "EventSelection/interface/CommonPlots_lt.h"
#include "EventSelection/interface/EventSelections.h"
#include "EventSelection/interface/TransverseMass.h"


#include "TDirectory.h"

class TauFakeRate_mm: public BaseSelector {
public:
  explicit TauFakeRate_mm(const ParameterSet& config, const TH1* skimCounters);
  virtual ~TauFakeRate_mm() {}

  /// Books histograms
  virtual void book(TDirectory *dir) override;
  /// Sets up branches for reading the TTree
  virtual void setupBranches(BranchManager& branchManager) override;
  /// Called for each event
  virtual void process(Long64_t entry) override;

private:
  // Input parameters
  const DirectionalCut<int> cfg_DileptonAbsCharge;
  const double cfg_DileptonInvariantMass;
  const double cfg_DileptonInvariantMassPlus;
  const double cfg_DileptonInvariantMassMinus;
  const DirectionalCut<int> cfg_DileptonDeltaR;

  // Common plots
  CommonPlots_lt fCommonPlots;

  // Event selection classes and event counters (in same order like they are applied)
  Count cAllEvents;
  Count cTrigger;
  METFilterSelection fMETFilterSelection;
  Count cVertexSelection;
  ElectronSelection fElectronSelection;
  MuonSelection fMuonSelection;
  Count cLeptonNCounter;
  Count cLeptonOSCounter;
  Count cLeptonMassCounter;
  TauSelection fTauSelection;
  TauSelection fLooseTauSelection;
  // Count cTauNCounter;
  JetSelection fJetSelection;
  BJetSelection fBJetSelection;
  Count cBTaggingSFCounter;
  METSelection fMETSelection;
  Count cTauSFCounter;
  Count cFakeTauSFCounter;
  Count cPreselected;
  Count cSelected;
    
  int getTauSrcBit(const Event& event, const TauSelection::Data& tauData);
  bool getIsFakeTau(const Event& event, const TauSelection::Data& tauData);
  void doLooseTaus(const Event& event, const TauSelection::Data& tauData, const bool bOnZMass);
  void doTightTaus(const Event& event, const TauSelection::Data& looseTauData, const bool bOnZMass);

  // Non-common histograms
  WrappedTH1 *hTauPt_num;
  WrappedTH1 *hTauPt_num_dm0;
  WrappedTH1 *hTauPt_num_dm1;
  WrappedTH1 *hTauPt_num_dm10;
  WrappedTH1 *hTauPt_den;
  WrappedTH1 *hTauPt_den_dm0;
  WrappedTH1 *hTauPt_den_dm1;
  WrappedTH1 *hTauPt_den_dm10;
  WrappedTH1 *hTauPt_num_g;
  WrappedTH1 *hTauPt_num_g_dm0;
  WrappedTH1 *hTauPt_num_g_dm1;
  WrappedTH1 *hTauPt_num_g_dm10;
  WrappedTH1 *hTauPt_den_g;
  WrappedTH1 *hTauPt_den_g_dm0;
  WrappedTH1 *hTauPt_den_g_dm1;
  WrappedTH1 *hTauPt_den_g_dm10;
  //
  WrappedTH1 *hTauEta_num;
  WrappedTH1 *hTauEta_num_dm0;
  WrappedTH1 *hTauEta_num_dm1;
  WrappedTH1 *hTauEta_num_dm10;
  WrappedTH1 *hTauEta_den;
  WrappedTH1 *hTauEta_den_dm0;
  WrappedTH1 *hTauEta_den_dm1;
  WrappedTH1 *hTauEta_den_dm10;
  WrappedTH1 *hTauEta_num_g;
  WrappedTH1 *hTauEta_num_g_dm0;
  WrappedTH1 *hTauEta_num_g_dm1;
  WrappedTH1 *hTauEta_num_g_dm10;
  WrappedTH1 *hTauEta_den_g;
  WrappedTH1 *hTauEta_den_g_dm0;
  WrappedTH1 *hTauEta_den_g_dm1;
  WrappedTH1 *hTauEta_den_g_dm10;

  // Barrel ( |eta| < 1.5)
  WrappedTH1 *hTauPt_num_barrel;
  WrappedTH1 *hTauPt_num_dm0_barrel;
  WrappedTH1 *hTauPt_num_dm1_barrel;
  WrappedTH1 *hTauPt_num_dm10_barrel;
  WrappedTH1 *hTauPt_den_barrel;
  WrappedTH1 *hTauPt_den_dm0_barrel;
  WrappedTH1 *hTauPt_den_dm1_barrel;
  WrappedTH1 *hTauPt_den_dm10_barrel;
  WrappedTH1 *hTauPt_num_g_barrel;
  WrappedTH1 *hTauPt_num_g_dm0_barrel;
  WrappedTH1 *hTauPt_num_g_dm1_barrel;
  WrappedTH1 *hTauPt_num_g_dm10_barrel;
  WrappedTH1 *hTauPt_den_g_barrel;
  WrappedTH1 *hTauPt_den_g_dm0_barrel;
  WrappedTH1 *hTauPt_den_g_dm1_barrel;
  WrappedTH1 *hTauPt_den_g_dm10_barrel;

  // Endcap ( |eta| >= 1.5)
  WrappedTH1 *hTauPt_num_endcap;
  WrappedTH1 *hTauPt_num_dm0_endcap;
  WrappedTH1 *hTauPt_num_dm1_endcap;
  WrappedTH1 *hTauPt_num_dm10_endcap;
  WrappedTH1 *hTauPt_den_endcap;
  WrappedTH1 *hTauPt_den_dm0_endcap;
  WrappedTH1 *hTauPt_den_dm1_endcap;
  WrappedTH1 *hTauPt_den_dm10_endcap;
  WrappedTH1 *hTauPt_num_g_endcap;
  WrappedTH1 *hTauPt_num_g_dm0_endcap;
  WrappedTH1 *hTauPt_num_g_dm1_endcap;
  WrappedTH1 *hTauPt_num_g_dm10_endcap;
  WrappedTH1 *hTauPt_den_g_endcap;
  WrappedTH1 *hTauPt_den_g_dm0_endcap;
  WrappedTH1 *hTauPt_den_g_dm1_endcap;
  WrappedTH1 *hTauPt_den_g_dm10_endcap;

  WrappedTH1 *hLeptonN_AfterLeptonSelection;
  WrappedTH1 *hLeptonPt_AfterLeptonSelection;
  WrappedTH1 *hLeptonEta_AfterLeptonSelection;
  WrappedTH1 *hDileptonPt_AfterLeptonSelection;
  WrappedTH1 *hDileptonEta_AfterLeptonSelection;
  WrappedTH1 *hDileptonMass_AfterLeptonSelection;
  WrappedTH1 *hDileptonCharge_AfterLeptonSelection;
  WrappedTH1 *hDileptonDeltaR_AfterLeptonSelection;
  
  WrappedTH1 *hDileptonMass_AfterTauSelection;

  WrappedTH1 *hNTau_Preselections;
  WrappedTH1 *hTauSrc_Preselections;
  WrappedTH1 *hTauSrcDM0_Preselections;
  WrappedTH1 *hTauSrcDM1_Preselections;
  WrappedTH1 *hTauSrcDM10_Preselections;
  WrappedTH1 *hNJet_Preselections;
  WrappedTH1 *hNBjet_Preselections;
  WrappedTH1 *hMet_Preselections;
  WrappedTH1 *hLt_Preselections;
  WrappedTH1 *hDileptonPt_Preselections;
  WrappedTH1 *hDileptonEta_Preselections;
  WrappedTH1 *hDileptonMass_Preselections;
  WrappedTH1 *hDileptonMassOnZ_Preselections;
  WrappedTH1 *hDileptonMassOffZ_Preselections;
  WrappedTH1 *hDileptonCharge_Preselections;
  WrappedTH1 *hDileptonDeltaR_Preselections;

  WrappedTH1 *hNTau_LooseTau;
  WrappedTH1 *hTauSrc_LooseTau;
  WrappedTH1 *hTauSrcDM0_LooseTau;
  WrappedTH1 *hTauSrcDM1_LooseTau;
  WrappedTH1 *hTauSrcDM10_LooseTau;
  WrappedTH1 *hNJet_LooseTau;
  WrappedTH1 *hNBjet_LooseTau;
  WrappedTH1 *hMet_LooseTau;
  WrappedTH1 *hLt_LooseTau;
  WrappedTH1 *hDileptonPt_LooseTau;
  WrappedTH1 *hDileptonEta_LooseTau;
  WrappedTH1 *hDileptonMass_LooseTau;

  WrappedTH1 *hNTau_TightTau;
  WrappedTH1 *hTauSrc_TightTau;
  WrappedTH1 *hTauSrcDM0_TightTau;
  WrappedTH1 *hTauSrcDM1_TightTau;
  WrappedTH1 *hTauSrcDM10_TightTau;
  WrappedTH1 *hNJet_TightTau;
  WrappedTH1 *hNBjet_TightTau;
  WrappedTH1 *hMet_TightTau;
  WrappedTH1 *hLt_TightTau;
  WrappedTH1 *hDileptonPt_TightTau;
  WrappedTH1 *hDileptonEta_TightTau;
  WrappedTH1 *hDileptonMass_TightTau;

};

#include "Framework/interface/SelectorFactory.h"
REGISTER_SELECTOR(TauFakeRate_mm);

TauFakeRate_mm::TauFakeRate_mm(const ParameterSet& config, const TH1* skimCounters)
  : BaseSelector(config, skimCounters),
    cfg_DileptonAbsCharge(config, "TauFakeRateMeasurement.dileptonAbsChargeCut"),
    cfg_DileptonInvariantMass(config.getParameter<double>("TauFakeRateMeasurement.dileptonInvariantMass")),
    cfg_DileptonInvariantMassPlus(config.getParameter<double>("TauFakeRateMeasurement.dileptonInvariantMassPlus")),
    cfg_DileptonInvariantMassMinus(config.getParameter<double>("TauFakeRateMeasurement.dileptonInvariantMassMinus")),      
    cfg_DileptonDeltaR(config, "TauFakeRateMeasurement.dileptonDeltaRCut"),
    // fCommonPlots(config.getParameter<ParameterSet>("CommonPlots"), CommonPlots::kHplus2hwAnalysisWithTop, fHistoWrapper),
    fCommonPlots(config.getParameter<ParameterSet>("CommonPlots"), CommonPlots_lt::kTauFakeRateMeasurement, fHistoWrapper),
    cAllEvents(fEventCounter.addCounter("All")),
    cTrigger(fEventCounter.addCounter("Trg")),
    fMETFilterSelection(config.getParameter<ParameterSet>("METFilter"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cVertexSelection(fEventCounter.addCounter("PV")),
    fElectronSelection(config.getParameter<ParameterSet>("ElectronSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, " Veto"),
    fMuonSelection(config.getParameter<ParameterSet>("MuonSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cLeptonNCounter(fEventCounter.addCounter("#mu#mu")),
    cLeptonOSCounter(fEventCounter.addCounter("OS #mu")),
    cLeptonMassCounter(fEventCounter.addCounter("m_{#mu#mu}")),
    // /fTauSelection(config.getParameter<ParameterSet>("TauSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fTauSelection(config.getParameter<ParameterSet>("TauSelection")),
    fLooseTauSelection(config.getParameter<ParameterSet>("LooseTauSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    //fLooseTauSelection(config.getParameter<ParameterSet>("LooseTauSelection")),
    // cTauNCounter(fEventCounter.addCounter("#geq 1 loose #tau_{h}")),
    fJetSelection(config.getParameter<ParameterSet>("JetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fBJetSelection(config.getParameter<ParameterSet>("BJetSelection")), //fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cBTaggingSFCounter(fEventCounter.addCounter("b-tag SF")),
    fMETSelection(config.getParameter<ParameterSet>("METSelection")),
    cTauSFCounter(fEventCounter.addCounter("#tau_{h} SF")),
    cFakeTauSFCounter(fEventCounter.addCounter("Fake #tau_{h} SF")),
    cPreselected(fEventCounter.addCounter("Preselections")),
    cSelected(fEventCounter.addCounter("Selections"))
{ }


void TauFakeRate_mm::book(TDirectory *dir) {

  if (0) std::cout << "=== TauFakeRate_mm::book()" << std::endl;
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
  const int   htN    = fCommonPlots.getHtBinSettings().bins();
  const float htMin  = fCommonPlots.getHtBinSettings().min();
  const float htMax  = fCommonPlots.getHtBinSettings().max();
  const int dRN      = fCommonPlots.getDeltaRBinSettings().bins();
  const double dRMin = fCommonPlots.getDeltaRBinSettings().min();
  const double dRMax = fCommonPlots.getDeltaRBinSettings().max();
  
  // Book non-common histograms 
  double bin[8] = {20,25,30,35,40,50,60,120};

  hTauPt_num        =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num"       , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_dm0    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm0"   , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_dm1    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm1"   , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_dm10   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm10"  , ";p_{T} (GeV)", 7, bin);
  hTauPt_den        =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den"       , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_dm0    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm0"   , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_dm1    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm1"   , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_dm10   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm10"  , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g      =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g"     , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm0  =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm0" , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm1  =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm1" , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm10 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm10", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g      =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g"     , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm0  =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm0" , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm1  =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm1" , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm10 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm10", ";p_{T} (GeV)", 7, bin);
  //
  hTauEta_num        =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num"       ,  ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_dm0    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm0"   ,  ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_dm1    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm1"   ,  ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_dm10   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_dm10"  ,  ";#eta", etaN, etaMin, etaMax);
  hTauEta_den        =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den"       ,  ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm0    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm0"   ,  ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm1    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm1"   ,  ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_dm10   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_dm10"  ,  ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_g      =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_g"     ,  ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_g_dm0  =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_g_dm0" ,  ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_g_dm1  =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_g_dm1" ,  ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_g_dm10 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_g_dm10",  ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_g      =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_g"     ,  ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_g_dm0  =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_g_dm0" ,  ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_g_dm1  =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_g_dm1" ,  ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_g_dm10 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_g_dm10",  ";#eta", etaN, etaMin, etaMax);

  // Barrel ( |eta| < 1.5)
  hTauPt_num_barrel        =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_barrel"       , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_dm0_barrel    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm0_barrel"   , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_dm1_barrel    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm1_barrel"   , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_dm10_barrel   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm10_barrel"  , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_barrel        =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_barrel"       , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_dm0_barrel    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm0_barrel"   , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_dm1_barrel    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm1_barrel"   , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_dm10_barrel   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm10_barrel"  , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_barrel      =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_barrel"     , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm0_barrel  =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm0_barrel" , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm1_barrel  =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm1_barrel" , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm10_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm10_barrel", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_barrel      =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_barrel"     , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm0_barrel  =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm0_barrel" , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm1_barrel  =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm1_barrel" , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm10_barrel =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm10_barrel", ";p_{T} (GeV)", 7, bin);

  // Endcap ( |eta| >= 1.5)
  hTauPt_num_endcap        =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_endcap"       , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_dm0_endcap    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm0_endcap"   , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_dm1_endcap    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm1_endcap"   , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_dm10_endcap   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_dm10_endcap"  , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_endcap        =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_endcap"       , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_dm0_endcap    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm0_endcap"   , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_dm1_endcap    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm1_endcap"   , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_dm10_endcap   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_dm10_endcap"  , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_endcap      =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_endcap"     , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm0_endcap  =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm0_endcap" , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm1_endcap  =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm1_endcap" , ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_dm10_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_dm10_endcap", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_endcap      =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_endcap"     , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm0_endcap  =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm0_endcap" , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm1_endcap  =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm1_endcap" , ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_dm10_endcap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_dm10_endcap", ";p_{T} (GeV)", 7, bin);

  hLeptonN_AfterLeptonSelection        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LeptonN_AfterLeptonSelection"       , ";lepton multiplicity", nN  , nMin  , nMax );
  hLeptonPt_AfterLeptonSelection       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LeptonPt_AfterLeptonSelection"      , ";p_{T} (GeV)", ptN , ptMin , ptMax );
  hLeptonEta_AfterLeptonSelection      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "LeptonEta_AfterLeptonSelection"     , ";#eta"       , etaN, etaMin, etaMax);
  hDileptonPt_AfterLeptonSelection     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonPt_AfterLeptonSelection"    , ";p_{T} (GeV)", ptN , ptMin , ptMax );
  hDileptonEta_AfterLeptonSelection    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonEta_AfterLeptonSelection"   , ";#eta"       , etaN, etaMin, etaMax);
  hDileptonMass_AfterLeptonSelection   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_AfterLeptonSelection"  , ";m_{ll} (GeV)", mN , mMin  , mMax  );
  hDileptonCharge_AfterLeptonSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonCharge_AfterLeptonSelection", ";Q_{ll} (GeV)", 6, -3, +3);
  hDileptonDeltaR_AfterLeptonSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonDeltaR_AfterLeptonSelection", ";#DeltaR_{ll}", dRN, dRMin , dRMax );

  hDileptonMass_AfterTauSelection    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_AfterTauSelection"  , ";m_{ll} (GeV)", mN, mMin, mMax);

  hTauSrc_Preselections           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrc_Preselections"      , ";#tau source", 80, 0.0, 80.0);
  hTauSrcDM0_Preselections        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrcDM0_Preselections"   , ";#tau source", 80, 0.0, 80.0);
  hTauSrcDM1_Preselections        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrcDM1_Preselections"   , ";#tau source", 80, 0.0, 80.0);
  hTauSrcDM10_Preselections       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrcDM10_Preselections"  , ";#tau source", 80, 0.0, 80.0);
  hNTau_Preselections             = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_Preselections"        , ";loose #tau-jet multiplicity", nN  , nMin  , nMax  );
  hNJet_Preselections             = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NJet_Preselections"        , ";jet multiplicity"  , nN  , nMin  , nMax  );
  hNBjet_Preselections            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NBjet_Preselections"       , ";b-jet multiplicity", nN  , nMin  , nMax  );
  hMet_Preselections              = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "Met_Preselections"         , ";E_{T}^{miss} (GeV)", metN, metMin, metMax);
  hLt_Preselections               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "Lt_Preselections"          , ";L_{T} (GeV)" , htN, htMin, htMax);
  hDileptonPt_Preselections       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonPt_Preselections"  , ";p_{T} (GeV)" , ptN , ptMin , ptMax );
  hDileptonEta_Preselections      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonEta_Preselections" , ";#eta"        , etaN, etaMin, etaMax);
  hDileptonMass_Preselections     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_Preselections", ";m_{ll} (GeV)", mN  , mMin  , mMax  );
  hDileptonMassOnZ_Preselections  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMassOnZ_Preselections" , ";m_{ll} (GeV)", mN  , mMin  , mMax);
  hDileptonMassOffZ_Preselections = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMassOffZ_Preselections", ";m_{ll} (GeV)", mN  , mMin  , mMax);
  hDileptonCharge_Preselections   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonCharge_Preselections"  , ";Q_{ll} (GeV)", 6, -3, +3);
  hDileptonDeltaR_Preselections   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonDeltaR_Preselections"  , ";#DeltaR_{ll}", dRN , dRMin , dRMax );

  hTauSrc_LooseTau         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrc_LooseTau"        , ";#tau source", 80, 0.0, 80.0);
  hTauSrcDM0_LooseTau      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrcDM0_LooseTau"     , ";#tau source", 80, 0.0, 80.0);
  hTauSrcDM1_LooseTau      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrcDM1_LooseTau"     , ";#tau source", 80, 0.0, 80.0);
  hTauSrcDM10_LooseTau     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrcDM10_LooseTau"    , ";#tau source", 80, 0.0, 80.0);
  hNTau_LooseTau           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_LooseTau"          , ";loose #tau-jet multiplicity", nN  , nMin  , nMax  );
  hNJet_LooseTau           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NJet_LooseTau"          , ";jet multiplicity"  , nN  , nMin  , nMax  );
  hNBjet_LooseTau          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NBjet_LooseTau"         , ";b-jet multiplicity", nN  , nMin  , nMax  );
  hMet_LooseTau            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "Met_LooseTau"           , ";E_{T}^{miss} (GeV)", metN, metMin, metMax);
  hLt_LooseTau             = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "Lt_LooseTau"            , ";L_{T} (GeV)" , htN, htMin, htMax);
  hDileptonPt_LooseTau     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonPt_LooseTau"    , ";p_{T} (GeV)" , ptN , ptMin , ptMax );
  hDileptonEta_LooseTau    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonEta_LooseTau"   , ";#eta"        , etaN, etaMin, etaMax);
  hDileptonMass_LooseTau   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_LooseTau"  , ";m_{ll} (GeV)", mN  , mMin  , mMax  );

  hTauSrc_TightTau         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrc_TightTau"        , ";#tau source", 80, 0.0, 80.0);
  hTauSrcDM0_TightTau      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrcDM0_TightTau"     , ";#tau source", 80, 0.0, 80.0);
  hTauSrcDM1_TightTau      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrcDM1_TightTau"     , ";#tau source", 80, 0.0, 80.0);
  hTauSrcDM10_TightTau     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauSrcDM10_TightTau"    , ";#tau source", 80, 0.0, 80.0);
  hNTau_TightTau           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NTau_TightTau"          , ";loose #tau-jet multiplicity", nN  , nMin  , nMax  );
  hNJet_TightTau           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NJet_TightTau"          , ";jet multiplicity"  , nN  , nMin  , nMax  );
  hNBjet_TightTau          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "NBjet_TightTau"         , ";b-jet multiplicity", nN  , nMin  , nMax  );
  hMet_TightTau            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "Met_TightTau"           , ";E_{T}^{miss} (GeV)", metN, metMin, metMax);
  hLt_TightTau             = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "Lt_TightTau"            , ";L_{T} (GeV)" , htN, htMin, htMax);
  hDileptonPt_TightTau     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonPt_TightTau"    , ";p_{T} (GeV)" , ptN , ptMin , ptMax );
  hDileptonEta_TightTau    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonEta_TightTau"   , ";#eta"        , etaN, etaMin, etaMax);
  hDileptonMass_TightTau   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_TightTau"  , ";m_{ll} (GeV)", mN  , mMin  , mMax  );

  return;
}


void TauFakeRate_mm::setupBranches(BranchManager& branchManager) {
  fEvent.setupBranches(branchManager);
  return;
}


void TauFakeRate_mm::process(Long64_t entry) {

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
  fCommonPlots.fillControlPlotsAfterTrigger(fEvent);


  //================================================================================================   
  // MET filters (to remove events with spurious sources of fake MET)
  //================================================================================================   
  if (0) std::cout << "=== MET Filter" << std::endl;
  const METFilterSelection::Data metFilterData = fMETFilterSelection.analyze(fEvent);
  if (!metFilterData.passedSelection()) return;
  fCommonPlots.fillControlPlotsAfterMETFilter(fEvent);  


  //================================================================================================   
  // Primarty Vertex (Check that a PV exists)
  //================================================================================================   
  if (0) std::cout << "=== Vertices" << std::endl;
  if (nVertices < 1) return;
  cVertexSelection.increment();
  fCommonPlots.fillControlPlotsAtVertexSelection(fEvent);

  
  //================================================================================================   
  // Electron Selection
  //================================================================================================   
  if (0) std::cout << "=== Electron veto" << std::endl;
  const ElectronSelection::Data eData = fElectronSelection.analyze(fEvent);
  if (eData.hasIdentifiedElectrons()) return;

  //================================================================================================
  // Muon Selection
  //================================================================================================
  if (0) std::cout << "=== Muon Selection" << std::endl;
  const MuonSelection::Data muData = fMuonSelection.analyze(fEvent);
  if(!(muData.hasIdentifiedMuons())) return;

  // Require exactly 2 muons
  if (muData.getSelectedMuons().size() != 2) return; // note: remember to disable trigger-matching option if using a single muon trigger
  cLeptonNCounter.increment(); 

  // Calculate variables for dilepton system
  double dilepton_pt      = ( muData.getSelectedMuons()[0].p4() + muData.getSelectedMuons()[1].p4() ).pt();
  double dilepton_eta     = ( muData.getSelectedMuons()[0].p4() + muData.getSelectedMuons()[1].p4() ).eta();
  double dilepton_invMass = ( muData.getSelectedMuons()[0].p4() + muData.getSelectedMuons()[1].p4() ).M();
  double dilepton_dR      = ROOT::Math::VectorUtil::DeltaR(muData.getSelectedMuons()[0].p4(), muData.getSelectedMuons()[1].p4());
  int dilepton_charge     = muData.getSelectedMuons()[0].charge() + muData.getSelectedMuons()[1].charge();
  bool bPassAbsCharge     = false;
  bool bOnZMass           = false; 
  bool bOnZMassAbove      = true; // yes, defualt is "true"
  bool bOnZMassBelow      = true; // yes, defualt is "true"
  bool bPassDeltaR        = false;

  // Apply muon trigger and ID scale factor (SF)
  // https://gitlab.cern.ch/HPlus/HiggsAnalysis/blob/master/NtupleAnalysis/python/parameters/scaleFactors.py
  // https://gitlab.cern.ch/HPlus/HiggsAnalysis/blob/master/NtupleAnalysis/data/TriggerEfficiency/muonPAGEff.json
  if (fEvent.isMC()) 
    {
      fEventWeight.multiplyWeight(muData.getMuonTriggerSF());
      fEventWeight.multiplyWeight(muData.getMuonIDSF());
    }
  
  // Require the muons to have opposite sign (OS)
  bPassAbsCharge = cfg_DileptonAbsCharge.passedCut( abs(dilepton_charge) );
  if (!bPassAbsCharge) return;
  cLeptonOSCounter.increment(); 
  if (0) std::cout << "=== Muon Selection: OS requirement" << std::endl;

  // Apply on-Z mass requirement for dilepton pair
  if (dilepton_invMass < ( cfg_DileptonInvariantMass - cfg_DileptonInvariantMassMinus ) ) bOnZMassBelow = false;
  if (dilepton_invMass > ( cfg_DileptonInvariantMass + cfg_DileptonInvariantMassPlus ) )  bOnZMassAbove = false;
  bOnZMass = bOnZMassAbove * bOnZMassBelow;
  if (bOnZMass) cLeptonMassCounter.increment();

  // Apply dR cut for dilepton system?
  bPassDeltaR = cfg_DileptonDeltaR.passedCut(dilepton_dR);
  if(!bPassDeltaR) return;

  // Fill histos
  hDileptonPt_AfterLeptonSelection->Fill(dilepton_pt);
  hDileptonEta_AfterLeptonSelection->Fill(dilepton_eta);
  hDileptonMass_AfterLeptonSelection->Fill(dilepton_invMass);
  hDileptonCharge_AfterLeptonSelection->Fill(dilepton_charge);
  hDileptonDeltaR_AfterLeptonSelection->Fill(dilepton_dR);
  hLeptonN_AfterLeptonSelection ->Fill(muData.getSelectedMuons().size());
  for(unsigned int i=0; i< muData.getSelectedMuons().size(); i++)
    {
      hLeptonPt_AfterLeptonSelection ->Fill(muData.getSelectedMuons()[i].pt());
      hLeptonEta_AfterLeptonSelection->Fill(muData.getSelectedMuons()[i].eta());
    }

  
  //================================================================================================   
  // Tau Selection
  //================================================================================================   
  if (0) std::cout << "=== Tau Selection (Loose)" << std::endl;
  const TauSelection::Data looseTauData = fLooseTauSelection.analyzeLoose(fEvent);

  if (!looseTauData.hasIdentifiedTaus() ) return;
  // if (looseTauData.getSelectedTaus().size() != 1) return; // needed?
  if (0) std::cout << "=== Tau Selection: Has at least 1 loose tau(s)" << std::endl;
  // cTauNCounter.increment(); 

  // Tau-related variables
  int looseTauSrcBit = getTauSrcBit(fEvent, looseTauData);
  bool looseTauFake  = getIsFakeTau(fEvent, looseTauData);
  
  // Fill histos ( also sets value for boolean bIsGenuineTau
  fCommonPlots.fillControlPlotsAfterTauSelection(fEvent, looseTauData); // uses: OBbIsGenuineTau = data.getSelectedTaus()[0].isGenuineTau();
  hDileptonMass_AfterTauSelection->Fill(dilepton_invMass);

  // Redefine what the "bGenuineTau" boolean is. N.B: Affects the genuineTau plots filled by control plots!
  fCommonPlots.setGenuineTauStatus(!looseTauFake); // CommonPlots.cc (L1292) bIsGenuineTau = data.getSelectedTaus()[0].isGenuineTau();


  //================================================================================================
  // Jet selection
  //================================================================================================
  if (0) std::cout << "=== Jet selection" << std::endl;
  const JetSelection::Data jetData = fJetSelection.analyze(fEvent, looseTauData.getSelectedTau());
  if (!jetData.passedSelection()) return;  
    
 
  //================================================================================================  
  // BJet selection
  //================================================================================================
  if (0) std::cout << "=== BJet selection" << std::endl;
  const BJetSelection::Data bjetData = fBJetSelection.analyze(fEvent, jetData);
  if (!bjetData.passedSelection()) return;
  fCommonPlots.fillControlPlotsAfterBjetSelection(fEvent, bjetData);
  if (fEvent.isMC()) 
    {
      fEventWeight.multiplyWeight(bjetData.getBTaggingScaleFactorEventWeight());
    }
  fCommonPlots.fillControlPlotsAfterBtagSF(fEvent, jetData, bjetData);
  cBTaggingSFCounter.increment();
  

  //================================================================================================
  // MET selection
  //================================================================================================
  if (0) std::cout << "=== MET selection" << std::endl;
  // const METSelection::Data metData = fMETSelection.silentAnalyze(fEvent, nVertices);
  const METSelection::Data metData = fMETSelection.analyze(fEvent, nVertices);
  if (!metData.passedSelection()) return;


  //================================================================================================
  // Preselections
  //================================================================================================
  if (0) std::cout << "=== Preselections" << std::endl;
  cPreselected.increment();
  fCommonPlots.fillControlPlotsAfterPreselections(fEvent, jetData, bjetData, metData);
  
  // Get decay mode (DM)
  unsigned int looseTau_dm = looseTauData.getSelectedTau().decayMode();
  if (fEvent.isMC() )
    {
      hTauSrc_Preselections->Fill(looseTauSrcBit);
      if (looseTau_dm == 0)  hTauSrcDM0_Preselections->Fill(looseTauSrcBit);
      if (looseTau_dm == 1)  hTauSrcDM1_Preselections->Fill(looseTauSrcBit);
      if (looseTau_dm == 10) hTauSrcDM10_Preselections->Fill(looseTauSrcBit);
    }

  // Fill other histos
  hNTau_Preselections->Fill( looseTauData.getSelectedTaus().size() );
  hNJet_Preselections ->Fill( jetData.getSelectedJets().size() );
  hNBjet_Preselections->Fill( bjetData.getSelectedBJets().size() );
  hDileptonPt_Preselections->Fill(dilepton_pt);
  hDileptonEta_Preselections->Fill(dilepton_eta);
  hDileptonMass_Preselections->Fill(dilepton_invMass);
  if (bOnZMass) hDileptonMassOnZ_Preselections ->Fill(dilepton_invMass);
  else hDileptonMassOffZ_Preselections->Fill(dilepton_invMass);
  hDileptonCharge_Preselections->Fill(dilepton_charge);
  hDileptonDeltaR_Preselections->Fill(dilepton_dR);
  hMet_Preselections->Fill(metData.getMET().R());
  hLt_Preselections->Fill(dilepton_pt + looseTauData.getSelectedTau().pt());


  //================================================================================================
  // "Tight" Tau
  //================================================================================================
  const TauSelection::Data tauData = fTauSelection.analyze(fEvent);
  if (tauData.hasIdentifiedTaus()) 
    {

      // Apply "tight" tau trigger and ID scale factor (SF)
      if (fEvent.isMC()) 
	{
	  fEventWeight.multiplyWeight(tauData.getTauIDSF());
	  cTauSFCounter.increment(); 

	  fEventWeight.multiplyWeight(tauData.getTauMisIDSF());
	  cFakeTauSFCounter.increment();
	}

      // Do rest of event selection
      doTightTaus(fEvent, tauData, bOnZMass);
      fCommonPlots.fillControlPlotsAfterSelections(fEvent);
      cSelected.increment();

      int tauSrcBit = getTauSrcBit(fEvent, looseTauData);
      unsigned int tau_dm = looseTauData.getSelectedTau().decayMode();

      // Fill histograms
      if (fEvent.isMC() )
	{
	  hTauSrc_TightTau->Fill(tauSrcBit);
	  if (tau_dm == 0)  hTauSrcDM0_TightTau->Fill(tauSrcBit);
	  if (tau_dm == 1)  hTauSrcDM1_TightTau->Fill(tauSrcBit);
	  if (tau_dm == 10) hTauSrcDM10_TightTau->Fill(tauSrcBit);
	}

      // Fill other histos
      hNTau_TightTau->Fill( tauData.getSelectedTaus().size() );
      hNJet_TightTau ->Fill( jetData.getSelectedJets().size() );
      hNBjet_TightTau->Fill( bjetData.getSelectedBJets().size() );
      hDileptonPt_TightTau->Fill(dilepton_pt);
      hDileptonEta_TightTau->Fill(dilepton_eta);
      hDileptonMass_TightTau->Fill(dilepton_invMass);
      hMet_TightTau->Fill(metData.getMET().R());
      hLt_TightTau->Fill(dilepton_pt + tauData.getSelectedTau().pt());
      
      // Undo "tight" tau trigger and ID scale factor (SF)
      if (fEvent.isMC())
	{ 
	  fEventWeight.multiplyWeight(1.0/tauData.getTauIDSF());
	  fEventWeight.multiplyWeight(1.0/tauData.getTauMisIDSF());
	}
    }


  //================================================================================================
  // "Loose" Tau
  //================================================================================================
  if (looseTauData.hasIdentifiedTaus()) 
    {
      
      // Apply "loose" tau trigger and ID scale factor (SF)
      if (fEvent.isMC())
	{
	  fEventWeight.multiplyWeight(looseTauData.getTauIDSF());
	  fEventWeight.multiplyWeight(looseTauData.getTauMisIDSF());
	}

      // Do rest of event selection
      doLooseTaus(fEvent, looseTauData, bOnZMass);

      int looseTauSrcBit = getTauSrcBit(fEvent, looseTauData);
      unsigned int looseTau_dm = looseTauData.getSelectedTau().decayMode();

      // Fill histograms
      if (fEvent.isMC())
	{
	  hTauSrc_LooseTau->Fill(looseTauSrcBit);
	  if (looseTau_dm == 0)  hTauSrcDM0_LooseTau->Fill(looseTauSrcBit);
	  if (looseTau_dm == 1)  hTauSrcDM1_LooseTau->Fill(looseTauSrcBit);
	  if (looseTau_dm == 10) hTauSrcDM10_LooseTau->Fill(looseTauSrcBit);
	}
      
      // Fill other histos
      hNTau_LooseTau->Fill( looseTauData.getSelectedTaus().size() );
      hNJet_LooseTau ->Fill( jetData.getSelectedJets().size() );
      hNBjet_LooseTau->Fill( bjetData.getSelectedBJets().size() );
      hDileptonPt_LooseTau->Fill(dilepton_pt);
      hDileptonEta_LooseTau->Fill(dilepton_eta);
      hDileptonMass_LooseTau->Fill(dilepton_invMass);
      hMet_LooseTau->Fill(metData.getMET().R());
      hLt_LooseTau->Fill(dilepton_pt + looseTauData.getSelectedTau().pt());
    }

  fEventSaver.save();

  return;
}


void TauFakeRate_mm::doLooseTaus(const Event& event, const TauSelection::Data& tauData, const bool bOnZMass)
{
  
  if (!bOnZMass) return;

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
	  // inclusive
	  hTauPt_den_g ->Fill( tau_pt );
	  hTauEta_den_g->Fill( tau_eta ); 
	  if (tau_inBarrel) hTauPt_den_g_barrel->Fill( tau_pt );
	  else hTauPt_den_g_endcap->Fill( tau_pt );

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
	}// (!bFakeTau)

      
      // inclusive 
      // std::cout << "Filling hTauPt_den with pt = " << tau_pt << " GeV. DM = " << tau_dm << std::endl;
      hTauPt_den ->Fill( tau_pt );
      hTauEta_den->Fill( tau_eta );
      if (tau_inBarrel) hTauPt_den_barrel->Fill( tau_pt );
      else hTauPt_den_endcap->Fill( tau_pt );
    
      // decay mode (DM) = 0
      if (tau_dm==0) 
	{
	  hTauPt_den_dm0 ->Fill( tau_pt  );
	  hTauEta_den_dm0->Fill( tau_eta );
	  if (tau_inBarrel) hTauPt_den_dm0_barrel ->Fill( tau_pt  );
	  else hTauPt_den_dm0_endcap ->Fill( tau_pt  );
	}
      
      // decay mode (DM) = 1
      if (tau_dm==1) 
	{
	  hTauPt_den_dm1 ->Fill( tau_pt  );
	  hTauEta_den_dm1->Fill( tau_eta );
	  if (tau_inBarrel) hTauPt_den_dm1_barrel ->Fill( tau_pt  );
	  else hTauPt_den_dm1_endcap ->Fill( tau_pt  );
	}
      
      // decay mode (DM) = 10
      if (tau_dm==10) 
	{
	  hTauPt_den_dm10 ->Fill( tau_pt  );
	  hTauEta_den_dm10->Fill( tau_eta );
	  if (tau_inBarrel) hTauPt_den_dm10_barrel ->Fill( tau_pt  );
	  else hTauPt_den_dm10_endcap ->Fill( tau_pt  );
	}
    }
  
  return;
}


void TauFakeRate_mm::doTightTaus(const Event& event, const TauSelection::Data& tauData, const bool bOnZMass) 
{
  
  if (!bOnZMass) return;

  // For-loop: All selected taus
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
	  bGenuineTau = tauData.getSelectedTaus()[i].isGenuineTau();
	  bEleToTau   = tauData.getSelectedTaus()[i].isElectronToTau();
	  bMuonToTau  = tauData.getSelectedTaus()[i].isMuonToTau();
	  bFakeTau    = (!bGenuineTau && !(bEleToTau || bMuonToTau) );
	}

      // Only if MC  and selected tau is genuine. (if not genuine tau, reject the events)
      //if (event.isMC() && tauData.getSelectedTaus()[i].isGenuineTau()) 
      if (bIsMC && !bFakeTau)
	{
	  // inclusive
	  hTauPt_num_g ->Fill( tau_pt ); 
	  hTauEta_num_g->Fill( tau_eta ); 
	  if (tau_inBarrel) hTauPt_num_g_barrel->Fill( tau_pt );
	  else hTauPt_num_g_endcap->Fill( tau_pt );

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
	}// (!bFakeTau)

      // inclusive
      hTauPt_num ->Fill( tau_pt ); 
      hTauEta_num->Fill( tau_eta );
      if (tau_inBarrel) hTauPt_num_barrel->Fill( tau_pt );
      else hTauPt_num_endcap->Fill( tau_pt );

      
      // decay mode (DM) = 0
      if (tau_dm==0) 
	{
	  hTauPt_num_dm0 ->Fill( tau_pt  );
	  hTauEta_num_dm0->Fill( tau_eta );

	  if (tau_inBarrel) hTauPt_num_dm0_barrel ->Fill( tau_pt  );
	  else hTauPt_num_dm0_endcap ->Fill( tau_pt  );
	}
      
      // decay mode (DM) = 1
      if (tau_dm==1) 
	{
	  hTauPt_num_dm1 ->Fill( tau_pt  );
	  hTauEta_num_dm1->Fill( tau_eta );
	  if (tau_inBarrel) hTauPt_num_dm1_barrel ->Fill( tau_pt  );
	  else hTauPt_num_dm1_endcap ->Fill( tau_pt  );
	}
      
      // decay mode (DM) = 10
      if (tau_dm==10) 
	{
	  hTauPt_num_dm10 ->Fill( tau_pt  );
	  hTauEta_num_dm10->Fill( tau_eta );
	  if (tau_inBarrel) hTauPt_num_dm10_barrel ->Fill( tau_pt  );
	  else hTauPt_num_dm10_endcap ->Fill( tau_pt  );
	}
    }

  return;
}


int TauFakeRate_mm::getTauSrcBit(const Event& event, const TauSelection::Data& tauData)
{
  
  // Definitions
  bool bEleToTau     = false; // ele  -> tau fakes only
  bool bMuonToTau    = false; // muon -> tau fakes only
  // bool bJetToTau     = false; // jet  -> tau fakes only
  bool bLightQToTau  = false; // light quark-> tau fakes only (u,d,s)
  bool bHeavyQToTau  = false; // heavy quark-> tau fakes only (c,b)
  bool bGluonToTau   = false; // jet  -> tau fakes only
  bool bUnknownToTau = false; // unknown -> tau fakes only

  if ( event.isMC() )
    {
      bEleToTau     = tauData.getSelectedTau().isElectronToTau();
      bMuonToTau    = tauData.getSelectedTau().isMuonToTau();
      // bJetToTau     = tauData.getSelectedTau().isJetToTau();
      bLightQToTau  = tauData.getSelectedTau().isQuarkToTau(1) || tauData.getSelectedTau().isQuarkToTau(2) || tauData.getSelectedTau().isQuarkToTau(3);
      bHeavyQToTau  = tauData.getSelectedTau().isQuarkToTau(4) || tauData.getSelectedTau().isQuarkToTau(5);
      bGluonToTau   = tauData.getSelectedTau().isGluonToTau();
      bUnknownToTau = tauData.getSelectedTau().isUnknownTauDecay();
    }

  // Define a bit to store the source of fake. (tauSrcBit = 0 for genuineTau)
  int tauSrcBit =  pow(2,0)*bEleToTau + pow(2,1)*bMuonToTau + pow(2,2)*bLightQToTau + pow(2,3)*bHeavyQToTau + pow(2,4)*bGluonToTau + pow(2,5)*bUnknownToTau;
  return tauSrcBit;
}


bool TauFakeRate_mm::getIsFakeTau(const Event& event, const TauSelection::Data& tauData) 
{
  
  // Definitions
  bool bEleToTau   = false;
  bool bMuonToTau  = false;
  bool bGenuineTau = false;

  if ( event.isMC() )
    {
      bEleToTau     = tauData.getSelectedTau().isElectronToTau();
      bMuonToTau    = tauData.getSelectedTau().isMuonToTau();
      bGenuineTau   = tauData.getSelectedTau().isGenuineTau();
    }

  bool bFakeTau = ( event.isMC() && !bGenuineTau && !(bEleToTau || bMuonToTau) );

  return bFakeTau;
}
