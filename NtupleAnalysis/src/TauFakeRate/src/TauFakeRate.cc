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
  Count cMuonOSCounter;
  Count cMuonMassCounter;
  TauSelection fTauSelection;
  TauSelection fLooseTauSelection;
  Count cLooseTauNCounter;
  // Count cTauNCounter;
  Count cTauSFCounter;
  Count cFakeTauSFCounter;
  // Count cTauDMCounter;
  JetSelection fJetSelection;
  BJetSelection fBJetSelection;
  Count cBTaggingSFCounter;
  METSelection fMETSelection;
  // FatJetSelection fFatJetSelection;
  Count cSelected;
    
  void doLooseTaus(const Event& event, const TauSelection::Data& tauData, const int nVertices);
  void doTightTaus(const Event& event, const TauSelection::Data& looseTauData, const int nVertices);

  // Non-common histograms
  WrappedTH1 *hTauPt_num_1pr;
  WrappedTH1 *hTauPt_num_2pr;
  WrappedTH1 *hTauPt_num_3pr;
  WrappedTH1 *hTauPt_den_1pr;
  WrappedTH1 *hTauPt_den_2pr;
  WrappedTH1 *hTauPt_den_3pr;

  WrappedTH1 *hTauPt_num_g_1pr;
  WrappedTH1 *hTauPt_num_g_2pr;
  WrappedTH1 *hTauPt_num_g_3pr;
  WrappedTH1 *hTauPt_den_g_1pr;
  WrappedTH1 *hTauPt_den_g_2pr;
  WrappedTH1 *hTauPt_den_g_3pr;

  WrappedTH1 *hTauEta_num_1pr;
  WrappedTH1 *hTauEta_num_2pr;
  WrappedTH1 *hTauEta_num_3pr;
  WrappedTH1 *hTauEta_den_1pr;
  WrappedTH1 *hTauEta_den_2pr;
  WrappedTH1 *hTauEta_den_3pr;

  WrappedTH1 *hMuonN;
  WrappedTH1 *hMuonPt;
  WrappedTH1 *hMuonEta;

  WrappedTH1 *hNTau;
  WrappedTH1 *hTauPt;
  WrappedTH1 *hTauEta;

  WrappedTH1 *hNJet;
  WrappedTH1 *hMET;

  WrappedTH1 *hDileptonMass_BeforeLeptonSelection;
  WrappedTH1 *hDileptonMass_AfterLeptonSelection;
  WrappedTH1 *hDileptonMass_AfterTauSelection;
  WrappedTH1 *hDileptonMass_AfterJetSelection;
  WrappedTH1 *hDileptonMass_AfterBJetSelection;
  WrappedTH1 *hDileptonMass_AfterMetSelection;
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
    cMuonOSCounter(fEventCounter.addCounter("#mu OS")),
    cMuonMassCounter(fEventCounter.addCounter("m_{#mu#mu} window")),
    fTauSelection(config.getParameter<ParameterSet>("TauSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fLooseTauSelection(config.getParameter<ParameterSet>("LooseTauSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cLooseTauNCounter(fEventCounter.addCounter("has loose #tau's")),
    // cTauNCounter(fEventCounter.addCounter("#tau N")),
    cTauSFCounter(fEventCounter.addCounter("#tau SF")),
    cFakeTauSFCounter(fEventCounter.addCounter("Fake #tau SF")),
    // cTauDMCounter(fEventCounter.addCounter("#tau DM")),
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

  hTauPt_num_1pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_1pr", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_2pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_2pr", "; p_{T} (GeV)", 7, bin);
  hTauPt_num_3pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_3pr", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_1pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_1pr", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_2pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_2pr", "; p_{T} (GeV)", 7, bin);
  hTauPt_den_3pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_3pr", "; p_{T} (GeV)", 7, bin);

  hTauPt_num_g_1pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_1pr", ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_2pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_2pr", ";p_{T} (GeV)", 7, bin);
  hTauPt_num_g_3pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_num_g_3pr", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_1pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_1pr", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_2pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_2pr", ";p_{T} (GeV)", 7, bin);
  hTauPt_den_g_3pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt_den_g_3pr", ";p_{T} (GeV)", 7, bin);

  hTauEta_num_1pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_1pr", ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_2pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_2pr", ";#eta", etaN, etaMin, etaMax);
  hTauEta_num_3pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_num_3pr", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_1pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_1pr", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_2pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_2pr", ";#eta", etaN, etaMin, etaMax);
  hTauEta_den_3pr =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta_den_3pr", ";#eta", etaN, etaMin, etaMax);

  hMuonN   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "muN"  , ";#mu multiplicity", nN  , nMin  , nMax );
  hMuonPt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "muPt" , ";p_{T} (GeV)"    , ptN , ptMin , ptMax );
  hMuonEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "muEta", ";#eta"           , etaN, etaMin, etaMax);

  // hTauN    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauN"  , ";#tau-jet multiplicity", nN  , nMin  , nMax  );
  // hTauPt   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt" , ";p_{T} (GeV)"          , ptN , ptMin , ptMax );
  // hTauEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta", ";#eta"                 , etaN, etaMin, etaMax);

  hMET  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "MET" , "; E_{T}^{miss} (GeV)", metN, metMin, metMax);
  hNJet = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "nJet", "; jet multiplicity"  , nN  , nMin  , nMax  );

  hDileptonMass_BeforeLeptonSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_BeforeLeptonSelection", ";m_{ll} (GeV)", mN, mMin, mMax);
  hDileptonMass_AfterLeptonSelection  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_AfterLeptonSelection" , ";m_{ll} (GeV)", mN, mMin, mMax);
  hDileptonMass_AfterTauSelection     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_AfterTauSelection"    , ";m_{ll} (GeV)", mN, mMin, mMax);
  hDileptonMass_AfterJetSelection     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_AfterJetSelection"    , ";m_{ll} (GeV)", mN, mMin, mMax);
  hDileptonMass_AfterBJetSelection    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_AfterBJetSelection"   , ";m_{ll} (GeV)", mN, mMin, mMax);
  hDileptonMass_AfterMetSelection     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_AfterMetSelection"    , ";m_{ll} (GeV)", mN, mMin, mMax);
  hDileptonMass_AfterAllSelections    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "DileptonMass_AfterAllSelections"   , ";m_{ll} (GeV)", mN, mMin, mMax);

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
  const double cfg_massWindow = 20.0; // 15.0;

  // Temporary Trigger SF
  if (fEvent.isMC()) 
    {
      if ( 26 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() <  30) fEventWeight.multiplyWeight(0.9664);
      else if ( 30 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() <  40) fEventWeight.multiplyWeight(0.9781);
      else if ( 40 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() <  50) fEventWeight.multiplyWeight(0.9819);
      else if ( 50 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() <  60) fEventWeight.multiplyWeight(0.9822);
      else if ( 60 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() <  80) fEventWeight.multiplyWeight(0.9804);
      else if ( 80 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 120) fEventWeight.multiplyWeight(0.9780);
      else if (120 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 200) fEventWeight.multiplyWeight(0.9752);
      else if (200 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 500) fEventWeight.multiplyWeight(0.9704);
      else {} //std::cout << "=== Muon Selection: No SF support for muons with Pt = " << muData.getSelectedMuons()[0].pt() << " GeV and eta = " << muData.getSelectedMuons()[0].eta() << std::endl;
    }

  // Temporary muon ID SF 
  if (fEvent.isMC()) {
    if (20 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 25) fEventWeight.multiplyWeight(0.991);
    else if (25 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 30) fEventWeight.multiplyWeight(0.982);
    else if (30 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 40) fEventWeight.multiplyWeight(0.981);
    else if (40 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 50) fEventWeight.multiplyWeight(0.981);
    else if (50 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 60) fEventWeight.multiplyWeight(0.978);
    else if (60 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 120) fEventWeight.multiplyWeight(0.986);
    else if (120 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 200) fEventWeight.multiplyWeight(1.024);
    else {}
  }


  // Require the muons to have opposite sign (OS)
  if( muData.getSelectedMuons()[0].charge() == muData.getSelectedMuons()[1].charge()) return;
  cMuonOSCounter.increment(); 
  if (0) std::cout << "=== Muon Selection: OS requirement" << std::endl;

  // Apply on-Z mass requirement for dilepton pair
  hDileptonMass_BeforeLeptonSelection->Fill(dilepton_invMass);
  if (dilepton_invMass < ( 91.1876 - cfg_massWindow ) ) return; // fixme. define in run.py 
  if (dilepton_invMass > ( 91.1876 + cfg_massWindow ) ) return; // fixme. define in run.py
  cMuonMassCounter.increment(); 

  // Apply dR cut for dilepton system?
  //  dilepton_dR = ROOT::Math::VectorUtil::DeltaR(muData.getSelectedMuons()[0].p4(),muData.getSelectedMuons()[1].p4());
  //  if(dilepton_dR < 0.3) return;


  // Fill histos
  fCommonPlots.fillControlPlotsAtMuonSelection(fEvent, muData);
  hDileptonMass_AfterLeptonSelection->Fill(dilepton_invMass);
  hMuonN ->Fill(muData.getSelectedMuons().size());
  // For-loop: All Selected muons
  for(unsigned int i=0; i< muData.getSelectedMuons().size(); i++)
    {
      hMuonPt ->Fill(muData.getSelectedMuons()[i].pt());
      hMuonEta->Fill(muData.getSelectedMuons()[i].eta());
    }
  

  //================================================================================================   
  // Tau Selection
  //================================================================================================   
  if (0) std::cout << "=== Tau Selection" << std::endl;
  const TauSelection::Data tauData = fTauSelection.analyze(fEvent);
  const TauSelection::Data looseTauData = fLooseTauSelection.analyzeLoose(fEvent);

  if (!looseTauData.hasIdentifiedTaus() ) return;
  if (0) std::cout << "=== Tau Selection: Has Identified taus" << std::endl;
  cLooseTauNCounter.increment(); 

  /*
  // -> Require # of taus
  if(looseTauData.getSelectedTaus().size() != 1) return;
  cTauNCounter.increment(); 
  if (0) std::cout << "=== Tau Selection: Exactly one selected tau" << std::endl;
  */
  
  /*
  // Determine if genuine or fake tau
  bool bGenuineTau      = false; // 
  bool bFakeTauElectron = false; // electron -> tau fakes only
  bool bFakeTauMuon     = false; // muon     -> tau fakes only
  bool bFakeTauJets     = false; // jet      -> tau fakes only
  if ( fEvent.isMC() )
    {
      bGenuineTau      = tauData.isGenuineTau();
      bFakeTauElectron = tauData.getSelectedTau().isElectronToTau();
      bFakeTauMuon     = tauData.getSelectedTau().isMuonToTau();
    }
  bFakeTauJets = (!bGenuineTau && !(bFakeTauElectron || bFakeTauMuon) );
  int isGenuineTau = int(!(bFakeTauJets));

  // Overwrite default boolean bIsGenuineTau = (data.getSelectedTaus()[0].isGenuineTau() || data.getSelectedTaus()[1].isGenuineTau()) 
  // [NOTE: Filled when  "fillControlPlotsAfterTauSelection" is called()]
  if (0) fCommonPlots.setGenuineTauStatus(isGenuineTau);
  // if (0) std::cout << "=== Tau Selection: isGenuineTau = " << isGenuineTau << std::endl;
  // cTauDMCounter.increment(); 
  */

  // Fill histos ( also sets value for boolean bIsGenuineTau
  hDileptonMass_AfterTauSelection->Fill(dilepton_invMass);
  fCommonPlots.fillControlPlotsAfterTauSelection(fEvent, looseTauData);


  //================================================================================================
  // Jet selection
  //================================================================================================
  if (0) std::cout << "=== Jet selection" << std::endl;
  const JetSelection::Data jetData = fJetSelection.analyze(fEvent, looseTauData.getSelectedTau());
  if (!jetData.passedSelection()) return;  
    
  // Fill histos
  fCommonPlots.fillControlPlotsAtJetSelection(fEvent, jetData);
  hDileptonMass_AfterJetSelection->Fill(dilepton_invMass);
  hNJet->Fill( jetData.getSelectedJets().size() );

  
  //================================================================================================  
  // BJet selection
  //================================================================================================
  if (0) std::cout << "=== BJet selection" << std::endl;
  const BJetSelection::Data bjetData = fBJetSelection.analyze(fEvent, jetData);
  if (!bjetData.passedSelection()) return;

  // Fill histos
  fCommonPlots.fillControlPlotsAtBtagging(fEvent, bjetData);

  //================================================================================================  
  // BJet SF  
  //================================================================================================
  if (0) std::cout << "=== BJet SF" << std::endl;
  if (fEvent.isMC()) 
    {
      fEventWeight.multiplyWeight(bjetData.getBTaggingScaleFactorEventWeight());
    }
  cBTaggingSFCounter.increment();

  // Fill histos
  fCommonPlots.fillControlPlotsAfterBtagSF(fEvent, jetData, bjetData);
  hDileptonMass_AfterBJetSelection->Fill(dilepton_invMass);

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
  hMET->Fill(metData.getMET().R());


  //================================================================================================
  // All Selections
  //================================================================================================
  if (0) std::cout << "=== All Selections" << std::endl;
  cSelected.increment();

  // Fill histos
  fCommonPlots.fillControlPlotsAfterAllSelections(fEvent);
  hDileptonMass_AfterAllSelections->Fill(dilepton_invMass);

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
      // Only if MC  and selected tau is genuine (not fake)
      if (event.isMC() && tauData.getSelectedTaus()[i].isGenuineTau()) 
	{
	  if (tauData.getSelectedTaus()[i].decayMode()==0)  hTauPt_den_g_1pr->Fill( tauData.getSelectedTaus()[i].pt() );
	  if (tauData.getSelectedTaus()[i].decayMode()==1)  hTauPt_den_g_2pr->Fill( tauData.getSelectedTaus()[i].pt() );
	  if (tauData.getSelectedTaus()[i].decayMode()==10) hTauPt_den_g_3pr->Fill( tauData.getSelectedTaus()[i].pt() );
	}
      
      // 1-prong decays; decay mode (DM) = 0
      if (tauData.getSelectedTaus()[i].decayMode()==0) 
      {
	hTauPt_den_1pr ->Fill( tauData.getSelectedTaus()[i].pt()  );
	hTauEta_den_1pr->Fill( tauData.getSelectedTaus()[i].eta() );
      }

      // 2-prong decays; decay mode (DM) = 1
      if (tauData.getSelectedTaus()[i].decayMode()==1) 
	{
	  hTauPt_den_2pr ->Fill( tauData.getSelectedTaus()[i].pt()  );
	  hTauEta_den_2pr->Fill( tauData.getSelectedTaus()[i].eta() );
	}
      
      // 3-prong decays; decay mode (DM) = 10
      if (tauData.getSelectedTaus()[i].decayMode()==10) 
	{
	  hTauPt_den_3pr ->Fill( tauData.getSelectedTaus()[i].pt()  );
	  hTauEta_den_3pr->Fill( tauData.getSelectedTaus()[i].eta() );
	}
    }
  
  return;
}


void TauFakeRate::doTightTaus(const Event& event, const TauSelection::Data& tightTauData, const int nVertices) {

  // For-loop: All selected taus
  for(unsigned int i=0; i<tightTauData.getSelectedTaus().size(); i++)
    {
      
      // Only if MC  and selected tau is genuine (not fake)
      if (event.isMC() && tightTauData.getSelectedTaus()[i].isGenuineTau()) 
	{
	  if (tightTauData.getSelectedTaus()[i].decayMode()==0)  hTauPt_num_g_1pr->Fill( tightTauData.getSelectedTaus()[i].pt() );
	  if (tightTauData.getSelectedTaus()[i].decayMode()==1)  hTauPt_num_g_2pr->Fill( tightTauData.getSelectedTaus()[i].pt() );
	  if (tightTauData.getSelectedTaus()[i].decayMode()==10) hTauPt_num_g_3pr->Fill( tightTauData.getSelectedTaus()[i].pt() );
	}
      
      // 1-prong decays; decay mode (DM) = 0
      if (tightTauData.getSelectedTaus()[i].decayMode()==0) 
	{
	  hTauPt_num_1pr ->Fill( tightTauData.getSelectedTaus()[i].pt()  );
	  hTauEta_num_1pr->Fill( tightTauData.getSelectedTaus()[i].eta() );
	}
      
      // 2-prong decays; decay mode (DM) = 1
      if (tightTauData.getSelectedTaus()[i].decayMode()==1) 
	{
	  hTauPt_num_2pr ->Fill( tightTauData.getSelectedTaus()[i].pt()  );
	  hTauEta_num_2pr->Fill( tightTauData.getSelectedTaus()[i].eta() );
	}
      
      // 3-prong decays; decay mode (DM) = 10
      if (tightTauData.getSelectedTaus()[i].decayMode()==10) 
	{
	  hTauPt_num_3pr ->Fill( tightTauData.getSelectedTaus()[i].pt()  );
	  hTauEta_num_3pr->Fill( tightTauData.getSelectedTaus()[i].eta() );
	}
      
    }
  
  return;
}
