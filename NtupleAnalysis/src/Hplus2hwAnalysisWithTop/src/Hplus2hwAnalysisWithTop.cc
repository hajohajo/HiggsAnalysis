// -*- c++ -*-
#include "Framework/interface/BaseSelector.h"
#include "Framework/interface/makeTH.h"

#include "EventSelection/interface/CommonPlots.h"
#include "EventSelection/interface/EventSelections.h"
#include "EventSelection/interface/TransverseMass.h"


#include "TDirectory.h"

class Hplus2hwAnalysisWithTop: public BaseSelector {
public:
  explicit Hplus2hwAnalysisWithTop(const ParameterSet& config, const TH1* skimCounters);
  virtual ~Hplus2hwAnalysisWithTop() {}

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
  Count cTauNCounter;
  // Count cTauOSCounter;
  Count cTauSFCounter;
  Count cFakeTauSFCounter;
  Count cTauDMCounter;
  JetSelection fJetSelection;
  AngularCutsCollinear fAngularCutsCollinear;
  BJetSelection fBJetSelection;
  Count cBTaggingSFCounter;
  METSelection fMETSelection;
  AngularCutsBackToBack fAngularCutsBackToBack;
  TopSelectionMVA fTopSelection;
  Count cTopTaggingSFCounter;
  // FatJetSelection fFatJetSelection;
  Count cSelected;
    
  // Non-common histograms
  // WrappedTH1 *hTauTauMass_AfterTauSelection;
  // WrappedTH1 *hTauTauMass_AfterJetSelection;
  // WrappedTH1 *hTauTauMass_AfterBjetSelection;
  // WrappedTH1 *hTauTauMass_AfterMetSelection;
  // WrappedTH1 *hTauTauMass_AfterTopSelection;
  // WrappedTH1 *hTauTauMass_AfterAllSelections;

  // WrappedTH1 *hTransverseMass_AfterTauSelection;
  WrappedTH1 *hTransverseMass_AfterJetSelection;
  WrappedTH1 *hTransverseMass_AfterBjetSelection;
  WrappedTH1 *hTransverseMass_AfterMetSelection;
  WrappedTH1 *hTransverseMass_AfterTopSelection;
  WrappedTH1 *hTransverseMass_AfterAllSelections;

};

#include "Framework/interface/SelectorFactory.h"
REGISTER_SELECTOR(Hplus2hwAnalysisWithTop);

Hplus2hwAnalysisWithTop::Hplus2hwAnalysisWithTop(const ParameterSet& config, const TH1* skimCounters)
  : BaseSelector(config, skimCounters),
    fCommonPlots(config.getParameter<ParameterSet>("CommonPlots"), CommonPlots::kHplus2hwAnalysisWithTop, fHistoWrapper),
    cAllEvents(fEventCounter.addCounter("all events")),
    cTrigger(fEventCounter.addCounter("passed trigger")),
    fMETFilterSelection(config.getParameter<ParameterSet>("METFilter"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cVertexSelection(fEventCounter.addCounter("passed PV")),
    fElectronSelection(config.getParameter<ParameterSet>("ElectronSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    fMuonSelection(config.getParameter<ParameterSet>("MuonSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fTauSelection(config.getParameter<ParameterSet>("TauSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cTauNCounter(fEventCounter.addCounter("#tau N")),
    // cTauOSCounter(fEventCounter.addCounter("#tau OS")),
    cTauSFCounter(fEventCounter.addCounter("#tau SF")),
    cFakeTauSFCounter(fEventCounter.addCounter("Fake #tau SF")),
    cTauDMCounter(fEventCounter.addCounter("#tau DM")),
    fJetSelection(config.getParameter<ParameterSet>("JetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fAngularCutsCollinear(config.getParameter<ParameterSet>("AngularCutsCollinear"), fEventCounter, fHistoWrapper, &fCommonPlots, ""), // fixme: keep counter?
    // fAngularCutsCollinear(config.getParameter<ParameterSet>("AngularCutsCollinear")), // fixme: keep counter?
    fBJetSelection(config.getParameter<ParameterSet>("BJetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cBTaggingSFCounter(fEventCounter.addCounter("b-tag SF")),
    // fMETSelection(config.getParameter<ParameterSet>("METSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fMETSelection(config.getParameter<ParameterSet>("METSelection")),
    fAngularCutsBackToBack(config.getParameter<ParameterSet>("AngularCutsBackToBack"), fEventCounter, fHistoWrapper, &fCommonPlots, ""), // fixme: keep counter?
    // fAngularCutsBackToBack(config.getParameter<ParameterSet>("AngularCutsBackToBack")), // fixme: keep counter?
    fTopSelection(config.getParameter<ParameterSet>("TopSelectionMVA"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cTopTaggingSFCounter(fEventCounter.addCounter("top-tag SF")),
    // fFatJetSelection(config.getParameter<ParameterSet>("FatJetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    cSelected(fEventCounter.addCounter("Selected Events"))
{ }


void Hplus2hwAnalysisWithTop::book(TDirectory *dir) {

  if (0) std::cout << "=== Hplus2hwAnalysisWithTop::book()" << std::endl;
  // Book common plots histograms
  fCommonPlots.book(dir, isData());

  // Book histograms in event selection classes
  fMETFilterSelection.bookHistograms(dir);
  fElectronSelection.bookHistograms(dir);
  fMuonSelection.bookHistograms(dir);
  fTauSelection.bookHistograms(dir);
  fJetSelection.bookHistograms(dir);
  fAngularCutsCollinear.bookHistograms(dir);
  fBJetSelection.bookHistograms(dir);
  fMETSelection.bookHistograms(dir);
  fAngularCutsBackToBack.bookHistograms(dir);
  fTopSelection.bookHistograms(dir);
  // fFatJetSelection.bookHistograms(dir);

  // Get binning from cfg file
  // const int   ptN    = fCommonPlots.getPtBinSettings().bins();
  // const float ptMin  = fCommonPlots.getPtBinSettings().min();
  // const float ptMax  = fCommonPlots.getPtBinSettings().max();
  // const int   etaN   = fCommonPlots.getEtaBinSettings().bins();
  // const float etaMin = fCommonPlots.getEtaBinSettings().min();
  // const float etaMax = fCommonPlots.getEtaBinSettings().max();
  const int   mtN    = fCommonPlots.getMtBinSettings().bins();
  const float mtMin  = fCommonPlots.getMtBinSettings().min();
  const float mtMax  = fCommonPlots.getMtBinSettings().max();
  
  // Book non-common histograms
  // hTauTauMass_AfterTauSelection  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauTauMass_AfterTauSelection" , ";m_{vis} (GeV/c^{2})", mtN, mtMin, mtMax);
  // hTauTauMass_AfterJetSelection  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauTauMass_AfterJetSelection" , ";m_{vis} (GeV/c^{2})", mtN, mtMin, mtMax);
  // hTauTauMass_AfterBjetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauTauMass_AfterBjetSelection", ";m_{vis} (GeV/c^{2})", mtN, mtMin, mtMax);
  // hTauTauMass_AfterMetSelection  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauTauMass_AfterMetSelection" , ";m_{vis} (GeV/c^{2})", mtN, mtMin, mtMax);
  // hTauTauMass_AfterTopSelection  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauTauMass_AfterTopSelection" , ";m_{vis} (GeV/c^{2})", mtN, mtMin, mtMax);
  // hTauTauMass_AfterAllSelections = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TauTauMass_AfterAllSelections", ";m_{vis} (GeV/c^{2})", mtN, mtMin, mtMax);

  // hTransverseMass_AfterTauSelection  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TransverseMass_AfterTauSelection" , ";m_{T} (GeV/c^{2})", mtN, mtMin, mtMax);
  hTransverseMass_AfterJetSelection  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TransverseMass_AfterJetSelection" , ";m_{T} (GeV/c^{2})", mtN, mtMin, mtMax);
  hTransverseMass_AfterBjetSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TransverseMass_AfterBjetSelection", ";m_{T} (GeV/c^{2})", mtN, mtMin, mtMax);
  hTransverseMass_AfterMetSelection  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TransverseMass_AfterMetSelection" , ";m_{T} (GeV/c^{2})", mtN, mtMin, mtMax);
  hTransverseMass_AfterTopSelection  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TransverseMass_AfterTopSelection" , ";m_{T} (GeV/c^{2})", mtN, mtMin, mtMax);
  hTransverseMass_AfterAllSelections = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TransverseMass_AfterAllSelections", ";m_{T} (GeV/c^{2})", mtN, mtMin, mtMax);

  return;
}


void Hplus2hwAnalysisWithTop::setupBranches(BranchManager& branchManager) {
  fEvent.setupBranches(branchManager);
  return;
}


void Hplus2hwAnalysisWithTop::process(Long64_t entry) {

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

  // Fill histos
  fCommonPlots.fillControlPlotsAtVertexSelection(fEvent);

  
  //================================================================================================   
  // 4) Electron veto
  //================================================================================================   
  if (0) std::cout << "=== Electron veto" << std::endl;
  const ElectronSelection::Data eData = fElectronSelection.analyze(fEvent);
  if (eData.hasIdentifiedElectrons()) return;

  // Fill histos
  fCommonPlots.fillControlPlotsAtElectronSelection(fEvent, eData);


  //================================================================================================
  // 5) Muon Selection
  //================================================================================================
  if (0) std::cout << "=== Muon Selection" << std::endl;
  const MuonSelection::Data muData = fMuonSelection.analyze(fEvent);
  if(!(muData.hasIdentifiedMuons())) return;
  if (muData.getSelectedMuons().size() != 1);

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

  // Fill histos
  fCommonPlots.fillControlPlotsAtMuonSelection(fEvent, muData);


  //================================================================================================
  // - MET selection
  //================================================================================================
  if (0) std::cout << "=== MET selection" << std::endl;
  const METSelection::Data metData = fMETSelection.silentAnalyze(fEvent, nVertices);
  // const METSelection::Data metData = fMETSelection.analyze(fEvent, nVertices);
  // if (!metData.passedSelection()) return;

  // Fill Histos
  fCommonPlots.fillControlPlotsAtMETSelection(fEvent, metData);
  // hTauTauMass_AfterMetSelection->Fill(mVis);
  // hTransverseMass_AfterMetSelection ->Fill(mT);


  //================================================================================================   
  // 6) Tau Selection
  //================================================================================================   
  if (0) std::cout << "=== Tau Selection" << std::endl;
  const TauSelection::Data tauData = fTauSelection.analyze(fEvent);

  if (!tauData.hasIdentifiedTaus() ) return;
  if (0) std::cout << "=== Tau Selection: Has Identified taus" << std::endl;
    
  // -> Require # of taus
  if(tauData.getSelectedTaus().size() != 1) return; // fixme. needed? check signal analysis
  cTauNCounter.increment(); 
  if (0) std::cout << "=== Tau Selection: 1 selected taus" << std::endl;

  // -> Require 2 taus to have opposite sign (OS)
  // if(tauData.getSelectedTaus()[0].charge() == tauData.getSelectedTaus()[1].charge()) return;
  // cTauOSCounter.increment(); 
  //if (0) std::cout << "=== Tau Selection: OS requirement" << std::endl;

  // if not genuine tau, reject the events (fake tau events are taken into account in QCDandFakeTau measurement)
  bool bGenuineTau      = false; // 
  bool bFakeTauElectron = false; // electron -> tau fakes only
  bool bFakeTauMuon     = false; // muon     -> tau fakes only
  bool bFakeTauJets     = false; // jet      -> tau fakes only
  if ( fEvent.isMC() )
    {
      bGenuineTau      = tauData.isGenuineTau();
      bFakeTauElectron = tauData.getSelectedTau().isElectronToTau();
      bFakeTauMuon     = tauData.getSelectedTau().isMuonToTau();
      // bGenuineLdgTau    = ( tauData.getSelectedTaus()[0].isGenuineTau() );
      // bGenuineSubldgTau = ( tauData.getSelectedTaus()[1].isGenuineTau() );
    }
  bFakeTauJets = (!bGenuineTau && !(bFakeTauElectron || bFakeTauMuon) );
  int isGenuineTau = int(!(bFakeTauJets));

  // Overwrite default boolean bIsGenuineTau = (data.getSelectedTaus()[0].isGenuineTau() || data.getSelectedTaus()[1].isGenuineTau()) 
  // [NOTE: Filled when  "fillControlPlotsAfterTauSelection" is called()]
  if (0) fCommonPlots.setGenuineTauStatus(isGenuineTau);
  // if (0) std::cout << "=== Tau Selection: isGenuineTau = " << isGenuineTau << " (LdgTau = " << bGenuineLdgTau << ", SubldTau = " << bGenuineSubldgTau << ")" << std::endl;
  if (0) std::cout << "=== Tau Selection: isGenuineTau = " << isGenuineTau << std::endl;
  cTauDMCounter.increment(); 

  // Tau ID SF
  if (fEvent.isMC()) 
    {
      fEventWeight.multiplyWeight(tauData.getTauIDSF());
      cTauSFCounter.increment(); 
    }

  // Tau misID SF
  if (fEvent.isMC()) 
    {
      fEventWeight.multiplyWeight(tauData.getTauMisIDSF());
      cFakeTauSFCounter.increment();
    }
  if (0) std::cout << "=== Tau Selection: SF applied" << std::endl;

  // Fill histos ( also sets value for boolean bIsGenuineTau
  fCommonPlots.fillControlPlotsAfterTauSelection(fEvent, tauData);


  //================================================================================================
  // 7) Jet selection
  //================================================================================================
  if (0) std::cout << "=== Jet selection" << std::endl;
  const JetSelection::Data jetData = fJetSelection.analyze(fEvent, tauData.getSelectedTau());
  if (!jetData.passedSelection()) return;  

  // Calculate transverse mass
  double mT  = -1.0;
  if (jetData.getSelectedJets().size() > 1) 
    {
      mT = TransverseMass::reconstruct(tauData.getSelectedTau().p2(), muData.getSelectedMuons()[0].p2(), jetData.getSelectedJets()[0].p2(), jetData.getSelectedJets()[1].p2(), metData.getMET());
    }
  

  // Fill histos
  fCommonPlots.fillControlPlotsAtJetSelection(fEvent, jetData);
  // hTauTauMass_AfterJetSelection->Fill(mVis);
  hTransverseMass_AfterJetSelection->Fill(mT);


  //================================================================================================  
  // 8) BJet selection
  //================================================================================================
  if (0) std::cout << "=== BJet selection" << std::endl;
  const BJetSelection::Data bjetData = fBJetSelection.analyze(fEvent, jetData);
  if (!bjetData.passedSelection()) return;

  // Fill histos
  fCommonPlots.fillControlPlotsAtBtagging(fEvent, bjetData);

  //================================================================================================  
  // 9) BJet SF  
  //================================================================================================
  if (0) std::cout << "=== BJet SF" << std::endl;
  if (fEvent.isMC()) 
    {
      fEventWeight.multiplyWeight(bjetData.getBTaggingScaleFactorEventWeight());
    }
  cBTaggingSFCounter.increment();

  // Fill histos
  fCommonPlots.fillControlPlotsAfterBtagSF(fEvent, jetData, bjetData);
  // hTauTauMass_AfterBjetSelection->Fill(mVis);
  hTransverseMass_AfterBjetSelection->Fill(mT);


  //================================================================================================  
  // Collinear angular cuts
  //================================================================================================  
  if (0) std::cout << "=== Angular Cuts (coll)" << std::endl;
  const AngularCutsCollinear::Data collinearData = fAngularCutsCollinear.analyze(fEvent, tauData.getSelectedTau(), jetData, metData);
  // const AngularCutsCollinear::Data collinearData = fAngularCutsCollinear.analyze(fEvent, muData.getSelectedMuons()[0], tauData.getSelectedTaus(), jetData, bjetData, metData);
  // if (!collinearData.passedSelection()) return;

  //================================================================================================
  // Back-to-back angular cuts
  //================================================================================================
  // if (0) std::cout << "=== Angular Cuts (b2b)" << std::endl;
  const AngularCutsBackToBack::Data backToBackData = fAngularCutsBackToBack.analyze(fEvent, tauData.getSelectedTau(), jetData, metData);
  // const AngularCutsBackToBack::Data backToBackData = fAngularCutsBackToBack.analyze(fEvent, muData.getSelectedMuons()[0], tauData.getSelectedTaus(), jetData, bjetData, metData);
  // if (!backToBackData.passedSelection()) return;

  /*
  //================================================================================================
  // - MET selection
  //================================================================================================
  if (0) std::cout << "=== MET selection" << std::endl;
  // const METSelection::Data metData = fMETSelection.analyze(fEvent, nVertices);
  if (!metData.passedSelection()) return;

  // Fill Histos
  fCommonPlots.fillControlPlotsAtMETSelection(fEvent, metData);
  hTauTauMass_AfterMetSelection->Fill(mVis);
  hTransverseMass_AfterMetSelection ->Fill(mT);
  */


  //================================================================================================
  // 10) Top selection
  //================================================================================================
  /*
    if (0) std::cout << "=== Top (MVA) selection" << std::endl;
  const TopSelectionMVA::Data topData = fTopSelection.analyze(fEvent, jetData, bjetData);

  // Fill histos after StandardSelections
   fCommonPlots.fillControlPlotsAfterStandardSelections(fEvent, jetData, bjetData, metData, topData);

  // Apply top selection
  if (!topData.passedSelection()) return;  

  // Apply top-tag SF
  if (fEvent.isMC()) 
    {
      fEventWeight.multiplyWeight(topData.getTopTaggingScaleFactorEventWeight());
    }
  cTopTaggingSFCounter.increment();

  // Fill histos
  // hTauTauMass_AfterTopSelection->Fill(mVis);
  hTransverseMass_AfterTopSelection ->Fill(mT);
  */

  //================================================================================================
  // All Selections
  //================================================================================================
  if (0) std::cout << "=== All Selections" << std::endl;
  cSelected.increment();

  // Fill histos after AllSelections: (After top-selections and top-tag SF)
  fCommonPlots.fillControlPlotsAfterAllSelections(fEvent);

  // Fill histos
  // hTauTauMass_AfterAllSelections->Fill(mVis);
  hTransverseMass_AfterAllSelections->Fill(mT);

  fEventSaver.save();

  return;
}
