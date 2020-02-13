// -*- c++ -*-
#include "Framework/interface/BaseSelector.h"
#include "Framework/interface/makeTH.h"

#include "EventSelection/interface/CommonPlots.h"
#include "EventSelection/interface/EventSelections.h"
#include "EventSelection/interface/TransverseMass.h"


#include "TDirectory.h"

class Hplus2hwAnalysis_background: public BaseSelector {
public:
  explicit Hplus2hwAnalysis_background(const ParameterSet& config, const TH1* skimCounters);
  virtual ~Hplus2hwAnalysis_background();

  /// Books histograms
  virtual void book(TDirectory *dir) override;
  /// Sets up branches for reading the TTree
  virtual void setupBranches(BranchManager& branchManager) override;
  /// Called for each event
  virtual void process(Long64_t entry) override;

private:
  // Input parameters

  /// Common plots
  CommonPlots fCommonPlots;

  // Event selection classes and event counters (in same order like they are applied)
  Count cAllEvents;

  METFilterSelection fMETFilterSelection;

  MuonSelection fMuonSelection;

  TauSelection fLooseTauSelection;
  TauSelection fTauSelection;

  Count cOverTwoTausCounter;

  Count cTauIDSFCounter;
  Count cFakeTauSFCounter;

  ElectronSelection fElectronSelection;

  JetSelection fJetSelection;

  BJetSelection fBJetSelection;

  METSelection fMETSelection;

  Count cSelected;

  void doBaselineAnalysis(const Event& event, const TauSelection::Data& tauData, const int nVertices);
  void doSignalAnalysis(const Event& event, const TauSelection::Data& tauData, const int nVertices);

  // Non-common histograms
  HistoSplitter::SplittedTripletTH1s hMtInvertedTauAfterStdSelections;

  HistoSplitter::SplittedTripletTH1s hBaselineTauTransverseMass;

  WrappedTH1 *hTauPt;
  WrappedTH1 *hMuonPt;
  WrappedTH1 *hNJet;

  WrappedTH1 *hTauEta;
  WrappedTH1 *hMuonEta;
  WrappedTH1 *hMET;

  WrappedTH1 *hMuonPt_afterMuonSelection;
  WrappedTH1 *hMuonEta_afterMuonSelection;

  WrappedTH1 *hNTau;
  WrappedTH1 *hTauCharge;

  WrappedTH1 *hTransverseMass;
  WrappedTH1 *hTransverseMass_genuine;
};


#include "Framework/interface/SelectorFactory.h"
REGISTER_SELECTOR(Hplus2hwAnalysis_background);


Hplus2hwAnalysis_background::Hplus2hwAnalysis_background(const ParameterSet& config, const TH1* skimCounters)
: BaseSelector(config, skimCounters),
  fCommonPlots(config.getParameter<ParameterSet>("CommonPlots"), CommonPlots::kQCDMeasurement_muon, fHistoWrapper),
  cAllEvents(fEventCounter.addCounter("All events")),
  fMETFilterSelection(config.getParameter<ParameterSet>("METFilter"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
  fMuonSelection(config.getParameter<ParameterSet>("MuonSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
  fLooseTauSelection(config.getParameter<ParameterSet>("LooseTauSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
  fTauSelection(config.getParameter<ParameterSet>("TauSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
  cOverTwoTausCounter(fEventCounter.addCounter("Over two selected tau leptons")),
  cTauIDSFCounter(fEventCounter.addCounter("Tau ID SF")),
  cFakeTauSFCounter(fEventCounter.addCounter("Fake tau SF")),
  fElectronSelection(config.getParameter<ParameterSet>("ElectronSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
  fJetSelection(config.getParameter<ParameterSet>("JetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
  fBJetSelection(config.getParameter<ParameterSet>("BJetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
  fMETSelection(config.getParameter<ParameterSet>("METSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
  cSelected(fEventCounter.addCounter("Selected events"))
{ }



Hplus2hwAnalysis_background::~Hplus2hwAnalysis_background() {
  fCommonPlots.getHistoSplitter().deleteHistograms(hMtInvertedTauAfterStdSelections);
  fCommonPlots.getHistoSplitter().deleteHistograms(hBaselineTauTransverseMass);
}

void Hplus2hwAnalysis_background::book(TDirectory *dir) {

  // Book common plots histograms
  fCommonPlots.book(dir, isData());

  // Book histograms in event selection classes
  fMETFilterSelection.bookHistograms(dir);


  fElectronSelection.bookHistograms(dir);
  fTauSelection.bookHistograms(dir);
  fLooseTauSelection.bookHistograms(dir);
  fMuonSelection.bookHistograms(dir);
  fJetSelection.bookHistograms(dir);
  fMETSelection.bookHistograms(dir);
  fBJetSelection.bookHistograms(dir);

  // Book non-common histograms



  hTauPt =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauPt", "Tau pT", 40, 0, 400);
  hMuonPt =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "muPt", "Muon pT", 40, 0, 400);
  hMET =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "MET", "MET", 40, 0, 400);

  hTauEta =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauEta", "Tau eta", 50, -2.5, 2.5);
  hMuonEta =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "muEta", "Muon eta", 50, -2.5, 2.5);

  hMuonEta_afterMuonSelection =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "muEta_afterMuonSelection", "Muon eta after muon selection", 50, -2.5, 2.5);
  hMuonPt_afterMuonSelection = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "muPt_afterMuonSelection", "Muon pT adter uon selection", 40, 0, 400);

  hNJet = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "nJet", "# of jets", 10, 0, 10);
  hNTau =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "nTau", "# of Selected taus", 10, 0, 10);
  hTauCharge = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "tauCharge", "charge of taus", 6,-3, 3);

  hTransverseMass = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TransverseMass", "TransverseMass", 200, 0, 2000);
  hTransverseMass_genuine = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, dir, "TransverseMass_genuine", "TransverseMass_genuine", 200, 0, 2000);

  return;
}


void Hplus2hwAnalysis_background::setupBranches(BranchManager& branchManager) {
  fEvent.setupBranches(branchManager);
  return;
}

void Hplus2hwAnalysis_background::process(Long64_t entry) {

  ////////////
  // Initialize
  ////////////

  fCommonPlots.initialize();

  cAllEvents.increment();


  ////////////
  // Apply Trigger
  ////////////

  if (!(fEvent.passTriggerDecision()))
    return;

  int nVertices = fEvent.vertexInfo().value();

  ////////////
  // Primarty Vertex (Check that a PV exists)
  ////////////

  if (nVertices < 1)
    return;

  ////////////
  // MET filters (to remove events with spurious sources of fake MET)
  ////////////

  const METFilterSelection::Data metFilterData = fMETFilterSelection.analyze(fEvent);
  if (!metFilterData.passedSelection()) 
    return;

  ////////////
  // Muon
  ////////////

  const MuonSelection::Data muData = fMuonSelection.analyze(fEvent);
  if(!(muData.hasIdentifiedMuons()))
    return;

  if (muData.getSelectedMuons().size() != 1)
    return;

  ////////////
  // Dummy Trigger SF for first check
  ////////////

  if (fEvent.isMC()) {
    if (26 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 30) fEventWeight.multiplyWeight(0.9664);
    if (30 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 40) fEventWeight.multiplyWeight(0.9781);
    if (40 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 50) fEventWeight.multiplyWeight(0.9819);
    if (50 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 60) fEventWeight.multiplyWeight(0.9822);
    if (60 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 80) fEventWeight.multiplyWeight(0.9804);
    if (80 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 120) fEventWeight.multiplyWeight(0.9780);
    if (120 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 200) fEventWeight.multiplyWeight(0.9752);
    if (200 <= muData.getSelectedMuons()[0].pt() && muData.getSelectedMuons()[0].pt() < 500) fEventWeight.multiplyWeight(0.9704);

  }

  hMuonPt_afterMuonSelection->Fill(muData.getSelectedMuons()[0].pt());
  hMuonEta_afterMuonSelection->Fill(muData.getSelectedMuons()[0].eta());


  ////////////
  // Tau
  ////////////

  const TauSelection::Data looseTauData = fLooseTauSelection.analyzeLoose(fEvent);
  const TauSelection::Data tauData = fTauSelection.analyze(fEvent);

  std::vector<float> myFactorisationInfo;

  if (!looseTauData.hasIdentifiedTaus())
    return;

  if(looseTauData.getSelectedTaus().size() != 2)
    return;

  if(looseTauData.getSelectedTaus()[0].charge() == looseTauData.getSelectedTaus()[1].charge())
    return;

//  if (fEvent.isMC() && looseTauData.getSelectedTaus()[0].isGenuineTau() && looseTauData.getSelectedTaus()[1].isGenuineTau()) { //if genuine tau, reject the event
//    std::cout << "two taus prompt" << "\n";
//    return;
//  }

/// -------------------------------------

  if (tauData.hasIdentifiedTaus() && tauData.getSelectedTaus().size()==2) {
    if (fEvent.isMC() && tauData.getSelectedTaus()[0].isGenuineTau() && tauData.getSelectedTaus()[1].isGenuineTau()) { //if genuine tau, reject the event
      return;
    }
  }
/*
  drMuTau1 = ROOT::Math::VectorUtil::DeltaR(muData.getSelectedMuons()[0].p4(),looseTauData.getSelectedTaus()[0].p4());
  if(drMuTau1 < 0.5) {
    return;
  }

  drMuTau2 = ROOT::Math::VectorUtil::DeltaR(muData.getSelectedMuons()[0].p4(),looseTauData.getSelectedTaus()[1].p4());
  if(drMuTau2 < 0.5) {
    return;
  }
*/

  if (tauData.getAntiIsolatedTaus().size() == 2 && tauData.getSelectedTaus().size()==0) { // both anti tight Iso

//    std::cout << "two or more anti tight iso" << "\n";

    myFactorisationInfo.push_back(tauData.getAntiIsolatedTaus()[0].pt());

    if(tauData.getAntiIsolatedTaus()[0].decayMode()==0) {
      myFactorisationInfo.push_back(1);
    }
    if(tauData.getAntiIsolatedTaus()[0].decayMode()==1) {
      myFactorisationInfo.push_back(2);
    }
    if(tauData.getAntiIsolatedTaus()[0].decayMode()==10) {
      myFactorisationInfo.push_back(3);
    }
    if(tauData.getAntiIsolatedTaus()[0].decayMode()>1 && tauData.getAntiIsolatedTaus()[0].decayMode()<10) {
      return;
    }

    myFactorisationInfo.push_back(tauData.getAntiIsolatedTaus()[1].pt());

    if(tauData.getAntiIsolatedTaus()[1].decayMode()==0) {
      myFactorisationInfo.push_back(1);
    }
    if(tauData.getAntiIsolatedTaus()[1].decayMode()==1) {
      myFactorisationInfo.push_back(2);
    }
    if(tauData.getAntiIsolatedTaus()[1].decayMode()==10) {
      myFactorisationInfo.push_back(3);
    }
    if(tauData.getAntiIsolatedTaus()[1].decayMode()>1 && tauData.getAntiIsolatedTaus()[1].decayMode()<10) {
      return;
    }


    fCommonPlots.setFactorisationBinForEvent(myFactorisationInfo);

//    cOverTwoTausCounter.increment();

    fCommonPlots.fillControlPlotsAfterTauSelection(fEvent, tauData);


  } else if (tauData.getAntiIsolatedTaus().size() >= 1 && tauData.getSelectedTaus().size()==1){  //only one tau is fake

//    std::cout << "one is anti tight iso" << "\n";

    myFactorisationInfo.push_back(tauData.getAntiIsolatedTaus()[0].pt());

    if(tauData.getAntiIsolatedTaus()[0].decayMode()==0) {
      myFactorisationInfo.push_back(1);
    }
    if(tauData.getAntiIsolatedTaus()[0].decayMode()==1) {
      myFactorisationInfo.push_back(2);
    }
    if(tauData.getAntiIsolatedTaus()[0].decayMode()==10) {
      myFactorisationInfo.push_back(3);
    }
    if(tauData.getAntiIsolatedTaus()[0].decayMode()>1 && tauData.getAntiIsolatedTaus()[0].decayMode()<10) {
      return;
    }

    myFactorisationInfo.push_back(0);
    myFactorisationInfo.push_back(0);


    fCommonPlots.setFactorisationBinForEvent(myFactorisationInfo);

    fCommonPlots.fillControlPlotsAfterTauSelection(fEvent, tauData);

  } else if (tauData.getAntiIsolatedTaus().size() == 0 || tauData.getAntiIsolatedTaus().size() >= 3) {

//    std::cout << "two or more pass tight iso" << "\n";

    return;

  } else {

//    this is reached when there are 2 anti-isolated but also identified taus
//    std::cout << "This should not be reached" << "\n";

    return;

  }//

/// -------------------------------------



  ////////////
  // Electron veto (Fully hadronic + orthogonality)
  ////////////

  const ElectronSelection::Data eData = fElectronSelection.analyze(fEvent);
  if (eData.hasIdentifiedElectrons())
    return;




  ////////////
  // Jet selection
  ////////////

//  const JetSelection::Data jetData = fJetSelection.analyze(fEvent, looseTauData.getSelectedTau());
  const JetSelection::Data jetData = fJetSelection.analyze(fEvent, looseTauData.getSelectedTaus()[0],looseTauData.getSelectedTaus()[1]);
  if (!jetData.passedSelection())
    return;


  ////////////
  // BJet selection
  ////////////

  const BJetSelection::Data bjetData = fBJetSelection.analyze(fEvent, jetData);

  if (!bjetData.passedSelection())
    return;

  ////////////
  // BJet SF
  ////////////

  if (fEvent.isMC()) {
    fEventWeight.multiplyWeight(bjetData.getBTaggingScaleFactorEventWeight());
  }
//  cBTaggingSFCounter.increment();
//  fCommonPlots.fillControlPlotsAfterBtagSF(fEvent,jetData,bjetData);


  ////////////
  // MET selection
  ////////////

  const METSelection::Data METData = fMETSelection.analyze(fEvent, nVertices);
  if (!METData.passedSelection())
    return;

//  cMETSelection.increment();

  ////////////
  // All cuts passed
  ////////////


  for(unsigned int i=0; i<2; i++){
    hTauPt->Fill(looseTauData.getSelectedTaus()[i].pt());
    hTauEta->Fill(looseTauData.getSelectedTaus()[i].eta());
    hTauCharge->Fill(looseTauData.getSelectedTaus()[i].charge());
  }

  hMuonPt->Fill(muData.getSelectedMuons()[0].pt());
  hMuonEta->Fill(muData.getSelectedMuons()[0].eta());

  hNTau->Fill(looseTauData.getSelectedTaus().size());

  hMET->Fill(METData.getMET().R());

//  hNJet->Fill(jetData.getNumberOfSelectedJets());

  double myTransverseMass = TransverseMass::reconstruct(looseTauData.getSelectedTaus()[0],looseTauData.getSelectedTaus()[1],muData.getSelectedMuons()[0], METData.getMET());


  hTransverseMass->Fill(myTransverseMass);

  if (fEvent.isMC()) {
    if (tauData.getAntiIsolatedTaus().size() == 2) {
      if (tauData.getAntiIsolatedTaus()[0].isGenuineTau() && tauData.getAntiIsolatedTaus()[1].isGenuineTau()) {
        myTransverseMass = TransverseMass::reconstruct(tauData.getAntiIsolatedTaus()[0],tauData.getAntiIsolatedTaus()[1],muData.getSelectedMuons()[0], METData.getMET());
        hTransverseMass_genuine->Fill(myTransverseMass);
      } else if (tauData.getAntiIsolatedTaus()[0].isGenuineTau()) {
        myTransverseMass = TransverseMass::reconstruct(looseTauData.getSelectedTaus()[0],tauData.getAntiIsolatedTaus()[0],muData.getSelectedMuons()[0], METData.getMET());
        hTransverseMass_genuine->Fill(myTransverseMass);
      } else if (tauData.getAntiIsolatedTaus()[1].isGenuineTau()) {
        myTransverseMass = TransverseMass::reconstruct(looseTauData.getSelectedTaus()[0],tauData.getAntiIsolatedTaus()[1],muData.getSelectedMuons()[0], METData.getMET());
        hTransverseMass_genuine->Fill(myTransverseMass);
      }
    } else if (tauData.getAntiIsolatedTaus()[0].isGenuineTau()) {
      myTransverseMass = TransverseMass::reconstruct(looseTauData.getSelectedTaus()[0],tauData.getAntiIsolatedTaus()[0],muData.getSelectedMuons()[0], METData.getMET());
      hTransverseMass_genuine->Fill(myTransverseMass);
    }
  }

  cSelected.increment();

  ////////////
  // Fill final plots
  ////////////

//  std::cout << "DEBUG: before all selections";

  fCommonPlots.fillControlPlotsAfterAllSelections(fEvent);
//  fCommonPlots.getHistoSplitter().fillShapeHistogramTriplet(hNormalizationInvertedTauAfterStdSelections, isGenuineTau, METvalue);
//  fCommonPlots.getHistoSplitter().fillShapeHistogramTriplet(hMtInvertedTauAfterStdSelections, tauData.getAntiIsolatedTauIsGenuineTau(), myTransverseMass);

  ////////////
  // Finalize
  ////////////

//  std::cout << "DEBUG: after all selections";

  fEventSaver.save();


}

void Hplus2hwAnalysis_background::doBaselineAnalysis(const Event& event, const TauSelection::Data& tauData, const int nVertices) {



}


void Hplus2hwAnalysis_background::doSignalAnalysis(const Event& event, const TauSelection::Data& tauData, const int nVertices) {



}
