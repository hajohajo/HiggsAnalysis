// -*- c++ -*-
#include "Framework/interface/BaseSelector.h"
#include "Framework/interface/makeTH.h"

#include "EventSelection/interface/CommonPlots.h"
#include "EventSelection/interface/EventSelections.h"

#include "TDirectory.h"

class Hplus2hwAnalysis: public BaseSelector {
public:
  explicit Hplus2hwAnalysis(const ParameterSet& config, const TH1* skimCounters);
  virtual ~Hplus2hwAnalysis() {}

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
  Count cTauOSCounter;
  Count cTauSFCounter;
  Count cFakeTauSFCounter;
  // Count cTauDMCounter;
  JetSelection fJetSelection;
  BJetSelection fBJetSelection;
  Count cBTaggingSFCounter;
  METSelection fMETSelection;
  Count cTopCleaningCounter;
  Count cTopTaggingSFCounter;
  TopSelectionBDT fTopSelection;
  // FatJetSelection fFatJetSelection;
  Count cSelected;
    
  // Non-common histograms
  WrappedTH1 *hMuon_Pt;
  WrappedTH1 *hMuon_Eta;

};

#include "Framework/interface/SelectorFactory.h"
REGISTER_SELECTOR(Hplus2hwAnalysis);

Hplus2hwAnalysis::Hplus2hwAnalysis(const ParameterSet& config, const TH1* skimCounters)
  : BaseSelector(config, skimCounters),
    fCommonPlots(config.getParameter<ParameterSet>("CommonPlots"), CommonPlots::kHplus2tbAnalysis, fHistoWrapper), // CommonPlots::kHplus2hwAnalysis
    cAllEvents(fEventCounter.addCounter("all events")),
    cTrigger(fEventCounter.addCounter("passed trigger")),
    fMETFilterSelection(config.getParameter<ParameterSet>("METFilter"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cVertexSelection(fEventCounter.addCounter("passed PV")),
    fElectronSelection(config.getParameter<ParameterSet>("ElectronSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    fMuonSelection(config.getParameter<ParameterSet>("MuonSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fTauSelection(config.getParameter<ParameterSet>("TauSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cTauNCounter(fEventCounter.addCounter("#tau N")),
    cTauOSCounter(fEventCounter.addCounter("#tau OS")),
    cTauSFCounter(fEventCounter.addCounter("#tau SF")),
    cFakeTauSFCounter(fEventCounter.addCounter("Fake #tau SF")),
    // cTauDMCounter(fEventCounter.addCounter("#tau DM")),
    fJetSelection(config.getParameter<ParameterSet>("JetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fBJetSelection(config.getParameter<ParameterSet>("BJetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cBTaggingSFCounter(fEventCounter.addCounter("b-tag SF")),
    fMETSelection(config.getParameter<ParameterSet>("METSelection")), // no subcounter in main counter
    cTopCleaningCounter(fEventCounter.addCounter("top cleaning")),
    cTopTaggingSFCounter(fEventCounter.addCounter("top-tag SF")),
    fTopSelection(config.getParameter<ParameterSet>("TopSelectionBDT"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    // fFatJetSelection(config.getParameter<ParameterSet>("FatJetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    cSelected(fEventCounter.addCounter("Selected Events"))
{ }


void Hplus2hwAnalysis::book(TDirectory *dir) {

  if (0) std::cout << "=== Hplus2hwAnalysis::book()" << std::endl;
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
  // fFatJetSelection.bookHistograms(dir);

  // Get binning from cfg file
  const int  nPtBins      = fCommonPlots.getPtBinSettings().bins();
  const float fPtMin      = fCommonPlots.getPtBinSettings().min();
  const float fPtMax      = fCommonPlots.getPtBinSettings().max();

  const int  nEtaBins     = fCommonPlots.getEtaBinSettings().bins();
  const float fEtaMin     = fCommonPlots.getEtaBinSettings().min();
  const float fEtaMax     = fCommonPlots.getEtaBinSettings().max();

  // Book non-common histograms
  hMuon_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "Muon_Pt" , ";p_{T} (GeV/c)", nPtBins , fPtMin , fPtMax);
  hMuon_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, dir, "Muon_Eta", ";#eta"         , nEtaBins, fEtaMin, fEtaMax);

  return;
}


void Hplus2hwAnalysis::setupBranches(BranchManager& branchManager) {
  fEvent.setupBranches(branchManager);
  return;
}


void Hplus2hwAnalysis::process(Long64_t entry) {

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
  // 4) Electron veto
  //================================================================================================   
  if (0) std::cout << "=== Electron veto" << std::endl;
  const ElectronSelection::Data eData = fElectronSelection.analyze(fEvent);
  if (eData.hasIdentifiedElectrons()) return;


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
  hMuon_Pt ->Fill(muData.getSelectedMuons()[0].pt());
  hMuon_Eta->Fill(muData.getSelectedMuons()[0].eta());  

  //================================================================================================   
  // 6) Tau Selection
  //================================================================================================   
  if (0) std::cout << "=== Tau Selection" << std::endl;
  const TauSelection::Data tauData = fTauSelection.analyze(fEvent);
  if (!tauData.hasIdentifiedTaus() ) return;
  if (0) std::cout << "=== Tau Selection: Has Identified taus" << std::endl;
    
  // Require 2 taus
  if(tauData.getSelectedTaus().size() != 2) return;
  cTauNCounter.increment(); 
  if (0) std::cout << "=== Tau Selection: 2 selected taus" << std::endl;

  // Require 2 taus to have opposite sign (OS)
  if(tauData.getSelectedTaus()[0].charge() == tauData.getSelectedTaus()[1].charge()) return;
  cTauOSCounter.increment(); 
  if (0) std::cout << "=== Tau Selection: OS requirement" << std::endl;

  // Uncomment to only perform Genuine #tau wit analysis (with fake #tau from data)
  // if (fEvent.isMC() && !tauData.getSelectedTaus()[0].isGenuineTau()) return;
  // if (fEvent.isMC() && !tauData.getSelectedTaus()[1].isGenuineTau()) return;
  int isGenuineTau = false; 
  if (fEvent.isMC()) isGenuineTau = tauData.getSelectedTaus()[0].isGenuineTau() && tauData.getSelectedTaus()[1].isGenuineTau();
  if (0) std::cout << "=== Tau Selection: isGenuineTau = " << isGenuineTau << std::endl;

  // Fake rates currently only available for three decay modes: 0, 1, 10 (only needed for the fake rate)
  // if(tauData.getSelectedTaus()[0].decayMode()>1 && tauData.getSelectedTaus()[0].decayMode()<10) return;
  // if(tauData.getSelectedTaus()[1].decayMode()>1 && tauData.getSelectedTaus()[1].decayMode()<10) return;
  // cTauDMCounter.increment(); 

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
  fCommonPlots.fillControlPlotsAfterTauSelection(fEvent, tauData);
  if (0) std::cout << "=== Tau Selection: SF applied" << std::endl;


  //================================================================================================
  // 7) Jet selection
  //================================================================================================
  if (0) std::cout << "=== Jet selection" << std::endl;
  // const JetSelection::Data jetData = fJetSelection.analyzeWithoutTau(fEvent);
  const JetSelection::Data jetData = fJetSelection.analyze(fEvent, tauData.getSelectedTau()); // fixme (only 1 tau is returned)
  if (!jetData.passedSelection()) return;
  fCommonPlots.fillControlPlotsAfterTopologicalSelections(fEvent, true);

 
  //================================================================================================  
  // 8) BJet selection
  //================================================================================================
  if (0) std::cout << "=== BJet selection" << std::endl;
  const BJetSelection::Data bjetData = fBJetSelection.analyze(fEvent, jetData);
  if (!bjetData.passedSelection()) return;
  // fCommonPlots.fillControlPlotsAfterBJetSelection(fEvent, bjetData);


  //================================================================================================  
  // 9) BJet SF  
  //================================================================================================
  if (0) std::cout << "=== BJet SF" << std::endl;
  if (fEvent.isMC()) 
    {
      fEventWeight.multiplyWeight(bjetData.getBTaggingScaleFactorEventWeight());
    }
  cBTaggingSFCounter.increment();
  // int isGenuineB = bjetData.isGenuineB();


  //================================================================================================
  // - MET selection
  //================================================================================================
  if (0) std::cout << "=== MET selection" << std::endl;
  const METSelection::Data METData = fMETSelection.analyze(fEvent, nVertices);
  // const METSelection::Data METData = fMETSelection.silentAnalyze(fEvent, nVertices);
  if (!METData.passedSelection()) return;
  // fixme - counter?


  //================================================================================================
  // 10) Top selection
  //================================================================================================
  if (0) std::cout << "=== Top (BDT) selection" << std::endl;
  const TopSelectionBDT::Data topData = fTopSelection.analyze(fEvent, jetData, bjetData);
  if (topData.getAllCleanedTopsSize() < 1) return;
  cTopCleaningCounter.increment();

  // Apply top-tag SF
  if (fEvent.isMC()) 
    {
      fEventWeight.multiplyWeight(topData.getTopTaggingScaleFactorEventWeight());
    }
  cTopTaggingSFCounter.increment();

  // Fill histos after StandardSelections: Require any two tops with BDT > -1.0 and presence of free b-jet (not taken up by any of the two best (in MVA) tops)
  fCommonPlots.fillControlPlotsAfterStandardSelections(fEvent, jetData, bjetData, METData, QuarkGluonLikelihoodRatio::Data(), topData, isGenuineTau);
  if (!topData.passedSelection()) return;  
  if (0) std::cout << "\nentry = " << entry << ", topData.getTopMVA() = " << topData.getTopMVA() << std::endl;

  //================================================================================================
  // All Selections
  // https://github.com/mlotti/HplusHW/blob/master/NtupleAnalysis/src/Hplus2hwAnalysis/src/Hplus2hwAnalysis.cc#L339
  //================================================================================================
  if (0) std::cout << "=== All Selections" << std::endl;
  cSelected.increment();

  // Fill histos after AllSelections: (After top-selections and top-tag SF)
  fCommonPlots.fillControlPlotsAfterAllSelections(fEvent, isGenuineTau);
  fEventSaver.save();

  return;
}
