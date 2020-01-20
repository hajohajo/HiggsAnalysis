// -*- c++ -*-
#include "EventSelection/interface/ElectronSelection.h"

#include "Framework/interface/ParameterSet.h"
#include "EventSelection/interface/CommonPlots.h"
#include "EventSelection/interface/CommonPlots_ttm.h"
#include "DataFormat/interface/Electron.h"
#include "DataFormat/interface/Event.h"
#include "Framework/interface/HistoWrapper.h"
#include "DataFormat/interface/Electron.h"
//#include "Framework/interface/makeTH.h"

ElectronSelection::Data::Data()
: fHighestSelectedElectronPt(0.0),
  fHighestSelectedElectronEta(0.0),
  fElectronIDSF(1.0),
  fElectronTriggerSF(1.0),
  fHLTElectronCharge(1) { }


ElectronSelection::Data::~Data() { }

ElectronSelection::ElectronSelection(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& postfix)
: BaseSelection(eventCounter, histoWrapper, commonPlots, postfix),
  cfg_ElectronPtCut(config.getParameter<float>("electronPtCut")),
  cfg_ElectronEtaCut(config.getParameter<float>("electronEtaCut")),
  cfg_ElectronMVACut(config.getParameter<string>("electronMVACut")),
  fRelIsoCut(-1.0),
  fMiniIsoCut(-1.0),
  fVetoMode(false),
  fMiniIsol(false),
  fElectronMVA(false),
  fElectronIDSFReader(config.getParameterOptional<ParameterSet>("electronIDSF")),
  fElectronTriggerSFReader(config.getParameterOptional<ParameterSet>("electronTriggerSF")),
  // Event counter for passing selection
  cPassedElectronSelection(fEventCounter.addCounter("passed e selection ("+postfix+")")),
  // Sub counters
  cSubAll(fEventCounter.addSubCounter("e selection ("+postfix+")", "All events")),
  cSubPassedIsPresent(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed is present")),
  cSubPassedTriggerMatching(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed trigger matching")),
  cSubPassedPt(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed pt cut")),
  cSubPassedEta(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed eta cut")),
  cSubPassedID(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed ID")),
  cSubPassedIsolation(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed isolation")),
  cSubPassedSelection(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed selection")),
  cSubPassedVeto(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed veto"))
{
  initialize(config, postfix);
}

ElectronSelection::ElectronSelection(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots_ttm* commonPlots, const std::string& postfix)
: BaseSelection(eventCounter, histoWrapper, commonPlots, postfix),
  cfg_ElectronPtCut(config.getParameter<float>("electronPtCut")),
  cfg_ElectronEtaCut(config.getParameter<float>("electronEtaCut")),
  cfg_ElectronMVACut(config.getParameter<string>("electronMVACut")),
  fRelIsoCut(-1.0),
  fMiniIsoCut(-1.0),
  fVetoMode(false),
  fMiniIsol(false),
  fElectronMVA(false),
  fElectronIDSFReader(config.getParameterOptional<ParameterSet>("electronIDSF")),
  fElectronTriggerSFReader(config.getParameterOptional<ParameterSet>("electronTriggerSF")),
  // Event counter for passing selection
  cPassedElectronSelection(fEventCounter.addCounter("passed e selection ("+postfix+")")),
  // Sub counters
  cSubAll(fEventCounter.addSubCounter("e selection ("+postfix+")", "All events")),
  cSubPassedIsPresent(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed is present")),
  cSubPassedTriggerMatching(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed trigger matching")),
  cSubPassedPt(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed pt cut")),
  cSubPassedEta(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed eta cut")),
  cSubPassedID(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed ID")),
  cSubPassedIsolation(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed isolation")),
  cSubPassedSelection(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed selection")),
  cSubPassedVeto(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed veto"))
{
  initialize(config, postfix);
}

ElectronSelection::ElectronSelection(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, std::nullptr_t, const std::string& postfix)
: BaseSelection(eventCounter, histoWrapper, nullptr, postfix),
  cfg_ElectronPtCut(config.getParameter<float>("electronPtCut")),
  cfg_ElectronEtaCut(config.getParameter<float>("electronEtaCut")),
  cfg_ElectronMVACut(config.getParameter<string>("electronMVACut")),
  fRelIsoCut(-1.0),
  fMiniIsoCut(-1.0),
  fVetoMode(false),
  fMiniIsol(false),
  fElectronMVA(false),
  fElectronIDSFReader(config.getParameterOptional<ParameterSet>("electronIDSF")),
  fElectronTriggerSFReader(config.getParameterOptional<ParameterSet>("electronTriggerSF")),
  // Event counter for passing selection
  cPassedElectronSelection(fEventCounter.addCounter("passed e selection ("+postfix+")")),
  // Sub counters
  cSubAll(fEventCounter.addSubCounter("e selection ("+postfix+")", "All events")),
  cSubPassedIsPresent(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed is present")),
  cSubPassedTriggerMatching(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed trigger matching")),
  cSubPassedPt(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed pt cut")),
  cSubPassedEta(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed eta cut")),
  cSubPassedID(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed ID")),
  cSubPassedIsolation(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed isolation")),
  cSubPassedSelection(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed selection")),
  cSubPassedVeto(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed veto"))
{
  initialize(config, postfix);
}


ElectronSelection::ElectronSelection(const ParameterSet& config, const std::string& postfix)
: BaseSelection(),
  cfg_ElectronPtCut(config.getParameter<float>("electronPtCut")),
  cfg_ElectronEtaCut(config.getParameter<float>("electronEtaCut")),
  cfg_ElectronMVACut(config.getParameter<string>("electronMVACut")),
  fRelIsoCut(-1.0),
  fMiniIsoCut(-1.0),
  fVetoMode(false),
  fMiniIsol(false),
  fElectronMVA(false),
  fElectronIDSFReader(config.getParameterOptional<ParameterSet>("electronIDSF")),
  fElectronTriggerSFReader(config.getParameterOptional<ParameterSet>("electronTriggerSF")),
  // Event counter for passing selection
  cPassedElectronSelection(fEventCounter.addCounter("passed e selection ("+postfix+")")),
  // Sub counters
  cSubAll(fEventCounter.addSubCounter("e selection ("+postfix+")", "All events")),
  cSubPassedIsPresent(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed is present")),
  cSubPassedTriggerMatching(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed trigger matching")),
  cSubPassedPt(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed pt cut")),
  cSubPassedEta(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed eta cut")),
  cSubPassedID(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed ID")),
  cSubPassedIsolation(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed isolation")),
  cSubPassedSelection(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed selection")),
  cSubPassedVeto(fEventCounter.addSubCounter("e selection ("+postfix+")", "Passed veto"))
{
  initialize(config, postfix);
  bookHistograms(new TDirectory());
}

ElectronSelection::~ElectronSelection() {
  delete hTriggerMatchDeltaR;
  delete hElectronNAll;
  delete hElectronPtAll;
  delete hElectronEtaAll;
  delete hElectronRelIsoAll;
  delete hElectronMiniIsoAll;

  delete hElectronNPassed;
  delete hElectronPtPassed;
  delete hElectronEtaPassed;
  delete hElectronRelIsoPassed;
  delete hElectronMiniIsoPassed;

  delete hPtResolution;
  delete hEtaResolution;
  delete hPhiResolution;

  delete hIsolPtBefore;
  delete hIsolEtaBefore;
  delete hIsolVtxBefore;
  delete hIsolRelIsoBefore;
  delete hIsolMiniIsoBefore;

  delete hIsolPtAfter;
  delete hIsolEtaAfter;
  delete hIsolVtxAfter;
  delete hIsolRelIsoAfter;
  delete hIsolMiniIsoAfter;
}

void ElectronSelection::initialize(const ParameterSet& config, const std::string& postfix) {
  if(config.getParameterOptional<bool>("applyTriggerMatching")) cfg_ApplyTriggerMatching = config.getParameter<bool>("applyTriggerMatching");
  else cfg_ApplyTriggerMatching = false;
  if(config.getParameterOptional<float>("triggerMatchingCone")) cfg_TriggerMatchingCone = config.getParameter<float>("triggerMatchingCone");

  if (postfix.find("veto") != std::string::npos || postfix.find("Veto") != std::string::npos)
    {
    fVetoMode = true;
    }
 
 std::string isolString = config.getParameter<std::string>("electronIsolation");
  if (isolString == "veto" || isolString == "Veto") {
    fRelIsoCut  = 0.15; // Loose iso sync'ed with MIT
    fMiniIsoCut = 0.4;  // from Brown/MIT sync
  } 
  else if (isolString == "tight" || isolString == "Tight") {
    fRelIsoCut  = 0.10; // Based on 2012 cut based isolation
    fMiniIsoCut = 0.10; // arbitrary value selected
  } 
  else
    {
      throw hplus::Exception("config") << "Invalid electronIsolation option '" << isolString << "'! Options: 'veto', 'tight'";
    } 

  std::string isolTypeString = config.getParameter<std::string>("electronIsolType");
  if (isolTypeString == "default")  fMiniIsol = false;
  else if (isolTypeString == "mini") fMiniIsol = true;
  else
   {
     throw hplus::Exception("config") << "Invalid electronIsolType option '" << isolTypeString << "'! Options: 'default', 'mini'";
   }
  
  std::string idTypeString = config.getParameter<std::string>("electronIDType");
  if (idTypeString == "default") fElectronMVA = false;
  else if (idTypeString == "MVA") fElectronMVA = true;
  else
    {
      throw hplus::Exception("config") << "Invalid electronIDType option '" << idTypeString << "'! Options: 'default', 'MVA'";
    }
  
  return;
}

void ElectronSelection::bookHistograms(TDirectory* dir) {
  TDirectory* subdir  = fHistoWrapper.mkdir(HistoLevel::kDebug, dir, "eSelection_"+sPostfix);

  // Electrons before any cuts
  hElectronNAll       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "electronNAll", ";e multiplicity;Occur / %.0f", 20, 0, 20.0);
  hElectronPtAll      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "electronPtAll", ";p_{T} (GeV/c);Occur / %.0f", 100, 0, 1000.0);
  hElectronEtaAll     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "electronEtaAll", ";#eta;Occur / %.2f", 50, -2.5, 2.5);
  hElectronRelIsoAll  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "electronRelIsoAll" , ";relative isolation;Occur / %.0f", 1000, 0.0, 200.0);
  hElectronMiniIsoAll = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "electronMiniIsoAll", ";relative mini-isolation;Occur / %.0f", 1000, 0.0, 200.0);
 
  // Electrons after all cuts
  hElectronNPassed       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "electronNPassed", ";e multiplicity;Occur / %.0f", 20, 0, 20.0);
  hElectronPtPassed      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "electronPtPassed", ";p_{T} (GeV/c);Occur / %.0f", 100, 0.0, 1000.0);
  hElectronEtaPassed     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "electronEtaPassed", ";#eta;Occur / %.2f", 50, -2.5, 2.5);
  hElectronRelIsoPassed  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "electronRelIsoPassed", ";relative isolation;Occur / %.2f", 1000, 0.0, 200.0);
  hElectronMiniIsoPassed = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "electronMiniIsoPassed", ";relative mini-isolation;Occur / %.2f", 1000, 0.0, 200.0);

  hTriggerMatchDeltaR   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "triggerMatchDeltaR"  , "Trigger match #DeltaR;#DeltaR", 60, 0, 3.);

  // Resolutions
  hPtResolution  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "ptResolution" , ";(p_{T}^{reco} - p_{T}^{gen})/p_{T}^{reco};Occur / %.2f", 400, -2.0, 2.0);
  hEtaResolution = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "etaResolution", ";(#eta^{reco} - #eta^{gen})/#eta^{reco};Occur / %.2f"   , 400, -2.0, 2.0);
  hPhiResolution = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "phiResolution", ";(#phi^{reco} - #phi^{gen})/#phi^{reco};Occur / %.2f"   , 400, -2.0, 2.0);

  // Isolation efficiency
  hIsolPtBefore      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "IsolPtBefore", ";p_{T} (GeV/c);Occur / %.0f", 100, 0.0, 1000.0);
  hIsolEtaBefore     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "IsolEtaBefore", ";#eta;Occur / %.2f", 50, -2.5, 2.5);
  hIsolVtxBefore     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "IsolVtxBefore", ";Number of Vertices;Occur / %.2f", 150, 0, 150.0);
  hIsolRelIsoBefore  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "IsolRelIsoBefore", ";relative isolation;Occur / %.2f", 1000, 0.0, 200.0);
  hIsolMiniIsoBefore = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "IsolMiniIsoBefore", ";relative mini-isolation;Occur / %.2f", 1000, 0.0, 200.0);

  hIsolPtAfter      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "IsolPtAfter", ";p_{T} (GeV/c);Occur / %.0f", 50, 0.0, 500.0);
  hIsolEtaAfter     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "IsolEtaAfter", ";#eta;Occur / %.2f", 50, -2.5, 2.5);
  hIsolVtxAfter     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "IsolVtxAfter", ";Number of Vertices;Occur / %.0f", 150, 0, 150);
  hIsolRelIsoAfter  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "IsolRelIsoAfter", ";relative isolation;Occur / %.2f", 1000, 0.0, 200.0);
  hIsolMiniIsoAfter = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "IsolMiniIsoAfter", ";relative mini-isolation;Occur / %.2f", 1000, 0.0, 200.0);
}

ElectronSelection::Data ElectronSelection::silentAnalyze(const Event& event) {
  ensureSilentAnalyzeAllowed(event.eventID());
  // Disable histogram filling and counter
  disableHistogramsAndCounters();
  Data myData = privateAnalyze(event);
  enableHistogramsAndCounters();
  return myData;
}

ElectronSelection::Data ElectronSelection::analyze(const Event& event) {
  ensureAnalyzeAllowed(event.eventID());
  ElectronSelection::Data data = privateAnalyze(event);
  // Send data to CommonPlots
  if (fCommonPlotsIsEnabled())
    {
      fCommonPlots->fillControlPlotsAtElectronSelection(event, data);
    }
  if (fCommonPlotsIsEnabled_ttm())
    {
      fCommonPlots_ttm->fillControlPlotsAtElectronSelection(event, data);
    }
  // Return data
  return data;
}

ElectronSelection::Data ElectronSelection::analyzeLoose(const Event& event) {
  ensureAnalyzeAllowed(event.eventID());
  ElectronSelection::Data data = privateAnalyzeLoose(event);
  // Send data to CommonPlots
//  if (fCommonPlotsIsEnabled())
//    {
//      fCommonPlots->fillControlPlotsAtElectronSelection(event, data);
//    }

  // Return data
  return data;
}


ElectronSelection::Data ElectronSelection::privateAnalyze(const Event& event) {
  Data output;
  cSubAll.increment();
  bool passedIsPresent = false;
  bool passedTrgMatch  = false;
  bool passedPt        = false;
  bool passedEta       = false;
  bool passedID        = false;
  bool passedIsol      = false;
  bool passedSelection = false;
  bool passedVeto      = false;

  // Cache vector of trigger ele 4-momenta
  std::vector<math::LorentzVectorT<double>> myTriggerElectronMomenta;
  if (cfg_ApplyTriggerMatching) 
    {
      // For-loop: All trigger electrons
      for (HLTElectron p: event.triggerElectrons())
	{
	  myTriggerElectronMomenta.push_back(p.p4());
	}
    }
  

  // For-loop: All muons to find trg match
  if (cfg_ApplyTriggerMatching) {
    for(Electron electron: event.electrons()) {
      // Apply trigger matching
      if (this->passTrgMatching(electron, myTriggerElectronMomenta)) {
        passedTrgMatch = true;
//        output.fHLTElectronCharge = electron.charge();
      }
    }
  }

  // For-loop: All electrons
  for(Electron electron: event.electrons()) 
    {
      passedIsPresent = true;
      
      // Apply trigger matching
//      if (cfg_ApplyTriggerMatching)
//	{
//	  if (!this->passTrgMatching(electron, myTriggerElectronMomenta)) continue;
//	}
      // Designate as trigger-matched
//      output.fHLTElectronCharge = electron.charge();
      if (cfg_ApplyTriggerMatching) {
        if(!passedTrgMatch) continue;
      }

//      passedTrgMatch = true;
      
      // Fill histograms before any cuts
      hElectronPtAll->Fill(electron.pt());
      hElectronEtaAll->Fill(electron.eta());
      hElectronRelIsoAll->Fill(electron.effAreaIsoDeltaBeta()); // electron.relIsoDeltaBeta()
      hElectronMiniIsoAll->Fill(electron.effAreaMiniIso());
      
      // Debug?
//      if (0) cout << "pt = " << electron.pt() << ", eta = " << electron.eta() << ", cut-ID = " << electron.electronIDDiscriminator() << endl;

      //=== Apply cut on pt    
      if (electron.pt() < cfg_ElectronPtCut) continue;
      passedPt = true;
      
      //=== Apply cut on eta
      if (std::fabs(electron.eta()) > cfg_ElectronEtaCut) continue;
      passedEta = true;
      
      // Determine if Cut-based ID passed
      bool passedCutBasedID = electron.electronIDDiscriminator();
      bool passedMVA        = false;
      bool passedIDCut      = false;
      // Determine if MVA ID passed
      if (fElectronMVA) passedMVA   = getMVADecision(electron, cfg_ElectronMVACut); 
      if (fElectronMVA) passedIDCut = passedMVA;
      else passedIDCut = passedCutBasedID;
      //=== Apply cut on ID (Cut-based or MVA)
      if (!passedIDCut) continue;
      passedID = true;
      
      // Fill histograms before isolation cut
      hIsolPtBefore->Fill(electron.pt());
      hIsolEtaBefore->Fill(electron.eta());
      hIsolRelIsoBefore->Fill(electron.effAreaIsoDeltaBeta());
      hIsolMiniIsoBefore->Fill(electron.effAreaMiniIso());
      if (fCommonPlotsIsEnabled())
	{
	  hIsolVtxBefore->Fill(fCommonPlots->nVertices());
	}
      
      // Determine Relative and Mini Isolation booleans
      bool passedRelIso  = (electron.effAreaIsoDeltaBeta() < fRelIsoCut);
      bool passedMiniIso = (electron.effAreaMiniIso() < fMiniIsoCut);
      bool passedIsolCut = false;
      if (fMiniIsol) passedIsolCut =  passedMiniIso;
      else passedIsolCut =  passedRelIso;

      //=== Apply cut on electron isolation
//      if (!passedIsolCut) continue;
//      passedIsol = true;
      if (!passedIsolCut) {
//      std::cout << "DEGUG: anti iso muon" << "\n";
        output.fAntiIsolatedElectrons.push_back(electron);
        passedIsol = false;
      } else {
        passedIsol = true;
        output.fSelectedElectrons.push_back(electron);
      }

      // Fill histograms after isolation cut
      hIsolPtAfter->Fill(electron.pt());
      hIsolEtaAfter->Fill(electron.eta());
      hIsolRelIsoAfter->Fill(electron.effAreaIsoDeltaBeta());
      hIsolMiniIsoAfter->Fill(electron.effAreaMiniIso());
      if (fCommonPlotsIsEnabled()) 
	{
	  hIsolVtxAfter->Fill(fCommonPlots->nVertices());
	}
      
      // Fill histograms after all cuts
      hElectronPtPassed->Fill(electron.pt());
      hElectronEtaPassed->Fill(electron.eta());
      hElectronRelIsoPassed->Fill(electron.effAreaIsoDeltaBeta()); // electron.relIsoDeltaBeta()
      hElectronMiniIsoPassed->Fill(electron.effAreaMiniIso());
      
      // Save the highest pt electron
      if (electron.pt() > output.fHighestSelectedElectronPt) 
	{
	  output.fHighestSelectedElectronPt = electron.pt();
	  output.fHighestSelectedElectronEta = electron.eta();
	}
      
      // Save all electrons surviving the cuts
//      output.fSelectedElectrons.push_back(electron);
      
      // Fill resolution histograms
      if (event.isMC()) 
	{
	  hPtResolution->Fill((electron.pt() - electron.MCelectron()->pt()) / electron.pt());
	  hEtaResolution->Fill((electron.eta() - electron.MCelectron()->eta()) / electron.eta());
	  hPhiResolution->Fill((electron.phi() - electron.MCelectron()->phi()) / electron.phi());
	}
      
    }//for-loop: electrons
  
  
  //sort electrons, needed comparisons defined in Electron.h
  std::sort(output.fSelectedElectrons.begin(), output.fSelectedElectrons.end());
  std::sort(output.fAntiIsolatedElectrons.begin(), output.fAntiIsolatedElectrons.end());

  // Assign booleans
  passedSelection = (output.fSelectedElectrons.size() > 0);
  passedVeto      = (output.fSelectedElectrons.size() == 0); 


  // Set electron ID SF value to data object
  if (event.isMC()) {
    if (output.hasIdentifiedElectrons()) {
      output.fElectronIDSF = fElectronIDSFReader.getScaleFactorValue(output.getSelectedElectrons()[0].eta());
    }
  }

  // Set electron trigger SF value to data object
  if (event.isMC()) {
    if (output.hasIdentifiedElectrons()) {
      output.fElectronTriggerSF = fElectronTriggerSFReader.getScaleFactorValue(output.getSelectedElectrons()[0].pt());
    }
  }

  // Fill histos
  hElectronNAll->Fill(event.electrons().size());
  hElectronNPassed->Fill(output.fSelectedElectrons.size());

  // Fill sub-counters
  if (passedIsPresent) cSubPassedIsPresent.increment();
  if (passedTrgMatch) cSubPassedTriggerMatching.increment();
  if (passedPt) cSubPassedPt.increment();
  if (passedEta) cSubPassedEta.increment();
  if (passedID) cSubPassedID.increment();
  if (passedIsol) cSubPassedIsolation.increment();
  if (passedSelection) cSubPassedSelection.increment();
  if (passedVeto) cSubPassedVeto.increment();
  if (fVetoMode) 
    {
      if (passedVeto) cPassedElectronSelection.increment();
    }
  else
    {
      if (passedSelection) cPassedElectronSelection.increment();
    }

  // Return data object
  return output;
}

ElectronSelection::Data ElectronSelection::privateAnalyzeLoose(const Event& event) {
  Data output;
  cSubAll.increment();
  bool passedIsPresent = false;
  bool passedTrgMatch  = false;
  bool passedPt        = false;
  bool passedEta       = false;
  bool passedID        = false;
  bool passedIsol      = false;
  bool passedSelection = false;
  bool passedVeto      = false;

  // Cache vector of trigger ele 4-momenta
  std::vector<math::LorentzVectorT<double>> myTriggerElectronMomenta;
//  if (cfg_ApplyTriggerMatching) 
//    {
      // For-loop: All trigger electrons
//      for (HLTElectron p: event.triggerElectrons())
//	{
//	  myTriggerElectronMomenta.push_back(p.p4());
//	}
//    }
  
  
  // For-loop: All electrons
  for(Electron electron: event.electrons()) 
    {
      passedIsPresent = true;
      
      // Apply trigger matching
//      if (cfg_ApplyTriggerMatching)
//	{
//	  if (!this->passTrgMatching(electron, myTriggerElectronMomenta)) continue;
//	}
      // Designate as trigger-matched
//      output.fHLTElectronCharge = electron.charge();
      passedTrgMatch = true;
      
      // Fill histograms before any cuts
      hElectronPtAll->Fill(electron.pt());
      hElectronEtaAll->Fill(electron.eta());
      hElectronRelIsoAll->Fill(electron.effAreaIsoDeltaBeta()); // electron.relIsoDeltaBeta()
      hElectronMiniIsoAll->Fill(electron.effAreaMiniIso());
      
      // Debug?
//      if (0) cout << "pt = " << electron.pt() << ", eta = " << electron.eta() << ", cut-ID = " << electron.electronIDDiscriminator() << endl;

      //=== Apply cut on pt    
      if (electron.pt() < cfg_ElectronPtCut) continue;
      passedPt = true;
      
      //=== Apply cut on eta
      if (std::fabs(electron.eta()) > cfg_ElectronEtaCut) continue;
      passedEta = true;
      
      // Determine if Cut-based ID passed
      bool passedCutBasedID = electron.electronIDDiscriminator();
      bool passedMVA        = false;
      bool passedIDCut      = false;
      // Determine if MVA ID passed
      if (fElectronMVA) passedMVA   = getMVADecision(electron, cfg_ElectronMVACut); 
      if (fElectronMVA) passedIDCut = passedMVA;
      else passedIDCut = passedCutBasedID;
      //=== Apply cut on ID (Cut-based or MVA)
      if (!passedIDCut) continue;
      passedID = true;
      
      // Fill histograms before isolation cut
      hIsolPtBefore->Fill(electron.pt());
      hIsolEtaBefore->Fill(electron.eta());
      hIsolRelIsoBefore->Fill(electron.effAreaIsoDeltaBeta());
      hIsolMiniIsoBefore->Fill(electron.effAreaMiniIso());
      if (fCommonPlotsIsEnabled())
	{
	  hIsolVtxBefore->Fill(fCommonPlots->nVertices());
	}
      
      // Determine Relative and Mini Isolation booleans
      bool passedRelIso  = (electron.effAreaIsoDeltaBeta() < fRelIsoCut);
      bool passedMiniIso = (electron.effAreaMiniIso() < fMiniIsoCut);
      bool passedIsolCut = false;
      if (fMiniIsol) passedIsolCut =  passedMiniIso;
      else passedIsolCut =  passedRelIso;

      //=== Apply cut on electron isolation
//      if (!passedIsolCut) continue;
      passedIsol = true;
      
      // Fill histograms after isolation cut
      hIsolPtAfter->Fill(electron.pt());
      hIsolEtaAfter->Fill(electron.eta());
      hIsolRelIsoAfter->Fill(electron.effAreaIsoDeltaBeta());
      hIsolMiniIsoAfter->Fill(electron.effAreaMiniIso());
      if (fCommonPlotsIsEnabled()) 
	{
	  hIsolVtxAfter->Fill(fCommonPlots->nVertices());
	}
      
      // Fill histograms after all cuts
      hElectronPtPassed->Fill(electron.pt());
      hElectronEtaPassed->Fill(electron.eta());
      hElectronRelIsoPassed->Fill(electron.effAreaIsoDeltaBeta()); // electron.relIsoDeltaBeta()
      hElectronMiniIsoPassed->Fill(electron.effAreaMiniIso());
      
      // Save the highest pt electron
      if (electron.pt() > output.fHighestSelectedElectronPt) 
	{
	  output.fHighestSelectedElectronPt = electron.pt();
	  output.fHighestSelectedElectronEta = electron.eta();
	}
      
      // Save all electrons surviving the cuts
      output.fSelectedElectrons.push_back(electron);
      
      // Fill resolution histograms
      if (event.isMC()) 
	{
	  hPtResolution->Fill((electron.pt() - electron.MCelectron()->pt()) / electron.pt());
	  hEtaResolution->Fill((electron.eta() - electron.MCelectron()->eta()) / electron.eta());
	  hPhiResolution->Fill((electron.phi() - electron.MCelectron()->phi()) / electron.phi());
	}
      
    }//for-loop: electrons
  
  
  //sort electrons, needed comparisons defined in Electron.h
  std::sort(output.fSelectedElectrons.begin(), output.fSelectedElectrons.end());

  // Assign booleans
  passedSelection = (output.fSelectedElectrons.size() > 0);
  passedVeto      = (output.fSelectedElectrons.size() == 0); 


  // Set electron ID SF value to data object
  if (event.isMC()) {
    if (output.hasIdentifiedElectrons()) {
      output.fElectronIDSF = fElectronIDSFReader.getScaleFactorValue(output.getSelectedElectrons()[0].eta());
    }
  }

  // Set electron trigger SF value to data object
  if (event.isMC()) {
    if (output.hasIdentifiedElectrons()) {
      output.fElectronTriggerSF = fElectronTriggerSFReader.getScaleFactorValue(output.getSelectedElectrons()[0].pt());
    }
  }

  // Fill histos
  hElectronNAll->Fill(event.electrons().size());
  hElectronNPassed->Fill(output.fSelectedElectrons.size());

  // Fill sub-counters
  if (passedIsPresent) cSubPassedIsPresent.increment();
  if (passedTrgMatch) cSubPassedTriggerMatching.increment();
  if (passedPt) cSubPassedPt.increment();
  if (passedEta) cSubPassedEta.increment();
  if (passedID) cSubPassedID.increment();
  if (passedIsol) cSubPassedIsolation.increment();
  if (passedSelection) cSubPassedSelection.increment();
  if (passedVeto) cSubPassedVeto.increment();
  if (fVetoMode) 
    {
      if (passedVeto) cPassedElectronSelection.increment();
    }
  else
    {
      if (passedSelection) cPassedElectronSelection.increment();
    }

  // Return data object
  return output;
}

bool ElectronSelection::getMVADecision(const Electron& ele, const std::string mvaCut){

  if (mvaCut == "loose" || mvaCut == "Loose")
    {
      double AbsEta = std::abs(ele.eta());
      
      if (AbsEta<=0.8 && ele.MVA()>=-0.041){
	return true;
      }
      if (AbsEta>0.8 && AbsEta<1.479 && ele.MVA()>=0.383){
	return true;
      }
      if (AbsEta>=1.479 && ele.MVA()>=-0.515){
	return true;
      }
    }
  else
    {
      throw hplus::Exception("config") << "Invalid electronMVACut option '" << mvaCut << "'! Options: 'Loose'";
    }
  return false;
}

bool ElectronSelection::passTrgMatching(const Electron& electron, std::vector<math::LorentzVectorT<double>>& trgElectrons) const {
  if (!cfg_ApplyTriggerMatching) return true;

  double myMinDeltaR = 9999.0;

  // For-loop: Trigger electrons
  for (auto& p: trgElectrons) 
    {
      double myDeltaR = ROOT::Math::VectorUtil::DeltaR(p, electron.p4());
      myMinDeltaR = std::min(myMinDeltaR, myDeltaR);
    }

  // Fill histos
  hTriggerMatchDeltaR->Fill(myMinDeltaR);
  return (myMinDeltaR < cfg_TriggerMatchingCone);
}
