// -*- c++ -*-
#include "EventSelection/interface/AngularCuts.h"

#include "Framework/interface/ParameterSet.h"
#include "EventSelection/interface/CommonPlots.h"
#include "DataFormat/interface/Event.h"
#include "Framework/interface/HistoWrapper.h"
#include "Framework/interface/Exception.h"

#include "Math/VectorUtil.h"

#include <sstream>
#include <cmath>

AngularCutsBase::Data::Data() 
: bPassedSelection(false) { }

AngularCutsBase::Data::~Data() { }

bool AngularCutsBase::Data::passedSelectionOnJet(size_t n) const {
//   if (n >= nMaxJets)
//     throw hplus::Exception("assert") << "AngularCuts::Data: Requested passedSelectionOnJet on jet " << n << ", but maximum is " << nMaxJets;
  if (n >= fPassedCutStatus.size())
    return bPassedSelection;
  return fPassedCutStatus[n];
}

double AngularCutsBase::Data::getDeltaPhiJetMET(size_t n) const {
//   if (n >= nMaxJets)
//     throw hplus::Exception("assert") << "AngularCuts::Data: Requested DeltaPhi(Jet,MET) on jet " << n << ", but maximum is " << nMaxJets;
  if (n >= fDeltaPhiJetMET.size())
    return -1.0;
  return fDeltaPhiJetMET[n];
}

double AngularCutsBase::Data::get1DCutVariable(size_t n) const {
//   if (n >= nMaxJets)
//     throw hplus::Exception("assert") << "AngularCuts::Data: Requested get1DCutVariable on jet " << n << ", but maximum is " << nMaxJets;
  if (n >= f1DCutVariables.size())
    return -1.0;
  return f1DCutVariables[n];
}

AngularCutsBase::AngularCutsBase(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& prefix, const AngularCutsBase::AngularCutsType type, const std::string& postfix)
: BaseSelection(eventCounter, histoWrapper, commonPlots, prefix+postfix),
  // Input parameters
  nMaxJets(4),
  nConsideredJets(static_cast<size_t>(config.getParameter<int>("nConsideredJets"))),
  bEnableOptimizationPlots(config.getParameter<bool>("enableOptimizationPlots")),
  sPrefix(prefix),
  fType(type),
  cPassedAngularCuts(fEventCounter.addCounter("passed angular cuts "+prefix+" ("+postfix+")")),
  cSubAllEvents(fEventCounter.addSubCounter("angular cuts "+prefix+" ("+postfix+")", "All events"))
{
  initialize(config, postfix);
}

AngularCutsBase::AngularCutsBase(const ParameterSet& config, const AngularCutsBase::AngularCutsType type)
: BaseSelection(),
  // Input parameters
  nMaxJets(4),
  nConsideredJets(static_cast<size_t>(config.getParameter<int>("nConsideredJets"))),
  bEnableOptimizationPlots(config.getParameter<bool>("enableOptimizationPlots")),
  sPrefix(""),
  fType(type),
  cPassedAngularCuts(fEventCounter.addCounter("passed angular cuts")),
  cSubAllEvents(fEventCounter.addSubCounter("angular cuts", "All events"))
{
  initialize(config, "");
  bookHistograms(new TDirectory());
}

AngularCutsBase::~AngularCutsBase() { 
  for (size_t i = 0; i < hOptimizationPlots.size(); ++i) {
    if (hOptimizationPlots[i] != nullptr)
      delete hOptimizationPlots[i];
  }
  hOptimizationPlots.clear();
}

void AngularCutsBase::initialize(const ParameterSet& config, const std::string& postfix) {
  // Check validity of parameters
  if (static_cast<size_t>(nConsideredJets) > nMaxJets)
    throw hplus::Exception("config") << "AngularCuts: Requested cuts on " << nConsideredJets << " jets, but maximum is " << nMaxJets;
  if (fType == kUndefined)
    throw hplus::Exception("assert") << "AngularCutsBase type has not been properly defined!";
  // Obtain cut points from config
  size_t maxIndex = std::min(nMaxJets, nConsideredJets);
  std::stringstream s;
  for (size_t i = 0; i < maxIndex; ++i) {
    s.str("");
    s << "cutValueJet" << i+1;
    fCutValue.push_back(static_cast<double>(config.getParameter<float>(s.str())));
  }
  // Subcounters
  for (size_t i = 0; i < nMaxJets; ++i) {
    s.str("");
    s << "Passed cut on jet " << i+1;
    cSubPassedCuts.push_back(fEventCounter.addSubCounter("angular cuts "+sPrefix+" ("+postfix+")", s.str()));
  }
}

void AngularCutsBase::bookHistograms(TDirectory* dir) {
  TDirectory* subdir = fHistoWrapper.mkdir(HistoLevel::kDebug, dir, "AngularCuts_"+sPostfix);
  if (bEnableOptimizationPlots) {
    std::stringstream sName;
    std::stringstream sLabel;
    for (size_t i = 0; i < nMaxJets; ++i) {
      sName.str("");
      sLabel.str("");
      sName << "AngularCuts_" << sPrefix << "_jet" << i+1;
      sLabel << "AngularCuts_" << sPrefix << "_jet" << i+1 << ";#Delta#phi(#tau,MET), ^{o};#Delta#phi(jet_{" << i+1 << "},MET), ^{o}";
      hOptimizationPlots.push_back(fHistoWrapper.makeTH<TH2F>(HistoLevel::kDebug, subdir, sName.str().c_str(), sLabel.str().c_str(),
                                                             18,0.,180.,18,0.,180.));
    }
  }
}

AngularCutsBase::Data AngularCutsBase::silentAnalyze(const Event& event, const Tau& tau, const JetSelection::Data& jetData, const METSelection::Data& metData) {
  ensureSilentAnalyzeAllowed(event.eventID());
  // Disable histogram filling and counter
  disableHistogramsAndCounters();
  Data myData = privateAnalyze(tau, jetData, metData);
  enableHistogramsAndCounters();
  return myData;
}

AngularCutsBase::Data AngularCutsBase::analyze(const Event& event, const Tau& tau, const JetSelection::Data& jetData, const METSelection::Data& metData) {
  ensureAnalyzeAllowed(event.eventID());
  AngularCutsBase::Data data = privateAnalyze(tau, jetData, metData);
  // Send data to CommonPlots
  if (fCommonPlots != nullptr) {
    if (fType == kCollinear)
      {
	// std::cout << "=== AngularCutBase::analyze() fType = " << fType << "(kCollinear)" << std::endl;
	fCommonPlots->fillControlPlotsAtAngularCutsCollinear(event, data);
      }
    else if (fType == kBackToBack)
      {
	// std::cout << "=== AngularCutBase::analyze() fType = " << fType << "(kBackToBack)" << std::endl;
	fCommonPlots->fillControlPlotsAtAngularCutsBackToBack(event, data);
      }
    else
      {
	throw hplus::Exception("assert") << "AngularCuts::Data: Unexpected fType = " << fType << " [options: kCollinear, kBackToBack]. This should not be reached";
      }
  }
  // Return data
  return data;
}

AngularCutsBase::Data AngularCutsBase::analyze(const Event& event, const Muon& muon, const std::vector<Tau>& taus, const JetSelection::Data& jetData, const BJetSelection::Data& bjetData, const METSelection::Data& metData) {
  ensureAnalyzeAllowed(event.eventID());
  AngularCutsBase::Data data = privateAnalyze(muon, taus, jetData, bjetData, metData);

  // Send data to CommonPlots
  if (fCommonPlots != nullptr) {
    if (fType == kCollinear)
      {
	// std::cout << "=== AngularCutBase::analyze() fType = " << fType << "(kCollinear)" << std::endl;
	fCommonPlots->fillControlPlotsAtAngularCutsCollinear(event, data);
      }
    else if (fType == kBackToBack)
      {
	// std::cout << "=== AngularCutBase::analyze() fType = " << fType << "(kBackToBack)" << std::endl;
	fCommonPlots->fillControlPlotsAtAngularCutsBackToBack(event, data);
      }
    else
      {
	throw hplus::Exception("assert") << "AngularCuts::Data: Unexpected fType = " << fType << " [options: kCollinear, kBackToBack]. This should not be reached";
      }
  }

  return data;
}

AngularCutsBase::Data AngularCutsBase::privateAnalyze(const Tau& tau, const JetSelection::Data& jetData, const METSelection::Data& metData) {
  Data output;
  cSubAllEvents.increment();
  output.bPassedSelection = true;
  // Calculate delta phi between MET and tau
  output.fDeltaPhiTauMET = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(tau.p4(), metData.getMET())*57.29578);
  size_t maxIndex = std::min(nMaxJets, jetData.getSelectedJets().size());
  for (size_t i = 0; i < maxIndex; ++i) {
    // Calculate delta phi between MET and jet
    double dphi = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(jetData.getSelectedJets()[i].p4(), metData.getMET())*57.29578);
    output.fDeltaPhiJetMET.push_back(dphi);
    if (bEnableOptimizationPlots)
      hOptimizationPlots[i]->Fill(output.fDeltaPhiTauMET, dphi);
    // Obtain cut value and status (i.e. calculate cut variable for all existing jets)
    bool passedStatus = true;
    if (fType == kCollinear)
      {
      passedStatus = doCollinearCuts(output.fDeltaPhiTauMET, dphi, fCutValue[i], output.f1DCutVariables);
      }
    else if (fType == kBackToBack)
      {
      passedStatus = doBackToBackCuts(output.fDeltaPhiTauMET, dphi, fCutValue[i], output.f1DCutVariables);
      }
    else
      {
	throw hplus::Exception("assert") << "AngularCuts::Data: Unexpected fType = " << fType << " [options: kCollinear, kBackToBack]. This should not be reached";
      }

    // Make cut by updating passed status
    if (i < nConsideredJets) 
      {
	output.bPassedSelection = output.bPassedSelection && passedStatus;

	if (output.bPassedSelection) cSubPassedCuts[i].increment();
      }
    output.fPassedCutStatus.push_back(output.bPassedSelection);
  }

  // Obtain minimum value
  output.fMinimumCutValue = 999.0;
  for (size_t i = 0; i < maxIndex; ++i) 
    {
      if (i < nConsideredJets)
	{
	  if (output.f1DCutVariables[i] >= 0.0 && output.f1DCutVariables[i] < output.fMinimumCutValue)
	    {
	      output.fMinimumCutValue = output.f1DCutVariables[i];
	    }
	}
    }
  
  // Fill main counter if passed
  if (output.bPassedSelection) cPassedAngularCuts.increment();

  // Return data object
  return output;
}


AngularCutsBase::Data AngularCutsBase::privateAnalyze(const Muon& muon, const std::vector<Tau>& taus, const JetSelection::Data& jetData, const BJetSelection::Data& bjetData, const METSelection::Data& metData) {
  Data output;
  cSubAllEvents.increment();
  output.bPassedSelection = true;

  // Sanity check
  if (taus.size() < 2) 
    {
      throw hplus::Exception("assert") << "AngularCuts::privateAnalyze(): Overloaded method called that required at least 2 selected taus, but only " << taus.size() << " are found.";
    }

  // Calculate delta phi between MET and tau
  output.fDeltaPhiLdgBJetMET = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(bjetData.getSelectedBJets()[0].p4(), metData.getMET())*57.29578);
  output.fDeltaPhiLdgJetMET = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(jetData.getSelectedJets()[0].p4(), metData.getMET())*57.29578);
  output.fDeltaPhiSubldgJetMET = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(jetData.getSelectedJets()[1].p4(), metData.getMET())*57.29578);
  output.fDeltaPhiMuonMET = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(muon.p4(), metData.getMET())*57.29578);
  output.fDeltaPhiTauMET = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(taus[0].p4() + taus[1].p4(), metData.getMET())*57.29578); // treat 2 taus as a single 4-vector 
  output.fDeltaPhiLdgTauMET = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(taus[0].p4(), metData.getMET())*57.29578);
  output.fDeltaPhiSubldgTauMET = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(taus[1].p4(), metData.getMET())*57.29578);
  output.fDeltaPhiLdgTauMuon = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(taus[0].p4(), muon.p4())*57.29578);
  output.fDeltaPhiSubldgTauMuon = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(taus[1].p4(), muon.p4())*57.29578);
  output.fDeltaPhiLdgTauLdgJet = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(taus[0].p4(), jetData.getSelectedJets()[0].p4())*57.29578);
  output.fDeltaPhiSubldgTauLdgJet = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(taus[1].p4(), jetData.getSelectedJets()[0].p4())*57.29578);
  output.fDeltaPhiTaus = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(taus[0].p4(), taus[1].p4())*57.29578);

  if (0)
    {
      std::cout << "output.fDeltaPhiLdgBJetMET = " << output.fDeltaPhiLdgBJetMET << std::endl;
      std::cout << "output.fDeltaPhiLdgJetMET = " << output.fDeltaPhiLdgJetMET << std::endl;
      std::cout << "output.fDeltaPhiSubldgJetMET = " << output.fDeltaPhiSubldgJetMET << std::endl;
      std::cout << "output.fDeltaPhiMuonMET = " << output.fDeltaPhiMuonMET << std::endl;
      std::cout << "output.fDeltaPhiTauMET = " << output.fDeltaPhiTauMET << std::endl;
      std::cout << "output.fDeltaPhiLdgTauMET = " << output.fDeltaPhiLdgTauMET << std::endl;
      std::cout << "output.fDeltaPhiSubldgTauMET = " << output.fDeltaPhiSubldgTauMET << std::endl;
      std::cout << "output.fDeltaPhiLdgTauMuon = " << output.fDeltaPhiLdgTauMuon << std::endl;
      std::cout << "output.fDeltaPhiSubldgTauMuon = " << output.fDeltaPhiSubldgTauMuon << std::endl;
      std::cout << "output.fDeltaPhiLdgTauLdgJet = " << output.fDeltaPhiLdgTauLdgJet << std::endl;
      std::cout << "output.fDeltaPhiSubldgTauLdgJet = " << output.fDeltaPhiSubldgTauLdgJet << std::endl;
      std::cout << "output.fDeltaPhiTaus = " << output.fDeltaPhiTaus << std::endl;
    }


  // FIXME: Implement code that refects HToHW needs
  size_t maxIndex = std::min(nMaxJets, jetData.getSelectedJets().size());
  for (size_t i = 0; i < maxIndex; ++i) {
    // Calculate delta phi between MET and jet
    double dphi = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(jetData.getSelectedJets()[i].p4(), metData.getMET())*57.29578);
    output.fDeltaPhiJetMET.push_back(dphi);

    // Fill optimisation histo?
    if (bEnableOptimizationPlots) hOptimizationPlots[i]->Fill(output.fDeltaPhiTauMET, dphi);

    // Obtain cut value and status (i.e. calculate cut variable for all existing jets)
    bool passedStatus = true;

    // Determine status
    if (fType == kCollinear)
      {
       // The method doCollinear only differs from doBackToBack by interchanging the x and y axes so that the corner of the 2d plots
       // changes meaning (back-to-back becomes collinear)
	passedStatus = doCollinearCuts(output.fDeltaPhiTauMET, dphi,  fCutValue[i], output.f1DCutVariables);
      }
    else if (fType == kBackToBack)
      {
	passedStatus = doBackToBackCuts(output.fDeltaPhiTauMET, dphi, fCutValue[i], output.f1DCutVariables);
      }
    else
      {
	// fixme - introduce a new category type for HToHW: tau-tau, tau-muon, bjet-MET
      }    

    // Make cut by updating passed status
    if (i < nConsideredJets) 
      {
	output.bPassedSelection = output.bPassedSelection && passedStatus;
	if (output.bPassedSelection) cSubPassedCuts[i].increment();
      }
    output.fPassedCutStatus.push_back(output.bPassedSelection);
  }
  
  // Obtain minimum value
  output.fMinimumCutValue = 999.0;
  for (size_t i = 0; i < maxIndex; ++i) 
    {
      if (i < nConsideredJets)
	{
	  if (output.f1DCutVariables[i] >= 0.0 && output.f1DCutVariables[i] < output.fMinimumCutValue) {
	    output.fMinimumCutValue = output.f1DCutVariables[i];
	  }
	}
    }

  // Fill main counter if passed
  if (output.bPassedSelection) cPassedAngularCuts.increment();

  return output;
}

bool AngularCutsBase::doCollinearCuts(const double deltaPhiTauMET, const double deltaPhiJetMET, double cutValue, std::vector< double >& results) {
  double x = deltaPhiJetMET - 180.0;
  double y = deltaPhiTauMET;
  double value = std::sqrt(x*x + y*y);
  results.push_back(value);
  return value > cutValue;
}

bool AngularCutsBase::doBackToBackCuts(const double deltaPhiTauMET, const double deltaPhiJetMET, double cutValue, std::vector< double >& results) {
  double x = deltaPhiTauMET - 180.0;
  double y = deltaPhiJetMET;
  double value = std::sqrt(x*x + y*y);
  results.push_back(value);
  return value > cutValue;
}

AngularCutsCollinear::AngularCutsCollinear(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& postfix)
: AngularCutsBase(config, eventCounter, histoWrapper, commonPlots, "Collinear", AngularCutsBase::kCollinear, postfix) { }

AngularCutsCollinear::~AngularCutsCollinear() { }

AngularCutsBackToBack::AngularCutsBackToBack(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& postfix)
: AngularCutsBase(config, eventCounter, histoWrapper, commonPlots, "BackToBack", AngularCutsBase::kBackToBack, postfix) { }

AngularCutsBackToBack::~AngularCutsBackToBack() { }
