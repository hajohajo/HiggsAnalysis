// -*- cut++ -*-
#include "EventSelection/interface/HplusSelectionNN.h"

#include "Framework/interface/ParameterSet.h"
#include "EventSelection/interface/CommonPlots.h"
#include "DataFormat/interface/Event.h"
#include "Framework/interface/HistoWrapper.h"
#include "Framework/interface/Exception.h"

#include "Tools/interface/MCTools.h"

#include "Math/VectorUtil.h"

HplusSelectionNN::Data::Data()
:
  bPassedSelection(false),
  bHasTwoTopsAndFreeB(false),
  nAllCleanedTops(),
  fMVAmax1(-1.0),
  fMVAmax2(-1.0),
  fTrijet1_p4(),
  fTrijet2_p4(),
  fTrijet1BJet(),
  fTrijet1Jet1(),
  fTrijet1Jet2(),
  fTrijet2BJet(),
  fTrijet2Jet1(),
  fTrijet2Jet2(),
  fTrijet1Dijet_p4(),
  fTrijet2Dijet_p4(),
  fLdgTetrajet_p4(),
  fTetrajetBJet(),
  fTopTaggingScaleFactorEventWeight(1.0)
{ }

HplusSelectionNN::Data::~Data() { }


HplusSelectionNN::HplusSelectionNN(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& postfix)
  : BaseSelection(eventCounter, histoWrapper, commonPlots, postfix),
    // Input parameters
    cfg_TopMVACut(config, "TopMVACut"),
    cfg_AnyTopMVACut(config, "AnyTopMVACut"),
    cfg_FreeBjetsCut(config, "FreeBjetsCut"),
    fTopTagSFCalculator(config)
{
  initialize(config);
}

HplusSelectionNN::HplusSelectionNN(const ParameterSet& config)
: BaseSelection(),
  // Input parameters
  cfg_TopMVACut(config, "TopMVACut"),
  cfg_AnyTopMVACut(config, "AnyTopMVACut"),
  cfg_FreeBjetsCut(config, "FreeBjetsCut"),
  fTopTagSFCalculator(config)
{
  initialize(config);
  bookHistograms(new TDirectory());

}

HplusSelectionNN::~HplusSelectionNN() {
  
  // Histograms
  //delete hTopMultiplicity_AllCandidates;
}

void HplusSelectionNN::initialize(const ParameterSet& config) {
  
  return;
}

void HplusSelectionNN::bookHistograms(TDirectory* dir) {
  // Fixed binning  
  /*
  const int nPtBins   = 100;  //2 * fCommonPlots->getPtBinSettings().bins();
  const double fPtMin = 0;    //2 * fCommonPlots->getPtBinSettings().min();
  const double fPtMax = 1000; //2 * fCommonPlots->getPtBinSettings().max();

  const int  nEtaBins = 50;   //fCommonPlots->getEtaBinSettings().bins();
  const float fEtaMin = -5.0; //fCommonPlots->getEtaBinSettings().min();
  const float fEtaMax = 5.0;  //fCommonPlots->getEtaBinSettings().max();

  const int nDRBins   = 100;  //fCommonPlots->getDeltaRBinSettings().bins();
  const double fDRMin = 0;    //fCommonPlots->getDeltaRBinSettings().min();
  const double fDRMax = 10;   //fCommonPlots->getDeltaRBinSettings().max();

  const int  nBDiscBins = 120; //fCommonPlots->getBJetDiscBinSettings().bins();
  const float fBDiscMin = 0;   //fCommonPlots->getBJetDiscBinSettings().min();
  const float fBDiscMax = 1.2; //fCommonPlots->getBJetDiscBinSettings().max();

  const int nWMassBins  = 200; //fCommonPlots->getWMassBinSettings().bins();
  const float fWMassMin = 0;   //fCommonPlots->getWMassBinSettings().min();
  const float fWMassMax = 1000;//fCommonPlots->getWMassBinSettings().max();

  const int nTopMassBins  = 200; //fCommonPlots->getTopMassBinSettings().bins();
  const float fTopMassMin = 0;   //fCommonPlots->getTopMassBinSettings().min();
  const float fTopMassMax = 1000;//fCommonPlots->getTopMassBinSettings().max();
  */
  // Histograms
  TDirectory* subdir = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "hplus2tbSelectionMVA_"    + sPostfix);

  // All Candidates
  //hTopMVA_AllCandidates           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMVA_AllCandidates",";top candidate MVA", 40, -1.0, 1.0) ; 
  fTopTagSFCalculator.bookHistograms(subdir, fHistoWrapper);
  return;
}


HplusSelectionNN::Data HplusSelectionNN::silentAnalyze(const Event& event, const JetSelection::Data& jetData, const BJetSelection::Data& bjetData, const TopSelectionNN::Data& topData) {
  ensureSilentAnalyzeAllowed(event.eventID());

  // Disable histogram filling and counter
  disableHistogramsAndCounters();
  Data myData = privateAnalyze(event, jetData.getSelectedJets(), bjetData.getSelectedBJets(), topData);
  enableHistogramsAndCounters();
  return myData;
}


HplusSelectionNN::Data HplusSelectionNN::analyze(const Event& event, const JetSelection::Data& jetData, const BJetSelection::Data& bjetData, const TopSelectionNN::Data& topData) {
  ensureAnalyzeAllowed(event.eventID());

  // Ready to analyze
  HplusSelectionNN::Data data = privateAnalyze(event, jetData.getSelectedJets(), bjetData.getSelectedBJets(), topData);

  // Send data to CommonPlots
  // if (fCommonPlots != nullptr) fCommonPlots->fillControlPlotsAtHplusSelectionNN(event, data);// fixme - implement
  return data;
}


HplusSelectionNN::Data HplusSelectionNN::privateAnalyze(const Event& event, const std::vector<Jet> selectedJets, const std::vector<Jet> selectedBjets, const TopSelectionNN::Data& topData) {
  Data output;

  std::vector<Jet> jets  = selectedJets;
  std::vector<Jet> bjets = selectedBjets;

  std::vector<Jet>	allTopsJet1	 = topData.getAllTopsJet1();
  std::vector<Jet>	allTopsJet2	 = topData.getAllTopsJet2();
  std::vector<Jet>	allTopsBJet	 = topData.getAllTopsBJet();  
  std::vector<float>	allTopsMVA      = topData.getAllTopsMVA();
  std::vector<bool>	allTopsIsGenuine = topData.getAllTopsIsGenuine();
  std::vector<bool>	allTopsIsTagged  = topData.getAllTopsIsTagged();

  std::vector<Jet>	allCleanedTopsJet1	= topData.getAllCleanedTopsJet1();
  std::vector<Jet>	allCleanedTopsJet2	= topData.getAllCleanedTopsJet2();
  std::vector<Jet>	allCleanedTopsBJet	= topData.getAllCleanedTopsBJet();
  std::vector<float>	allCleanedTopsMVA      = topData.getAllCleanedTopsMVA();
  std::vector<bool>	allCleanedTopsIsGenuine = topData.getAllCleanedTopsIsGenuine();
  std::vector<bool>	allCleanedTopsIsTagged  = topData.getAllCleanedTopsIsTagged();
  
  //Define variables
  Jet tetrajetBjet;
  double tetrajetBjetPt_max = -999.99;

  output.fTopTaggingScaleFactorEventWeight = topData.getTopTaggingScaleFactorEventWeight();
  //================================================================================================  
  // Hplus2tb Selection
  //================================================================================================  
  if (0) std::cout << "=== HplusSelectionNN:: " << std::endl;
  
  //Check if there are free bjets after reconstructing selected cleaned tops with highest MVA.
  int nbjets_Top1 = 0, nbjets_Top2 = 0; 
  
  if (topData.getAllCleanedTopsSize() > 0) nbjets_Top1 = isBJet(allCleanedTopsJet1.at(0), bjets) + isBJet(allCleanedTopsJet2.at(0), bjets) + isBJet(allCleanedTopsBJet.at(0), bjets);  
  if (topData.getAllCleanedTopsSize() > 1) nbjets_Top2 = isBJet(allCleanedTopsJet1.at(1), bjets) + isBJet(allCleanedTopsJet2.at(1), bjets) + isBJet(allCleanedTopsBJet.at(1), bjets);

  int n_freeBjets1 = (int)bjets.size() - nbjets_Top1;
  int n_freeBjets2 = (int)jets.size() - nbjets_Top2;

  //If the number of free bjets is less that fFreeBjetsCuts, skip trijet and cross clean the new candidates
  if ( !cfg_FreeBjetsCut.passedCut(n_freeBjets1) || !cfg_FreeBjetsCut.passedCut(n_freeBjets2)){
    //Clear old collections
    allCleanedTopsJet1.clear(); allCleanedTopsJet2.clear(); allCleanedTopsBJet.clear();
    allCleanedTopsMVA.clear(); allCleanedTopsIsGenuine.clear(); allCleanedTopsIsTagged.clear();   
    
    std::vector< math::XYZTLorentzVector> allCleanedTops_P4;
    for (size_t i = 0; i < allTopsBJet.size(); i++){
      int nbjets_i = isBJet(allTopsJet1.at(i), bjets) + isBJet(allTopsJet2.at(i), bjets) + isBJet(allTopsBJet.at(i), bjets);
      bool isCrossCleaned = TopIsCrossCleaned(i, allTopsJet1, allTopsJet2, allTopsBJet) && ((size_t)nbjets_i != bjets.size());
      if (!isCrossCleaned) continue;      
      
      //Store new coss-cleaned candidates.
      allCleanedTopsJet1.push_back(allTopsJet1.at(i));
      allCleanedTopsJet2.push_back(allTopsJet2.at(i));
      allCleanedTopsBJet.push_back(allTopsBJet.at(i));
      allCleanedTopsMVA.push_back(allTopsMVA.at(i));
      allCleanedTopsIsGenuine.push_back(allTopsIsGenuine.at(i));
      allCleanedTopsIsTagged.push_back(allTopsIsTagged.at(i));
      math::XYZTLorentzVector top_p4;
      top_p4 = allTopsJet1.at(i).p4() + allTopsJet2.at(i).p4() + allTopsBJet.at(i).p4();
      allCleanedTops_P4.push_back( top_p4 );
    }
    // Re-calculate and store top tagging SF weight and it's uncertainty
    if (event.isMC())
      {
	output.fTopTaggingScaleFactorEventWeight = fTopTagSFCalculator.calculateSF(allCleanedTops_P4, allCleanedTopsMVA, allCleanedTopsIsTagged, allCleanedTopsIsGenuine);    
      }
  } //if (nbjets_Top1 == bjets.size() || nbjets_Top2 == bjets.size()){
  
  
  if (0) std::cout << "=== HplusSelectionNN:: Find the tetrajet b-jet" << std::endl;

  //For-loop: All selected b jets.
  for (auto& bjet: bjets)
    {
      // Skip if tetrajet bjet pT is greater that this pt
      if (tetrajetBjetPt_max > bjet.pt()) continue;

       // Check if this bjet is matched with one of the jets assigned to the best TWO tops
      if (isMatchedJet(bjet, allCleanedTopsJet1, allCleanedTopsJet2, allCleanedTopsBJet, 0)) continue;
      if (isMatchedJet(bjet, allCleanedTopsJet1, allCleanedTopsJet2, allCleanedTopsBJet, 1)) continue;

      // TetrajetBjet: Save variables
      tetrajetBjetPt_max = bjet.pt();
      tetrajetBjet       = bjet;
    }


  //Fill output data
  if (0) std::cout << "=== HplusSelectionNN:: Fill output" << std::endl;

  output.nAllCleanedTops = allCleanedTopsMVA.size();
  
  bool bPass_LdgMVA     = false;   
  bool bPass_SubldgMVA  = false;

  if (allCleanedTopsMVA.size() == 1)
    {
      bPass_LdgMVA    = cfg_TopMVACut.passedCut( allCleanedTopsMVA.at(0) );

      // Leading-in-MVA top
      output.fMVAmax1           = allCleanedTopsMVA.at(0);
      output.fTrijet1Jet1       = allCleanedTopsJet1.at(0);
      output.fTrijet1Jet2       = allCleanedTopsJet2.at(0);
      output.fTrijet1BJet       = allCleanedTopsBJet.at(0);
      output.fTrijet1Dijet_p4   = allCleanedTopsJet1.at(0).p4() + allCleanedTopsJet2.at(0).p4();
      output.fTrijet1_p4        = allCleanedTopsJet1.at(0).p4() + allCleanedTopsJet2.at(0).p4() + allCleanedTopsBJet.at(0).p4();
    }
  else if (allCleanedTopsMVA.size() > 1)
    {
      bPass_LdgMVA     = cfg_TopMVACut.passedCut( allCleanedTopsMVA.at(0) );
      bPass_SubldgMVA  = cfg_TopMVACut.passedCut( allCleanedTopsMVA.at(1) );
      // Leading-in-MVA top
      output.fMVAmax1           = allCleanedTopsMVA.at(0);
      output.fTrijet1Jet1       = allCleanedTopsJet1.at(0);
      output.fTrijet1Jet2       = allCleanedTopsJet2.at(0);
      output.fTrijet1BJet       = allCleanedTopsBJet.at(0);
      output.fTrijet1Dijet_p4   = allCleanedTopsJet1.at(0).p4() + allCleanedTopsJet2.at(0).p4();
      output.fTrijet1_p4        = allCleanedTopsJet1.at(0).p4() + allCleanedTopsJet2.at(0).p4() + allCleanedTopsBJet.at(0).p4();
      // Subleading-in-MVA top
      output.fMVAmax2           = allCleanedTopsMVA.at(1);
      output.fTrijet2Jet1       = allCleanedTopsJet1.at(1);
      output.fTrijet2Jet2       = allCleanedTopsJet2.at(1);
      output.fTrijet2BJet       = allCleanedTopsBJet.at(1);
      output.fTrijet2Dijet_p4   = allCleanedTopsJet1.at(1).p4() + allCleanedTopsJet2.at(1).p4();
      output.fTrijet2_p4        = allCleanedTopsJet1.at(1).p4() + allCleanedTopsJet2.at(1).p4() + allCleanedTopsBJet.at(1).p4();
    }

  if (tetrajetBjetPt_max > 0)
    {
      output.fTetrajetBJet = tetrajetBjet;
      
      if (allCleanedTopsMVA.size() == 1) output.fLdgTetrajet_p4 =  output.fTrijet1_p4 + output.fTetrajetBJet.p4();
      else if (allCleanedTopsMVA.size() > 1)
	{
	  bool bLdgMVAisLdgPt = output.fTrijet1_p4.Pt() > output.fTrijet2_p4.Pt();
	  output.fLdgTetrajet_p4 = (bLdgMVAisLdgPt) ? (output.fTrijet1_p4 + output.fTetrajetBJet.p4()) : (output.fTrijet2_p4 + output.fTetrajetBJet.p4());
	}
    }
  
  //Booleans
  bool bPass_FreeBjet   = (tetrajetBjetPt_max > 0);
  bool bPass_BothMVA    = bPass_LdgMVA * bPass_SubldgMVA;
  bool bPass_AnyTwoTops = false; // at least TWO tops with MVA > -1.0 
  if (allCleanedTopsMVA.size() > 1) bPass_AnyTwoTops = cfg_AnyTopMVACut.passedCut(allCleanedTopsMVA.at(1) );

  output.bPassedSelection = bPass_FreeBjet * bPass_BothMVA;
  output.bHasTwoTopsAndFreeB = bPass_AnyTwoTops * bPass_FreeBjet;

  return output;
}

bool HplusSelectionNN::areSameJets(const Jet& jet1, const Jet& jet2) {
  float dR = ROOT::Math::VectorUtil::DeltaR(jet1.p4(), jet2.p4());
  float dR_match = 0.1;
  if (dR <= dR_match) return true;
  else return false;
}

bool HplusSelectionNN::isBJet(const Jet& jet, const std::vector<Jet>& bjets) {
  for (auto bjet: bjets)
    {
      if (areSameJets(jet, bjet)) return true;
    }
  return false;
}


bool HplusSelectionNN::isMatchedJet(const Jet& jet, std::vector<Jet> Jet1, std::vector<Jet> Jet2, std::vector<Jet> BJet, const unsigned int index) {
  // Sanity check: If index is bigger than the number of tops that top does not exist => not possible to be matched
  if ( index+1 > BJet.size()) return false;

  std::vector<Jet> jets;
  jets.push_back(Jet1.at(index));
  jets.push_back(Jet2.at(index));
  jets.push_back(BJet.at(index));
  
  for (auto iJet: jets)
    {
      if (areSameJets(jet, iJet)) return true;
    }
  return false;
}

bool HplusSelectionNN::TopIsCrossCleaned( int Index, std::vector<Jet> jet1, std::vector<Jet> jet2, std::vector<Jet> bjet){
  // Description:
  // Used to find the cross-cleaned trijet multiplicity. The function takes as input the index of a trijet and the total trijet collection and returns 
  // false true if If at least one of the subjets of the trijets are used by the trijets with Higher MVA value (Higher MVA value -> smaller index: Trijets sorted in MVA value) 
  // Returns true: If cross-cleaned trijet 
  
  // Assume sorted in MVA score top candidats
  if (Index > 0)
    {

      // For-loop: All top candidates with higher MVA score than this one
      for (size_t i=0; i<(size_t)Index; i++)
	{
	  // Skip top candidates with same jets as Trijets with higher MVA value
	  bool bMatchedJet1 = isMatchedJet(jet1.at(Index), jet1, jet2, bjet, i);
	  bool bMatchedJet2 = isMatchedJet(jet2.at(Index), jet1, jet2, bjet, i);
	  bool bMatchedBJet = isMatchedJet(bjet.at(Index), jet1, jet2, bjet, i);
	  bool sharedJets   = (bMatchedJet1 || bMatchedJet2 || bMatchedBJet);

	  // Candidate is not cross-cleaned if a shared jet is found and the candidate this is shared with is itself cross-cleaned
	  if (sharedJets && TopIsCrossCleaned(i, jet1, jet2, bjet)) return false;
	}
    }
  
  return true;
}
