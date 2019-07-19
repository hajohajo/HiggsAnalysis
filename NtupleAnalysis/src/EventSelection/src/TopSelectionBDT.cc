// -*- cut++ -*-
#include "EventSelection/interface/TopSelectionBDT.h"

#include "Framework/interface/ParameterSet.h"
#include "EventSelection/interface/CommonPlots.h"
#include "DataFormat/interface/Event.h"
#include "Framework/interface/HistoWrapper.h"
#include "Framework/interface/Exception.h"

#include "Tools/interface/MCTools.h"

#include "Math/VectorUtil.h"

TopSelectionBDT::Data::Data()
:
  bPassedBDTGCut(false),
  bPassedNTopsCut(false),
  bPassedSelection(false),
  fTopBDTG(-1.0),
  fTopJet1(),
  fTopJet2(),
  fTopBJet(),
  fTopDijet_p4(),
  fTop_p4(), 
  fSelectedTopsJet1(),
  fSelectedTopsJet2(),
  fSelectedTopsBJet(),
  fSelectedTopsBDTG(),
  fNotSelectedTopsJet1(),
  fNotSelectedTopsJet2(),
  fNotSelectedTopsBJet(),
  fNotSelectedTopsBDTG(),
  fAllTopsJet1(),
  fAllTopsJet2(),
  fAllTopsBJet(),
  fAllTopsBDTG(),
  fSelectedCleanedTopsJet1(),
  fSelectedCleanedTopsJet2(),
  fSelectedCleanedTopsBJet(),
  fSelectedCleanedTopsBDTG(),
  fAllCleanedTopsJet1(),
  fAllCleanedTopsJet2(),
  fAllCleanedTopsBJet(),
  fAllCleanedTopsBDTG(),
  fTopTaggingScaleFactorEventWeight(1.0)
{ }

TopSelectionBDT::Data::~Data() { }


TopSelectionBDT::TopSelectionBDT(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& postfix)
  : BaseSelection(eventCounter, histoWrapper, commonPlots, postfix),
    // Input parameters
    cfg_NumberOfTopsCut(config, "NumberOfTopsCut"),
    cfg_TopBDTGCut(config     , "TopBDTGCut"),
    cfg_TopMassLowCut(config  , "TopMassLowCut"),
    cfg_TopMassUppCut(config  , "TopMassUppCut"),
    cfg_CSV_bDiscCut(config   , "CSV_bDiscCut"),
    // Event counter for passing selection
    cPassedTopSelectionBDT(fEventCounter.addCounter("passed top selection ("+postfix+")")),
    // Sub counters
    cSubAll(fEventCounter.addSubCounter("top selection ("+postfix+")", "All")),
    cSubPassedJetsCut(fEventCounter.addSubCounter("top selection ("+postfix+")", "#geq 3 jets")),
    cSubPassedBjetsCut(fEventCounter.addSubCounter("top selection ("+postfix+")", "#geq 1 b jets")),
    cSubPassedBDTGCut(fEventCounter.addSubCounter("top selection ("+postfix+")", "BDT")),
    cSubPassedNTopsCut(fEventCounter.addSubCounter("top selection ("+postfix+")", "# tops")),
    // Top candidates
    cTopsAll(fEventCounter.addSubCounter("top candidates ("+postfix+")", "All")),
    cTopsPassTopMassLowCut(fEventCounter.addSubCounter("top candidates ("+postfix+")", "top mass (low")),
    cTopsPassTopMassUppCut(fEventCounter.addSubCounter("top candidates ("+postfix+")", "top mass (upp)")),
    cTopsPassBDiscCut(fEventCounter.addSubCounter("top candidates ("+postfix+")", "b-disc")),
    cTopsPassBDTCut(fEventCounter.addSubCounter("top candidates ("+postfix+")", "BDTG")),
    cTopsPassCrossCleanCut(fEventCounter.addSubCounter("top candidates ("+postfix+")", "cross-clean")),
    fTopTagSFCalculator(config)
{
  initialize(config);
}

TopSelectionBDT::TopSelectionBDT(const ParameterSet& config)
: BaseSelection(),
  // Input parameters
  cfg_NumberOfTopsCut(config, "NumberOfTopsCut"),
  cfg_TopBDTGCut(config     , "TopBDTGCut"),
  cfg_TopMassLowCut(config  , "TopMassLowCut"),
  cfg_TopMassUppCut(config  , "TopMassUppCut"),
  cfg_CSV_bDiscCut(config   , "CSV_bDiscCut"),
  // Event counter for passing selection
  cPassedTopSelectionBDT(fEventCounter.addCounter("top")),
  // Sub counters
  cSubAll(fEventCounter.addSubCounter("top selection", "All")),
  cSubPassedJetsCut(fEventCounter.addSubCounter("top selection", "#geq 3 jets")),
  cSubPassedBjetsCut(fEventCounter.addSubCounter("top selection", "#geq 1 b jets")),
  cSubPassedBDTGCut(fEventCounter.addSubCounter("top selection", "BDTG")),
  cSubPassedNTopsCut(fEventCounter.addSubCounter("top selection", "# tops")),
  // Top candidates
  cTopsAll(fEventCounter.addSubCounter("top candidates", "All")),
  cTopsPassTopMassLowCut(fEventCounter.addSubCounter("top candidates", "top mass (low)")),
  cTopsPassTopMassUppCut(fEventCounter.addSubCounter("top candidates", "top mass (upp)")),
  cTopsPassBDiscCut(fEventCounter.addSubCounter("top candidates", "b-disc")),
  cTopsPassBDTCut(fEventCounter.addSubCounter("top candidates", "BDTG")),
  cTopsPassCrossCleanCut(fEventCounter.addSubCounter("top candidates", "cross-clean")),
  fTopTagSFCalculator(config)
{
  initialize(config);
  bookHistograms(new TDirectory());

}

TopSelectionBDT::~TopSelectionBDT() {
  
  // Histograms
  delete hTopMultiplicity_AllCandidates;
  delete hTopBDT_AllCandidates;
  delete hTopMass_AllCandidates;
  delete hTopPt_AllCandidates;

  delete hTopBDT_SelectedCandidates;
  delete hTopMass_SelectedCandidates;
  delete hTopPt_SelectedCandidates;
  delete hTopMultiplicity_SelectedCandidates;
  delete hTopBDT_SelectedCleanedCandidates;
  delete hTopMass_SelectedCleanedCandidates;
  delete hTopPt_SelectedCleanedCandidates;
  delete hTopMultiplicity_SelectedCleanedCandidates;

  delete hTopBDT_NotSelectedCandidates;
  delete hTopMass_NotSelectedCandidates;
  delete hTopPt_NotSelectedCandidates;
  delete hTopMultiplicity_NotSelectedCandidates;
  delete hTopBDT_AllCleanedCandidates;
  delete hTopMass_AllCleanedCandidates;
  delete hTopPt_AllCleanedCandidates;
  delete hTopMultiplicity_AllCleanedCandidates;

  delete hTopPt;
  delete hTopMass;
  delete hTopJet1Pt;
  delete hTopJet1Eta;
  delete hTopJet1BDisc;
  delete hTopJet2Pt;
  delete hTopJet2Eta;
  delete hTopJet2BDisc;
  delete hTopBJetPt;
  delete hTopBJetEta;
  delete hTopBJetBDisc;
  delete hTopDiJetPt;
  delete hTopDiJetEta;
  delete hTopDiJetMass;
  delete hTopDijetDeltaR;
  delete hTopMassWMassRatio;

  
  // TMVA reader
  delete reader;
  
}

void TopSelectionBDT::initialize(const ParameterSet& config) {
  
  // Load TMVA library
  TMVA::Tools::Instance();
    
  // Create the reader
  reader = new TMVA::Reader( "!Color:Silent" );
  
  // Add variables
  reader->AddVariable( "TrijetPtDR",              &TrijetPtDR              );
  reader->AddVariable( "TrijetDijetPtDR",         &TrijetDijetPtDR         );
  reader->AddVariable( "TrijetBjetMass",          &TrijetBjetMass          );
  reader->AddVariable( "TrijetLdgJetBDisc",       &TrijetLdgJetBDisc       );
  reader->AddVariable( "TrijetSubldgJetBDisc",    &TrijetSubldgJetBDisc    );
  reader->AddVariable( "TrijetBJetLdgJetMass",    &TrijetBJetLdgJetMass    );
  reader->AddVariable( "TrijetBJetSubldgJetMass", &TrijetBJetSubldgJetMass );
  reader->AddVariable( "TrijetMass",              &TrijetMass              );
  reader->AddVariable( "TrijetDijetMass",         &TrijetDijetMass         );
  reader->AddVariable( "TrijetBJetBDisc",         &TrijetBJetBDisc         );
  reader->AddVariable( "TrijetSoftDrop_n2",       &TrijetSoftDrop_n2       );
  reader->AddVariable( "TrijetLdgJetCvsL",        &TrijetLdgJetCvsL        );
  reader->AddVariable( "TrijetSubldgJetCvsL",     &TrijetSubldgJetCvsL     );
  reader->AddVariable( "TrijetLdgJetPtD",         &TrijetLdgJetPtD         );
  reader->AddVariable( "TrijetSubldgJetPtD",      &TrijetSubldgJetPtD      );
  reader->AddVariable( "TrijetLdgJetAxis2",       &TrijetLdgJetAxis2       );
  reader->AddVariable( "TrijetSubldgJetAxis2",    &TrijetSubldgJetAxis2    );
  reader->AddVariable( "TrijetLdgJetMult",        &TrijetLdgJetMult        );
  reader->AddVariable( "TrijetSubldgJetMult",     &TrijetSubldgJetMult     );

  // Construct the relative path toe the xml file with the BDT weights
  std::string relPath  = "../../../data/TopTaggerWeights/";
  std::string fileName = config.getParameter<std::string>("WeightFile");
  std::string fullPath = relPath + fileName;
  // std::cout << "Opening BDT weight file " << fullPath << std::endl;
  reader->BookMVA("BTDG method", fullPath);

  return;
}

void TopSelectionBDT::bookHistograms(TDirectory* dir) {
  // Fixed binning  
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

  // Histograms
  TDirectory* subdir = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "topSelectionBDT_"    + sPostfix);

  // All Candidates
  hTopBDT_AllCandidates           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopBDT_AllCandidates",";top candidate BDT", 40, -1.0, 1.0) ; 
  hTopMass_AllCandidates          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMass_AllCandidates" ,";top candidate M (GeVc^{-2})", nTopMassBins, fTopMassMin, fTopMassMax);
  hTopPt_AllCandidates            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopPt_AllCandidates", ";top candidate p_{T} (GeV/c)", nPtBins, fPtMin, fPtMax);
  hTopMultiplicity_AllCandidates  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMultiplicity_AllCandidates", ";top candidate multiplicity", 400, 0.0, 400.0);

  // Cross-cleaned Candidates
  hTopBDT_AllCleanedCandidates           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopBDT_AllCleanedCandidates",";top candidate BDT", 40, -1.0, 1.0) ; 
  hTopMass_AllCleanedCandidates          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMass_AllCleanedCandidates" ,";top candidate M (GeVc^{-2})", nTopMassBins, fTopMassMin, fTopMassMax);
  hTopPt_AllCleanedCandidates            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopPt_AllCleanedCandidates", ";top candidate p_{T} (GeV/c)", nPtBins, fPtMin, fPtMax);
  hTopMultiplicity_AllCleanedCandidates  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMultiplicity_AllCleanedCandidates", ";top candidate multiplicity", 10, 0.0, 10.0);

  // Selected Candidates
  hTopBDT_SelectedCandidates           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopBDT_SelectedCandidates",";top candidate BDT", 40, -1.0, 1.0) ; 
  hTopMass_SelectedCandidates          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMass_SelectedCandidates" ,";top candidate M (GeVc^{-2})", nTopMassBins, fTopMassMin, fTopMassMax);
  hTopPt_SelectedCandidates            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopPt_SelectedCandidates", ";top candidate p_{T} (GeV/c)", nPtBins, fPtMin, fPtMax);
  hTopMultiplicity_SelectedCandidates  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMultiplicity_SelectedCandidates", ";top candidate multiplicity", 50, 0.0, 50.0);

  // Selected & Cleaned Candidates
  hTopBDT_SelectedCleanedCandidates           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopBDT_SelectedCleanedCandidates",";top candidate BDT", 40, -1.0, 1.0) ; 
  hTopMass_SelectedCleanedCandidates          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMass_SelectedCleanedCandidates" ,";top candidate M (GeVc^{-2})", nTopMassBins, fTopMassMin, fTopMassMax);
  hTopPt_SelectedCleanedCandidates            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopPt_SelectedCleanedCandidates", ";top candidate p_{T} (GeV/c)", nPtBins, fPtMin, fPtMax);
  hTopMultiplicity_SelectedCleanedCandidates  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMultiplicity_SelectedCleanedCandidates", ";top candidate multiplicity", 10, 0.0, 10.0);

  // Candidates that were rejected (due to BDTG score)
  hTopBDT_NotSelectedCandidates           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopBDT_NotSelectedCandidates",";top candidate BDT", 40, -1.0, 1.0) ; 
  hTopMass_NotSelectedCandidates          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMass_NotSelectedCandidates" ,";top candidate M (GeVc^{-2})", nTopMassBins, fTopMassMin, fTopMassMax);
  hTopPt_NotSelectedCandidates            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopPt_NotSelectedCandidates", ";top candidate p_{T} (GeV/c)", nPtBins, fPtMin, fPtMax);
  hTopMultiplicity_NotSelectedCandidates  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMultiplicity_NotSelectedCandidates", ";top candidate multiplicity", 50, 0.0, 50.0);

  // Selected Top
  hTopPt             = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Top_Pt"        , ";p_{T} (GeV/c)", nPtBins     , fPtMin     , fPtMax);
  hTopMass           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Top_Mass"      , ";M (GeV/c^{2})", nTopMassBins, fTopMassMin, fTopMassMax);
  hTopJet1Pt         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Top_Jet1Pt"    , ";p_{T} (GeV/c)", nPtBins     , fPtMin     , fPtMax);
  hTopJet1Eta        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Top_Jet1Eta"   , ";#eta"         , nEtaBins    , fEtaMin    , fEtaMax);
  hTopJet1BDisc      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Top_Jet1BDisc" , ";b-tag discr." , nBDiscBins  , fBDiscMin  , fBDiscMax);
  hTopJet2Pt         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Top_Jet2Pt"    , ";p_{T} (GeV/c)", nPtBins     , fPtMin     , fPtMax);
  hTopJet2Eta        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Top_Jet2Eta"   , ";#eta"         , nEtaBins    , fEtaMin    , fEtaMax);
  hTopJet2BDisc      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Top_Jet2BDisc" , ";b-tag discr." , nBDiscBins  , fBDiscMin  , fBDiscMax);
  hTopBJetPt         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Top_BJetPt"    , ";p_{T} (GeV/c)", nPtBins     , fPtMin     , fPtMax);
  hTopBJetEta        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Top_BJetEta"   , ";#eta"         , nEtaBins    , fEtaMin    , fEtaMax);
  hTopBJetBDisc      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Top_BJetBDisc" , ";b-tag discr." , nBDiscBins  , fBDiscMin  , fBDiscMax);
  hTopDiJetPt        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Top_DiJetPt"   , ";p_{T} (GeV/c)", nPtBins     , fPtMin     , fPtMax);
  hTopDiJetEta       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Top_DiJetEta"  , ";#eta"         , nEtaBins    , fEtaMin    , fEtaMax);
  hTopDiJetMass      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Top_DiJetMass" , ";M (GeV/c^{2})", nWMassBins  , fWMassMin  , fWMassMax);
  hTopDijetDeltaR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Top_DijetDR"   ,  ";#Delta R(j_{1},j_{2})", 2*nDRBins, fDRMin, fDRMax);
  hTopMassWMassRatio = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Top_WMassRatio", ";R_{32}", 100 , 0.0, 10.0);

  fTopTagSFCalculator.bookHistograms(subdir, fHistoWrapper);
  return;
}


TopSelectionBDT::Data TopSelectionBDT::silentAnalyze(const Event& event, const JetSelection::Data& jetData, const BJetSelection::Data& bjetData) {
  ensureSilentAnalyzeAllowed(event.eventID());

  // Disable histogram filling and counter
  disableHistogramsAndCounters();
  Data myData = privateAnalyze(event, jetData.getSelectedJets(), bjetData.getSelectedBJets());
  enableHistogramsAndCounters();
  return myData;
}


TopSelectionBDT::Data TopSelectionBDT::analyze(const Event& event, const JetSelection::Data& jetData, const BJetSelection::Data& bjetData) {
  ensureAnalyzeAllowed(event.eventID());

  // Ready to analyze
  TopSelectionBDT::Data data = privateAnalyze(event, jetData.getSelectedJets(), bjetData.getSelectedBJets());

  // Send data to CommonPlots
  // if (fCommonPlots != nullptr) fCommonPlots->fillControlPlotsAtTopSelectionBDT(event, data);// fixme - implement
  return data;
}


TopSelectionBDT::Data TopSelectionBDT::privateAnalyze(const Event& event, const std::vector<Jet> selectedJets, const std::vector<Jet> selectedBjets) {
  Data output;
  cSubAll.increment();

  // Initialise variables
  output.fJetsUsedAsBJets = selectedBjets;
  output.bPassedSelection = true;// iro

  // Sanity check
  if (selectedJets.size() < 3) return output;
  cSubPassedJetsCut.increment();

  if (selectedBjets.size() < 1) return output;
  cSubPassedBjetsCut.increment();

  // Only resize if their size exceeds max allowed value
  std::vector<Jet> jets  = selectedJets;
  std::vector<Jet> bjets = selectedBjets;
  
  //================================================================================================  
  // Top Candidates
  //================================================================================================  
  if (0) std::cout << "=== TopSelectionBDT:: Obtaining all top candidates" << std::endl;
  TrijetSelection fAllTops;
  TrijetSelection fSelectedTops;
  TrijetSelection fSelectedCleanedTops;
  TrijetSelection fNotSelectedTops;
  TrijetSelection fAllCleanedTops;

  // For-loop: All b-jets
  for (auto& bjet: bjets)
    {
      int index1 = 0;

      // For-loop: All jets
      for (auto& jet1: jets)
	{
	  index1++;
	  int index2 = 0;

	  // Skip if jet1 is same as bjet
	  if (areSameJets(jet1, bjet)) continue;

	  // For-loop: All jets
	  for (auto& jet2: jets)
	    {
	      index2++;

	      // Do not consider duplicate compinations
	      if (index2 < index1) continue;

	      // Skip if jet2 is same as jet1, or jet2 same as bjet
	      if (areSameJets(jet2,  jet1) || areSameJets(jet2,  bjet)) continue;

	      // Increment top candidate counter
	      cTopsAll.increment();

	      // Get 4-momentum of top (trijet) and W (dijet)
	      math::XYZTLorentzVector top_p4, w_p4;
	      top_p4 = bjet.p4() + jet1.p4() + jet2.p4();
	      w_p4   = jet1.p4() + jet2.p4();
	      
	      // Skip trijet combinations which do not fulfil the top invariant mass (lower threshold)
	      if (!cfg_TopMassLowCut.passedCut(top_p4.M())) continue;
	      cTopsPassTopMassLowCut.increment();

	      // Skip trijet combinations which do not fulfil the top invariant mass (upper threshold)
	      if (!cfg_TopMassUppCut.passedCut(top_p4.M())) continue;
	      cTopsPassTopMassUppCut.increment();

	      // Skip trijet combinations which do not fulfil bjet_CSV threshold
	      if (!cfg_CSV_bDiscCut.passedCut(bjet.bjetDiscriminator())) continue;
	      cTopsPassBDiscCut.increment();

	      // Calculate variables
	      double dr_sd = ROOT::Math::VectorUtil::DeltaR( jet1.p4(), jet2.p4());
	      double softDrop_n2 = min(jet2.pt(), jet1.pt()) / ( (jet2.pt() + jet1.pt()) * dr_sd * dr_sd);

	      // Calculate our BDT discriminating variables for BDTG use
	      TrijetPtDR              = top_p4.Pt() * ROOT::Math::VectorUtil::DeltaR( w_p4  , bjet.p4() );
	      TrijetDijetPtDR         = w_p4.Pt() * ROOT::Math::VectorUtil::DeltaR( jet1.p4() , jet2.p4() );
	      TrijetBjetMass          = bjet.p4().M();
	      TrijetLdgJetBDisc       = jet1.bjetDiscriminator();
	      TrijetSubldgJetBDisc    = jet2.bjetDiscriminator();
	      TrijetBJetLdgJetMass    = (bjet.p4() + jet1.p4()).M();
	      TrijetBJetSubldgJetMass = (bjet.p4() + jet2.p4()).M();
	      TrijetMass              = top_p4.M();
	      TrijetDijetMass         = w_p4.M();
	      TrijetBJetBDisc         = bjet.bjetDiscriminator();
	      TrijetSoftDrop_n2       = softDrop_n2;
	      TrijetLdgJetCvsL        = jet1.pfCombinedCvsLJetTags();
	      TrijetSubldgJetCvsL     = jet2.pfCombinedCvsLJetTags();
	      TrijetLdgJetPtD         = jet1.QGTaggerAK4PFCHSptD();
	      TrijetSubldgJetPtD      = jet2.QGTaggerAK4PFCHSptD();
	      TrijetLdgJetAxis2       = jet1.QGTaggerAK4PFCHSaxis2();
	      TrijetSubldgJetAxis2    = jet2.QGTaggerAK4PFCHSaxis2();
	      TrijetLdgJetMult        = jet1.QGTaggerAK4PFCHSmult();
	      TrijetSubldgJetMult     = jet2.QGTaggerAK4PFCHSmult();

	      // Evaluate the MVA discriminator value
	      float BDTGoutput = reader->EvaluateMVA("BTDG method");
	      // std::cout << "MVA = " << MVAoutput << ", Pt = " << top_p4.pt() << ", M = " << TrijetMass << std::endl;

	      // Fill top candidate BDT values
	      hTopBDT_AllCandidates -> Fill(BDTGoutput);
	      hTopMass_AllCandidates-> Fill(top_p4.M());
	      hTopPt_AllCandidates  -> Fill(top_p4.pt());

	      // Save top candidates
	      fAllTops.BDTG.push_back(BDTGoutput);
	      fAllTops.TrijetP4.push_back(top_p4);
	      fAllTops.DijetP4.push_back(w_p4);
	      fAllTops.Jet1.push_back(getLeadingSubleadingJet(jet1, jet2, "leading"));
	      fAllTops.Jet2.push_back(getLeadingSubleadingJet(jet1, jet2, "subleading"));
	      fAllTops.BJet.push_back(bjet);
	      fAllTops.isGenuine.push_back(false);
	      fAllTops.isTagged.push_back(cfg_TopBDTGCut.passedCut(BDTGoutput));

	      // Get top candidates above BDTG cut
	      if (cfg_TopBDTGCut.passedCut(BDTGoutput))
		{
		  cTopsPassBDTCut.increment();
		  fSelectedTops.BDTG.push_back(BDTGoutput);
		  fSelectedTops.TrijetP4.push_back(top_p4);
		  fSelectedTops.DijetP4.push_back(w_p4);
		  fSelectedTops.Jet1.push_back(getLeadingSubleadingJet(jet1, jet2, "leading"));
		  fSelectedTops.Jet2.push_back(getLeadingSubleadingJet(jet1, jet2, "subleading"));
		  fSelectedTops.BJet.push_back(bjet);
		  fSelectedTops.isGenuine.push_back(false);
		  fSelectedTops.isTagged.push_back(true);  // fixme: which BDTG cut? ldg, or subldg, or?
		}
	      else
		{
		  // Get top candidates failing BDTG cut
		  fNotSelectedTops.BDTG.push_back(BDTGoutput);
		  fNotSelectedTops.TrijetP4.push_back(top_p4);
		  fNotSelectedTops.DijetP4.push_back(w_p4);
		  fNotSelectedTops.Jet1.push_back(getLeadingSubleadingJet(jet1, jet2, "leading"));
		  fNotSelectedTops.Jet2.push_back(getLeadingSubleadingJet(jet1, jet2, "subleading"));
		  fNotSelectedTops.BJet.push_back(bjet);
		  fNotSelectedTops.isGenuine.push_back(false); // fixme
		  fNotSelectedTops.isTagged.push_back(false); // fixme: which BDTG cut? ldg, or subldg, or?
		}
		
	    }// For-loop: All jets
	}// For-loop: All jets
    }// For-loop: All b-jets
  hTopMultiplicity_AllCandidates->Fill(fAllTops.BDTG.size());

  //================================================================================================  
  // Sort top candidates in descending BDTG values
  //================================================================================================    
  if (0) std::cout << "=== TopSelectionBDT::Sort top candidates in BDTG" << std::endl;
  fAllTops         = SortInBDTGvalue(fAllTops);
  fSelectedTops    = SortInBDTGvalue(fSelectedTops);
  fNotSelectedTops = SortInBDTGvalue(fNotSelectedTops);
  
  //Debug: Check cleaned top candidates
  bool printCleanedTops = false;

  if (printCleanedTops)
    {
      // std::cout << "\n" << std::endl;
      // std::cout << "Entry: "<< cSubAll.value() << std::endl;
      // std::cout << std::string(10*10, '=') << std::endl;
      // std::cout << std::setw(12) << "Jet1 indx "  << std::setw(12) << "Jet2 indx"<< std::setw(12) << "Bjet indx"   << std::setw(12) << "BDT"   << std::setw(12) << "is cleaned"  << std::endl;
      // std::cout << std::string(10*10, '=') << std::endl;
      std::cout << "\nEntry: "<< cSubAll.value() << ", selected tops = " << fSelectedTops.BDTG.size() << ", jets = " << selectedJets.size() << ", bjets = " << selectedBjets.size() << std::endl;
    }

  // For-loop: All (BDTG sorted) selected tops
  for (size_t i = 0; i < fSelectedTops.BDTG.size(); i++)
    {
      // Fill Histos (selected candidates)
      hTopBDT_SelectedCandidates -> Fill(fSelectedTops.BDTG.at(i) );
      hTopMass_SelectedCandidates-> Fill(fSelectedTops.TrijetP4.at(i).M());
      hTopPt_SelectedCandidates  -> Fill(fSelectedTops.TrijetP4.at(i).pt());
      
      // Is this top cross-cleaned (b-jet availability check, object sharing with higher BDT value)
      bool isCrossCleaned = TopIsCrossCleaned(i, fSelectedTops);

      if (printCleanedTops)
	{
	  // std::cout << std::setw(12) << fSelectedTops.Jet1.at(i).index() << std::setw(12) << fSelectedTops.Jet2.at(i).index() << std::setw(12) << fSelectedTops.BJet.at(i).index() << std::setw(12) << fSelectedTops.BDTG.at(i) << std::setw(12) << isCrossCleaned <<  std::endl;
	  std::cout << i+1 << ") jet1 = " << fSelectedTops.Jet1.at(i).index() << ", jet2 = " << fSelectedTops.Jet2.at(i).index() 
		    << ", bjet = " << fSelectedTops.BJet.at(i).index() << ", BDTG = " << fSelectedTops.BDTG.at(i)
		    << ", X-Cleaned = " << isCrossCleaned <<  std::endl;
	}

      

      // If top candidate is not cross-cleaned skip to the next one
      if (!isCrossCleaned) continue;
      cTopsPassCrossCleanCut.increment();

      // Save the top selected && cross-cleaned top candidate (sorted in BDTG)
      fSelectedCleanedTops.Jet1.push_back( fSelectedTops.Jet1.at(i) );
      fSelectedCleanedTops.Jet2.push_back( fSelectedTops.Jet2.at(i) );
      fSelectedCleanedTops.BJet.push_back( fSelectedTops.BJet.at(i) );
      fSelectedCleanedTops.BDTG.push_back( fSelectedTops.BDTG.at(i) );
      fSelectedCleanedTops.TrijetP4.push_back( fSelectedTops.TrijetP4.at(i) );
      fSelectedCleanedTops.DijetP4.push_back( fSelectedTops.DijetP4.at(i) );
      fSelectedCleanedTops.isGenuine.push_back( fSelectedTops.isGenuine.at(i) ); // fixme: always false, overwite later.
      fSelectedCleanedTops.isTagged.push_back( fSelectedTops.isTagged.at(i) );   // fixme: which BDTG cut? ldg, or subldg, or?

      // Fill Histos (selected, cross-cleaned candidates)
      hTopBDT_SelectedCleanedCandidates -> Fill(fSelectedTops.BDTG.at(i) );
      hTopMass_SelectedCleanedCandidates-> Fill(fSelectedTops.TrijetP4.at(i).M());
      hTopPt_SelectedCleanedCandidates  -> Fill(fSelectedTops.TrijetP4.at(i).pt());
    }


  // For-loop: All candidates
  for (size_t i = 0; i < fAllTops.BDTG.size(); i++)
    {      
      // Is this top cross-cleaned (b-jet availability check, object sharing with higher BDT value)
      bool isCrossCleaned = TopIsCrossCleaned(i, fAllTops);
      printCleanedTops = false;
      
      if (printCleanedTops){
	if (i==0){
	  std::cout << "\n" << std::endl;
	  std::cout<<"Entry: "<<cSubAll.value()<<std::endl;
	  std::cout << std::string(10*10, '=') << std::endl;
	  std::cout << std::setw(12) << "Jet1 indx "  << std::setw(12) << "Jet2 indx"<< std::setw(12) << "Bjet indx"   << std::setw(12) << "BDT"   << std::setw(12) << "is cleaned"  << std::endl;
	  std::cout << std::string(10*10, '=') << std::endl;	
	}
	std::cout << std::setw(12) << fAllTops.Jet1.at(i).index()      << std::setw(12)   << fAllTops.Jet2.at(i).index()   << std::setw(12)   << fAllTops.BJet.at(i).index()
		  << std::setw(12) << fAllTops.BDTG.at(i)               << std::setw(12)   << isCrossCleaned
		  <<  std::endl;
	
      }
      if (!isCrossCleaned) continue;

      // Save the top selected && cross-cleaned top candidate
      fAllCleanedTops.Jet1.push_back( fAllTops.Jet1.at(i) );
      fAllCleanedTops.Jet2.push_back( fAllTops.Jet2.at(i) );
      fAllCleanedTops.BJet.push_back( fAllTops.BJet.at(i) );
      fAllCleanedTops.BDTG.push_back( fAllTops.BDTG.at(i) );
      fAllCleanedTops.TrijetP4.push_back( fAllTops.TrijetP4.at(i) );
      fAllCleanedTops.DijetP4.push_back( fAllTops.DijetP4.at(i) );
      fAllCleanedTops.isGenuine.push_back( fAllTops.isGenuine.at(i) ); // fixme: always false, overwite later.
      fAllCleanedTops.isTagged.push_back( fAllTops.isTagged.at(i) );   // fixme: which BDTG cut? ldg, or subldg, or?

      // Fill Histos (selected, cross-cleaned candidates)
      hTopBDT_AllCleanedCandidates -> Fill(fAllTops.BDTG.at(i) );
      hTopMass_AllCleanedCandidates-> Fill(fAllTops.TrijetP4.at(i).M());
      hTopPt_AllCleanedCandidates  -> Fill(fAllTops.TrijetP4.at(i).pt());
    }


  // For-loop: Failed candidates
  for (size_t i = 0; i < fNotSelectedTops.BDTG.size(); i++)
    {
      // Fill Histos (selected candidates)
      hTopBDT_NotSelectedCandidates -> Fill(fNotSelectedTops.BDTG.at(i) );
      hTopMass_NotSelectedCandidates-> Fill(fNotSelectedTops.TrijetP4.at(i).M());
      hTopPt_NotSelectedCandidates  -> Fill(fNotSelectedTops.TrijetP4.at(i).pt());
    }

  // Fill multiplicities
  hTopMultiplicity_SelectedCandidates->Fill(fSelectedTops.BDTG.size());
  hTopMultiplicity_SelectedCleanedCandidates->Fill(fSelectedCleanedTops.BDTG.size());
  hTopMultiplicity_NotSelectedCandidates->Fill(fNotSelectedTops.BDTG.size());
  hTopMultiplicity_AllCleanedCandidates->Fill(fAllCleanedTops.BDTG.size());

  //================================================================================================
  // Find top decay products
  //================================================================================================
  std::vector<genParticle> GenTops, GenTops_BQuark, GenTops_SubldgQuark, GenTops_LdgQuark;
  std::vector <Jet>        MCtrue_LdgJet, MCtrue_SubldgJet, MCtrue_Bjet, MC_BJets;
  std::vector <double>     dRminB;
  const double twoSigmaDpt = 0.32, dRcut = 0.4;

  if (event.isMC())
    {
    GenTops = GetGenParticles(event.genparticles().getGenParticles(), 6);
    }
  
  // For-loop: All top quarks
  for (auto& top: GenTops){
    std::vector<genParticle> quarks, bquark;
    // For-loop: Top quark daughters (Nested)
    for (size_t i=0; i<top.daughters().size(); i++)
      {
        genParticle dau = event.genparticles().getGenParticles()[top.daughters().at(i)];
        // B-Quark, W-Boson
        if (std::abs(dau.pdgId()) ==  5) bquark.push_back(dau);
	if (std::abs(dau.pdgId()) == 24) quarks = GetWpartons(dau, event);
      }//Top-quark daughters

    // Skip top if b-quark is not found (i.e. top decays to W and c)
    if (bquark.size() < 1) continue;

    // Skip top if at least one W decays leptonically
    if (quarks.size() < 2) continue;

    // Fill vectors for b-quarks, leading and subleading quarks coming from tops
    GenTops_BQuark.push_back(bquark.at(0));
    GenTops_LdgQuark.push_back(getLeadingSubleadingParton(quarks.at(0),    quarks.at(1), "leading"));
    GenTops_SubldgQuark.push_back(getLeadingSubleadingParton(quarks.at(0), quarks.at(1), "subleading"));
  }
  
  // Skip matcing if top does not decay to b
  bool doMatching = (GenTops_BQuark.size() == GenTops.size());	    

  //================================================================================================
  // Find truth-matched top candidates
  //================================================================================================
  if (doMatching)
    {
      // For-loop: All top quarks (for b-jet matching)
      for (size_t i=0; i<GenTops.size(); i++)
        {
          genParticle BQuark = GenTops_BQuark.at(i);
          Jet mcMatched_BJet;
          double dRmin  = 99999.9;

          // For-loop: All selected jets
          for (auto& bjet: jets)
            {
              double dR  = ROOT::Math::VectorUtil::DeltaR( bjet.p4(), BQuark.p4());
              double dPtOverPt = std::abs((bjet.pt() - BQuark.pt())/BQuark.pt());

	      // Find minimum dR
              if (dR > dRmin) continue;

              //Skip if dPtOverPt > twoSigmaDpt/Pt
              if (dPtOverPt > twoSigmaDpt) continue;

              // Store values
              dRmin  = dR;
              mcMatched_BJet = bjet;
            }// For-loop: All selected jets

          // Store match
          dRminB.push_back(dRmin);
          MC_BJets.push_back(mcMatched_BJet);
        }// For-loop: All top quarks    
      
      // For-loop: All top quarks (for dijet matching)
      for (size_t i=0; i<GenTops.size(); i++)
	{
	  genParticle top         = GenTops.at(i);
	  genParticle LdgQuark    = GenTops_LdgQuark.at(i);
	  genParticle SubldgQuark = GenTops_SubldgQuark.at(i);
	  Jet mcMatched_LdgJet, mcMatched_SubldgJet;	  
	  double dR1min, dR2min;
	  dR1min = dR2min = 99999.9;
	  
	  // For-loop: All selected jets
	  for (auto& jet: jets)
	    {
	      bool same = false;
	      
	      // For-loop: All top-quarks (Skip the jets that are matched to bquarks)
	      for (size_t k=0; k<GenTops.size(); k++)
		{
		  if (dRminB.at(k) > dRcut) continue;
		  if (areSameJets(jet,MC_BJets.at(k))) same = true;
		}// For-loop: All top-quarks
	      
	      if (same) continue;

	      // Find dR for the two jets in top-decay dijet
	      double dR1 = ROOT::Math::VectorUtil::DeltaR(jet.p4(), LdgQuark.p4());
	      double dR2 = ROOT::Math::VectorUtil::DeltaR(jet.p4(), SubldgQuark.p4());

	      // Calculate dPtOverPt for each jet in top-decay dijet
	      double dPtOverPt1 = std::abs((jet.pt() - LdgQuark.pt())/LdgQuark.pt());
	      double dPtOverPt2 = std::abs((jet.pt() - SubldgQuark.pt())/SubldgQuark.pt());
	      
	      // Find which of the two is the correct match                                                     
	      if (dR1 < dR2)
		{
		  // Is Jet1 closer in eta-phi AND has smaller pT difference?
		  if (dR1 < dR1min)
		    {
		      if (dPtOverPt1 < twoSigmaDpt)
			{
			  dR1min = dR1;
			  mcMatched_LdgJet = jet;
			}
		    }		  
		  // Is Jet2 closer in eta-phi AND has smaller pT difference?
		  else if (dR2 < dR2min)
		    {
		      if (dPtOverPt2 < twoSigmaDpt)
			{
			  dR2min  = dR2;
			  mcMatched_SubldgJet = jet;
			}
		    }
		}
	      else
		{
		  // Is Jet2 closer in eta-phi AND has smaller pT difference?
		  if (dR2 < dR2min)
		    {
		      if (dPtOverPt2 < twoSigmaDpt)
			{
			  dR2min  = dR2;
			  mcMatched_SubldgJet = jet;
			}
		    }		  
		  // Is Jet2 closer in eta-phi AND has smaller pT difference?
		  else if (dR1 < dR1min)
		    {
		      if  (dPtOverPt1 < twoSigmaDpt)
			{
			  dR1min  = dR1;
			  mcMatched_LdgJet = jet;
			}
		    }
		}
	    }//For-loop: All selected jets
	  
	  // Check if top-tagged object is genuine (fully truth-matched)
	  bool isGenuine = (dR1min<= dRcut && dR2min <= dRcut && dRminB.at(i) <= dRcut);

	  if (isGenuine)
	    {
	      // Store genuine top candidates
	      MCtrue_LdgJet.push_back(mcMatched_LdgJet);
	      MCtrue_SubldgJet.push_back(mcMatched_SubldgJet);
	      MCtrue_Bjet.push_back(MC_BJets.at(i));
	    }
	}
    }

  //================================================================================================
  // Calculate the event top-tagging Scale Factor (SF)
  //================================================================================================
  if (event.isMC())
    {

      // First of all, update the mc-truth info
      for (size_t i = 0; i < fAllCleanedTops.BDTG.size(); i++)
	{      
	  bool isFullyMatched = IsGenuineTop(fAllCleanedTops.Jet1.at(i), fAllCleanedTops.Jet2.at(i), fAllCleanedTops.BJet.at(i), MCtrue_LdgJet, MCtrue_SubldgJet, MCtrue_Bjet);
	  fAllCleanedTops.isGenuine.at(i) = isFullyMatched;
	}

      // Calculate and store b-jet scale factor weight and it's uncertainty
      output.fTopTaggingScaleFactorEventWeight = fTopTagSFCalculator.calculateSF(fAllCleanedTops.TrijetP4, fAllCleanedTops.BDTG, fAllCleanedTops.isTagged, fAllCleanedTops.isGenuine);
    }

    
  //================================================================================================
  // Fill output data
  //================================================================================================
  if (0) std::cout << "=== TopSelectionBDT:: Fill output" << std::endl;
  if (fAllCleanedTops.BDTG.size() > 0) 
    {
      
      output.bPassedBDTGCut = cfg_TopBDTGCut.passedCut( fAllCleanedTops.BDTG.at(0) );
      // Leading-in-BDTG top
      output.fTopBDTG      = fAllCleanedTops.BDTG.at(0);
      output.fTopJet1     = fAllCleanedTops.Jet1.at(0);
      output.fTopJet2     = fAllCleanedTops.Jet2.at(0);
      output.fTopBJet     = fAllCleanedTops.BJet.at(0);
      output.fTopDijet_p4 = fAllCleanedTops.Jet1.at(0).p4() + fAllCleanedTops.Jet2.at(0).p4();
      output.fTop_p4      = fAllCleanedTops.TrijetP4.at(0);      
    }

  // Fill in remaining data
  for (size_t i = 0; i < fSelectedTops.BDTG.size(); i++)
    {
      output.fSelectedTopsJet1.push_back(fSelectedTops.Jet1.at(i));
      output.fSelectedTopsJet2.push_back(fSelectedTops.Jet2.at(i));
      output.fSelectedTopsBJet.push_back(fSelectedTops.BJet.at(i));
      output.fSelectedTopsBDTG.push_back(fSelectedTops.BDTG.at(i));
    }

  for (size_t i = 0; i < fNotSelectedTops.BDTG.size(); i++)
    {
      output.fNotSelectedTopsJet1.push_back(fNotSelectedTops.Jet1.at(i));
      output.fNotSelectedTopsJet2.push_back(fNotSelectedTops.Jet2.at(i));
      output.fNotSelectedTopsBJet.push_back(fNotSelectedTops.BJet.at(i));
      output.fNotSelectedTopsBDTG.push_back(fNotSelectedTops.BDTG.at(i));
    }

  for (size_t i = 0; i < fAllTops.BDTG.size(); i++)
    {
      output.fAllTopsJet1.push_back(fAllTops.Jet1.at(i));
      output.fAllTopsJet2.push_back(fAllTops.Jet2.at(i));
      output.fAllTopsBJet.push_back(fAllTops.BJet.at(i));
      output.fAllTopsBDTG.push_back(fAllTops.BDTG.at(i));
    }

  for (size_t i = 0; i < fSelectedCleanedTops.BDTG.size(); i++)
    {
      output.fSelectedCleanedTopsJet1.push_back(fSelectedCleanedTops.Jet1.at(i));
      output.fSelectedCleanedTopsJet2.push_back(fSelectedCleanedTops.Jet2.at(i));
      output.fSelectedCleanedTopsBJet.push_back(fSelectedCleanedTops.BJet.at(i));
      output.fSelectedCleanedTopsBDTG.push_back(fSelectedCleanedTops.BDTG.at(i));
    }

  for (size_t i = 0; i < fAllCleanedTops.BDTG.size(); i++)
    {
      output.fAllCleanedTopsJet1.push_back(fAllCleanedTops.Jet1.at(i));
      output.fAllCleanedTopsJet2.push_back(fAllCleanedTops.Jet2.at(i));
      output.fAllCleanedTopsBJet.push_back(fAllCleanedTops.BJet.at(i));
      output.fAllCleanedTopsBDTG.push_back(fAllCleanedTops.BDTG.at(i));
    }

   //===============================================================================================
  // Increment counters
  //================================================================================================
  output.bPassedNTopsCut  = cfg_NumberOfTopsCut.passedCut(output.fSelectedCleanedTopsBDTG.size());
  output.bPassedSelection = output.bPassedBDTGCut * output.bPassedNTopsCut;

  if (0) std::cout << "=== TopSelectionBDT:: Increment counters" << std::endl;
  //=== Apply cut on leading BDTG top (at least 1 top with desired BDTG score)
  if (output.bPassedBDTGCut) cSubPassedBDTGCut.increment();
  
  //=== Apply cut on number of tops with BDTG score above BDTG score threshold
  if (output.bPassedNTopsCut) cSubPassedNTopsCut.increment();

  //================================================================================================
  // Fill histograms
  //================================================================================================
  if (0) std::cout << "=== TopSelectionBDT:: Filling histograms" << std::endl;
  if (output.bPassedSelection) cPassedTopSelectionBDT.increment();

  // Get the top candidate with highest BDTG score
  TrijetSelection myTops = fAllCleanedTops;
  SelectedTrijets top;
  if (myTops.Jet1.size() > 0)
    {
      top.Jet1      = getLeadingSubleadingJet(myTops.Jet1.at(0), myTops.Jet2.at(0), "leading");
      top.Jet2      = getLeadingSubleadingJet(myTops.Jet1.at(0), myTops.Jet2.at(0), "subleading");
      top.BJet      = myTops.BJet.at(0);
      top.DijetP4   = myTops.DijetP4.at(0);
      top.TrijetP4  = myTops.TrijetP4.at(0);
      top.BDTG       = myTops.BDTG.at(0);
      top.isGenuine = myTops.isGenuine.at(0);
      top.isTagged  = myTops.isTagged.at(0);
      
      double dijetMass = (top.Jet1.p4() +  top.Jet2.p4()).M();
      hTopPt          -> Fill(top.TrijetP4.Pt());
      hTopMass        -> Fill(top.TrijetP4.M());
      hTopJet1Pt      -> Fill(top.Jet1.pt());
      hTopJet1Eta     -> Fill(top.Jet1.eta());
      hTopJet1BDisc   -> Fill(top.Jet1.bjetDiscriminator());
      hTopJet2Pt      -> Fill(top.Jet2.pt());
      hTopJet2Eta     -> Fill(top.Jet2.eta());
      hTopJet2BDisc   -> Fill(top.Jet2.bjetDiscriminator());
      hTopBJetPt      -> Fill(top.BJet.pt());
      hTopBJetEta     -> Fill(top.BJet.eta());
      hTopBJetBDisc   -> Fill(top.BJet.bjetDiscriminator());
      hTopDiJetPt     -> Fill(top.DijetP4.Pt());
      hTopDiJetEta    -> Fill(top.DijetP4.Eta());
      hTopDiJetMass   -> Fill(top.DijetP4.M());
      hTopMassWMassRatio -> Fill(top.TrijetP4.M()/dijetMass);
      hTopDijetDeltaR -> Fill(ROOT::Math::VectorUtil::DeltaR(top.Jet1.p4(), top.Jet2.p4()));
      // double Top_Rapidity    = 0.5*log((top.TrijetP4.E() + top.TrijetP4.Pz())/(top.TrijetP4.E() - top.TrijetP4.Pz()));
    }
      
  return output;
}

bool TopSelectionBDT::areSameJets(const Jet& jet1, const Jet& jet2) {
  float dR = ROOT::Math::VectorUtil::DeltaR(jet1.p4(), jet2.p4());
  float dR_match = 0.1;
  if (dR <= dR_match) return true;
  else return false;
}

bool TopSelectionBDT::isBJet(const Jet& jet, const std::vector<Jet>& bjets) {
  for (auto bjet: bjets)
    {
      if (areSameJets(jet, bjet)) return true;
    }
  return false;
}

Jet TopSelectionBDT::getLeadingSubleadingJet(const Jet& jet0, const Jet& jet1, string selectedJet){
  if (selectedJet != "leading" && selectedJet!="subleading") std::cout<<"WARNING! Unknown option "<<selectedJet<<". Function getLeadingSubleadingJet returns leading Jet"<<std::endl;
  Jet leadingJet, subleadingJet;
  if (jet0.pt() > jet1.pt()){                                                                                                   
    leadingJet    = jet0;                  
    subleadingJet = jet1;      
  }           
  else{                         
    leadingJet    = jet1;                                          
    subleadingJet = jet0;
  }
  if (selectedJet == "subleading") return subleadingJet;
  return leadingJet;
}

genParticle TopSelectionBDT::getLeadingSubleadingParton(const genParticle& quark0, const genParticle& quark1, string selectedParton){
  if (selectedParton != "leading" && selectedParton!="subleading") std::cout<<"WARNING! Unknown option "<<selectedParton<<". Function getLeadingSubleadingParton returns leading Parton"<<std::endl;
  genParticle leadingParton, subleadingParton;
  if (quark0.pt() > quark1.pt()){                                                                                                   
    leadingParton    = quark0;                  
    subleadingParton = quark1;      
  }           
  else{                         
    leadingParton    = quark1;                                          
    subleadingParton = quark0;
  }
  if (selectedParton == "subleading") return subleadingParton;
  return leadingParton;
}


bool TopSelectionBDT::isMatchedJet(const Jet& jet, const TrijetSelection& myTops, const unsigned int index) {
  // Sanity check: If index is bigger than the number of tops that top does not exist => not possible to be matched
  if ( index+1 > myTops.BDTG.size()) return false;

  std::vector<Jet> jets;
  jets.push_back(myTops.Jet1.at(index));
  jets.push_back(myTops.Jet2.at(index));
  jets.push_back(myTops.BJet.at(index));
  
  for (auto iJet: jets)
    {
      if (areSameJets(jet, iJet)) return true;
    }
  return false;
}

TrijetSelection TopSelectionBDT::SortInBDTGvalue(TrijetSelection TopCand){
  //  Description
  //  Takes as input a collection of Top Candidates and returns the collection sorted in BDT(BDTG) value
  size_t size = TopCand.BDTG.size();
  if (size < 1) return TopCand;

  for (size_t i=0; i<size-1; i++)
    {
      for  (size_t j=i+1; j<size; j++)
	{
	  Jet Jet1_i   = TopCand.Jet1.at(i);
	  Jet Jet2_i   = TopCand.Jet2.at(i);
	  Jet BJet_i   = TopCand.BJet.at(i);
	  double mva_i = TopCand.BDTG.at(i);
	  math::XYZTLorentzVector TrijetP4_i = TopCand.TrijetP4.at(i);
	  math::XYZTLorentzVector DijetP4_i  = TopCand.DijetP4.at(i);
	  bool isGenuine_i = TopCand.isGenuine.at(i);
          bool isTagged_i  = TopCand.isTagged.at(i);

	  Jet Jet1_j   = TopCand.Jet1.at(j);
	  Jet Jet2_j   = TopCand.Jet2.at(j);
	  Jet BJet_j   = TopCand.BJet.at(j);
	  double mva_j = TopCand.BDTG.at(j);
	  math::XYZTLorentzVector TrijetP4_j = TopCand.TrijetP4.at(j);
	  math::XYZTLorentzVector DijetP4_j  = TopCand.DijetP4.at(j);
	  bool isGenuine_j = TopCand.isGenuine.at(j);
          bool isTagged_j  = TopCand.isTagged.at(j);

	  if (mva_i >= mva_j) continue;
	  TopCand.Jet1.at(i) = Jet1_j;
	  TopCand.Jet2.at(i) = Jet2_j;
	  TopCand.BJet.at(i) = BJet_j;
	  TopCand.BDTG.at(i)  = mva_j;
	  TopCand.TrijetP4.at(i)  = TrijetP4_j;
	  TopCand.DijetP4.at(i)   = DijetP4_j;
	  TopCand.isGenuine.at(i) = isGenuine_j;
          TopCand.isTagged.at(i)  = isTagged_j;

	  TopCand.Jet1.at(j)      = Jet1_i;
	  TopCand.Jet2.at(j)      = Jet2_i;
	  TopCand.BJet.at(j)      = BJet_i;
	  TopCand.BDTG.at(j)       = mva_i;
	  TopCand.TrijetP4.at(j)  = TrijetP4_i;
	  TopCand.DijetP4.at(j)   = DijetP4_i;
	  TopCand.isGenuine.at(j) = isGenuine_i;
          TopCand.isTagged.at(j)  = isTagged_i;
	}
    }
  return TopCand;
}


SelectedTrijets TopSelectionBDT::GetSelectedTopCandidate(TrijetSelection TopCand, int index){
  SelectedTrijets trijet;
  trijet.Jet1      = TopCand.Jet1.at(index);
  trijet.Jet2      = TopCand.Jet2.at(index);
  trijet.BJet      = TopCand.BJet.at(index);
  trijet.BDTG      = TopCand.BDTG.at(index);
  trijet.TrijetP4  = TopCand.TrijetP4.at(index);
  trijet.DijetP4   = TopCand.DijetP4.at(index);
  trijet.isGenuine = TopCand.isGenuine.at(index);
  trijet.isTagged  = TopCand.isTagged.at(index);
  return trijet;
}

bool TopSelectionBDT::TopIsCrossCleaned(int Index, TrijetSelection TopCand){
  // Description:
  // Used to find the cross-cleaned trijet multiplicity. The function takes as input the index of a trijet and the total trijet collection and returns 
  // false true if If at least one of the subjets of the trijets are used by the trijets with Higher BDT value (Higher BDT value -> smaller index: Trijets sorted in BDT value) 
  // Returns true: If cross-cleaned trijet 
  
  // Assume sorted in BDTG score top candidats
  if (Index > 0)
    {

      // For-loop: All top candidates with higher BDTG score than this one
      for (size_t i=0; i<(size_t)Index; i++)
	{
	  // Skip top candidates with same jets as Trijets with higher BDDT value
	  bool bMatchedJet1 = isMatchedJet(TopCand.Jet1.at(Index), TopCand, i);
	  bool bMatchedJet2 = isMatchedJet(TopCand.Jet2.at(Index), TopCand, i);
	  bool bMatchedBJet = isMatchedJet(TopCand.BJet.at(Index), TopCand, i);
	  bool sharedJets   = (bMatchedJet1 || bMatchedJet2 || bMatchedBJet);

	  // Candidate is not cross-cleaned if a shared jet is found and the candidate this is shared with is itself cross-cleaned
	  if (sharedJets && TopIsCrossCleaned(i, TopCand)) return false;
	}
    }
  
  return true;
}

//Get all gen particles by pdgId
vector<genParticle> TopSelectionBDT::GetGenParticles(const vector<genParticle> genParticles, const int pdgId)
{
  std::vector<genParticle> particles;
  // For-loop: All genParticles  
  for (auto& p: genParticles){
    // Find last copy of a given particle
    if (!p.isLastCopy()) continue;
    // Consider only particles
    if (std::abs(p.pdgId()) != pdgId) continue;
    // Save this particle
    particles.push_back(p);
  }
  return particles;
}

//  Get the last copy of a particle.
const genParticle TopSelectionBDT::GetLastCopy(const vector<genParticle> genParticles, const genParticle &p){

  int gen_pdgId = p.pdgId();
  for (size_t i=0; i<p.daughters().size(); i++){
    const genParticle genDau = genParticles[p.daughters().at(i)];
    int genDau_pdgId   = genDau.pdgId();

    if (gen_pdgId == genDau_pdgId)  return GetLastCopy(genParticles, genDau);
  }
  return p;
}


vector<genParticle> TopSelectionBDT::GetWpartons( genParticle daughter, const Event& event){
  vector<genParticle> quarks;

  // Get the last copy
  genParticle W = GetLastCopy(event.genparticles().getGenParticles(), daughter);
  // For-loop: W-boson daughters
  for (size_t idau = 0; idau < W.daughters().size(); idau++)
    {
      // Find the decay products of W-boson
      int Wdau_index   = W.daughters().at(idau);
      genParticle Wdau = event.genparticles().getGenParticles()[Wdau_index];
      // Consider only quarks as decaying products
      if (std::abs(Wdau.pdgId()) > 5) continue;
      // Save daughter
      quarks.push_back(Wdau);
    }//W-boson daughters
  
  return quarks;
}

bool TopSelectionBDT::IsGenuineTop(const Jet& trijetJet1, const Jet& trijetJet2, const Jet& trijetBJet,
				   const std::vector<Jet>& MCtrue_LdgJet,  const std::vector<Jet>& MCtrue_SubldgJet, const std::vector<Jet>& MCtrue_Bjet){
  
  for (size_t k=0; k<MCtrue_Bjet.size(); k++){
    bool same1 = areSameJets(trijetJet1, MCtrue_LdgJet.at(k))       && areSameJets(trijetJet2, MCtrue_SubldgJet.at(k)) && areSameJets(trijetBJet,  MCtrue_Bjet.at(k));
    bool same2 = areSameJets(trijetJet1, MCtrue_SubldgJet.at(k))    && areSameJets(trijetJet2, MCtrue_LdgJet.at(k))    && areSameJets(trijetBJet,  MCtrue_Bjet.at(k));
    if (same1 || same2) return true;
  }
  return false;
}
