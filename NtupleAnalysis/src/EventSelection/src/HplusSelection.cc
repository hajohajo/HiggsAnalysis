// -*- cut++ -*-
#include "EventSelection/interface/HplusSelection.h"

#include "Framework/interface/ParameterSet.h"
#include "EventSelection/interface/CommonPlots.h"
#include "DataFormat/interface/Event.h"
#include "Framework/interface/HistoWrapper.h"
#include "Framework/interface/Exception.h"

#include "Tools/interface/MCTools.h"

#include "Math/VectorUtil.h"

HplusSelection::Data::Data()
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

HplusSelection::Data::~Data() { }


HplusSelection::HplusSelection(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& postfix)
  : BaseSelection(eventCounter, histoWrapper, commonPlots, postfix),
    // Input parameters
    cfg_TopBDTGCut(config, "TopBDTGCut"),
    cfg_AnyTopBDTGCut(config, "AnyTopBDTGCut"),
    cfg_FreeBjetsCut(config, "FreeBjetsCut"),
    fTopTagSFCalculator(config)
{
  initialize(config);
}

HplusSelection::HplusSelection(const ParameterSet& config)
: BaseSelection(),
  // Input parameters
  cfg_TopBDTGCut(config, "TopBDTGCut"),
  cfg_AnyTopBDTGCut(config, "AnyTopBDTGCut"),
  cfg_FreeBjetsCut(config, "FreeBjetsCut"),
  fTopTagSFCalculator(config)
{
  initialize(config);
  bookHistograms(new TDirectory());

}

HplusSelection::~HplusSelection() {
  
  // Histograms
  delete hTopBDT_AllCandidates;
  delete hTopMultiplicity_AllCandidates;
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
  delete hTetrajetBJetPt;
  delete hTetrajetBJetEta;
  delete hTetrajetBJetBDisc;
  delete hTetrajetPt;
  delete hTetrajetMass;
  delete hTetrajetEta;
  delete hLdgTrijetPt;
  delete hLdgTrijetMass;
  delete hLdgTrijetJet1Pt;
  delete hLdgTrijetJet1Eta;
  delete hLdgTrijetJet1BDisc;
  delete hLdgTrijetJet2Pt;
  delete hLdgTrijetJet2Eta;
  delete hLdgTrijetJet2BDisc;
  delete hLdgTrijetBJetPt;
  delete hLdgTrijetBJetEta;
  delete hLdgTrijetBJetBDisc;
  delete hLdgTrijetDiJetPt;
  delete hLdgTrijetDiJetEta;
  delete hLdgTrijetDiJetMass;
  delete hLdgTrijetDijetDeltaR;
  delete hLdgTrijetTopMassWMassRatio;
  delete hLdgTrijet_DeltaR_Trijet_TetrajetBjet;
  delete hLdgTrijet_DeltaEta_Trijet_TetrajetBjet;
  delete hLdgTrijet_DeltaPhi_Trijet_TetrajetBjet;
  delete hLdgTrijet_DeltaY_Trijet_TetrajetBjet;
  delete hSubldgTrijetPt;
  delete hSubldgTrijetMass;
  delete hSubldgTrijetJet1Pt;
  delete hSubldgTrijetJet1Eta;
  delete hSubldgTrijetJet1BDisc;
  delete hSubldgTrijetJet2Pt;
  delete hSubldgTrijetJet2Eta;
  delete hSubldgTrijetJet2BDisc;
  delete hSubldgTrijetBJetPt;
  delete hSubldgTrijetBJetEta;
  delete hSubldgTrijetBJetBDisc;
  delete hSubldgTrijetDiJetPt;
  delete hSubldgTrijetDiJetEta;
  delete hSubldgTrijetDiJetMass;
  delete hSubldgTrijetDijetDeltaR;
  delete hSubldgTrijetTopMassWMassRatio;
  delete hSubldgTrijet_DeltaR_Trijet_TetrajetBjet;
  delete hSubldgTrijet_DeltaEta_Trijet_TetrajetBjet;
  delete hSubldgTrijet_DeltaPhi_Trijet_TetrajetBjet;
  delete hSubldgTrijet_DeltaY_Trijet_TetrajetBjet;
  delete hDeltaR_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet;
  delete hDeltaEta_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet;
  delete hDeltaPhi_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet;
  delete hDeltaY_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet;
}

void HplusSelection::initialize(const ParameterSet& config) {
  
  return;
}

void HplusSelection::bookHistograms(TDirectory* dir) {
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

  const int nDEtaBins   = fCommonPlots->getDeltaEtaBinSettings().bins();
  const double fDEtaMin = fCommonPlots->getDeltaEtaBinSettings().min();
  const double fDEtaMax = fCommonPlots->getDeltaEtaBinSettings().max();

  const int nDPhiBins   = fCommonPlots->getDeltaPhiBinSettings().bins();
  const double fDPhiMin = fCommonPlots->getDeltaPhiBinSettings().min();
  const double fDPhiMax = fCommonPlots->getDeltaPhiBinSettings().max();

  const int  nBDiscBins = 120; //fCommonPlots->getBJetDiscBinSettings().bins();
  const float fBDiscMin = 0;   //fCommonPlots->getBJetDiscBinSettings().min();
  const float fBDiscMax = 1.2; //fCommonPlots->getBJetDiscBinSettings().max();

  const int nWMassBins  = 200; //fCommonPlots->getWMassBinSettings().bins();
  const float fWMassMin = 0;   //fCommonPlots->getWMassBinSettings().min();
  const float fWMassMax = 1000;//fCommonPlots->getWMassBinSettings().max();

  const int nTopMassBins  = 200; //fCommonPlots->getTopMassBinSettings().bins();
  const float fTopMassMin = 0;   //fCommonPlots->getTopMassBinSettings().min();
  const float fTopMassMax = 1000;//fCommonPlots->getTopMassBinSettings().max();

  const int nInvMassBins  = fCommonPlots->getInvMassBinSettings().bins();
  const float fInvMassMin = fCommonPlots->getInvMassBinSettings().min();
  const float fInvMassMax = fCommonPlots->getInvMassBinSettings().max();


  // Histograms
  TDirectory* subdir = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "hplus2tbSelectionBDT_"    + sPostfix);

  // All Candidates
  hTopBDT_AllCandidates           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopBDT_AllCandidates",";top candidate BDT", 40, -1.0, 1.0) ;
  hTopMass_AllCandidates          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMass_AllCandidates" ,";top candidate M (GeVc^{-2})", nTopMassBins, fTopMassMin, fTopMassMax);
  hTopPt_AllCandidates            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopPt_AllCandidates", ";top candidate p_{T} (GeV/c)", nPtBins, fPtMin, fPtMax);
  hTopMultiplicity_AllCandidates  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMultiplicity_AllCandidates", ";top candidate multiplicity", 400, 0.0, 400.0);

  hTopBDT_SelectedCandidates           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopBDT_SelectedCandidates",";top candidate BDT", 40, -1.0, 1.0) ;
  hTopMass_SelectedCandidates          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMass_SelectedCandidates" ,";top candidate M (GeVc^{-2})", nTopMassBins, fTopMassMin, fTopMassMax);
  hTopPt_SelectedCandidates            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopPt_SelectedCandidates", ";top candidate p_{T} (GeV/c)", nPtBins, fPtMin, fPtMax);
  hTopMultiplicity_SelectedCandidates  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMultiplicity_SelectedCandidates", ";top candidate multiplicity", 50, 0.0, 50.0);
  hTopBDT_SelectedCleanedCandidates           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopBDT_SelectedCleanedCandidates",";top candidate BDT", 40, -1.0, 1.0) ;
  hTopMass_SelectedCleanedCandidates          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMass_SelectedCleanedCandidates" ,";top candidate M (GeVc^{-2})", nTopMassBins, fTopMassMin, fTopMassMax);
  hTopPt_SelectedCleanedCandidates            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopPt_SelectedCleanedCandidates", ";top candidate p_{T} (GeV/c)", nPtBins, fPtMin, fPtMax);
  hTopMultiplicity_SelectedCleanedCandidates  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMultiplicity_SelectedCleanedCandidates", ";top candidate multiplicity", 10, 0.0, 10.0);

  hTopBDT_NotSelectedCandidates           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopBDT_NotSelectedCandidates",";top candidate BDT", 40, -1.0, 1.0) ;
  hTopMass_NotSelectedCandidates          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMass_NotSelectedCandidates" ,";top candidate M (GeVc^{-2})", nTopMassBins, fTopMassMin, fTopMassMax);
  hTopPt_NotSelectedCandidates            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopPt_NotSelectedCandidates", ";top candidate p_{T} (GeV/c)", nPtBins, fPtMin, fPtMax);
  hTopMultiplicity_NotSelectedCandidates  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMultiplicity_NotSelectedCandidates", ";top candidate multiplicity", 50, 0.0, 50.0);
  hTopBDT_AllCleanedCandidates           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopBDT_AllCleanedCandidates",";top candidate BDT", 40, -1.0, 1.0) ;
  hTopMass_AllCleanedCandidates          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMass_AllCleanedCandidates" ,";top candidate M (GeVc^{-2})", nTopMassBins, fTopMassMin, fTopMassMax);
  hTopPt_AllCleanedCandidates            = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopPt_AllCleanedCandidates", ";top candidate p_{T} (GeV/c)", nPtBins, fPtMin, fPtMax);
  hTopMultiplicity_AllCleanedCandidates  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TopMultiplicity_AllCleanedCandidates", ";top candidate multiplicity", 10, 0.0, 10.0);

  // Ldg in pt free b-jet
  hTetrajetBJetPt     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TetrajetBJetPt"   ,";p_{T} (GeV/c)"      , nPtBins     , fPtMin     , fPtMax);
  hTetrajetBJetEta    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TetrajetBJetEta"  ,";#eta"               , nEtaBins    , fEtaMin    , fEtaMax);
  hTetrajetBJetBDisc  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TetrajetBJetBDisc",";b-tag discriminator", nBDiscBins  , fBDiscMin  , fBDiscMax);

  // Tetrajet
  hTetrajetPt         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TetrajetPt"  , ";p_{T} (GeV/c)"     , nPtBins     , fPtMin     , fPtMax);
  hTetrajetMass       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TetrajetMass", ";M (GeV/c^{2})"     , nInvMassBins, fInvMassMin, fInvMassMax);
  hTetrajetEta        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "TetrajetEta" , ";#eta"              , nEtaBins    , fEtaMin    , fEtaMax);

  // Leading in pt top
  hLdgTrijetPt          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijetPt"       ,";p_{T} (GeV/c)", nPtBins     , fPtMin     , fPtMax);
  hLdgTrijetMass        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijetMass"     ,";M (GeV/c^{2})", nTopMassBins, fTopMassMin, fTopMassMax);
  hLdgTrijetJet1Pt      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijetJet1Pt"   ,";p_{T} (GeV/c)", nPtBins     , fPtMin     , fPtMax);
  hLdgTrijetJet1Eta     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijetJet1Eta"  ,";#eta"         , nEtaBins    , fEtaMin    , fEtaMax);
  hLdgTrijetJet1BDisc   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijetJet1BDisc",";b-tag discr." , nBDiscBins  , fBDiscMin  , fBDiscMax);
  hLdgTrijetJet2Pt      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijetJet2Pt"   ,";p_{T} (GeV/c)", nPtBins     , fPtMin     , fPtMax);
  hLdgTrijetJet2Eta     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijetJet2Eta"  ,";#eta"         , nEtaBins    , fEtaMin    , fEtaMax);
  hLdgTrijetJet2BDisc   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijetJet2BDisc",";b-tag discr." , nBDiscBins  , fBDiscMin  , fBDiscMax);
  hLdgTrijetBJetPt      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijetBJetPt"   ,";p_{T} (GeV/c)", nPtBins     , fPtMin     , fPtMax);
  hLdgTrijetBJetEta     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijetBJetEta"  ,";#eta"         , nEtaBins    , fEtaMin    , fEtaMax);
  hLdgTrijetBJetBDisc   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijetBJetBDisc",";b-tag discr." , nBDiscBins  , fBDiscMin  , fBDiscMax);
  hLdgTrijetDiJetPt     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijetDiJetPt"  ,";p_{T} (GeV/c)", nPtBins     , fPtMin     , fPtMax);
  hLdgTrijetDiJetEta    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijetDiJetEta" ,";#eta"         , nEtaBins    , fEtaMin    , fEtaMax);
  hLdgTrijetDiJetMass   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijetDiJetMass",";M (GeV/c^{2})", nWMassBins  , fWMassMin  , fWMassMax);
  hLdgTrijetDijetDeltaR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijetDijetDR"  , ";#Delta R(j_{1},j_{2})"  , 2*nDRBins     , fDRMin     , fDRMax);
  hLdgTrijetTopMassWMassRatio    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijetTopMassWMassRatio"   ,";R_{32}", 100 , 0.0, 10.0);
  hLdgTrijet_DeltaR_Trijet_TetrajetBjet   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijet_DeltaR_Trijet_TetrajetBjet" , ";#Delta R(top, b_{free})" , nDRBins, fDRMin, fDRMax);
  hLdgTrijet_DeltaEta_Trijet_TetrajetBjet = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijet_DeltaEta_Trijet_TetrajetBjet", ";#Delta#eta(top, b_{free})", nDEtaBins, fDEtaMin, fDEtaMax);
  hLdgTrijet_DeltaPhi_Trijet_TetrajetBjet = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijet_DeltaPhi_Trijet_TetrajetBjet", ";#Delta#phi(top, b_{free})", nDPhiBins, fDPhiMin, fDPhiMax);
  hLdgTrijet_DeltaY_Trijet_TetrajetBjet   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "LdgTrijet_DeltaY_Trijet_TetrajetBjet", ";#Delta Y(top, b_{free})", nDRBins, fDRMin, fDRMax);

  // Sub-Leading in pt top
  hSubldgTrijetPt          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijetPt"       ,";p_{T} (GeV/c)", nPtBins     , fPtMin     , fPtMax);
  hSubldgTrijetMass        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijetMass"     ,";M (GeV/c^{2})", nTopMassBins, fTopMassMin, fTopMassMax);
  hSubldgTrijetJet1Pt      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijetJet1Pt"   ,";p_{T} (GeV/c)", nPtBins     , fPtMin     , fPtMax);
  hSubldgTrijetJet1Eta     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijetJet1Eta"  ,";#eta"         , nEtaBins    , fEtaMin    , fEtaMax);
  hSubldgTrijetJet1BDisc   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijetJet1BDisc",";b-tag discr." , nBDiscBins  , fBDiscMin  , fBDiscMax);
  hSubldgTrijetJet2Pt      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijetJet2Pt"   ,";p_{T} (GeV/c)", nPtBins     , fPtMin     , fPtMax);
  hSubldgTrijetJet2Eta     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijetJet2Eta"  ,";#eta"         , nEtaBins    , fEtaMin    , fEtaMax);
  hSubldgTrijetJet2BDisc   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijetJet2BDisc",";b-tag discr." , nBDiscBins  , fBDiscMin  , fBDiscMax);
  hSubldgTrijetBJetPt      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijetBJetPt"   ,";p_{T} (GeV/c)", nPtBins     , fPtMin     , fPtMax);
  hSubldgTrijetBJetEta     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijetBJetEta"  ,";#eta"         , nEtaBins    , fEtaMin    , fEtaMax);
  hSubldgTrijetBJetBDisc   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijetBJetBDisc",";b-tag discr." , nBDiscBins  , fBDiscMin  , fBDiscMax);
  hSubldgTrijetDiJetPt     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijetDiJetPt"  ,";p_{T} (GeV/c)", nPtBins     , fPtMin     , fPtMax);
  hSubldgTrijetDiJetEta    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijetDiJetEta" ,";#eta"         , nEtaBins    , fEtaMin    , fEtaMax);
  hSubldgTrijetDiJetMass   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijetDiJetMass",";M (GeV/c^{2})", nWMassBins  , fWMassMin  , fWMassMax);
  hSubldgTrijetDijetDeltaR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SbldgTrijetDijetDR"  , ";#Delta R(j_{1},j_{2})"  , nDRBins     , fDRMin     , fDRMax);
  hSubldgTrijetTopMassWMassRatio = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijetTopMassWMassRatio",";R_{32}", 100 , 0.0, 10.0);
  hSubldgTrijet_DeltaR_Trijet_TetrajetBjet   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijet_DeltaR_Trijet_TetrajetBjet", ";#Delta R(top, b_{free})", nDRBins, fDRMin, fDRMax);
  hSubldgTrijet_DeltaEta_Trijet_TetrajetBjet = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijet_DeltaEta_Trijet_TetrajetBjet", ";#Delta#eta(top, b_{free})", nDEtaBins, fDEtaMin, fDEtaMax);
  hSubldgTrijet_DeltaPhi_Trijet_TetrajetBjet = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijet_DeltaPhi_Trijet_TetrajetBjet", ";#Delta#phi(top, b_{free})", nDPhiBins, fDPhiMin, fDPhiMax);
  hSubldgTrijet_DeltaY_Trijet_TetrajetBjet   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "SubldgTrijet_DeltaY_Trijet_TetrajetBjet", ";#Delta Y(top, b_{free})", nDRBins, fDRMin, fDRMax);

  // Histograms (2D)
  hDeltaR_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet   = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, subdir, "DeltaR_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet",
                                                                                             ";#Delta R (Trijet_{Ldg}, b_{free}^{ldg});#Delta R (Trijet_{Sbldg}, b_{free}^{ldg})",
                                                                                             nDRBins     , fDRMin     , 6., nDRBins     , fDRMin     , 6.);
  hDeltaEta_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, subdir, "DeltaEta_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet",
                                                                                             ";#Delta #eta (Trijet_{Ldg}, b_{free}^{ldg}) ;#Delta #eta (Trijet_{Sbldg}, b_{free}^{ldg})",
                                                                                             nDEtaBins, fDEtaMin, 6., nDEtaBins, fDEtaMin, 6.);
  hDeltaPhi_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, subdir, "DeltaPhi_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet",
                                                                                             ";#Delta #phi (Trijet_{Ldg}, b_{free}^{ldg});#Delta #phi (Trijet_{Sbldg}, b_{free}^{ldg})",
                                                                                             nDPhiBins , fDPhiMin , fDPhiMax, nDPhiBins , fDPhiMin , fDPhiMax);
  hDeltaY_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet   = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, subdir, "DeltaY_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet",
                                                                                             ";#Delta Y (Trijet_{Ldg}, b_{free}^{ldg});#Delta Y (Trijet_{Sbldg}, b_{free}^{ldg})",
                                                                                             nDRBins     , fDRMin     , 6., nDRBins     , fDRMin     , 6.);

  fTopTagSFCalculator.bookHistograms(subdir, fHistoWrapper);
  return;
}


HplusSelection::Data HplusSelection::silentAnalyze(const Event& event, const JetSelection::Data& jetData, const BJetSelection::Data& bjetData, const TopSelectionBDT::Data& topData) {
  ensureSilentAnalyzeAllowed(event.eventID());

  // Disable histogram filling and counter
  disableHistogramsAndCounters();
  Data myData = privateAnalyze(event, jetData.getSelectedJets(), bjetData.getSelectedBJets(), topData);
  enableHistogramsAndCounters();
  return myData;
}


HplusSelection::Data HplusSelection::analyze(const Event& event, const JetSelection::Data& jetData, const BJetSelection::Data& bjetData, const TopSelectionBDT::Data& topData) {
  ensureAnalyzeAllowed(event.eventID());

  // Ready to analyze
  HplusSelection::Data data = privateAnalyze(event, jetData.getSelectedJets(), bjetData.getSelectedBJets(), topData);

  // Send data to CommonPlots
  // if (fCommonPlots != nullptr) fCommonPlots->fillControlPlotsAtHplusSelection(event, data);// fixme - implement
  return data;
}


HplusSelection::Data HplusSelection::privateAnalyze(const Event& event, const std::vector<Jet> selectedJets, const std::vector<Jet> selectedBjets, const TopSelectionBDT::Data& topData) {
  Data output;

  std::vector<Jet> jets  = selectedJets;
  std::vector<Jet> bjets = selectedBjets;

  std::vector<Jet>	allTopsJet1	 = topData.getAllTopsJet1();
  std::vector<Jet>	allTopsJet2	 = topData.getAllTopsJet2();
  std::vector<Jet>	allTopsBJet	 = topData.getAllTopsBJet();  
  std::vector<float>	allTopsBDTG      = topData.getAllTopsBDTG();
  std::vector<bool>	allTopsIsGenuine = topData.getAllTopsIsGenuine();
  std::vector<bool>	allTopsIsTagged  = topData.getAllTopsIsTagged();

  std::vector<Jet>	allCleanedTopsJet1	= topData.getAllCleanedTopsJet1();
  std::vector<Jet>	allCleanedTopsJet2	= topData.getAllCleanedTopsJet2();
  std::vector<Jet>	allCleanedTopsBJet	= topData.getAllCleanedTopsBJet();
  std::vector<float>	allCleanedTopsBDTG      = topData.getAllCleanedTopsBDTG();
  std::vector<bool>	allCleanedTopsIsGenuine = topData.getAllCleanedTopsIsGenuine();
  std::vector<bool>	allCleanedTopsIsTagged  = topData.getAllCleanedTopsIsTagged();
  
  //Define variables
  Jet tetrajetBjet;
  double tetrajetBjetPt_max = -999.99;

  output.fTopTaggingScaleFactorEventWeight = topData.getTopTaggingScaleFactorEventWeight();
  //================================================================================================  
  // Hplus2tb Selection
  //================================================================================================  
  if (0) std::cout << "=== HplusSelection:: " << std::endl;
  
  //Check if there are free bjets after reconstructing selected cleaned tops with highest BDTG.
  int nbjets_Top1 = 0, nbjets_Top2 = 0; 
  
  if (topData.getAllCleanedTopsSize() > 0) nbjets_Top1 = isBJet(allCleanedTopsJet1.at(0), bjets) + isBJet(allCleanedTopsJet2.at(0), bjets) + isBJet(allCleanedTopsBJet.at(0), bjets);  
  if (topData.getAllCleanedTopsSize() > 1) nbjets_Top2 = isBJet(allCleanedTopsJet1.at(1), bjets) + isBJet(allCleanedTopsJet2.at(1), bjets) + isBJet(allCleanedTopsBJet.at(1), bjets);

  int n_freeBjets1 = (int)bjets.size() - nbjets_Top1;
  int n_freeBjets2 = (int)jets.size() - nbjets_Top2;

  //If the number of free bjets is less that fFreeBjetsCuts, skip trijet and cross clean the new candidates
  if ( !cfg_FreeBjetsCut.passedCut(n_freeBjets1) || !cfg_FreeBjetsCut.passedCut(n_freeBjets2)){
    //Clear old collections
    allCleanedTopsJet1.clear(); allCleanedTopsJet2.clear(); allCleanedTopsBJet.clear();
    allCleanedTopsBDTG.clear(); allCleanedTopsIsGenuine.clear(); allCleanedTopsIsTagged.clear();   
    
    std::vector< math::XYZTLorentzVector> allCleanedTops_P4;
    for (size_t i = 0; i < allTopsBJet.size(); i++){
      int nbjets_i = isBJet(allTopsJet1.at(i), bjets) + isBJet(allTopsJet2.at(i), bjets) + isBJet(allTopsBJet.at(i), bjets);
      bool isCrossCleaned = TopIsCrossCleaned(i, allTopsJet1, allTopsJet2, allTopsBJet, bjets) && ((size_t)nbjets_i != bjets.size());
      if (!isCrossCleaned) continue;      
      
      //Store new coss-cleaned candidates.
      allCleanedTopsJet1.push_back(allTopsJet1.at(i));
      allCleanedTopsJet2.push_back(allTopsJet2.at(i));
      allCleanedTopsBJet.push_back(allTopsBJet.at(i));
      allCleanedTopsBDTG.push_back(allTopsBDTG.at(i));
      allCleanedTopsIsGenuine.push_back(allTopsIsGenuine.at(i));
      allCleanedTopsIsTagged.push_back(allTopsIsTagged.at(i));
      math::XYZTLorentzVector top_p4;
      top_p4 = allTopsJet1.at(i).p4() + allTopsJet2.at(i).p4() + allTopsBJet.at(i).p4();
      allCleanedTops_P4.push_back( top_p4 );
    }
    // Re-calculate and store top tagging SF weight and it's uncertainty
    if (event.isMC())
      {
	output.fTopTaggingScaleFactorEventWeight = fTopTagSFCalculator.calculateSF(allCleanedTops_P4, allCleanedTopsBDTG, allCleanedTopsIsTagged, allCleanedTopsIsGenuine);    
      }
  } //if (nbjets_Top1 == bjets.size() || nbjets_Top2 == bjets.size()){  
  
  if (0) std::cout << "=== HplusSelection:: Find the tetrajet b-jet" << std::endl;

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
  if (0) std::cout << "=== HplusSelection:: Fill output" << std::endl;

  output.nAllCleanedTops = allCleanedTopsBDTG.size();

  bool bPass_LdgMVA     = false;   
  bool bPass_SubldgMVA  = false;

  if (allCleanedTopsBDTG.size() == 1)
    {
      bPass_LdgMVA    = cfg_TopBDTGCut.passedCut( allCleanedTopsBDTG.at(0) );

      // Leading-in-BDTG top
      output.fMVAmax1           = allCleanedTopsBDTG.at(0);
      output.fTrijet1Jet1       = allCleanedTopsJet1.at(0);
      output.fTrijet1Jet2       = allCleanedTopsJet2.at(0);
      output.fTrijet1BJet       = allCleanedTopsBJet.at(0);
      output.fTrijet1Dijet_p4   = allCleanedTopsJet1.at(0).p4() + allCleanedTopsJet2.at(0).p4();
      output.fTrijet1_p4        = allCleanedTopsJet1.at(0).p4() + allCleanedTopsJet2.at(0).p4() + allCleanedTopsBJet.at(0).p4();
    }
  else if (allCleanedTopsBDTG.size() > 1)
    {
      bPass_LdgMVA     = cfg_TopBDTGCut.passedCut( allCleanedTopsBDTG.at(0) );
      bPass_SubldgMVA  = cfg_TopBDTGCut.passedCut( allCleanedTopsBDTG.at(1) );
      // Leading-in-BDTG top
      output.fMVAmax1           = allCleanedTopsBDTG.at(0);
      output.fTrijet1Jet1       = allCleanedTopsJet1.at(0);
      output.fTrijet1Jet2       = allCleanedTopsJet2.at(0);
      output.fTrijet1BJet       = allCleanedTopsBJet.at(0);
      output.fTrijet1Dijet_p4   = allCleanedTopsJet1.at(0).p4() + allCleanedTopsJet2.at(0).p4();
      output.fTrijet1_p4        = allCleanedTopsJet1.at(0).p4() + allCleanedTopsJet2.at(0).p4() + allCleanedTopsBJet.at(0).p4();
      // Subleading-in-MVA top
      output.fMVAmax2           = allCleanedTopsBDTG.at(1);
      output.fTrijet2Jet1       = allCleanedTopsJet1.at(1);
      output.fTrijet2Jet2       = allCleanedTopsJet2.at(1);
      output.fTrijet2BJet       = allCleanedTopsBJet.at(1);
      output.fTrijet2Dijet_p4   = allCleanedTopsJet1.at(1).p4() + allCleanedTopsJet2.at(1).p4();
      output.fTrijet2_p4        = allCleanedTopsJet1.at(1).p4() + allCleanedTopsJet2.at(1).p4() + allCleanedTopsBJet.at(1).p4();
    }

  if (tetrajetBjetPt_max > 0)
    {
      output.fTetrajetBJet = tetrajetBjet;
      
      if (allCleanedTopsBDTG.size() == 1) output.fLdgTetrajet_p4 =  output.fTrijet1_p4 + output.fTetrajetBJet.p4();
      else if (allCleanedTopsBDTG.size() > 1)
	{
	  bool bLdgBDTGisLdgPt = output.fTrijet1_p4.Pt() > output.fTrijet2_p4.Pt();
	  output.fLdgTetrajet_p4 = (bLdgBDTGisLdgPt) ? (output.fTrijet1_p4 + output.fTetrajetBJet.p4()) : (output.fTrijet2_p4 + output.fTetrajetBJet.p4());
	}
    }
  
  //Booleans
  bool bPass_FreeBjet   = (tetrajetBjetPt_max > 0);
  bool bPass_BothMVA    = bPass_LdgMVA * bPass_SubldgMVA;
  bool bPass_AnyTwoTops = false; // at least TWO tops with BDT > -1.0 
  if (allCleanedTopsBDTG.size() > 1) bPass_AnyTwoTops = cfg_AnyTopBDTGCut.passedCut(allCleanedTopsBDTG.at(1) );

  output.bPassedSelection = bPass_FreeBjet * bPass_BothMVA;
  output.bHasTwoTopsAndFreeB = bPass_AnyTwoTops * bPass_FreeBjet;


  //Fill Histograms
  if (0) std::cout << "=== HplusSelection:: Fill histograms" << std::endl;
  
  //Fill Histograms for all candidates
  int nSelected = 0;
  for (int i=0; i< (int)allTopsBDTG.size(); i++){
    hTopBDT_AllCandidates -> Fill(allTopsBDTG.at(i));    
    hTopMass_AllCandidates -> Fill( (allTopsJet1.at(i).p4() + allTopsJet2.at(i).p4() + allTopsBJet.at(i).p4()).M());
    hTopPt_AllCandidates -> Fill( (allTopsJet1.at(i).p4() + allTopsJet2.at(i).p4() + allTopsBJet.at(i).p4()).Pt());
    if (!cfg_TopBDTGCut.passedCut(allTopsBDTG.at(i))) continue;
    nSelected ++;
    hTopBDT_SelectedCandidates -> Fill(allTopsBDTG.at(i));
    hTopMass_SelectedCandidates -> Fill((allTopsJet1.at(i).p4() + allTopsJet2.at(i).p4() + allTopsBJet.at(i).p4()).M());
    hTopPt_SelectedCandidates -> Fill((allTopsJet1.at(i).p4() + allTopsJet2.at(i).p4() + allTopsBJet.at(i).p4()).Pt());
  }

  //Top multiplicity
  hTopMultiplicity_AllCandidates      -> Fill(allTopsBDTG.size());
  hTopMultiplicity_SelectedCandidates -> Fill(nSelected);

  int nCleanedSelected = 0;
  for (int i=0; i<(int)allCleanedTopsBDTG.size(); i++){    
    hTopBDT_AllCleanedCandidates  -> Fill(allCleanedTopsBDTG.at(i));
    hTopMass_AllCleanedCandidates -> Fill( (allCleanedTopsJet1.at(i).p4() + allCleanedTopsJet2.at(i).p4() + allCleanedTopsBJet.at(i).p4()).M());
    hTopPt_AllCleanedCandidates   -> Fill( (allCleanedTopsJet1.at(i).p4() + allCleanedTopsJet2.at(i).p4() + allCleanedTopsBJet.at(i).p4()).Pt());
    if (!cfg_TopBDTGCut.passedCut(allCleanedTopsBDTG.at(i))) continue;
    nCleanedSelected ++;
    hTopBDT_SelectedCleanedCandidates  -> Fill(allCleanedTopsBDTG.at(i));
    hTopMass_SelectedCleanedCandidates -> Fill( (allCleanedTopsJet1.at(i).p4() + allCleanedTopsJet2.at(i).p4() + allCleanedTopsBJet.at(i).p4()).M());
    hTopPt_SelectedCleanedCandidates   -> Fill( (allCleanedTopsJet1.at(i).p4() + allCleanedTopsJet2.at(i).p4() + allCleanedTopsBJet.at(i).p4()).Pt());
  }
  //Top multiplicity
  hTopMultiplicity_AllCleanedCandidates      -> Fill(allCleanedTopsBDTG.size());
  hTopMultiplicity_SelectedCleanedCandidates -> Fill(nCleanedSelected);
  
  //Pass selections
  if (!output.bPassedSelection) return output;
  
  //Ldg in BDTG top is ldg in pt
  if (output.fTrijet1_p4.Pt() > output.fTrijet2_p4.Pt()){
    //Ldg Top
    hLdgTrijetPt          -> Fill(output.fTrijet1_p4.Pt());
    hLdgTrijetMass        -> Fill(output.fTrijet1_p4.M());
    hLdgTrijetJet1Pt      -> Fill(output.fTrijet1Jet1.pt());
    hLdgTrijetJet1Eta     -> Fill(output.fTrijet1Jet1.eta());
    hLdgTrijetJet1BDisc   -> Fill(output.fTrijet1Jet1.bjetDiscriminator());
    hLdgTrijetJet2Pt      -> Fill(output.fTrijet1Jet2.pt());
    hLdgTrijetJet2Eta     -> Fill(output.fTrijet1Jet2.eta());
    hLdgTrijetJet2BDisc   -> Fill(output.fTrijet1Jet2.bjetDiscriminator());
    hLdgTrijetBJetPt      -> Fill(output.fTrijet1BJet.pt());
    hLdgTrijetBJetEta     -> Fill(output.fTrijet1BJet.eta());
    hLdgTrijetBJetBDisc   -> Fill(output.fTrijet1BJet.bjetDiscriminator());
    hLdgTrijetDiJetPt     -> Fill(output.fTrijet1Dijet_p4.Pt());
    hLdgTrijetDiJetEta    -> Fill(output.fTrijet1Dijet_p4.Eta());
    hLdgTrijetDiJetMass   -> Fill(output.fTrijet1Dijet_p4.M());
    hLdgTrijetDijetDeltaR -> Fill(ROOT::Math::VectorUtil::DeltaR(output.fTrijet1Jet1.p4(), output.fTrijet1Jet2.p4()));
    hLdgTrijetTopMassWMassRatio -> Fill(output.fTrijet1_p4.M()/output.fTrijet1Dijet_p4.M());
    hLdgTrijet_DeltaR_Trijet_TetrajetBjet -> Fill(ROOT::Math::VectorUtil::DeltaR(output.fTrijet1_p4, output.fTetrajetBJet.p4()));
    hLdgTrijet_DeltaEta_Trijet_TetrajetBjet -> Fill(std::abs(output.fTrijet1_p4.Eta() - output.fTetrajetBJet.eta()));
    hLdgTrijet_DeltaPhi_Trijet_TetrajetBjet -> Fill(std::abs(ROOT::Math::VectorUtil::DeltaPhi(output.fTrijet1_p4, output.fTetrajetBJet.p4())));

    double LdgTrijet_Rapidity    = 0.5*log((output.fTrijet1_p4.E() + output.fTrijet1_p4.Pz())/(output.fTrijet1_p4.E() - output.fTrijet1_p4.Pz()));
    double TetrajetBjet_Rapidity = 0.5*log((output.fTetrajetBJet.p4().E() + output.fTetrajetBJet.p4().Pz())/(output.fTetrajetBJet.p4().E() - output.fTetrajetBJet.p4().Pz()));
    hLdgTrijet_DeltaY_Trijet_TetrajetBjet   -> Fill(std::abs(LdgTrijet_Rapidity - TetrajetBjet_Rapidity));

   //Subldg Top
    hSubldgTrijetPt          -> Fill(output.fTrijet2_p4.Pt());
    hSubldgTrijetMass        -> Fill(output.fTrijet2_p4.M());
    hSubldgTrijetJet1Pt      -> Fill(output.fTrijet2Jet1.pt());
    hSubldgTrijetJet1Eta     -> Fill(output.fTrijet2Jet1.eta());
    hSubldgTrijetJet1BDisc   -> Fill(output.fTrijet2Jet1.bjetDiscriminator());
    hSubldgTrijetJet2Pt      -> Fill(output.fTrijet2Jet2.pt());
    hSubldgTrijetJet2Eta     -> Fill(output.fTrijet2Jet2.eta());
    hSubldgTrijetJet2BDisc   -> Fill(output.fTrijet2Jet2.bjetDiscriminator());
    hSubldgTrijetBJetPt      -> Fill(output.fTrijet2BJet.pt());
    hSubldgTrijetBJetEta     -> Fill(output.fTrijet2BJet.eta());
    hSubldgTrijetBJetBDisc   -> Fill(output.fTrijet2BJet.bjetDiscriminator());
    hSubldgTrijetDiJetPt     -> Fill(output.fTrijet2Dijet_p4.Pt());
    hSubldgTrijetDiJetEta    -> Fill(output.fTrijet2Dijet_p4.Eta());
    hSubldgTrijetDiJetMass   -> Fill(output.fTrijet2Dijet_p4.M());
    hSubldgTrijetDijetDeltaR -> Fill(ROOT::Math::VectorUtil::DeltaR(output.fTrijet2Jet1.p4(), output.fTrijet2Jet2.p4()));
    hSubldgTrijetTopMassWMassRatio -> Fill(output.fTrijet2_p4.M()/output.fTrijet2Dijet_p4.M());
    hSubldgTrijet_DeltaR_Trijet_TetrajetBjet   -> Fill(ROOT::Math::VectorUtil::DeltaR(output.fTrijet2_p4, output.fTetrajetBJet.p4()));
    hSubldgTrijet_DeltaEta_Trijet_TetrajetBjet -> Fill(std::abs(output.fTrijet2_p4.Eta() - output.fTetrajetBJet.eta()));
    hSubldgTrijet_DeltaPhi_Trijet_TetrajetBjet -> Fill(std::abs(ROOT::Math::VectorUtil::DeltaPhi(output.fTrijet2_p4, output.fTetrajetBJet.p4())));
    double SubldgTrijet_Rapidity    = 0.5*log((output.fTrijet2_p4.E() + output.fTrijet2_p4.Pz())/(output.fTrijet2_p4.E() - output.fTrijet2_p4.Pz()));
    hSubldgTrijet_DeltaY_Trijet_TetrajetBjet   -> Fill(std::abs(SubldgTrijet_Rapidity - TetrajetBjet_Rapidity));

    hDeltaR_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet   -> Fill(ROOT::Math::VectorUtil::DeltaR(output.fTrijet1_p4, output.fTetrajetBJet.p4()), 
									ROOT::Math::VectorUtil::DeltaR(output.fTrijet2_p4, output.fTetrajetBJet.p4()));
    hDeltaEta_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet -> Fill(std::abs(output.fTrijet1_p4.Eta() - output.fTetrajetBJet.eta()), 
									  std::abs(output.fTrijet2_p4.Eta() - output.fTetrajetBJet.eta()));
    hDeltaPhi_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet -> Fill(std::abs(ROOT::Math::VectorUtil::DeltaPhi(output.fTrijet1_p4, output.fTetrajetBJet.p4())),
									  std::abs(ROOT::Math::VectorUtil::DeltaPhi(output.fTrijet2_p4, output.fTetrajetBJet.p4())));
    hDeltaY_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet   -> Fill(std::abs(LdgTrijet_Rapidity - TetrajetBjet_Rapidity),
									  std::abs(SubldgTrijet_Rapidity - TetrajetBjet_Rapidity));
  }
  //Subldg in BDTG top is ldg in pt
  else{
   //Ldg Top
    hLdgTrijetPt          -> Fill(output.fTrijet2_p4.Pt());
    hLdgTrijetMass        -> Fill(output.fTrijet2_p4.M());
    hLdgTrijetJet1Pt      -> Fill(output.fTrijet2Jet1.pt());
    hLdgTrijetJet1Eta     -> Fill(output.fTrijet2Jet1.eta());
    hLdgTrijetJet1BDisc   -> Fill(output.fTrijet2Jet1.bjetDiscriminator());
    hLdgTrijetJet2Pt      -> Fill(output.fTrijet2Jet2.pt());
    hLdgTrijetJet2Eta     -> Fill(output.fTrijet2Jet2.eta());
    hLdgTrijetJet2BDisc   -> Fill(output.fTrijet2Jet2.bjetDiscriminator());
    hLdgTrijetBJetPt      -> Fill(output.fTrijet2BJet.pt());
    hLdgTrijetBJetEta     -> Fill(output.fTrijet2BJet.eta());
    hLdgTrijetBJetBDisc   -> Fill(output.fTrijet2BJet.bjetDiscriminator());
    hLdgTrijetDiJetPt     -> Fill(output.fTrijet2Dijet_p4.Pt());
    hLdgTrijetDiJetEta    -> Fill(output.fTrijet2Dijet_p4.Eta());
    hLdgTrijetDiJetMass   -> Fill(output.fTrijet2Dijet_p4.M());
    hLdgTrijetDijetDeltaR -> Fill(ROOT::Math::VectorUtil::DeltaR(output.fTrijet2Jet1.p4(), output.fTrijet2Jet2.p4()));
    hLdgTrijetTopMassWMassRatio -> Fill(output.fTrijet2_p4.M()/output.fTrijet2Dijet_p4.M());
    hLdgTrijet_DeltaR_Trijet_TetrajetBjet   -> Fill(ROOT::Math::VectorUtil::DeltaR(output.fTrijet2_p4, output.fTetrajetBJet.p4()));
    hLdgTrijet_DeltaEta_Trijet_TetrajetBjet -> Fill(std::abs(output.fTrijet2_p4.Eta() - output.fTetrajetBJet.eta()));
    hLdgTrijet_DeltaPhi_Trijet_TetrajetBjet -> Fill(std::abs(ROOT::Math::VectorUtil::DeltaPhi(output.fTrijet2_p4, output.fTetrajetBJet.p4())));
    double LdgTrijet_Rapidity    = 0.5*log((output.fTrijet2_p4.E() + output.fTrijet2_p4.Pz())/(output.fTrijet2_p4.E() - output.fTrijet2_p4.Pz()));
    double TetrajetBjet_Rapidity = 0.5*log((output.fTetrajetBJet.p4().E() + output.fTetrajetBJet.p4().Pz())/(output.fTetrajetBJet.p4().E() - output.fTetrajetBJet.p4().Pz()));
    hLdgTrijet_DeltaY_Trijet_TetrajetBjet   -> Fill(std::abs(LdgTrijet_Rapidity - TetrajetBjet_Rapidity));

    //Subldg Top
    hSubldgTrijetPt          -> Fill(output.fTrijet1_p4.Pt());
    hSubldgTrijetMass        -> Fill(output.fTrijet1_p4.M());
    hSubldgTrijetJet1Pt      -> Fill(output.fTrijet1Jet1.pt());
    hSubldgTrijetJet1Eta     -> Fill(output.fTrijet1Jet1.eta());
    hSubldgTrijetJet1BDisc   -> Fill(output.fTrijet1Jet1.bjetDiscriminator());
    hSubldgTrijetJet2Pt      -> Fill(output.fTrijet1Jet2.pt());
    hSubldgTrijetJet2Eta     -> Fill(output.fTrijet1Jet2.eta());
    hSubldgTrijetJet2BDisc   -> Fill(output.fTrijet1Jet2.bjetDiscriminator());
    hSubldgTrijetBJetPt      -> Fill(output.fTrijet1BJet.pt());
    hSubldgTrijetBJetEta     -> Fill(output.fTrijet1BJet.eta());
    hSubldgTrijetBJetBDisc   -> Fill(output.fTrijet1BJet.bjetDiscriminator());
    hSubldgTrijetDiJetPt     -> Fill(output.fTrijet1Dijet_p4.Pt());
    hSubldgTrijetDiJetEta    -> Fill(output.fTrijet1Dijet_p4.Eta());
    hSubldgTrijetDiJetMass   -> Fill(output.fTrijet1Dijet_p4.M());
    hSubldgTrijetDijetDeltaR -> Fill(ROOT::Math::VectorUtil::DeltaR(output.fTrijet1Jet1.p4(), output.fTrijet1Jet2.p4()));
    hSubldgTrijetTopMassWMassRatio -> Fill(output.fTrijet1_p4.M()/output.fTrijet1Dijet_p4.M());
    hSubldgTrijet_DeltaR_Trijet_TetrajetBjet -> Fill(ROOT::Math::VectorUtil::DeltaR(output.fTrijet1_p4, output.fTetrajetBJet.p4()));
    hSubldgTrijet_DeltaEta_Trijet_TetrajetBjet -> Fill(std::abs(output.fTrijet1_p4.Eta() - output.fTetrajetBJet.eta()));
    hSubldgTrijet_DeltaPhi_Trijet_TetrajetBjet -> Fill(std::abs(ROOT::Math::VectorUtil::DeltaPhi(output.fTrijet1_p4, output.fTetrajetBJet.p4())));

    double SubldgTrijet_Rapidity    = 0.5*log((output.fTrijet1_p4.E() + output.fTrijet1_p4.Pz())/(output.fTrijet1_p4.E() - output.fTrijet1_p4.Pz()));
    hSubldgTrijet_DeltaY_Trijet_TetrajetBjet   -> Fill(std::abs(SubldgTrijet_Rapidity - TetrajetBjet_Rapidity));

    hDeltaR_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet   -> Fill(ROOT::Math::VectorUtil::DeltaR(output.fTrijet2_p4, output.fTetrajetBJet.p4()), 
									ROOT::Math::VectorUtil::DeltaR(output.fTrijet1_p4, output.fTetrajetBJet.p4()));
    hDeltaEta_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet -> Fill(std::abs(output.fTrijet2_p4.Eta() - output.fTetrajetBJet.eta()), 
									  std::abs(output.fTrijet1_p4.Eta() - output.fTetrajetBJet.eta()));
    hDeltaPhi_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet -> Fill(std::abs(ROOT::Math::VectorUtil::DeltaPhi(output.fTrijet2_p4, output.fTetrajetBJet.p4())),
									  std::abs(ROOT::Math::VectorUtil::DeltaPhi(output.fTrijet1_p4, output.fTetrajetBJet.p4())));
    hDeltaY_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet   -> Fill(std::abs(LdgTrijet_Rapidity - TetrajetBjet_Rapidity),
									  std::abs(SubldgTrijet_Rapidity - TetrajetBjet_Rapidity));
  }
  
  //Tetrajet variables
  hTetrajetBJetPt    -> Fill(output.fTetrajetBJet.pt());
  hTetrajetBJetEta   -> Fill(output.fTetrajetBJet.eta());
  hTetrajetBJetBDisc -> Fill(output.fTetrajetBJet.bjetDiscriminator());
  hTetrajetPt        -> Fill(output.fLdgTetrajet_p4.Pt());
  hTetrajetMass      -> Fill(output.fLdgTetrajet_p4.M());
  hTetrajetEta       -> Fill(output.fLdgTetrajet_p4.Eta());

  return output;
}

bool HplusSelection::areSameJets(const Jet& jet1, const Jet& jet2) {
  float dR = ROOT::Math::VectorUtil::DeltaR(jet1.p4(), jet2.p4());
  float dR_match = 0.1;
  if (dR <= dR_match) return true;
  else return false;
}

bool HplusSelection::isBJet(const Jet& jet, const std::vector<Jet>& bjets) {
  for (auto bjet: bjets)
    {
      if (areSameJets(jet, bjet)) return true;
    }
  return false;
}


bool HplusSelection::isMatchedJet(const Jet& jet, std::vector<Jet> Jet1, std::vector<Jet> Jet2, std::vector<Jet> BJet, const unsigned int index) {
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

bool HplusSelection::TopIsCrossCleaned( int Index, std::vector<Jet> jet1, std::vector<Jet> jet2, std::vector<Jet> bjet, std::vector<Jet> bjets){
  // Description:
  // Used to find the cross-cleaned trijet multiplicity. The function takes as input the index of a trijet and the total trijet collection and returns 
  // false true if If at least one of the subjets of the trijets are used by the trijets with Higher BDT value (Higher BDT value -> smaller index: Trijets sorted in BDT value) 
  // Returns true: If cross-cleaned trijet 
  
  int nbjets_i = isBJet(jet1.at(Index), bjets) + isBJet(jet2.at(Index), bjets) + isBJet(bjet.at(Index), bjets); //fixme!
  int n_freeBjets = (int)bjets.size() - nbjets_i;
  if ( !cfg_FreeBjetsCut.passedCut(n_freeBjets)  ) return false;

  // Assume sorted in BDTG score top candidats
  if (Index > 0)
    {

      // For-loop: All top candidates with higher BDTG score than this one
      for (size_t i=0; i<(size_t)Index; i++)
	{
	  // Skip top candidates with same jets as Trijets with higher BDT value
	  bool bMatchedJet1 = isMatchedJet(jet1.at(Index), jet1, jet2, bjet, i);
	  bool bMatchedJet2 = isMatchedJet(jet2.at(Index), jet1, jet2, bjet, i);
	  bool bMatchedBJet = isMatchedJet(bjet.at(Index), jet1, jet2, bjet, i);
	  bool sharedJets   = (bMatchedJet1 || bMatchedJet2 || bMatchedBJet);

	  // Candidate is not cross-cleaned if a shared jet is found and the candidate this is shared with is itself cross-cleaned
	  if (sharedJets && TopIsCrossCleaned(i, jet1, jet2, bjet, bjets)) return false;
	}
    }
  
  return true;
}
