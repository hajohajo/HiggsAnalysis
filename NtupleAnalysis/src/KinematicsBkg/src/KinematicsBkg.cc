// -*- c++ -*-
#include "Framework/interface/BaseSelector.h"
#include "Framework/interface/makeTH.h"

// User
#include "Auxiliary/interface/Table.h"
#include "Auxiliary/interface/Tools.h"
#include "Tools/interface/MCTools.h"
#include "Tools/interface/DirectionalCut.h"
#include "DataFormat/interface/Particle.h"
#include "EventSelection/interface/CommonPlots.h"
#include "EventSelection/interface/EventSelections.h"

// ROOT
#include "TDirectory.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

struct PtComparator
{
  bool operator() (const genParticle p1, const genParticle p2) const { return ( p1.pt() > p2.pt() ); }
  bool operator() (const math::XYZTLorentzVector p1, const math::XYZTLorentzVector p2) const { return ( p1.pt() > p2.pt() ); }
};


class KinematicsBkg: public BaseSelector {
public:
  explicit KinematicsBkg(const ParameterSet& config, const TH1* skimCounters);
  virtual ~KinematicsBkg() {}

  /// Books histograms
  virtual void book(TDirectory *dir) override;
  /// Sets up branches for reading the TTree
  virtual void setupBranches(BranchManager& branchManager) override;
  /// Called for each event
  virtual void process(Long64_t entry) override;
  virtual double GetTransverseMass(const math::XYZVector tau1, const math::XYZVector tau2, const math::XYZVector muon, const math::XYZVector met);
  virtual double GetMt(const math::XYVector tau1, const math::XYVector tau2, const math::XYVector muon, const math::XYVector& met);
  virtual double GetEffectiveMass(const math::XYZTLorentzVector tau1, const math::XYZTLorentzVector tau2, const math::XYVector& met);
  virtual double GetCollinearMass(const math::XYZTLorentzVector tau1, const math::XYZTLorentzVector tau2, const math::XYVector& met);
  virtual vector<genParticle> GetAllPreviousCopies(const vector<genParticle> genParticles, genParticle genP);
  virtual vector<genParticle> GetGenParticles(const vector<genParticle> genParticles, double ptCut, double etaCut, const int pdgId, const bool isLastCopy=true, const bool hasNoDaughters=false);
  virtual vector<GenJet> GetGenJets(const GenJetCollection& genJets, std::vector<float> ptCut, std::vector<float> etaCut);
  virtual vector<GenJet> GetGenJets(vector<GenJet>& genJets, std::vector<float> ptCuts, std::vector<float> etaCuts, vector<genParticle> genParticlesToMatch);
  virtual TMatrixDSym ComputeMomentumTensor(std::vector<math::XYZTLorentzVector> jets, double r = 2.0); 
  virtual TMatrixDSym ComputeMomentumTensor2D(std::vector<math::XYZTLorentzVector> jets);
  virtual vector<float> GetMomentumTensorEigenValues(std::vector<math::XYZTLorentzVector> jets,
						     float &C,
						     float &D,
						     float &H2);
  virtual vector<float> GetMomentumTensorEigenValues2D(std::vector<math::XYZTLorentzVector> jets, 
						       float &Circularity);
  virtual vector<float> GetSphericityTensorEigenValues(std::vector<math::XYZTLorentzVector> jets,
						       float &y23, float &Sphericity, float &SphericityT, float &Aplanarity, float &Planarity, float &Y);
  virtual double GetAlphaT(std::vector<math::XYZTLorentzVector> jets,
			   float &HT,
			   float &JT,
			   float &MHT,
			   float &Centrality);

   
private:
  // Input parameters
  const double cfg_Verbose;
  const ParameterSet PSet_ElectronSelection;
  const double cfg_ElectronPtCut;
  const double cfg_ElectronEtaCut;
  const ParameterSet PSet_MuonSelection;
  const double cfg_MuonPtCut;
  const double cfg_MuonEtaCut;
  const ParameterSet PSet_TauSelection;
  const std::vector<float> cfg_TauPtCut;
  const std::vector<float> cfg_TauEtaCut;
  const ParameterSet PSet_JetSelection;
  const std::vector<float> cfg_JetPtCuts;
  const std::vector<float> cfg_JetEtaCuts;
  const DirectionalCut<int> cfg_JetNumberCut;
  const DirectionalCut<int> cfg_MuonNumberCut;
  const DirectionalCut<int> cfg_TaujetNumberCut;
  const DirectionalCut<float> cfg_HtCut;
  const ParameterSet PSet_BJetSelection;
  const std::vector<float> cfg_BJetPtCuts;
  const std::vector<float> cfg_BJetEtaCuts;
  const DirectionalCut<int> cfg_BJetNumberCut;
  const ParameterSet PSet_METSelection;
  const DirectionalCut<float> cfg_METCut;
  // METSelection PSet_METSelection;
  // TopSelection PSet_TopSelection;
  const HistogramSettings cfg_PtBinSetting;
  const HistogramSettings cfg_EtaBinSetting;
  const HistogramSettings cfg_PhiBinSetting;
  const HistogramSettings cfg_MassBinSetting;
  const HistogramSettings cfg_DeltaEtaBinSetting;
  const HistogramSettings cfg_DeltaPhiBinSetting;
  const HistogramSettings cfg_DeltaRBinSetting;
  
  Tools auxTools;
  
  // Event Counters
  Count cAllEvents;  
  Count cTrigger;
  Count cElectronVeto;
  Count cMuonVeto;
  Count cTauSelection;
  Count cJetSelection;
  Count cBJetSelection;
  Count cMETSelection;
  Count cSelected;

  // Event Variables
  WrappedTH1 *h_genMET_Et;
  WrappedTH1 *h_genHT_GenParticles;
  WrappedTH1 *h_genHT_GenJets;  
  WrappedTH1 *h_HplusMt_genJ;
  WrappedTH1 *h_HplusMt_genP;
  WrappedTH1 *h_HiggsMVis_genJ;
  WrappedTH1 *h_HiggsMVis_genP;
  WrappedTH1 *h_HiggsMEff_genJ;
  WrappedTH1 *h_HiggsMEff_genP;
  WrappedTH1 *h_HiggsMColl_genJ;
  WrappedTH1 *h_HiggsMColl_genP;

  // Event-Shape Variables
  WrappedTH1 *h_y23;
  WrappedTH1 *h_Sphericity;
  WrappedTH1 *h_SphericityT;
  WrappedTH1 *h_Y;
  WrappedTH1 *h_Aplanarity;
  WrappedTH1 *h_Planarity;
  WrappedTH1 *h_CParameter;
  WrappedTH1 *h_DParameter;
  WrappedTH1 *h_H2;
  WrappedTH1 *h_Circularity;
  WrappedTH1 *h_Centrality;
  WrappedTH1 *h_HT;
  WrappedTH1 *h_JT;
  WrappedTH1 *h_MHT;
  WrappedTH1 *h_AlphaT;

  // GenParticles

  // GenJets
  WrappedTH1 *h_GenJets_N;  
  WrappedTH1 *h_GenHadronicJets_N;  
  WrappedTH1 *h_GenTauJets_N;  
  WrappedTH1 *h_GenBJets_N;  
  WrappedTH1 *h_GenMuons_N;  
  WrappedTH1 *h_GenElectrons_N;  
  WrappedTH1 *h_GenJet1_Pt;
  WrappedTH1 *h_GenJet2_Pt;
  WrappedTH1 *h_GenJet3_Pt;
  WrappedTH1 *h_GenJet4_Pt;
  WrappedTH1 *h_GenJet5_Pt;
  WrappedTH1 *h_GenJet6_Pt;
  WrappedTH1 *h_GenJet1_Eta;
  WrappedTH1 *h_GenJet2_Eta;
  WrappedTH1 *h_GenJet3_Eta;
  WrappedTH1 *h_GenJet4_Eta;
  WrappedTH1 *h_GenJet5_Eta;
  WrappedTH1 *h_GenJet6_Eta;
  WrappedTH1 *h_MaxDiJetMass_Pt;
  WrappedTH1 *h_MaxDiJetMass_Eta;
  WrappedTH1 *h_MaxDiJetMass_Rap; 
  WrappedTH1 *h_MaxDiJetMass_Mass;
  WrappedTH1 *h_MaxDiJetMass_dR;
  WrappedTH1 *h_MaxDiJetMass_dRrap;
  WrappedTH1 *h_MaxDiJetMass_dEta;
  WrappedTH1 *h_MaxDiJetMass_dPhi;
  WrappedTH1 *h_MaxDiJetMass_dRap;

  // Taus
  WrappedTH1 *h_TauJet1_TauJet2_dEt;
  WrappedTH1 *h_TauJet1_TauJet2_dEta;
  WrappedTH1 *h_TauJet1_TauJet2_dPhi;
  WrappedTH1 *h_TauJet1_TauJet2_dR;
  WrappedTH1 *h_TauJet1_TauJet2_dQ;
  WrappedTH1 *h_TauJet1_MET_dPhi;
  WrappedTH1 *h_TauJet2_MET_dPhi;
  WrappedTH1 *h_TauJets_MET_dPhi;
  WrappedTH1 *h_TauJet1_BJet1_dR;
  WrappedTH1 *h_TauJet1_BJet1_dEta;
  WrappedTH1 *h_TauJet1_BJet1_dPhi;
  WrappedTH1 *h_TauJet2_BJet1_dR;
  WrappedTH1 *h_TauJet2_BJet1_dEta;
  WrappedTH1 *h_TauJet2_BJet1_dPhi;
  WrappedTH1 *h_TauJet1_Jet1_dR;
  WrappedTH1 *h_TauJet1_Jet1_dEta;
  WrappedTH1 *h_TauJet1_Jet1_dPhi;
  WrappedTH1 *h_TauJet2_Jet1_dR;
  WrappedTH1 *h_TauJet2_Jet1_dEta;
  WrappedTH1 *h_TauJet2_Jet1_dPhi;
  WrappedTH1 *h_TauJet1_Muon1_dR;
  WrappedTH1 *h_TauJet1_Muon1_dEta;
  WrappedTH1 *h_TauJet1_Muon1_dPhi;
  WrappedTH1 *h_TauJet2_Muon1_dR;
  WrappedTH1 *h_TauJet2_Muon1_dEta;
  WrappedTH1 *h_TauJet2_Muon1_dPhi;
  WrappedTH1 *h_Muon1_MET_dPhi;
  WrappedTH1 *h_Muon1_BJet1_dR;
  WrappedTH1 *h_Muon1_BJet1_dEta;
  WrappedTH1 *h_Muon1_BJet1_dPhi;
  WrappedTH1 *h_Muon1_Jet1_dR;
  WrappedTH1 *h_Muon1_Jet1_dEta;
  WrappedTH1 *h_Muon1_Jet1_dPhi;
  WrappedTH1 *h_Muon1_Jet2_dR;
  WrappedTH1 *h_Muon1_Jet2_dEta;
  WrappedTH1 *h_Muon1_Jet2_dPhi;
  WrappedTH1 *h_BJet1_MET_dPhi;
  WrappedTH1 *h_BJet1_Jet1_dR;
  WrappedTH1 *h_BJet1_Jet1_dEta;
  WrappedTH1 *h_BJet1_Jet1_dPhi;
  WrappedTH1 *h_BJet1_Jet2_dR;
  WrappedTH1 *h_BJet1_Jet2_dEta;
  WrappedTH1 *h_BJet1_Jet2_dPhi;
  WrappedTH1 *h_Jet1_Jet2_dR;
  WrappedTH1 *h_Jet1_Jet2_dEta;
  WrappedTH1 *h_Jet1_Jet2_dPhi;

  // TH2
  WrappedTH2 *h_S_Vs_Y;
  WrappedTH2 *h_Jet1Jet2_dEta_Vs_Jet3Jet4_dEta;
  WrappedTH2 *h_Jet1Jet2_dPhi_Vs_Jet3Jet4_dPhi;  
  WrappedTH2 *h_Jet1Jet2_dEta_Vs_Jet1Jet2_Mass;
  WrappedTH2 *h_Jet3Jet4_dEta_Vs_Jet3Jet4_Mass;
  WrappedTH2 *h_MaxDiJetMass_dEta_Vs_dPhi;
  WrappedTH2 *h_MaxDiJetMass_dRap_Vs_dPhi;  
  WrappedTH2 *h_TauJets_Pt_Vs_Pt;
  WrappedTH2 *h_TauJets_Eta_Vs_Eta;
  WrappedTH2 *h_TauJets_Phi_Vs_Phi;
  WrappedTH2 *h_TauJets_dEta_Vs_dPhi;
  WrappedTH2 *h_TauJets_Pt1_Vs_dR;
  WrappedTH2 *h_TauJets_Pt2_Vs_dR;
  WrappedTH2 *h_TauJets_Pt_Vs_MET;
  WrappedTH2 *h_TauJet1_Muon1_dPhi_Vs_TauJet2_Muon1_dPhi;
  WrappedTH2 *h_TauJet1_BJet1_dPhi_Vs_TauJet2_BJet1_dPhi;
  WrappedTH2 *h_TauJet1_MET_dPhi_Vs_TauJet2_MET_dPhi;
  WrappedTH2 *h_TauJet1_MET_dPhi_Vs_Jet1_MET_dPhi;
  WrappedTH2 *h_TauJet1_MET_dPhi_Vs_Jet2_MET_dPhi;
  WrappedTH2 *h_TauJet1_MET_dPhi_Vs_BJet1_MET_dPhi;
  WrappedTH2 *h_TauJet1_MET_dPhi_Vs_Muon1_MET_dPhi;
  WrappedTH2 *h_TauJet2_MET_dPhi_Vs_Jet1_MET_dPhi;
  WrappedTH2 *h_TauJet2_MET_dPhi_Vs_Jet2_MET_dPhi;
  WrappedTH2 *h_TauJet2_MET_dPhi_Vs_BJet1_MET_dPhi;
  WrappedTH2 *h_TauJet2_MET_dPhi_Vs_Muon1_MET_dPhi;

};

#include "Framework/interface/SelectorFactory.h"
REGISTER_SELECTOR(KinematicsBkg);

KinematicsBkg::KinematicsBkg(const ParameterSet& config, const TH1* skimCounters)
  : BaseSelector(config, skimCounters),
    cfg_Verbose(config.getParameter<bool>("verbose")),
    PSet_ElectronSelection(config.getParameter<ParameterSet>("ElectronSelection")),
    cfg_ElectronPtCut(config.getParameter<float>("ElectronSelection.electronPtCut")),  
    cfg_ElectronEtaCut(config.getParameter<float>("ElectronSelection.electronEtaCut")),
    PSet_MuonSelection(config.getParameter<ParameterSet>("MuonSelection")),
    cfg_MuonPtCut(config.getParameter<float>("MuonSelection.muonPtCut")),
    cfg_MuonEtaCut(config.getParameter<float>("MuonSelection.muonEtaCut")),
    PSet_TauSelection(config.getParameter<ParameterSet>("TauSelection")),
    cfg_TauPtCut(config.getParameter<std::vector<float>>("TauSelection.tauPtCut")),
    cfg_TauEtaCut(config.getParameter<std::vector<float>>("TauSelection.tauEtaCut")),
    PSet_JetSelection(config.getParameter<ParameterSet>("JetSelection")),
    cfg_JetPtCuts(config.getParameter<std::vector<float>>("JetSelection.jetPtCuts")),
    cfg_JetEtaCuts(config.getParameter<std::vector<float>>("JetSelection.jetEtaCuts")),
    cfg_JetNumberCut(config, "JetSelection.numberOfJetsCut"),
    cfg_MuonNumberCut(config, "MuonSelection.numberCut"), //fixme
    cfg_TaujetNumberCut(config, "TauSelection.numberOfJetsCut"), //fixme
    cfg_HtCut(config, "JetSelection.HTCut"),
    PSet_BJetSelection(config.getParameter<ParameterSet>("BJetSelection")),
    cfg_BJetPtCuts(config.getParameter<std::vector<float>>("BJetSelection.jetPtCuts")),
    cfg_BJetEtaCuts(config.getParameter<std::vector<float>>("BJetSelection.jetEtaCuts")),
    cfg_BJetNumberCut(config, "BJetSelection.numberOfBJetsCut"),
    PSet_METSelection(config.getParameter<ParameterSet>("METSelection")),
    cfg_METCut(config, "METSelection.METCut"),
    cfg_PtBinSetting(config.getParameter<ParameterSet>("CommonPlots.ptBins")),
    cfg_EtaBinSetting(config.getParameter<ParameterSet>("CommonPlots.etaBins")),
    cfg_PhiBinSetting(config.getParameter<ParameterSet>("CommonPlots.phiBins")),
    cfg_MassBinSetting(config.getParameter<ParameterSet>("CommonPlots.invMassBins")),
    cfg_DeltaEtaBinSetting(config.getParameter<ParameterSet>("CommonPlots.deltaEtaBins")),
    cfg_DeltaPhiBinSetting(config.getParameter<ParameterSet>("CommonPlots.deltaPhiBins")),
    cfg_DeltaRBinSetting(config.getParameter<ParameterSet>("CommonPlots.deltaRBins")),
    cAllEvents(fEventCounter.addCounter("All events")),
    cTrigger(fEventCounter.addCounter("Trigger")),
    cElectronVeto(fEventCounter.addCounter("e-veto")),
    cMuonVeto(fEventCounter.addCounter("#mu-veto")),
    cTauSelection(fEventCounter.addCounter("#tau-jets")),
    cJetSelection(fEventCounter.addCounter("Jets + H_{T}")),
    cBJetSelection(fEventCounter.addCounter("b-jets")),
    cMETSelection(fEventCounter.addCounter("MET")),
    cSelected(fEventCounter.addCounter("All Selections"))
{ }

void KinematicsBkg::book(TDirectory *dir) {

  // Binning
  const int nBinsPt   = 2*cfg_PtBinSetting.bins();
  const double minPt  = cfg_PtBinSetting.min();
  const double maxPt  = 2*cfg_PtBinSetting.max();

  const int nBinsEta  = 2*cfg_EtaBinSetting.bins();
  const double minEta = cfg_EtaBinSetting.min();
  const double maxEta = 2*cfg_EtaBinSetting.max();

  const int nBinsRap  = cfg_EtaBinSetting.bins();
  const double minRap = cfg_EtaBinSetting.min();
  const double maxRap = cfg_EtaBinSetting.max();

  const int nBinsPhi  = cfg_PhiBinSetting.bins();
  const double minPhi = cfg_PhiBinSetting.min();
  const double maxPhi = cfg_PhiBinSetting.max();

  const int nBinsM  = cfg_MassBinSetting.bins();
  const double minM = cfg_MassBinSetting.min();
  const double maxM = cfg_MassBinSetting.max();
  
  const int nBinsdEta  = 2*cfg_DeltaEtaBinSetting.bins();
  const double mindEta = cfg_DeltaEtaBinSetting.min();
  const double maxdEta = 2*cfg_DeltaEtaBinSetting.max();

  const int nBinsdRap  = 2*cfg_DeltaEtaBinSetting.bins();
  const double mindRap = cfg_DeltaEtaBinSetting.min();
  const double maxdRap = 2*cfg_DeltaEtaBinSetting.max();

  const int nBinsdPhi  = cfg_DeltaPhiBinSetting.bins();
  const double mindPhi = 0.0; // cfg_DeltaPhiBinSetting.min();
  const double maxdPhi = cfg_DeltaPhiBinSetting.max();

  const int nBinsdR  = 2*cfg_DeltaRBinSetting.bins();
  const double mindR = cfg_DeltaRBinSetting.min();
  const double maxdR = cfg_DeltaRBinSetting.max();

  TDirectory* th1 = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "TH1");
  TDirectory* th2 = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "TH2");
    
  // Event Variables
  h_genMET_Et         =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "genMET_Et"         , ";E_{T}^{miss} (GeV)",      200,    0.0, 1000.0);
  h_genHT_GenParticles=  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "genHT_GenParticles", ";GenP H_{T} (GeV)"  , nBinsM, minM, maxM   );
  h_genHT_GenJets     =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "genHT_GenJets"     , ";GenJ H_{T} (GeV)"  , nBinsM, minM, maxM   );
  h_HplusMt_genJ      =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "HplusMt_genJ"      , ";m_{T}(#tau_{h,1},#tau_{h,2}, #mu_{1}, E_{T}^{miss}) (GeV)" , nBinsM, minM, maxM);
  h_HplusMt_genP      =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "HplusMt_genP"      , ";m_{T}(#tau_{h,1},#tau_{h,2}, #mu_{1}, E_{T}^{miss}) (GeV)" , nBinsM, minM, maxM);
  h_HiggsMVis_genP    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "HiggsMVis_genP"    , ";m_{vis}(#tau_{h,1},#tau_{h,2}) (GeV)" , nBinsM, minM, maxM);
  h_HiggsMVis_genJ    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "HiggsMVis_genJ"    , ";m_{vis}(#tau_{h,1},#tau_{h,2}) (GeV)" , nBinsM, minM, maxM);
  h_HiggsMEff_genP    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "HiggsMEff_genP"    , ";m_{eff}(#tau_{h,1},#tau_{h,2}, E_{T}^{miss}) (GeV)" , nBinsM, minM, maxM);
  h_HiggsMEff_genJ    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "HiggsMEff_genJ"    , ";m_{eff}(#tau_{h,1},#tau_{h,2}, E_{T}^{miss}) (GeV)" , nBinsM, minM, maxM);
  h_HiggsMColl_genP   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "HiggsMColl_genP"   , ";m_{coll}(#tau_{h,1},#tau_{h,2}, E_{T}^{miss}) (GeV)" , nBinsM, minM, maxM);
  h_HiggsMColl_genJ   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "HiggsMColl_genJ"   , ";m_{coll}(#tau_{h,1},#tau_{h,2}, E_{T}^{miss}) (GeV)" , nBinsM, minM, maxM);

  // Event-Shape Variables
  h_y23         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "y23"        , ";y_{23}"        , 50, 0.0,    0.50);
  h_Sphericity  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Sphericity" , ";Sphericity"    , 20, 0.0,    1.00);
  h_SphericityT = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "SphericityT", ";Sphericity_{T}", 20, 0.0,    1.00);
  h_Y           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Y"          , ";Y"             , 50, 0.0,    0.50);
  h_S_Vs_Y      = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "S_Vs_Y"     , ";Sphericity;Y=#frac{#sqrt{3}}{2}x(Q1-Q2)", 100, 0.0, 1.0, 50, 0.0, 0.5);
  h_Aplanarity  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Aplanarity" , ";Aplanarity" , 25, 0.0, 0.5);
  h_Planarity   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Planarity"  , ";Planarity"  , 25, 0.0, 0.5);
  h_CParameter  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "CParameter" , ";C"          , 20, 0.0, 1.0);
  h_DParameter  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "DParameter" , ";D"          , 20, 0.0, 1.0);
  h_H2          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "H2"         , ";H_{2}"      , 20, 0.0, 1.0);
  h_Circularity = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Circularity", ";Circularity", 20, 0.0, 1.0);
  h_Centrality  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Centrality" , ";Centrality" , 20, 0.0, 1.0);
  h_HT          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "HT"         , ";H_{T} (GeV)", 250, 0.0, 2500.0);
  h_JT          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "JT"         , ";J_{T} (GeV)", 150, 0.0, 1500.0);
  h_MHT         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MHT"        , ";MHT (GeV)"  , 150, 0.0, 1500.0);
  h_AlphaT      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "AlphaT"     , ";#alpha_{T}" , 50, 0.0,    5.0);

  // GenParticles  
  
  // GenJets
  h_GenJets_N         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJets_N"         , ";genJet multiplicity"      , 20, 0.0, 20.0);
  h_GenHadronicJets_N = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenHadronicJets_N" , ";hadronic jet multiplicity", 10, 0.0, 10.0);
  h_GenTauJets_N      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenTauJets_N"      , ";#tau jet multiplicity"    , 10, 0.0, 10.0);
  h_GenBJets_N        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenBJets_N"        , ";b jet multiplicity"       , 10, 0.0, 10.0);
  h_GenMuons_N        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenMuons_N"        , ";#mu multiplicity"         , 10, 0.0, 10.0);
  h_GenElectrons_N    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenElectrons_N"    , ";e multiplicity"           , 10, 0.0, 10.0);
  h_GenJet1_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet1_Pt", ";p_{T} (GeV)", nBinsPt, minPt, maxPt);
  h_GenJet2_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet2_Pt", ";p_{T} (GeV)", nBinsPt, minPt, maxPt);
  h_GenJet3_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet3_Pt", ";p_{T} (GeV)", nBinsPt, minPt, maxPt);
  h_GenJet4_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet4_Pt", ";p_{T} (GeV)", nBinsPt, minPt, maxPt);
  h_GenJet5_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet5_Pt", ";p_{T} (GeV)", nBinsPt, minPt, maxPt);
  h_GenJet6_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet6_Pt", ";p_{T} (GeV)", nBinsPt, minPt, maxPt);
  h_GenJet1_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet1_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_GenJet2_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet2_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_GenJet3_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet3_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_GenJet4_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet4_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_GenJet5_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet5_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_GenJet6_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet6_Eta", ";#eta", nBinsEta, minEta, maxEta);

  h_MaxDiJetMass_Pt    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_Pt"   , ";p_{T} (GeV)"      , nBinsPt , minPt , maxPt);
  h_MaxDiJetMass_Eta   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_Eta"  , ";#eta"             , nBinsEta, minEta, maxEta);
  h_MaxDiJetMass_Rap   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_Rap"  , ";#omega"           , nBinsRap, minRap, maxRap);
  h_MaxDiJetMass_Mass  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_Mass" , ";M (GeV)"          ,      250,    0.0, +2500.0);  
  h_MaxDiJetMass_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_MaxDiJetMass_dRap  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_dRap" , ";#Delta#omega"     , nBinsdRap, mindRap, maxdRap);
  h_MaxDiJetMass_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_dPhi" , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);  
  h_MaxDiJetMass_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_dR"   , ";#DeltaR"          , nBinsdR, mindR, maxdR);  
  h_MaxDiJetMass_dRrap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_dRrap", ";#DeltaR_{#omega}" , nBinsdR, mindR, maxdR);  

  // Taujets
  h_TauJet1_TauJet2_dEt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet1_TauJet2_dEt" , ";#DeltaE_{T}(#tau_{h,1}, #tau_{h,2}) (GeV) "  ,2*nBinsPt , minPt  , maxPt);
  h_TauJet1_TauJet2_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet1_TauJet2_dEta", ";#Delta#eta(#tau_{h,1}, #tau_{h,2})"          , nBinsdEta, mindEta, maxdEta);
  h_TauJet1_TauJet2_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet1_TauJet2_dPhi", ";#Delta#phi(#tau_{h,1}, #tau_{h,2}) (rads)"   , nBinsdPhi, mindPhi, maxdPhi);
  h_TauJet1_TauJet2_dR   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet1_TauJet2_dR"  , ";#DeltaR(#tau_{h,1}, #tau_{h,2})"             , nBinsdR  , mindR  , maxdR);
  h_TauJet1_TauJet2_dQ   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet1_TauJet2_dQ"  , ";#DeltaQ(#tau_{h,1}, #tau_{h,2}) (e)"         ,       3  ,    0.0 , 3.0);
  h_TauJet1_MET_dPhi     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet1_MET_dPhi"    , ";#Delta#phi(#tau_{h,1}, E_{T}^{miss}) (rads)" , nBinsdPhi, mindPhi, maxdPhi);  
  h_TauJet2_MET_dPhi     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet2_MET_dPhi"    , ";#Delta#phi(#tau_{h,2}, E_{T}^{miss}) (rads)" , nBinsdPhi, mindPhi, maxdPhi);  
  h_TauJets_MET_dPhi     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJets_MET_dPhi"    , ";#Delta#phi(#tau_{h}'s, E_{T}^{miss}) (rads)" , nBinsdPhi, mindPhi, maxdPhi);  
  h_TauJet1_BJet1_dR     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet1_BJet1_dR"    , ";#DeltaR(#tau_{h,1},b_{1})"             , nBinsdR  , mindR , maxdR);
  h_TauJet1_BJet1_dEta   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet1_BJet1_dEta"  , ";#Delta#eta(#tau_{h,1},b_{1})"          , nBinsdEta, mindEta, maxdEta);
  h_TauJet1_BJet1_dPhi   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet1_BJet1_dPhi"  , ";#Delta#phi(#tau_{h,1},b_{1}) (rads)"   , nBinsdPhi, mindPhi, maxdPhi);  
  h_TauJet2_BJet1_dR     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet2_BJet1_dR"    , ";#DeltaR(#tau_{h,2},b_{1})"             , nBinsdR  , mindR  , maxdR);
  h_TauJet2_BJet1_dEta   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet2_BJet1_dEta"  , ";#Delta#eta(#tau_{h,2},b_{1})"          , nBinsdEta, mindEta, maxdEta); 
  h_TauJet2_BJet1_dPhi   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet2_BJet1_dPhi"  , ";#Delta#phi(#tau_{h,2},b_{1}) (rads)"   , nBinsdPhi, mindPhi, maxdPhi);  
  h_TauJet1_Jet1_dR      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet1_Jet1_dR"     , ";#DeltaR(#tau_{h,1},j_{1})"             , nBinsdR  , mindR  , maxdR);
  h_TauJet1_Jet1_dEta    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet1_Jet1_dEta"   , ";#Delta#eta(#tau_{h,1},j_{1})"          , nBinsdEta, mindEta, maxdEta); 
  h_TauJet1_Jet1_dPhi    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet1_Jet1_dPhi"   , ";#Delta#phi(#tau_{h,1},j_{1}) (rads)"   , nBinsdPhi, mindPhi, maxdPhi);  
  h_TauJet2_Jet1_dR      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet2_Jet1_dR"     , ";#DeltaR(#tau_{h,2},j_{1})"             , nBinsdR  , mindR  , maxdR);
  h_TauJet2_Jet1_dEta    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet2_Jet1_dEta"   , ";#Delta#eta(#tau_{h,2},j_{1})"          , nBinsdEta, mindEta, maxdEta); 
  h_TauJet2_Jet1_dPhi    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet2_Jet1_dPhi"   , ";#Delta#phi(#tau_{h,2},j_{1}) (rads)"   , nBinsdPhi, mindPhi, maxdPhi);  
  h_TauJet1_Muon1_dR     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet1_Muon1_dR"    , ";#DeltaR(#tau_{h,1},#mu_{1})"           , nBinsdR  , mindR  , maxdR);
  h_TauJet1_Muon1_dEta   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet1_Muon1_dEta"  , ";#Delta#eta(#tau_{h,1},#mu_{1})"        , nBinsdEta, mindEta, maxdEta); 
  h_TauJet1_Muon1_dPhi   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet1_Muon1_dPhi"  , ";#Delta#phi(#tau_{h,1},#mu_{1}) (rads)" , nBinsdPhi, mindPhi, maxdPhi);  
  h_TauJet2_Muon1_dR     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet2_Muon1_dR"    , ";#DeltaR(#tau_{h,2},#mu_{1})"           , nBinsdR  , mindR  , maxdR);
  h_TauJet2_Muon1_dEta   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet2_Muon1_dEta"  , ";#Delta#eta(#tau_{h,2},#mu_{1})"        , nBinsdEta, mindEta, maxdEta); 
  h_TauJet2_Muon1_dPhi   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "TauJet2_Muon1_dPhi"  , ";#Delta#phi(#tau_{h,2},#mu_{1}) (rads)" , nBinsdPhi, mindPhi, maxdPhi);
  h_Muon1_MET_dPhi       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Muon1_MET_dPhi"      , ";#Delta#phi(#mu_{1}, E_{T}^{miss}) (rads)", nBinsdPhi, mindPhi, maxdPhi);
  h_Muon1_BJet1_dR       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Muon1_BJet1_dR"      , ";#DeltaR(#mu_{1}, b_{1})"            , nBinsdR  , mindR  , maxdR);
  h_Muon1_BJet1_dEta     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Muon1_BJet1_dEta"    , ";#Delta#eta(#mu_{1}, b_{1})"         , nBinsdEta, mindEta, maxdEta); 
  h_Muon1_BJet1_dPhi     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Muon1_BJet1_dPhi"    , ";#Delta#phi(#mu_{1}, b_{1}) (rads)"  , nBinsdPhi, mindPhi, maxdPhi);  
  h_Muon1_Jet1_dR        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Muon1_Jet1_dR"       , ";#DeltaR(#mu_{1}, j_{1})"            , nBinsdR  , mindR  , maxdR);
  h_Muon1_Jet1_dEta      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Muon1_Jet1_dEta"     , ";#Delta#eta(#mu_{1}, j_{1})"         , nBinsdEta, mindEta, maxdEta); 
  h_Muon1_Jet1_dPhi      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Muon1_Jet1_dPhi"     , ";#Delta#phi(#mu_{1}, j_{1}) (rads)"  , nBinsdPhi, mindPhi, maxdPhi);  
  h_Muon1_Jet2_dR        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Muon1_Jet2_dR"       , ";#DeltaR(#mu_{1}, j_{2})"            , nBinsdR  , mindR  , maxdR);
  h_Muon1_Jet2_dEta      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Muon1_Jet2_dEta"     , ";#Delta#eta(#mu_{1}, j_{2})"         , nBinsdEta, mindEta, maxdEta); 
  h_Muon1_Jet2_dPhi      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Muon1_Jet2_dPhi"     , ";#Delta#phi(#mu_{1}, j_{2}) (rads)"  , nBinsdPhi, mindPhi, maxdPhi);  
  h_BJet1_MET_dPhi       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "BJet1_MET_dPhi"      , ";#Delta#phi(#b_{1}, E_{T}^{miss}) (rads)" , nBinsdPhi, mindPhi, maxdPhi);
  h_BJet1_Jet1_dR        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "BJet1_Jet1_dR"       , ";#DeltaR(#b_{1}, j_{1})"             , nBinsdR  , mindR  , maxdR);
  h_BJet1_Jet1_dEta      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "BJet1_Jet1_dEta"     , ";#Delta#eta(#b_{1}, j_{1})"          , nBinsdEta, mindEta, maxdEta); 
  h_BJet1_Jet1_dPhi      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "BJet1_Jet1_dPhi"     , ";#Delta#phi(#b_{1}, j_{1}) (rads)"   , nBinsdPhi, mindPhi, maxdPhi);  
  h_BJet1_Jet2_dR        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "BJet1_Jet2_dR"       , ";#DeltaR(#b_{1}, j_{2})"             , nBinsdR  , mindR  , maxdR);
  h_BJet1_Jet2_dEta      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "BJet1_Jet2_dEta"     , ";#Delta#eta(#b_{1}, j_{2})"          , nBinsdEta, mindEta, maxdEta); 
  h_BJet1_Jet2_dPhi      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "BJet1_Jet2_dPhi"     , ";#Delta#phi(#b_{1}, j_{2}) (rads)"   , nBinsdPhi, mindPhi, maxdPhi);  
  h_Jet1_Jet2_dR         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Jet1_Jet2_dR"        , ";#DeltaR(#j_{1}, j_{2})"             , nBinsdR  , mindR  , maxdR);
  h_Jet1_Jet2_dEta       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Jet1_Jet2_dEta"      , ";#Delta#eta(#j_{1}, j_{2})"          , nBinsdEta, mindEta, maxdEta); 
  h_Jet1_Jet2_dPhi       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "Jet1_Jet2_dPhi"      , ";#Delta#phi(#j_{1}, j_{2}) (rads)"   , nBinsdPhi, mindPhi, maxdPhi);  


  // TH2 
  h_Jet1Jet2_dEta_Vs_Jet3Jet4_dEta = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Jet1Jet2_dEta_Vs_Jet3Jet4_dEta", 
								";#Delta#eta(j_{1},j_{2});#Delta#eta(j_{3},j_{4})", nBinsdEta, mindEta, maxdEta, nBinsdEta, mindEta, maxdEta);
  h_Jet1Jet2_dPhi_Vs_Jet3Jet4_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Jet1Jet2_dPhi_Vs_Jet3Jet4_dPhi", 
								";#Delta#phi(j_{1},j_{2}) (rads);#Delta#phi(j_{3},j_{4}) (rads)", nBinsdPhi, mindPhi, maxdPhi, nBinsdPhi, mindPhi, maxdPhi);
  h_Jet1Jet2_dEta_Vs_Jet1Jet2_Mass = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Jet1Jet2_dEta_Vs_Jet1Jet2_Mass", 
								";#Delta#eta(j_{1},j_{2});M(j_{1},j_{2}) (GeV)", nBinsdEta, mindEta, maxdEta, nBinsM, minM, maxM);
  h_Jet3Jet4_dEta_Vs_Jet3Jet4_Mass = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Jet3Jet4_dEta_Vs_Jet3Jet4_Mass", 
								";#Delta#eta(j_{3},j_{4});M(j_{4},j_{4}) (GeV)", nBinsdEta, mindEta, maxdEta, nBinsM, minM, maxM);

  h_MaxDiJetMass_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "MaxDiJetMass_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)"  , nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_MaxDiJetMass_dRap_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "MaxDiJetMass_dRap_Vs_dPhi", ";#Delta#omega;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);

  h_TauJets_Pt_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJets_Pt_Vs_Pt"    , ";p_{T}^{#tau_{h,1}} (GeV);p_{T}^{#tau_{h,2}} (GeV)", nBinsPt, minPt, maxPt, nBinsPt, minPt, maxPt);
  h_TauJets_Pt_Vs_MET    = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJets_Pt_Vs_MET"   , ";p_{T}^{#tau_{h,1}} (GeV);E_{T}^{miss} (GeV)"      , nBinsPt, minPt, maxPt, 200, 0.0, 1000.0);
  h_TauJets_Eta_Vs_Eta   = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJets_Eta_Vs_Eta"  , ";#eta^{#tau_{h,1}};#eta^{#tau_{h,2}}"              , nBinsEta, minEta, maxEta, nBinsEta, minEta, maxEta);
  h_TauJets_Phi_Vs_Phi   = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJets_Phi_Vs_Phi"  , ";#phi^{#tau_{h,1}} (rads);#phi^{#tau_{h,2}} (rads)", nBinsPhi, minPhi, maxPhi, nBinsPhi, minPhi, maxPhi);
  h_TauJets_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJets_dEta_Vs_dPhi", ";#Delta#eta(#tau_{h,1}, #tau_{h,2});#Delta#phi(#tau_{h,1}, #tau_{h,2}) (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_TauJets_Pt1_Vs_dR    = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJets_Pt1_Vs_dR"   , ";p_{T}^{#tau_{h,1}} (GeV);#DeltaR(#tau_{h,1}, #tau_{h,2})", nBinsPt, minPt, maxPt, nBinsdR, mindR, maxdR);  
  h_TauJets_Pt2_Vs_dR    = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJets_Pt2_Vs_dR"   , ";p_{T}^{#tau_{h,2}} (GeV);#DeltaR(#tau_{h,1}, #tau_{h,2})", nBinsPt, minPt, maxPt, nBinsdR, mindR, maxdR);  

  h_TauJet1_Muon1_dPhi_Vs_TauJet2_Muon1_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJet1_Muon1_dPhi_Vs_TauJet2_Muon1_dPhi", ";#Delta#phi(#tau_{h,1},#mu_{1}) (rads);#Delta#phi(#tau_{h,2}, #mu_{1}) (rads)", nBinsPhi, minPhi, maxPhi, nBinsPhi, minPhi, maxPhi);

  h_TauJet1_BJet1_dPhi_Vs_TauJet2_BJet1_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJet1_BJet1_dPhi_Vs_TauJet2_BJet1_dPhi", ";#Delta#phi(#tau_{h,1}, b_{1}) (rads);#Delta#phi(#tau_{h,2}, b_{1}) (rads)", nBinsPhi, minPhi, maxPhi, nBinsPhi, minPhi, maxPhi);

  h_TauJet1_MET_dPhi_Vs_TauJet2_MET_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJet1_MET_dPhi_Vs_TauJet2_MET_dPhi", ";#Delta#phi(#tau_{h,1},E_{T}^{miss}) (rads);#Delta#phi(#tau_{h,2}, E_{T}^{miss}) (rads)", nBinsPhi, minPhi, maxPhi, nBinsPhi, minPhi, maxPhi);

  h_TauJet1_MET_dPhi_Vs_Jet1_MET_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJet1_MET_dPhi_Vs_Jet1_MET_dPhi", ";#Delta#phi(#tau_{h,1},E_{T}^{miss}) (rads);#Delta#phi(jet_{1}, E_{T}^{miss}) (rads)", nBinsPhi, minPhi, maxPhi, nBinsPhi, minPhi, maxPhi); 
   
  h_TauJet1_MET_dPhi_Vs_Jet2_MET_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJet1_MET_dPhi_Vs_Jet2_MET_dPhi", ";#Delta#phi(#tau_{h,1},E_{T}^{miss}) (rads);#Delta#phi(jet_{2}, E_{T}^{miss}) (rads)", nBinsPhi, minPhi, maxPhi, nBinsPhi, minPhi, maxPhi); 

  h_TauJet1_MET_dPhi_Vs_BJet1_MET_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJet1_MET_dPhi_Vs_BJet1_MET_dPhi", ";#Delta#phi(#tau_{h,1},E_{T}^{miss}) (rads);#Delta#phi(bjet_{1}, E_{T}^{miss}) (rads)", nBinsPhi, minPhi, maxPhi, nBinsPhi, minPhi, maxPhi); 

  h_TauJet1_MET_dPhi_Vs_Muon1_MET_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJet1_MET_dPhi_Vs_Muon1_MET_dPhi", ";#Delta#phi(#tau_{h,1},E_{T}^{miss}) (rads);#Delta#phi(#mu_{1}, E_{T}^{miss}) (rads)", nBinsPhi, minPhi, maxPhi, nBinsPhi, minPhi, maxPhi); 

  h_TauJet2_MET_dPhi_Vs_Jet1_MET_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJet2_MET_dPhi_Vs_Jet1_MET_dPhi", ";#Delta#phi(#tau_{h,2},E_{T}^{miss}) (rads);#Delta#phi(jet_{1}, E_{T}^{miss}) (rads)", nBinsPhi, minPhi, maxPhi, nBinsPhi, minPhi, maxPhi); 
   
  h_TauJet2_MET_dPhi_Vs_Jet2_MET_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJet2_MET_dPhi_Vs_Jet2_MET_dPhi", ";#Delta#phi(#tau_{h,1},E_{T}^{miss}) (rads);#Delta#phi(jet_{2}, E_{T}^{miss}) (rads)", nBinsPhi, minPhi, maxPhi, nBinsPhi, minPhi, maxPhi); 

  h_TauJet2_MET_dPhi_Vs_BJet1_MET_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJet2_MET_dPhi_Vs_BJet1_MET_dPhi", ";#Delta#phi(#tau_{h,1},E_{T}^{miss}) (rads);#Delta#phi(bjet_{1}, E_{T}^{miss}) (rads)", nBinsPhi, minPhi, maxPhi, nBinsPhi, minPhi, maxPhi); 

  h_TauJet2_MET_dPhi_Vs_Muon1_MET_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "TauJet2_MET_dPhi_Vs_Muon1_MET_dPhi", ";#Delta#phi(#tau_{h,1},E_{T}^{miss}) (rads);#Delta#phi(#mu_{1}, E_{T}^{miss}) (rads)", nBinsPhi, minPhi, maxPhi, nBinsPhi, minPhi, maxPhi); 


  return;
}

void KinematicsBkg::setupBranches(BranchManager& branchManager) {
  fEvent.setupBranches(branchManager);
}


void KinematicsBkg::process(Long64_t entry) {

  if ( !fEvent.isMC() ) return;
  
  // Create MCools object
  MCTools mcTools(fEvent);
  
  // Increment Counter
  cAllEvents.increment();

  // if (entry != 9081) return;
  if (cfg_Verbose) std::cout << "\n=== Event " << entry << std::endl;


  //================================================================================================
  // 1) Apply trigger
  //================================================================================================
  // if (cfg_Verbose) std::cout << "=== Trigger" << std::endl;
  // if ( !(fEvent.passTriggerDecision()) ) return;
  // cTrigger.increment();


  //================================================================================================
  // 2) MET filters (to remove events with spurious sources of fake MET)       
  //================================================================================================


  //================================================================================================
  // 3) Primarty Vertex (Check that a PV exists)
  //================================================================================================


  //================================================================================================
  // 4) Electron veto (fully hadronic + orthogonality)  
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== Electron veto" << std::endl;
  vector<genParticle> selectedElectrons = GetGenParticles(fEvent.genparticles().getGenParticles(), cfg_ElectronPtCut, cfg_ElectronEtaCut, 11, true, false);
  if (0)
    {
      std::cout << "\nnElectrons = " << selectedElectrons.size() << std::endl;
      for (auto& p: selectedElectrons) std::cout << "\tpT = " << p.pt() << " (GeV), eta = " << p.eta() << ", phi = " << p.phi() << " (rads), status = " << p.status() << std::endl;
    }
  if ( selectedElectrons.size() > 0 ) return;
  cElectronVeto.increment();


  //================================================================================================
  // 5) Muon selection
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== Muon selection" << std::endl;
  vector<genParticle> selectedMuons = GetGenParticles(fEvent.genparticles().getGenParticles(), cfg_MuonPtCut, cfg_MuonEtaCut, 13, true, false);
  if (0)
    {
      std::cout << "nMuons = " << selectedMuons.size() << std::endl;
      for (auto& p: selectedMuons) std::cout << "\tpT = " << p.pt() << " (GeV), eta = " << p.eta() << ", phi = " << p.phi() << " (rads), status = " << p.status() << std::endl;
    } 
  // if ( selectedMuons.size() != 1 ) return;
  if ( !cfg_MuonNumberCut.passedCut(selectedMuons.size() ) ) return;
  cMuonVeto.increment();

  //================================================================================================
  // 6) Tau selection
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== Tau selection" << std::endl;
  vector<genParticle> selectedTaus = GetGenParticles(fEvent.genparticles().getGenParticles(), cfg_TauPtCut.at(0), cfg_TauEtaCut.at(0), 15, true, false); // fixme - dirty fix
  for (auto& p: selectedTaus)
    {
      // std::cout << "\tpT = " << p.pt() << " (GeV), eta = " << p.eta() << ", phi = " << p.phi() << " (rads), status = " << p.status() << std::endl;
      for (unsigned int i=0; i < p.daughters().size(); i++) 
	{
	  int d_pdgId = fEvent.genparticles().getGenParticles()[p.daughters().at(i)].pdgId(); //fixme
	  if (0) cout << "tau daughter = " << d_pdgId << endl;

	  // Only consider hadronic tau decays
	  if( (abs(d_pdgId) == 11) || (abs(d_pdgId) == 12)  ||
	      (abs(d_pdgId) == 13) || (abs(d_pdgId) == 14) ) return;
	}
    }
  if (0) std::cout << "nTaus = " << selectedTaus.size() << std::endl;
  
  // Match tau-leptons with GenJets
  vector<GenJet> selectedJets    = GetGenJets(fEvent.genjets(), cfg_JetPtCuts, cfg_JetEtaCuts); // all hadronic jets
  vector<GenJet> selectedHadJets = GetGenJets(fEvent.genjets(), cfg_JetPtCuts, cfg_JetEtaCuts); // all hadronic jets (without taujets and bjets)
  vector<GenJet> selectedTauJets = GetGenJets(selectedHadJets, cfg_TauPtCut, cfg_TauEtaCut, selectedTaus);
  if ( !cfg_TaujetNumberCut.passedCut(selectedTauJets.size() ) ) return;
  cTauSelection.increment();

  if (cfg_Verbose) 
    {
      std::cout << "nTauJets = " << selectedTauJets.size() << std::endl;
      for (auto& p: selectedTauJets) std::cout << "\tpT = " << p.pt() << " (GeV), eta = " << p.eta() << ", phi = " << p.phi() << " (rads)" << std::endl;
      std::cout << "" << std::endl;
    }


  //===============================================================================================
  // 7) Jet Selection
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== Jet Selection" << std::endl;

  if (cfg_Verbose)
    {
      std::cout << "nJets = " << selectedJets.size() << std::endl;
      for (auto& p: selectedJets) std::cout << "\tpT = " << p.pt() << " (GeV), eta = " << p.eta() << ", phi = " << p.phi() << " (rads)" << std::endl;
    }
  if (!cfg_JetNumberCut.passedCut(selectedJets.size())) return;

  // HT Selection
  double genJ_HT = 0.0;
  std::vector<math::XYZTLorentzVector> selJets_p4;
  math::XYZTLorentzVector jet_p4;
  for(auto jet: selectedJets) 
    {
        jet_p4 = jet.p4();
	genJ_HT += jet.pt();
	selJets_p4.push_back( jet_p4 );
    }

  if ( !cfg_HtCut.passedCut(genJ_HT) ) return;
  cJetSelection.increment();


  //================================================================================================
  // 8) BJet Selection
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== BJet Selection" << std::endl;
  vector<genParticle> selectedBQuarks_All = GetGenParticles(fEvent.genparticles().getGenParticles(), 10, 3, 5, true, false);
  vector<genParticle> selectedBQuarks;
  for (auto& p: selectedBQuarks_All)
    {
      // Skip bquark from hard process (incoming bquark)
      genParticle m = fEvent.genparticles().getGenParticles()[p.mothers().at(0)]; 
      selectedBQuarks.push_back(p);
    }
  std::sort( selectedBQuarks.begin(), selectedBQuarks.end(), PtComparator() );  

  if (0)
    {
      std::cout << "\nnBQuarks = " << selectedBQuarks.size() << std::endl;
      for (auto& p: selectedBQuarks)
	{
	  mcTools.PrintGenParticle(p);
	}
      // for (auto& p: selectedBQuarks) std::cout << "\tpT = " << p.pt() << " (GeV), eta = " << p.eta() << ", phi = " << p.phi() << " (rads), status = " << p.status() << std::endl;
    }

  // Match b-quarks with GenJets
  vector<GenJet> selectedBJets = GetGenJets(selectedHadJets, cfg_BJetPtCuts, cfg_BJetEtaCuts, selectedBQuarks);
  if (cfg_Verbose) 
    {
      std::cout << "nBJets = " << selectedBJets.size() << std::endl;
      for (auto& p: selectedBJets) std::cout << "\tpT = " << p.pt() << " (GeV), eta = " << p.eta() << ", phi = " << p.phi() << " (rads)" << std::endl;
      std::cout << "" << std::endl;
    }
  if (!cfg_BJetNumberCut.passedCut(selectedBJets.size())) return;
  cBJetSelection.increment();


  //================================================================================================
  // 9) BJet SF
  //================================================================================================
  

  //================================================================================================
  // 10) MET selection 
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== MET Selection" << std::endl;
  if (!cfg_METCut.passedCut(fEvent.genMET().et())) return;
  if (0) std::cout << "=== MET = " << fEvent.genMET().et() << std::endl;      
  cMETSelection.increment();


  //================================================================================================
  // 11) Topology selection 
  //================================================================================================
  float C, D, H2;
  float Circularity;
  float y23, Sphericity, SphericityT, Aplanarity, Planarity, Y; // functions to return values when properly implemented
  float HT, JT, MHT, Centrality;
  vector<float> a = GetMomentumTensorEigenValues(selJets_p4, C, D, H2);
  vector<float> b = GetMomentumTensorEigenValues2D(selJets_p4, Circularity);
  vector<float> c = GetSphericityTensorEigenValues(selJets_p4, y23, Sphericity, SphericityT, Aplanarity, Planarity, Y);
  double alphaT   = GetAlphaT(selJets_p4, HT, JT, MHT, Centrality);
  
  if (cfg_Verbose)
    {
      std::cout << "\nalphaT = " << alphaT << ", MHT = " << MHT << ", Centrality = " << Centrality << ", y23 = " << y23 
		<< ", Sphericity = " << Sphericity << ", SphericityT = " << SphericityT << ", Aplanarity = " << Aplanarity
		<< ", Planarity = " << Planarity << ", Y = " << Y << std::endl;
    }   
   

  //================================================================================================
  // All Selections
  //================================================================================================
  // Clean the hadronic jets (many events whereby a hadronic jet is a muon)
  vector<GenJet> selectedHadJetsCleaned;
  for(auto j: selectedHadJets) 
    {
      double dEta  = abs( selectedMuons.at(0).eta() - j.eta() );
      double dPhi  = abs(ROOT::Math::VectorUtil::DeltaPhi( selectedMuons.at(0).p4(), j.p4() ));
      double dR    = ROOT::Math::VectorUtil::DeltaR(selectedMuons.at(0).p4(), j.p4() );
      bool bIsMuon = (dEta < 5e-2) && (dPhi < 5e-2);
      if (bIsMuon) continue;
      if (0) cout << "dEta = " << dEta << ", dPhi = " << dPhi << ", dR = " << dR << endl;
      selectedHadJetsCleaned.push_back(j);
    }

  if (cfg_Verbose) std::cout << "=== All Selections" << std::endl;
  cSelected.increment();
  if (cfg_Verbose) 
    {
      std::cout << "" << std::endl;
      std::cout << "nElectrons = " << selectedElectrons.size() << std::endl;
      std::cout << "nMuons     = " << selectedMuons.size() << std::endl;
      std::cout << "nTauJets   = " << selectedTauJets.size() << " (nTaus = " << selectedTaus.size() << ")" <<  std::endl;
      std::cout << "nJets      = " << selectedJets.size() << std::endl;
      std::cout << "nHadJets   = " << selectedHadJetsCleaned.size() << std::endl;
      std::cout << "nBJets     = " << selectedBJets.size() << " (nBQuarks = " << selectedBQuarks.size() << ")" << std::endl;
      std::cout << "MET        = " << fEvent.genMET().et() << " GeV" << std::endl;
    }
  

  // Definitions
  const math::XYVector tau1p = selectedTaus.at(0).p2();
  const math::XYVector tau2p = selectedTaus.at(1).p2();
  const math::XYVector tau1j = selectedTauJets.at(0).p2();
  const math::XYVector tau2j = selectedTauJets.at(1).p2();
  const math::XYVector muon  = selectedMuons.at(0).p2();
  const math::XYVector met   = (const math::XYVector) fEvent.genMET().p2();

  // Transverse Mass
  double mT_genP = GetMt(tau1p, tau2p, muon, met);
  double mT_genJ = GetMt(tau1j, tau2j, muon, met);

  // Visible mass
  double mVis_genP = ( selectedTaus.at(0).p4() + selectedTaus.at(1).p4() ).M();
  double mVis_genJ = ( selectedTauJets.at(0).p4() + selectedTauJets.at(1).p4() ).M();

  // Effective visible mass (Provides better discrimination against backgrounds than "visible mass")
  double mEff_genP = GetEffectiveMass((const math::XYZTLorentzVector) selectedTaus.at(0).p4()   , (const math::XYZTLorentzVector) selectedTaus.at(1).p4()   , met);
  double mEff_genJ = GetEffectiveMass((const math::XYZTLorentzVector) selectedTauJets.at(0).p4(), (const math::XYZTLorentzVector) selectedTauJets.at(1).p4(), met);

  // Collinear Approximation - mVis/sqrt(x1 x2) where x1 = pvis1/(pvis1 + pmiss1) is= momentum fraction carried away by visible tau decay products
  double mColl_genP = GetCollinearMass((const math::XYZTLorentzVector) selectedTaus.at(0).p4(), (const math::XYZTLorentzVector) selectedTaus.at(1).p4(), met);
  double mColl_genJ = GetCollinearMass((const math::XYZTLorentzVector) selectedTauJets.at(0).p4(), (const math::XYZTLorentzVector) selectedTauJets.at(1).p4(), met);

  if (0) cout << "mT_genJ = " << mT_genJ << ", mT_genP = " << mT_genP << ", mEff_genP = " << mEff_genP << ", mEff_genJ = " << mEff_genJ << ", mColl_genP = " << mColl_genP << ", mColl_genJ = " << mColl_genJ << std::endl;


  // Fill Histograms
  h_HplusMt_genJ   ->Fill( mT_genJ );
  h_HplusMt_genP   ->Fill( mT_genP );
  h_HiggsMVis_genJ ->Fill( mVis_genJ );
  h_HiggsMVis_genP ->Fill( mVis_genP );
  h_HiggsMEff_genP ->Fill( mEff_genP );
  h_HiggsMEff_genJ ->Fill( mEff_genJ );
  h_HiggsMColl_genP->Fill( mColl_genP );
  h_HiggsMColl_genJ->Fill( mColl_genJ );

  h_genMET_Et    ->Fill(fEvent.genMET().et()); 
  h_genHT_GenJets->Fill(genJ_HT);
  h_y23          ->Fill(y23);
  h_Sphericity   ->Fill(Sphericity);
  h_SphericityT  ->Fill(SphericityT);
  h_Y            ->Fill(Y);
  h_S_Vs_Y       ->Fill(Sphericity, Y);
  h_Aplanarity   ->Fill(Aplanarity);
  h_Planarity    ->Fill(Planarity);
  h_CParameter   ->Fill(C);
  h_DParameter   ->Fill(D);
  h_H2           ->Fill(H2);
  h_Circularity  ->Fill(Circularity);
  h_Centrality   ->Fill(Centrality);
  h_HT           ->Fill(HT);
  h_JT           ->Fill(JT);
  h_MHT          ->Fill(MHT);
  h_AlphaT       ->Fill(alphaT);


  ///////////////////////////////////////////////////////////////////////////
  // GenParticles
  ///////////////////////////////////////////////////////////////////////////
  // Define the table
  Table table("Evt | Index | PdgId | Status | Charge | Pt | Eta | Phi | E | Vertex (mm) | D0 (mm) | Lxy (mm) | Mom | Daughters", "Text"); //LaTeX or Text

  int row = 0;
  double genP_HT = 0.0;

  // For-loop: GenParticles
  for (auto& p: fEvent.genparticles().getGenParticles()) {

    // Particle properties
    short genP_index     = p.index();
    int genP_pdgId       = p.pdgId();
    int genP_status      = p.status();
    double genP_pt       = p.pt();
    double genP_eta      = p.eta();
    double genP_phi      = p.phi();
    double genP_energy   = p.e();
    int genP_charge      = p.charge();

    // // Associated genParticles
    std::vector<genParticle> genP_daughters;
    for (unsigned int i=0; i < p.daughters().size(); i++) genP_daughters.push_back(fEvent.genparticles().getGenParticles()[p.daughters().at(i)]);
    std::vector<unsigned int> genDaus_index;
    std::vector<unsigned int> genDaus_pdgId;
    
    for (unsigned int i=0; i < genP_daughters.size(); i++) 
      {
	if (genP_daughters.at(i).index() < 0) continue;
	genDaus_index.push_back(genP_daughters.at(i).index());
	genDaus_pdgId.push_back(genP_daughters.at(i).pdgId());
      }

    // Properties that need to be calculated
    bool bIsLastCopy = std::find(genDaus_pdgId.begin(), genDaus_pdgId.end(), genP_pdgId) == genDaus_pdgId.end();
    if (!bIsLastCopy) continue;

    if ( !mcTools.IsLepton(p.pdgId()) )  genP_HT += p.p4().Et();

    // Print genParticle properties or decay tree ?
    // mcTools.PrintGenParticle(p, row==0);
    // mcTools.PrintGenDaughters(p);

    // Add table rows
    table.AddRowColumn(row, auxTools.ToString(entry)           );
    table.AddRowColumn(row, auxTools.ToString(genP_index)      );
    table.AddRowColumn(row, auxTools.ToString(genP_pdgId)      );
    table.AddRowColumn(row, auxTools.ToString(genP_status)     );
    table.AddRowColumn(row, auxTools.ToString(genP_charge)     );
    table.AddRowColumn(row, auxTools.ToString(genP_pt , 3)     );
    table.AddRowColumn(row, auxTools.ToString(genP_eta, 4)     );
    table.AddRowColumn(row, auxTools.ToString(genP_phi, 3)     );
    table.AddRowColumn(row, auxTools.ToString(genP_energy, 3)  );
    if (genDaus_index.size() < 6)
      {
	table.AddRowColumn(row, auxTools.ConvertIntVectorToString(genDaus_index) );
      }
    else table.AddRowColumn(row, ".. Too many .." );
    row++;
	           
    // if ( mcTools.IsQuark(d.pdgId()) )
    // if ( mcTools.IsLepton(d.pdgId()) )
	
	
  }// for-loop: genParticles
  
  ///////////////////////////////////////////////////////////////////////////
  // GenJets (selected)
  ///////////////////////////////////////////////////////////////////////////
  if (cfg_Verbose) std::cout << "=== GenJets" << std::endl;
  std::vector<math::XYZTLorentzVector> v_dijet_p4;
  std::vector<double> v_dijet_masses;
  std::vector<double> v_dijet_dR;
  std::vector<double> v_dijet_dRrap;
  std::vector<double> v_dijet_dEta;
  std::vector<double> v_dijet_dPhi;
  std::vector<double> v_dijet_dRap;
  

  double maxDijetMass_mass;
  math::XYZTLorentzVector maxDijetMass_p4;
  int maxDijetMass_pos;
  double maxDijetMass_dR;
  double maxDijetMass_dRrap;
  double maxDijetMass_dEta;
  double maxDijetMass_dPhi;
  double maxDijetMass_dRap;
  double maxDijetMass_rapidity;
  
  int iJet = 0;
  if (selJets_p4.size() > 1) {

    // For-loop: All selected jets 
    for (size_t i=0; i < selJets_p4.size(); i++)
      {
	iJet++;
	double genJ_Pt  = selJets_p4.at(i).pt();
	double genJ_Eta = selJets_p4.at(i).eta();
	// double genJ_Rap = mcTools.GetRapidity(selJets_p4.at(i));
	
	if (iJet==1)
	  {
	    h_GenJet1_Pt -> Fill( genJ_Pt  );
	    h_GenJet1_Eta-> Fill( genJ_Eta );
	  }
	else if (iJet==2)
	  {
	    h_GenJet2_Pt -> Fill( genJ_Pt  );
	    h_GenJet2_Eta-> Fill( genJ_Eta );
	  }
	else if (iJet==3)
	  {
	    h_GenJet3_Pt -> Fill( genJ_Pt  );
	    h_GenJet3_Eta-> Fill( genJ_Eta );
	  }
	else if (iJet==4)
	  {
	    h_GenJet4_Pt -> Fill( genJ_Pt  );
	    h_GenJet4_Eta-> Fill( genJ_Eta );
	  }
	else if (iJet==5)
	  {
	    h_GenJet5_Pt -> Fill( genJ_Pt  );
	    h_GenJet5_Eta-> Fill( genJ_Eta );
	  }
	else if (iJet==6)
	  {
	    h_GenJet6_Pt -> Fill( genJ_Pt  );
	    h_GenJet6_Eta-> Fill( genJ_Eta );
	  }
	else{}
	
	// For-loop: All selected jets (nested)
	for (size_t j=i+1; j < selJets_p4.size(); j++)
	  {
	    math::XYZTLorentzVector p4_i = selJets_p4.at(i);
	    math::XYZTLorentzVector p4_j = selJets_p4.at(j);
	    math::XYZTLorentzVector p4   = p4_i + p4_j;
	    double rap_i = mcTools.GetRapidity(p4_i);
	    double rap_j = mcTools.GetRapidity(p4_j);
	    double dR    = ROOT::Math::VectorUtil::DeltaR(p4_i, p4_j);
	    double dRap  = abs(rap_i - rap_j);
	    double dEta  = abs(p4_i.eta() - p4_j.eta());
	    double dPhi  = abs(ROOT::Math::VectorUtil::DeltaPhi(p4_i, p4_j));
	    
	    v_dijet_p4.push_back( p4 );
	    v_dijet_masses.push_back( p4.mass() );
	    v_dijet_dR.push_back( dR );
	    v_dijet_dRrap.push_back( sqrt( pow(dRap, 2) + pow(dPhi, 2) ) ); 
	    v_dijet_dEta.push_back( dEta ); 
	    v_dijet_dRap.push_back( dRap );
	    v_dijet_dPhi.push_back( dPhi );

	  }
      }

    // MaxDiJet: DiJet combination with largest mass
    maxDijetMass_pos      = std::max_element(v_dijet_masses.begin(), v_dijet_masses.end()) - v_dijet_masses.begin();
    maxDijetMass_mass     = v_dijet_masses.at(maxDijetMass_pos);
    maxDijetMass_p4       = v_dijet_p4.at(maxDijetMass_pos);
    maxDijetMass_dR       = v_dijet_dR.at(maxDijetMass_pos);
    maxDijetMass_dRrap    = v_dijet_dRrap.at(maxDijetMass_pos);
    maxDijetMass_dEta     = v_dijet_dEta.at(maxDijetMass_pos);
    maxDijetMass_dPhi     = v_dijet_dPhi.at(maxDijetMass_pos);
    maxDijetMass_dRap     = v_dijet_dRap.at(maxDijetMass_pos);
    maxDijetMass_rapidity = mcTools.GetRapidity(maxDijetMass_p4);
  
  }// if (selJets_p4.size() > 1) {
    
 

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Fill histograms
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (cfg_Verbose) std::cout << "=== Fill histograms (Calculations)" << std::endl;
  double taujets_dEta   = abs( selectedTauJets.at(0).eta() - selectedTauJets.at(1).eta() );
  double taujets_dPhi   = abs(ROOT::Math::VectorUtil::DeltaPhi( selectedTauJets.at(0).p4(), selectedTauJets.at(1).p4() ));
  double taujets_dR     = ROOT::Math::VectorUtil::DeltaR( selectedTauJets.at(0).p4(), selectedTauJets.at(1).p4() );
  double tau1_MET_dPhi  = abs( auxTools.DeltaPhi( selectedTauJets.at(0).phi(), fEvent.genMET().Phi()) );
  double tau2_MET_dPhi  = abs( auxTools.DeltaPhi( selectedTauJets.at(1).phi(), fEvent.genMET().Phi()) );
  double bjet1_MET_dPhi = abs( auxTools.DeltaPhi( selectedBJets.at(0).phi(), fEvent.genMET().Phi()) );
  double muon1_MET_dPhi = abs( auxTools.DeltaPhi( selectedMuons.at(0).phi(), fEvent.genMET().Phi()) );
  double jet1_MET_dPhi  = +999.9;
  double jet2_MET_dPhi  = +999.9;
  double jet1_jet2_dEta = +999.9;
  double jet1_jet2_dPhi = +999.9;
  double jet1_jet2_dR   = +999.9;
  double tau1_muon1_dEta = abs(selectedTauJets.at(0).eta() - selectedMuons.at(0).eta() );
  double tau1_muon1_dPhi = abs(ROOT::Math::VectorUtil::DeltaPhi( selectedTauJets.at(0).p4(), selectedMuons.at(0).p4() ));
  double tau1_muon1_dR   = ROOT::Math::VectorUtil::DeltaR( selectedTauJets.at(0).p4(), selectedMuons.at(0).p4());
  double tau2_muon1_dEta = abs(selectedTauJets.at(1).eta() - selectedMuons.at(0).eta() );
  double tau2_muon1_dPhi = abs(ROOT::Math::VectorUtil::DeltaPhi( selectedTauJets.at(1).p4(), selectedMuons.at(0).p4() ));
  double tau2_muon1_dR   = ROOT::Math::VectorUtil::DeltaR( selectedTauJets.at(1).p4(), selectedMuons.at(0).p4());
  double tau1_bjet1_dEta = abs(selectedTauJets.at(0).eta() - selectedBJets.at(0).eta() );
  double tau1_bjet1_dPhi = abs(ROOT::Math::VectorUtil::DeltaPhi( selectedTauJets.at(0).p4(), selectedBJets.at(0).p4() ));
  double tau1_bjet1_dR   = ROOT::Math::VectorUtil::DeltaR( selectedTauJets.at(0).p4(), selectedBJets.at(0).p4());
  double tau2_bjet1_dEta = abs(selectedTauJets.at(1).eta() - selectedBJets.at(0).eta() );
  double tau2_bjet1_dPhi = abs(ROOT::Math::VectorUtil::DeltaPhi( selectedTauJets.at(1).p4(), selectedBJets.at(0).p4() ));
  double tau2_bjet1_dR   = ROOT::Math::VectorUtil::DeltaR( selectedTauJets.at(1).p4(), selectedBJets.at(0).p4());



  if (cfg_Verbose) std::cout << "=== Fill histograms (TH1)" << std::endl;
  h_genHT_GenParticles->Fill(genP_HT);
  h_genHT_GenJets     ->Fill(genJ_HT);
  h_genMET_Et         ->Fill(fEvent.genMET().et()); 
  h_GenJets_N         ->Fill(selectedJets.size());
  h_GenHadronicJets_N ->Fill(selectedHadJetsCleaned.size());
  h_GenTauJets_N      ->Fill(selectedTauJets.size());
  h_GenBJets_N        ->Fill(selectedBJets.size());
  h_GenMuons_N        ->Fill(selectedMuons.size());
  h_GenElectrons_N    ->Fill(selectedElectrons.size());
  h_MaxDiJetMass_Mass ->Fill( maxDijetMass_mass     );
  h_MaxDiJetMass_Pt   ->Fill( maxDijetMass_p4.pt()  );
  h_MaxDiJetMass_Eta  ->Fill( maxDijetMass_p4.eta() );
  h_MaxDiJetMass_Rap  ->Fill( maxDijetMass_rapidity );
  h_MaxDiJetMass_dR   ->Fill( maxDijetMass_dR       );
  h_MaxDiJetMass_dRrap->Fill( maxDijetMass_dRrap    );
  h_MaxDiJetMass_dEta ->Fill( maxDijetMass_dEta     );
  h_MaxDiJetMass_dPhi ->Fill( maxDijetMass_dPhi     );
  h_MaxDiJetMass_dRap ->Fill( maxDijetMass_dRap     );
 if (selectedHadJetsCleaned.size() > 0)
    {
      jet1_MET_dPhi = abs( auxTools.DeltaPhi( selectedHadJetsCleaned.at(0).phi(), fEvent.genMET().Phi()) );
    }
 if (selectedHadJetsCleaned.size() > 1) 
   {
     jet2_MET_dPhi  = abs( auxTools.DeltaPhi( selectedHadJetsCleaned.at(1).phi(), fEvent.genMET().Phi()) );
     jet1_jet2_dEta = abs( selectedHadJetsCleaned.at(0).eta() - selectedHadJetsCleaned.at(1).eta() );
     jet1_jet2_dPhi = abs( auxTools.DeltaPhi( selectedHadJetsCleaned.at(0).phi(), selectedHadJetsCleaned.at(1).phi()) );
     jet1_jet2_dR   = ROOT::Math::VectorUtil::DeltaR( selectedHadJetsCleaned.at(0).p4(), selectedHadJetsCleaned.at(1).p4());
   }
    

  if (cfg_Verbose) std::cout << "=== Fill histograms (TH2)" << std::endl;
  if (selJets_p4.size() > 3) 
    {
      double jet1_Eta   = selJets_p4.at(0).eta();
      double jet2_Eta   = selJets_p4.at(1).eta();
      double jet3_Eta   = selJets_p4.at(2).eta();
      double jet4_Eta   = selJets_p4.at(3).eta();
      double jet1_Phi   = selJets_p4.at(0).phi();
      double jet2_Phi   = selJets_p4.at(1).phi();
      double jet3_Phi   = selJets_p4.at(2).phi();
      double jet4_Phi   = selJets_p4.at(3).phi();
      double jet12_dEta = abs(jet1_Eta - jet2_Eta);
      double jet34_dEta = abs(jet3_Eta - jet4_Eta);
      double jet12_dPhi = abs(jet1_Phi - jet2_Phi);
      double jet34_dPhi = abs(jet3_Phi - jet4_Phi);
      double jet12_Mass = (selJets_p4.at(0) + selJets_p4.at(1)).mass(); 
      double jet34_Mass = (selJets_p4.at(2) + selJets_p4.at(3)).mass();

      h_Jet1Jet2_dEta_Vs_Jet3Jet4_dEta ->Fill(jet12_dEta, jet34_dEta);
      h_Jet1Jet2_dPhi_Vs_Jet3Jet4_dPhi ->Fill(jet12_dPhi, jet34_dPhi);
      h_Jet1Jet2_dEta_Vs_Jet1Jet2_Mass ->Fill(jet12_dEta, jet12_Mass);
      h_Jet3Jet4_dEta_Vs_Jet3Jet4_Mass ->Fill(jet34_dEta, jet34_Mass);
    }

  h_MaxDiJetMass_dEta_Vs_dPhi->Fill( maxDijetMass_dEta, maxDijetMass_dPhi );
  h_MaxDiJetMass_dRap_Vs_dPhi->Fill( maxDijetMass_dRap, maxDijetMass_dPhi );
  h_TauJets_Pt_Vs_Pt         ->Fill( selectedTauJets.at(0).pt() , selectedTauJets.at(1).pt() );
  h_TauJets_Pt_Vs_MET        ->Fill( (selectedTauJets.at(0).p4() + selectedTauJets.at(1).p4()).pt(), fEvent.genMET().et() );
  h_TauJets_Eta_Vs_Eta       ->Fill( selectedTauJets.at(0).eta(), selectedTauJets.at(1).eta() );
  h_TauJets_Phi_Vs_Phi       ->Fill( selectedTauJets.at(0).phi(), selectedTauJets.at(1).phi() );
  h_TauJets_dEta_Vs_dPhi     ->Fill( taujets_dEta, taujets_dPhi);
  h_TauJets_Pt1_Vs_dR        ->Fill( selectedTauJets.at(0).pt(), taujets_dR );
  h_TauJets_Pt2_Vs_dR        ->Fill( selectedTauJets.at(1).pt(), taujets_dR );  
  h_TauJet1_TauJet2_dEt      ->Fill( selectedTauJets.at(0).p4().Et() - selectedTauJets.at(1).p4().Et() );
  h_TauJet1_TauJet2_dEta     ->Fill( abs(selectedTauJets.at(0).eta() - selectedTauJets.at(1).eta()) );
  h_TauJet1_TauJet2_dPhi     ->Fill( taujets_dPhi);
  h_TauJet1_TauJet2_dR       ->Fill( taujets_dR);
  h_TauJet1_TauJet2_dQ       ->Fill( abs(selectedTaus.at(0).charge() - selectedTaus.at(1).charge()) );
  h_TauJet1_MET_dPhi         ->Fill( tau1_MET_dPhi);
  h_TauJet2_MET_dPhi         ->Fill( tau2_MET_dPhi);
  h_Jet1_Jet2_dEta           ->Fill( jet1_jet2_dEta);
  h_Jet1_Jet2_dPhi           ->Fill( jet1_jet2_dPhi);
  h_Jet1_Jet2_dR             ->Fill( jet1_jet2_dR);

  h_TauJet1_Muon1_dPhi_Vs_TauJet2_Muon1_dPhi->Fill(tau1_muon1_dPhi, tau2_muon1_dPhi);
  h_TauJet1_BJet1_dPhi_Vs_TauJet2_BJet1_dPhi->Fill(tau1_bjet1_dPhi, tau2_bjet1_dPhi);
  h_TauJet1_MET_dPhi_Vs_TauJet2_MET_dPhi->Fill(tau1_MET_dPhi, tau2_MET_dPhi);
  h_TauJet1_MET_dPhi_Vs_BJet1_MET_dPhi  ->Fill(tau1_MET_dPhi, bjet1_MET_dPhi);
  h_TauJet2_MET_dPhi_Vs_BJet1_MET_dPhi  ->Fill(tau2_MET_dPhi, bjet1_MET_dPhi);
  h_TauJet1_MET_dPhi_Vs_Muon1_MET_dPhi  ->Fill(tau1_MET_dPhi, muon1_MET_dPhi);
  h_TauJet2_MET_dPhi_Vs_Muon1_MET_dPhi  ->Fill(tau2_MET_dPhi, muon1_MET_dPhi);
 if (selectedHadJetsCleaned.size() > 0)
   {
     h_TauJet1_MET_dPhi_Vs_Jet1_MET_dPhi ->Fill(tau1_MET_dPhi, jet1_MET_dPhi);
     h_TauJet2_MET_dPhi_Vs_Jet1_MET_dPhi ->Fill(tau2_MET_dPhi, jet1_MET_dPhi);
   }
 if (selectedHadJetsCleaned.size() > 1)
   {
     h_TauJet1_MET_dPhi_Vs_Jet2_MET_dPhi ->Fill(tau1_MET_dPhi, jet2_MET_dPhi);
     h_TauJet2_MET_dPhi_Vs_Jet2_MET_dPhi ->Fill(tau2_MET_dPhi, jet2_MET_dPhi);
   }
 
  h_TauJets_MET_dPhi    ->Fill( abs( auxTools.DeltaPhi( (selectedTauJets.at(0).p4() + selectedTauJets.at(1).p4()).Phi(), fEvent.genMET().Phi())) );
  h_TauJet1_BJet1_dEta  ->Fill( tau1_bjet1_dEta );
  h_TauJet1_BJet1_dPhi  ->Fill( tau1_bjet1_dPhi );
  h_TauJet1_BJet1_dR    ->Fill( tau1_bjet1_dR   );
  h_TauJet2_BJet1_dEta  ->Fill( tau2_bjet1_dEta );
  h_TauJet2_BJet1_dPhi  ->Fill( tau2_bjet1_dPhi );
  h_TauJet2_BJet1_dR    ->Fill( tau2_bjet1_dR   );
  h_TauJet1_Muon1_dEta  ->Fill( tau1_muon1_dEta );
  h_TauJet1_Muon1_dPhi  ->Fill( tau1_muon1_dPhi );
  h_TauJet1_Muon1_dR    ->Fill( tau1_muon1_dR   );
  h_TauJet2_Muon1_dEta  ->Fill( tau2_muon1_dEta );
  h_TauJet2_Muon1_dPhi  ->Fill( tau2_muon1_dPhi );
  h_TauJet2_Muon1_dR    ->Fill( tau2_muon1_dR   );

  // Ensure a hadronic jet exists
  if (selectedHadJetsCleaned.size() > 0)
    {
      // std::cout << "selectedHadJetsCleaned.size() = " << selectedHadJetsCleaned.size() << std::endl;
      h_TauJet1_Jet1_dR     ->Fill(ROOT::Math::VectorUtil::DeltaR( selectedTauJets.at(0).p4(), selectedHadJetsCleaned.at(0).p4() ));
      h_TauJet1_Jet1_dEta   ->Fill(abs(selectedTauJets.at(1).eta() - selectedHadJetsCleaned.at(0).eta() ));
      h_TauJet1_Jet1_dPhi   ->Fill(abs(ROOT::Math::VectorUtil::DeltaPhi( selectedTauJets.at(0).p4(), selectedHadJetsCleaned.at(0).p4() )));
      h_TauJet2_Jet1_dR     ->Fill(ROOT::Math::VectorUtil::DeltaR( selectedTauJets.at(1).p4(), selectedHadJetsCleaned.at(0).p4() )); 
      h_TauJet2_Jet1_dEta   ->Fill(abs(selectedTauJets.at(1).eta() - selectedHadJetsCleaned.at(0).eta() ));
      h_TauJet2_Jet1_dPhi   ->Fill(abs(ROOT::Math::VectorUtil::DeltaPhi( selectedTauJets.at(1).p4(), selectedHadJetsCleaned.at(0).p4() )));
    }

  h_Muon1_MET_dPhi   ->Fill( abs( auxTools.DeltaPhi( selectedMuons.at(0).p4().Phi(), fEvent.genMET().Phi())) );
  h_BJet1_MET_dPhi   ->Fill( abs( auxTools.DeltaPhi( selectedBJets.at(0).p4().Phi(), fEvent.genMET().Phi())) );
  h_Muon1_BJet1_dEta ->Fill(abs(selectedMuons.at(0).eta() - selectedBJets.at(0).eta() ));
  h_Muon1_BJet1_dPhi ->Fill(abs(ROOT::Math::VectorUtil::DeltaPhi( selectedMuons.at(0).p4(), selectedBJets.at(0).p4() )));  
  h_Muon1_BJet1_dR   ->Fill(ROOT::Math::VectorUtil::DeltaR( selectedMuons.at(0).p4(), selectedBJets.at(0).p4() )); 

  if (selectedHadJetsCleaned.size() > 0)
    {
      h_Muon1_Jet1_dEta ->Fill(abs(selectedMuons.at(0).eta() - selectedHadJetsCleaned.at(0).eta() ));
      h_Muon1_Jet1_dPhi ->Fill(abs(ROOT::Math::VectorUtil::DeltaPhi( selectedMuons.at(0).p4(), selectedHadJetsCleaned.at(0).p4() )));
      h_Muon1_Jet1_dR   ->Fill(ROOT::Math::VectorUtil::DeltaR( selectedMuons.at(0).p4(), selectedHadJetsCleaned.at(0).p4() )); 

      h_BJet1_Jet1_dEta ->Fill(abs(selectedBJets.at(0).eta() - selectedHadJetsCleaned.at(0).eta() ));
      h_BJet1_Jet1_dPhi ->Fill(abs(ROOT::Math::VectorUtil::DeltaPhi( selectedBJets.at(0).p4(), selectedHadJetsCleaned.at(0).p4() )));  
      h_BJet1_Jet1_dR   ->Fill(ROOT::Math::VectorUtil::DeltaR( selectedBJets.at(0).p4(), selectedHadJetsCleaned.at(0).p4() )); 
    }
  if (selectedHadJetsCleaned.size() > 1)
    {
      h_Muon1_Jet2_dEta ->Fill(abs(selectedMuons.at(0).eta() - selectedHadJetsCleaned.at(1).eta() ));
      h_Muon1_Jet2_dPhi ->Fill(abs(ROOT::Math::VectorUtil::DeltaPhi( selectedMuons.at(0).p4(), selectedHadJetsCleaned.at(1).p4() )));  
      h_Muon1_Jet2_dR   ->Fill(ROOT::Math::VectorUtil::DeltaR( selectedMuons.at(0).p4(), selectedHadJetsCleaned.at(1).p4() )); 

      h_BJet1_Jet2_dEta ->Fill(abs(selectedBJets.at(0).eta() - selectedHadJetsCleaned.at(1).eta() ));
      h_BJet1_Jet2_dPhi ->Fill(abs(ROOT::Math::VectorUtil::DeltaPhi( selectedBJets.at(0).p4(), selectedHadJetsCleaned.at(1).p4() )));  
      h_BJet1_Jet2_dR   ->Fill(ROOT::Math::VectorUtil::DeltaR( selectedBJets.at(0).p4(), selectedHadJetsCleaned.at(1).p4() )); 
    }

  return;
}



vector<float> KinematicsBkg::GetMomentumTensorEigenValues(std::vector<math::XYZTLorentzVector> jets,
							  float &C,
							  float &D,
							  float &H2){

  // Tensor required for calculation of: Sphericity, Aplanarity, Planarity
  // Need all particles in event to calculate kinematic variables. Use all tracks (ch. particles) instead.
  // Links:
  // https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_24/doc/html/d5/d29/EventShapeVariables_8h_source.html
  // https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_24/doc/html/dd/d99/classEventShapeVariables.html

  // Get the Linear Momentum Tensor
  TMatrixDSym MomentumTensor = ComputeMomentumTensor(jets, 1.0);

  // Find the Momentum-Tensor EigenValues (Q1, Q2, Q3)
  TMatrixDSymEigen eigen(MomentumTensor);
  TVectorD eigenvals = eigen.GetEigenValues();

  // Store & Sort the eigenvalues
  vector<float> eigenvalues(3);
  eigenvalues.at(0) = eigenvals(0); // Q1
  eigenvalues.at(1) = eigenvals(1); // Q2
  eigenvalues.at(2) = eigenvals(2); // Q3
  sort( eigenvalues.begin(), eigenvalues.end(), std::greater<float>() );

  // Calculate the eigenvalues sum (Requirement: Q1 + Q2 + Q3 = 1)
  float eigenSum = std::accumulate(eigenvalues.begin(), eigenvalues.end(), 0);
  if ( (eigenSum - 1.0) > 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that Q1+Q2+Q3=1. Found that Q1+Q2+Q3 = " << eigenSum << ", instead.";
    }

  // Save the final eigenvalues
  float Q1 = eigenvalues.at(0);
  float Q2 = eigenvalues.at(1);
  float Q3 = eigenvalues.at(2);

  // Sanity check on eigenvalues: Q1 >= Q2 >= Q3 (Q1 >= 0)
  bool bQ1Zero = (Q1 >= 0.0);
  bool bQ1Q2   = (Q1 >= Q2);
  bool bQ2Q3   = (Q2 >= Q3);
  bool bInequality = bQ1Zero * bQ1Q2 * bQ2Q3;
  
  if ( !(bInequality) )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that eigenvalues are ordered as Q1 >= Q2 >= Q3 (Q1 >= 0). Q1 = " << Q1 << ", Q2 = " << Q2 << ", Q3 = " << Q3;
    }

  // Calculate the linear combinations C and D
  C  = 3*(Q1*Q2 + Q1*Q3 + Q2*Q3); // Used to measure the 3-jet structure. Vanishes for perfece 2-jet event. Related to the 2nd Fox-Wolfram Moment (H2)
  D  = 27*Q1*Q2*Q3; // Used to measure the 4-jet structure. Vanishes for a planar event
  H2 = 1-C; // The C-measure is related to the second Fox-Wolfram moment (see below), $C = 1 - H_2$.

  // C
  if ( abs(C-1.0) > 1 + 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that the quantity C satisfies the inequality: 0.0 <= C <= 1.0. Found that C = " << C;
    }

  // D
  if ( abs(C-1.0) > 1 + 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that the quantity C satisfies the inequality: 0.0 <= D <= 1.0. Found that D = " << D;
    }

  // 2nd Fox-Wolfram Moment
  if ( abs(H2-1.0) > 1 + 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that the 2nd Fox-Wolfram Moment (H2) satisfies the inequality: 0.0 <= H2 <= 1.0. Found that H2 = " << H2;
    }
  
  if (0)
    {

      Table vars("Variable | Value | Allowed Range | Definition", "Text"); //LaTeX or Text
      vars.AddRowColumn(0, "C");
      vars.AddRowColumn(0, auxTools.ToString(C) );
      vars.AddRowColumn(0, "0.0 <= C <= 1.0");
      vars.AddRowColumn(0, "C = 3 x (Q1Q2 + Q1Q3 + Q2Q3");
      //
      vars.AddRowColumn(1, "D");
      vars.AddRowColumn(1, auxTools.ToString(D) );
      vars.AddRowColumn(1, "0.0 <= D <= 1.0");
      vars.AddRowColumn(1, "D = 27 x Q1 x Q2 x Q3");
      //
      vars.AddRowColumn(2, "2nd F-W Moment");
      vars.AddRowColumn(2, auxTools.ToString(H2) );
      vars.AddRowColumn(2, "0.0 <= H2 <= 1.0 ");
      vars.AddRowColumn(2, "H2 = 1-C");
      vars.Print();
    }

  return eigenvalues;
}


vector<float> KinematicsBkg::GetMomentumTensorEigenValues2D(std::vector<math::XYZTLorentzVector> jets,
							 float &Circularity){

  // For Circularity, the momentum tensor is the 2¡Á2 submatrix of Mjk, normalized by the sum of pT instead by the sum of
  // This matrix has two eigenvalues Qi with 0 < Q1 < Q2. The following definition for the circularity C has been used:
  //  C = 2 ¡Á min (Q1,Q2) / (Q1 +Q2)
  // The circularity C is especially interesting for hadron colliders because it only uses the momentum values in x and y direction transverse
  // to the beam line. So C is a two dimensional event shape variable and is therefore independent from a boost along z. 
  // In addition, the normalization by the sum of the particle momenta makes C highly independent from energy calibration effects  (systematic uncertainty).
  // C takes small values for linear and high values for circular events. 
  
  // Find the Momentum-Tensor EigenValues (E1, E2)
  TMatrixDSym MomentumTensor = ComputeMomentumTensor2D(jets);
  TMatrixDSymEigen eigen(MomentumTensor);
  TVectorD eigenvals = eigen.GetEigenValues();

  // Store & Sort the eigenvalues
  vector<float> eigenvalues(2);
  eigenvalues.at(0) = eigenvals[0]; // Q1
  eigenvalues.at(1) = eigenvals[1]; // Q2
  sort( eigenvalues.begin(), eigenvalues.end(), std::greater<float>() );

  // Save the final eigenvalues
  float Q1 = eigenvalues.at(0);
  float Q2 = eigenvalues.at(1);

  // Sanity check on eigenvalues: (Q1 > Q2)
  if (Q1 == 0) return eigenvalues;

  bool bInequality = (Q1 > 0 && Q1 > Q2);
  if ( !(bInequality) )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that eigenvalues are ordered as Q1 >= Q2. Found Q1 = " << Q1 << ", Q2 " << Q2;
    }
  
  // Calculate circularity
  Circularity = 2*std::min(Q1, Q2)/(Q1+Q2); // is this definition correct?

  if (0)
    {      
      Table vars("Variable | Value | Allowed Range | Definition", "Text"); //LaTeX or Text
      vars.AddRowColumn(0, "Circularity");
      vars.AddRowColumn(0, auxTools.ToString(Circularity) );
      vars.AddRowColumn(0, "0.0 <= C <= 1.0 ");
      vars.AddRowColumn(0, "C = 2 ¡Á min (Q1,Q2)/(Q1+Q2)");
      vars.Print();
    }
  
  return eigenvalues;
}


vector<float> KinematicsBkg::GetSphericityTensorEigenValues(std::vector<math::XYZTLorentzVector> jets,
							 float &y23, float &Sphericity, float &SphericityT, float &Aplanarity, 
							 float &Planarity, float &Y){

  // C, D parameters
  // Need all particles in event to calculate kinematic variables. Use all tracks (ch. particles) instead.
  // Links:
  // http://home.fnal.gov/~mrenna/lutp0613man2/node234.html

  // Sanity check: at least 3 jets (for 3rd-jet resolution)
  if( (jets.size()) < 3 )
    {
      vector<float> zeros(3, 0);
      return zeros;
    }

  // Sort the jets by pT (leading jet first)
  std::sort( jets.begin(), jets.end(), PtComparator() );  

  // Get the Sphericity Tensor
  TMatrixDSym SphericityTensor = ComputeMomentumTensor(jets, 2.0);

  // Find the Momentum-Tensor EigenValues (Q1, Q2, Q3)
  TMatrixDSymEigen eigen(SphericityTensor);
  TVectorD eigenvals = eigen.GetEigenValues();

  // Store & Sort the eigenvalues
  vector<float> eigenvalues(3);
  eigenvalues.at(0) = eigenvals(0); // Q1
  eigenvalues.at(1) = eigenvals(1); // Q2
  eigenvalues.at(2) = eigenvals(2); // Q3
  sort( eigenvalues.begin(), eigenvalues.end(), std::greater<float>() );

  // Calculate the eigenvalues sum (Requirement: Q1 + Q2 + Q3 = 1)
  float eigenSum = std::accumulate(eigenvalues.begin(), eigenvalues.end(), 0);
  if ( (eigenSum - 1.0) > 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that Q1+Q2+Q3=1. Found that Q1+Q2+Q3 = " << eigenSum << ", instead.";
    }

  // Save the final eigenvalues
  float Q1 = eigenvalues.at(0);
  float Q2 = eigenvalues.at(1);
  float Q3 = eigenvalues.at(2);

  // Sanity check on eigenvalues: Q1 >= Q2 >= Q3 (Q1 >= 0)
  bool bQ1Zero = (Q1 >= 0.0);
  bool bQ1Q2   = (Q1 >= Q2);
  bool bQ2Q3   = (Q2 >= Q3);
  bool bInequality = bQ1Zero * bQ1Q2 * bQ2Q3;
  
  if ( !(bInequality) )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that eigenvalues are ordered as Q1 >= Q2 >= Q3 (Q1 >= 0)";
    }


  // Calculate the event-shape variables
  float pT3Squared  = pow(jets.at(2).Pt(), 2);
  float HT2Squared  = pow(jets.at(0).Pt() + jets.at(1).Pt(), 2);
  y23         = pT3Squared/HT2Squared;
  Sphericity  = -1.0;
  SphericityT = -1.0;
  Aplanarity  = -1.0;
  Planarity   = -1.0;
  Y = (sqrt(3.0)/2.0)*(Q2-Q3); // (Since Q1>Q2, then my Q1 corresponds to Q3 when Q's are reversed-ordered). Calculate the Y (for Y-S plane)

  // Check the value of the third-jet resolution
  if (abs(y23-0.25) > 0.25 + 1e-4)
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that y23 satisfies the inequality: 0.0 <= y23 <= 0.25. Found that y23 = " << y23;
    }
  
  // Calculate the Sphericity (0 <= S <= 1). S~0 for a 2-jet event, and S~1 for an isotropic one
  Sphericity = 1.5*(Q2 + Q3); 
  if ( abs(Sphericity-1.0) > 1 + 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that Sphericity (S) satisfies the inequality: 0.0 <= S <= 1.0. Found that S = " << Sphericity;
    }

  // Calculate the Sphericity (0 <= S <= 1). S~0 for a 2-jet event, and S~1 for an isotropic one
  SphericityT = 2.0*Q2/(Q1 + Q2); 
  if ( abs(SphericityT-1.0) > 1 + 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that Transverse Sphericity (ST) satisfies the inequality: 0.0 <= ST <= 1.0. Found that ST = " << SphericityT;
    }

  // Calculate the Aplanarity (0 <= A <= 0.5).  It measures the transverse momentum component out of the event plane
  // A~0 for a planar event, A~0.5 for an isotropic one
  Aplanarity = 1.5*(Q3);
  if ( abs(Aplanarity-0.5) > 0.5 + 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that Aplanarity (A) satisfies the inequality: 0.0 <= A <= 0.5";
    }

  // Calculate the Aplanarity (0 <= P <= 0.5)
  Planarity  = (2.0/3.0)*(Sphericity-2*Aplanarity);
  if ( abs(Planarity-0.5) > 0.5 + 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that Planarity (P) satisfies the inequality: 0.0 <= P <= 0.5";
    }

  if (0)
    {

      Table vars("Variable | Value | Allowed Range | Definition", "Text"); //LaTeX or Text
      vars.AddRowColumn(0, "y23");
      vars.AddRowColumn(0, auxTools.ToString(y23) );
      vars.AddRowColumn(0, "0.0 <= y23 <= 0.25");
      vars.AddRowColumn(0, "y23 = pow(jet3_Pt, 2) / pow(jet1_Pt + jet2_Pt, 2)" );
      //
      vars.AddRowColumn(1, "Sphericity");
      vars.AddRowColumn(1, auxTools.ToString(Sphericity) );
      vars.AddRowColumn(1, "0.0 <= S <= 1.0");
      vars.AddRowColumn(1, "S = 1.5 x (Q2 + Q3)");
      //
      vars.AddRowColumn(2, "Sphericity (T)");
      vars.AddRowColumn(2, auxTools.ToString(SphericityT) );
      vars.AddRowColumn(2, "0.0 <= S (T) <= 1.0");
      vars.AddRowColumn(2, "S (T) = (2 x Q2)/(Q1 + Q2)");
      //
      vars.AddRowColumn(3, "Aplanarity");
      vars.AddRowColumn(3, auxTools.ToString(Aplanarity) );
      vars.AddRowColumn(3, "0.0 <= A <= 0.5 ");
      vars.AddRowColumn(3, "A = 1.5 x Q3");
      //
      vars.AddRowColumn(4, "Planarity");
      vars.AddRowColumn(4, auxTools.ToString(Planarity) );
      vars.AddRowColumn(4, "0.0 <= P <= 0.5 ");
      vars.AddRowColumn(4, "P (2/3) x (S - 2A)");
      //
      vars.AddRowColumn(5, "Y");
      vars.AddRowColumn(5, auxTools.ToString(Y) );
      vars.AddRowColumn(5, "");
      vars.AddRowColumn(5, "Y = sqrt(3)/2 x (Q1 - Q2)");
      vars.Print();
    }
    
  return eigenvalues;
  
}

TMatrixDSym KinematicsBkg::ComputeMomentumTensor(std::vector<math::XYZTLorentzVector> jets, double r)
{

  // r = 2: Corresponds to sphericity tensor (Get: Sphericity, Aplanarity, Planarity, ..)
  // r = 1: Corresponds to linear measures (Get: C, D, Second Fox-Wolfram moment, ...)
  TMatrixDSym momentumTensor(3);
  momentumTensor.Zero();
  
  if (r!=1.0 && r!=2.0)
    {
      throw hplus::Exception("LogicError") << "Invalid value r-value in computing the Momentum Tensor (r=" << r << ").  Supported valued are r=2.0 and r=1.0.";
    }
  
  // Sanity Check
  if ( jets.size() < 2 )
    {
      return momentumTensor;
    }
  
  // Declare the Matrix normalisation (sum of momentum magnitutes to power r). That is: sum(|p|^{r})
  double normalisation = 0.0;
  double trace = 0.0;

  // For-loop: Jets
  for (auto& jet: jets){
    
    // Get the |p|^2 of the jet
    double p2 = pow(jet.P(), 2); // jet.P(); 
    
    // For r=2, use |p|^{2}, for r=1 use |p| as the momentum weight
    double pR = ( r == 2.0 ) ? p2 : TMath::Power(p2, 0.5*r);
    
    // For r=2, use |1|, for r=1 use   (|p|^{2})^{-0.5} = |p|^{2 (-1/2)} = |p|^{-1} = 1.0/|p|
    double pRminus2 = ( r == 2.0 ) ? 1.0 : TMath::Power(p2, 0.5*r - 1.0); // 
    
    // Add pR to the matrix normalisation factor
    normalisation += pR;
       
    // Fill the momentum (r=1) or  sphericity (r=2) tensor (Must be symmetric: Mij = Mji)
    momentumTensor(0,0) += pRminus2*jet.px()*jet.px(); // xx
    momentumTensor(0,1) += pRminus2*jet.px()*jet.py(); // xy
    momentumTensor(0,2) += pRminus2*jet.px()*jet.pz(); // xz
    
    momentumTensor(1,0) += pRminus2*jet.py()*jet.px(); // yx
    momentumTensor(1,1) += pRminus2*jet.py()*jet.py(); // yy
    momentumTensor(1,2) += pRminus2*jet.py()*jet.pz(); // yz
    
    momentumTensor(2,0) += pRminus2*jet.pz()*jet.px(); // zx
    momentumTensor(2,1) += pRminus2*jet.pz()*jet.py(); // zy
    momentumTensor(2,2) += pRminus2*jet.pz()*jet.pz(); // zz
    
  }// for (auto& jet: jets){

  // Normalise the tensors to have unit trace (Mxx + Myy + Mzz = 1)
  momentumTensor *= (1.0/normalisation);
  trace = momentumTensor(0,0) + momentumTensor(1,1) + momentumTensor(2,2);
  
  // Print the tensor
  if (0)
    {
      std::cout << "\nMomentum Tensor (r = " << r << "):" << std::endl;
      Table tensor(" |  | ", "Text"); //LaTeX or Text    
      tensor.AddRowColumn(0, auxTools.ToString( momentumTensor(0,0) ) );
      tensor.AddRowColumn(0, auxTools.ToString( momentumTensor(0,1) ) );
      tensor.AddRowColumn(0, auxTools.ToString( momentumTensor(0,2) ) );
      tensor.AddRowColumn(1, auxTools.ToString( momentumTensor(1,0) ) );
      tensor.AddRowColumn(1, auxTools.ToString( momentumTensor(1,1) ) );
      tensor.AddRowColumn(1, auxTools.ToString( momentumTensor(1,2) ) );
      tensor.AddRowColumn(2, auxTools.ToString( momentumTensor(2,0) ) );
      tensor.AddRowColumn(2, auxTools.ToString( momentumTensor(2,1) ) );
      tensor.AddRowColumn(2, auxTools.ToString( momentumTensor(2,2) ) );
      tensor.AddRowColumn(3, "");
      tensor.AddRowColumn(4, "Normalisation");
      tensor.AddRowColumn(4, auxTools.ToString(normalisation));
      tensor.AddRowColumn(5, "IsSymmetric");
      tensor.AddRowColumn(5, auxTools.ToString(momentumTensor.IsSymmetric()));
      tensor.AddRowColumn(6, "Determinant");
      tensor.AddRowColumn(6, auxTools.ToString(momentumTensor.Determinant()));
      tensor.AddRowColumn(7, "Trace");
      tensor.AddRowColumn(7, auxTools.ToString(trace));
      tensor.Print(false);
    }

  if ( abs(trace-1.0) > 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that the Momentum-Tensor (r = " << r << ") Trace is 1.0. Found that abs(trace-1) = " << abs(trace-1) << ", instead.";
    }

  return momentumTensor;
}



TMatrixDSym KinematicsBkg::ComputeMomentumTensor2D(std::vector<math::XYZTLorentzVector> jets)
{

  TMatrixDSym momentumTensor(3);
  momentumTensor.Zero();
  
  // Sanity Check
  if ( jets.size() < 2 )
    {
      return momentumTensor;
    }
  
  // Declare the Matrix normalisation (sum of momentum magnitutes to power r). That is: sum(|p|^{r})
  double normalisation = 0.0;
  double trace = 0.0;

  // For-loop: Jets
  for (auto& jet: jets){
    
    // Get the pT
    double pT = jet.Pt();
    
    // Add pT to the matrix normalisation factor
    normalisation += pow(pT,2);
       
    // Fill the two-dimensional momentum tensor
    momentumTensor(0,0) += jet.px()*jet.px(); // xx
    momentumTensor(0,1) += jet.px()*jet.py(); // xy
    momentumTensor(1,0) += jet.py()*jet.px(); // yx
    momentumTensor(1,1) += jet.py()*jet.py(); // yy
    
  }// for (auto& jet: jets){

  // Normalise tensor to get the normalised 2-d momentum tensor
  momentumTensor *= (1.0/normalisation);
  trace = momentumTensor(0,0) + momentumTensor(1,1);
  
  // Print the tensor
  if (0)
    {
      std::cout << "\nNormalied 2-D Momentum Tensor"  << std::endl;
      Table tensor(" |  | ", "Text"); //LaTeX or Text    
      tensor.AddRowColumn(0, auxTools.ToString( momentumTensor(0,0) ) );
      tensor.AddRowColumn(0, auxTools.ToString( momentumTensor(0,1) ) );
      tensor.AddRowColumn(1, auxTools.ToString( momentumTensor(1,0) ) );
      tensor.AddRowColumn(1, auxTools.ToString( momentumTensor(1,1) ) );
      tensor.AddRowColumn(2, "");
      tensor.AddRowColumn(3, "Normalisation");
      tensor.AddRowColumn(3, auxTools.ToString(normalisation));
      tensor.AddRowColumn(4, "IsSymmetric");
      tensor.AddRowColumn(4, auxTools.ToString(momentumTensor.IsSymmetric()));
      tensor.AddRowColumn(5, "Determinant");
      tensor.AddRowColumn(5, auxTools.ToString(momentumTensor.Determinant()));
      tensor.AddRowColumn(6, "Trace");
      tensor.AddRowColumn(6, auxTools.ToString(trace));
      tensor.Print(false);
    }

  if ( abs(trace-1.0) > 1e-4 )
    {
      throw hplus::Exception("LogicError") << "Failure of requirement that the 2D Momentum-Tensor Trace is 1.0. Found that abs(trace-1) = " << abs(trace-1) << ", instead.";
    }

  return momentumTensor;
}


double KinematicsBkg::GetAlphaT(std::vector<math::XYZTLorentzVector> jets,
			     float &HT,
			     float &JT,
			     float &MHT,
			     float &Centrality){

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// AlphaT:
  /// Calculates the AlphaT variable, defined as an N-jets. This definition reproduces the kinematics of a 
  // di-jet system by constructing two pseudo-jets, which balance one another in Ht. 
  // The two pseudo-jets are formed from the combination of the N objects that minimizes the quantity
  /// DeltaHt = |Ht_pseudoJet1 - Ht_pseudoJet2| of the pseudo-jets.                                             
  //
  // Detailed Explanation: 
  // The method "alphaT()" of this class takes as input all jets in the event and uses them to form 
  // two Pseudo-Jets to describe the event. 
  // If there are NJets in a given event this means there are 2^{NJets-1} combinations to do this.
  // The methods does exactly that and for the combination which minimises the quantity
  // DeltaHt = Ht_PseudoJet1 - Ht_PseudoJet2,
  // it calculates the quantity alphaT.
  // The method "alphaT()" employs a double loop to recreate all the possilbe jet combinations 
  // out of NJets, by the use of an NJets-binary system. For example, if NJets=5, the loop
  // indices ("k" outside, "l" inside) run both from "k"=0 to "k"=2^{4}=16 . The upper limit of 
  // the outside loop is given by the expression:
  // 1<<(NJets-1)) = shift the number 1 by (NJets-1) positions to the left. 
  // So, for NJets=5  (i.e. 1  --> 1 0 0 0 0 ) 
  // This is now the way we will represent grouping into 2 Pseudo-Jets. The 0's represent one group and the 1's the other.
  // So, for example 1 0 0 0 0 means 1 jet forms Pseudo-Jet1 and 4 jets form Pseudo-Jet2. 
  // Also, for example, 1 0 0 1 0 means 2 jets form Pseudo-Jet1 and 3 jets form Pseudo-Jet2.
  // The inside loop performs a bitwise right shift of index "k" by "l" positions and then
  // compares the resulting bit to 1. So, for "k"=0, all the resulting comparisons in the 
  // inside loop will result to 0, except the one with "l"=4.
  // This gives the first combination: 0 0 0 0 0   ( i.e. 0 jets form Pseudo-Jet1 and 5 jets form Pseudo-Jet2 )
  // For "k"=1 (00000001 in 8bit representation), the first comparison is 1, since k is shifted by zero positions 
  // and then compared to 1. The rest comparisons yield zero, since by shifting the bit by any position and comparing to 1 gives zero. 
  // Thus, for "k"=1 we have after the second loop: 0 0 0 0 1
  // In the same manner, we get for "k"=2 (00000001 in 8bit representation) we have after the second loop: 0 0 0 1 0
  //  To summarise, for NJets=5 we have 16 combinations:
  // For "k"=0  ( 00000000 in 8bit representation) we have after the second loop: 0 0 0 0 0
  // For "k"=1  ( 00000001 in 8bit representation) we have after the second loop: 0 0 0 0 1
  // For "k"=2  ( 00000001 in 8bit representation) we have after the second loop: 0 0 0 1 0
  // For "k"=3  ( 00000011 in 8bit representation) we have after the second loop: 0 0 0 1 1
  // For "k"=4  ( 00000100 in 8bit representation) we have after the second loop: 0 0 1 0 0
  // For "k"=5  ( 00000101 in 8bit representation) we have after the second loop: 0 0 1 0 1
  // For "k"=6  ( 00000110 in 8bit representation) we have after the second loop: 0 0 1 1 0
  // For "k"=7  ( 00000111 in 8bit representation) we have after the second loop: 0 0 1 1 1
  // For "k"=8  ( 00001000 in 8bit representation) we have after the second loop: 0 1 0 0 0
  // For "k"=9  ( 00001001 in 8bit representation) we have after the second loop: 0 1 0 0 1
  // For "k"=10 ( 00010000 in 8bit representation) we have after the second loop: 0 1 0 0 0
  // For "k"=11 ( 00010001 in 8bit representation) we have after the second loop: 0 1 0 0 1
  // For "k"=12 ( 00010010 in 8bit representation) we have after the second loop: 0 1 0 1 0
  // For "k"=13 ( 00010011 in 8bit representation) we have after the second loop: 0 1 0 1 1
  // For "k"=14 ( 00010100 in 8bit representation) we have after the second loop: 0 1 1 0 0
  // For "k"=15 ( 00010101 in 8bit representation) we have after the second loop: 0 1 1 0 1
  // For "k"=16 ( 00010110 in 8bit representation) we have after the second loop: 0 1 1 1 0
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Declaration of variables 
  unsigned nJets = jets.size();

  // Sanity Check
  if ( jets.size() < 2 )
    {
      return -1.0;
    }

  std::vector<float> vE, vEt, vPx, vPy, vPz;
  std::vector<bool> vPseudo_jet1;
  const bool bList = true;

  // For-Loop: All jets
  for (auto& jet: jets){
    vE.push_back( jet.E() );
    vEt.push_back( jet.Et() );
    vPx.push_back( jet.Px() );
    vPy.push_back( jet.Py() );
    vPz.push_back( jet.Pz() );
  }

  // Calculate sums
  float fSum_e  = accumulate( vE.begin() , vE.end() , 0.0 );
  float fSum_et = accumulate( vEt.begin(), vEt.end(), 0.0 );
  float fSum_px = accumulate( vPx.begin(), vPx.end(), 0.0 );
  float fSum_py = accumulate( vPy.begin(), vPy.end(), 0.0 );

  // Minimum Delta Et for two pseudo-jets
  float fMin_delta_sum_et = -1.0;

  // Iterate through different combinations
  for ( unsigned k=0; k < unsigned(1<<(nJets-1)); k++ ) { 
    float fDelta_sum_et = 0.0;
    std::vector<bool> jet;

    // Iterate through jets
    for ( unsigned l=0; l < vEt.size(); l++ ) { 
      /// Bitwise shift of "k" by "l" positions to the right and compare to 1 (&1)
      /// i.e.: fDelta_sum_et += vEt[l] * ( 1 - 2*0 );  if comparison is un-successful
      ///  or   fDelta_sum_et += vEt[l] * ( 1 - 2*1 );  if comparison is successful
      // in this way you add up all Et from PseudoJetsGroupA (belonging to 0's group) and subtract that from PseudoJetsGroupB (1's group)
      fDelta_sum_et += vEt[l] * ( 1 - 2 * (int(k>>l)&1) ); 
      if ( bList ) { jet.push_back( (int(k>>l)&1) == 0 ); } 
    }

    // Find configuration with minimum value of DeltaHt 
    if ( ( fabs(fDelta_sum_et) < fMin_delta_sum_et || fMin_delta_sum_et < 0.0 ) ) {
      fMin_delta_sum_et = fabs(fDelta_sum_et);
      if ( bList && jet.size() == vEt.size() ){vPseudo_jet1.resize(jet.size());}
    }

  }
    
  // Sanity check
  if ( fMin_delta_sum_et < 0.0 )
    { 
      throw hplus::Exception("LogicError") << "Minimum Delta(Sum_Et) is less than zero! fMin_delta_sum_et = " << fMin_delta_sum_et;
    }
  
  // Calculate Event-Shape Variables
  float dHT = fMin_delta_sum_et;
  HT  = fSum_et;
  JT  = fSum_et - vEt.at(0); // Ht without considering the Ldg Jet of the Event
  MHT = sqrt(pow(fSum_px,2) + pow(fSum_py,2));
  Centrality = fSum_et/fSum_e;
  float AlphaT = ( 0.5 * ( HT - dHT ) / sqrt( pow(HT,2) - pow(MHT,2) ) );

  if (0)
    {

      Table vars("Variable | Value | Definition", "Text"); //LaTeX or Text
      vars.AddRowColumn(0, "HT");
      vars.AddRowColumn(0, auxTools.ToString(HT) );
      vars.AddRowColumn(0, "HT = Sum(Jet_Et)");
      //
      vars.AddRowColumn(1, "JT");
      vars.AddRowColumn(1, auxTools.ToString(JT) );
      vars.AddRowColumn(1, "JT = Ht - Jet1_Et");
      //
      vars.AddRowColumn(2, "dHT");
      vars.AddRowColumn(2, auxTools.ToString(dHT) );
      vars.AddRowColumn(2, "DeltaHT = min[Delta(Pseudojet1_Et, Pseudojet2_Et)]");
      //
      vars.AddRowColumn(3, "MHT");
      vars.AddRowColumn(3, auxTools.ToString(MHT) );
      vars.AddRowColumn(3, "MHT = sqrt( pow(Sum(px), 2) + pow(Sum(py), 2))");
      //
      vars.AddRowColumn(4, "AlphaT");
      vars.AddRowColumn(4, auxTools.ToString(AlphaT) );
      vars.AddRowColumn(4, "AlphaT = 0.5 x (HT - dHT) /sqr( pow(HT, 2) - pow(MHT, 2))");
      vars.Print();
    }

  return AlphaT;
}

vector<GenJet> KinematicsBkg::GetGenJets(vector<GenJet>& genJets, std::vector<float> ptCuts, std::vector<float> etaCuts, vector<genParticle> genParticlesToMatch)
{
  /*
    Jet-Flavour Definitions (https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools)

    Algorithmic definition: (NOTE: Algorithmic definition is used by default for all b-tagging purposes)
    - Try to find the parton that most likely determines the properties of the jet and assign that flavour as the true flavour
    - Here, the ¡°final state¡± partons (after showering, radiation) are analyzed (within ¦¤R < 0.3 of the reconstructed jet axis). 
    Partons selected for the algorithmic definition are those partons that don't have other partons as daughters, 
    without any explicit requirement on their status (for Pythia6 these are status=2 partons).
    - Jets from radiation are matched with full efficiency
    -If there is a b/c within the jet cone: label as b/c
    -Otherwise: assign flavour of the hardest parton
   */

  // Definitions
  std::vector<GenJet> jets;
  unsigned int jet_index   = -1;
  unsigned int ptCut_index  = 0;
  unsigned int etaCut_index = 0;

  // For-loop: All genParticles
  for (auto& p: genParticlesToMatch) 
    {
      
      // Comparison variables
      double dR   = 1e6;
      double dPt  = 1e6;
      double dEta = 1e6;
      double dPhi = 1e6;
      unsigned int genJet_index = -1;

      // For-loop: All Generated Jets
      for (auto jet: genJets) 
	{
	  // Jet index (for pT and eta cuts)
	  jet_index++;
	  genJet_index++;
	  
	  dPt  = jet.pt() - p.pt();
	  dEta = jet.eta() - p.eta();
	  dPhi = jet.phi() - p.phi();
       	  dR   = ROOT::Math::VectorUtil::DeltaR(jet.p4(), p.p4());
      
	  // Fail to match
	  if (dR > 0.3) continue;
	  
	  // Apply cuts after the matching
	  const float ptCut  = ptCuts.at(ptCut_index);
	  const float etaCut = etaCuts.at(etaCut_index);
	  if (jet.pt() < ptCut) continue;
	  if (std::abs(jet.eta()) > etaCut) continue;
	  
	  // Save this particle
	  jets.push_back(jet);

	  // Remove from vector
	  if (0) std::cout << "Erasing jet with pt  = " << jet.pt() << " GeV"  << std::endl;
	  genJets.erase(genJets.begin()+genJet_index);

	  if (0) std::cout << "dR = " << dR << ": dPt = " << dPt << ", dEta = " << dEta << ", dPhi = " << dPhi << std::endl;

	  // Increment cut index only. Cannot be bigger than the size of the cut list provided
	  if (ptCut_index  < ptCuts.size()-1  ) ptCut_index++;
	  if (etaCut_index < etaCuts.size()-1  ) etaCut_index++;
	  break;
	}
    }

  // Debugging
  if (0) 
    {
      std::cout << "=== AFTER" << std::endl;
      for (auto jet: genJets) std::cout << "genJet.pt() = " << jet.pt() << std::endl;
      for (auto jet: jets) std::cout << "jet.pt() = " << jet.pt() << std::endl;
      std::cout << "jets.size() = " << jets.size() << ", genJets.size() = " << genJets.size() << std::endl;
    }
  return jets;
}

vector<GenJet> KinematicsBkg::GetGenJets(const GenJetCollection& genJets, std::vector<float> ptCuts, std::vector<float> etaCuts)
{
  std::vector<GenJet> jets;

  // Definitions
  unsigned int jet_index   = -1;
  unsigned int ptCut_index  = 0;
  unsigned int etaCut_index = 0;

  // For-loop: All Generated Jets
  for (auto jet: genJets) 
    {

      // Jet index (for pT and eta cuts)
      jet_index++;
      
      // Apply cuts
      const float ptCut  = ptCuts.at(ptCut_index);
      const float etaCut = etaCuts.at(etaCut_index);
      if (jet.pt() < ptCut) continue;
      if (std::abs(jet.eta()) > etaCut) continue;

      // Save this particle
      jets.push_back(jet);

      // Increment cut index only. Cannot be bigger than the size of the cut list provided
      if (ptCut_index  < ptCuts.size()-1  ) ptCut_index++;
      if (etaCut_index < etaCuts.size()-1  ) etaCut_index++;
    }
  return jets;
}


vector<genParticle> KinematicsBkg::GetAllPreviousCopies(const vector<genParticle> genParticles, genParticle genP)
{

  std::vector<short int> genMoms_index = genP.mothers();
  std::vector<genParticle> copies; 
  copies.push_back(genP); // Include self in copies vector!

  while (genMoms_index.size() > 0)
    {
      
      for (unsigned int i = 0; i < genMoms_index.size(); i++ )
	{
	  genParticle p = genParticles.at( genMoms_index.at(i) );

	  // If previous mom does not have same pdgId skip
	  if ( p.pdgId() != genP.pdgId() )  
	    {
	      genMoms_index = p.mothers();
	      continue;
	    }
	  
	  // Save mom as copy 
	  copies.push_back(p);
	  
	  // Replace moms
	  genMoms_index = p.mothers();
	}
    }

  // Sort vector with ascending index order
  std::reverse(copies.begin(),copies.end());


  if (0)
    {
      for (auto& p: copies) 
	{
	  if (abs(p.pdgId()) != 5) continue;
	  cout << "p.index() = " << p.index() << ", p.pdgId() = " << p.pdgId()  << endl;
	}
      cout << "copies.size() = " << copies.size() << endl;
    }

  return copies;
}


vector<genParticle> KinematicsBkg::GetGenParticles(const vector<genParticle> genParticles, double ptCut, double etaCut, const int pdgId, const bool isLastCopy, const bool hasNoDaughters)
  {
  
  std::vector<genParticle> particles;

  // For-loop: All genParticles
  for (auto& p: genParticles) 
    {

      // std::cout << "p.pdgId() = " << p.pdgId() << std::endl;

      // Find last copy of a given particle
      // if (isLastCopy) if (!p.isLastCopy()) continue; // crashes if branch missing from ntuple

      if (isLastCopy)
	{
	  if (abs(p.pdgId()) == 5){ if (p.status() != 23) continue;} // Consider only status=23 (outgoing) particles
	  else{ if (p.status() != 1 and p.status() != 2) continue;}
	}
      
      // Commonly enables for parton-based jet flavour definition
      if (hasNoDaughters) if (p.daughters().size() > 0) continue;

      // Consider only particles
      if (std::abs(p.pdgId()) != pdgId) continue;

      if (0) if (std::abs(p.pdgId()) == 13) std::cout << "p.pdgId() = " << p.pdgId() << ", p.pt() = " << p.pt() << ", p.status() = " << p.status() << ", p.mothers().at(0) = " << p.mothers().at(0) << std::endl;
      if (0) if (std::abs(p.pdgId()) == 15) std::cout << "p.pdgId() = " << p.pdgId() << ", p.pt() = " << p.pt() << ", p.status() = " << p.status() << ", p.mothers().at(0) = " << p.mothers().at(0) << std::endl;

      // Apply cuts
      if ( p.pt() < ptCut ) continue;
      if (std::abs(p.eta()) > etaCut) continue;
      
      // Save this particle
      particles.push_back(p);
    }

  // std::cout << "Returning " << particles.size() << " particles." << std::endl;
  return particles;
    }

double KinematicsBkg::GetTransverseMass(const math::XYZVector tau1, const math::XYZVector tau2, const math::XYZVector muon, const math::XYZVector met) {
  //
  // mTSquared = (ET_Taujet1 + ET_Taujet2 + ET_Muon1 + MET)**2 - (p_Taujet1 + p_Taujet2 + p_Muon1 + p_MET)**2
  //

  // Initialise Variables
  double mT        = -999.9;
  double mTSquared =    0.0;
  
  mTSquared = ( pow(sqrt(tau1.Mag2()) + sqrt(tau2.Mag2()) + sqrt(muon.Mag2()) + sqrt(met.Mag2()), 2) - (tau1 + tau2 + muon + met).Dot(tau1 + tau2 + muon + met) );
  mT = std::sqrt(mTSquared);
  
  return mT;
}

double KinematicsBkg::GetMt(const math::XYVector tau1, const math::XYVector tau2, const math::XYVector muon, const math::XYVector& met) {
  // Use scalar sums to get the transverse mass
  double metEt  = met.R();
  double tau1Et = tau1.r();
  double tau2Et = tau2.r();
  double muonEt = muon.r();
  double mT     =-999.9;
  double mTSq   =   0.0;
  double EtSq   = (metEt + tau1Et + tau2Et + muonEt) * (metEt + tau1Et + tau2Et + muonEt);
  double EtXSq  = (tau1.x() + tau2.x() + muon.x() + met.x()) * (tau1.x() + tau2.x() + muon.x() + met.x());
  double EtYSq  = (tau1.y() + tau2.y() + muon.y() + met.y()) * (tau1.y() + tau2.y() + muon.y() + met.y());
  mTSq = EtSq - (EtXSq + EtYSq);

  if (mTSq >= 0) mT = std::sqrt(mTSq);
  return mT;
}

double KinematicsBkg::GetEffectiveMass(const math::XYZTLorentzVector tau1, const math::XYZTLorentzVector tau2, const math::XYVector& met) {
  // https://root.cern.ch/root/html/ROOT__Math__LorentzVector_-p1PxPyPzE4D_double___.html#ROOT__Math__LorentzVector_-p1PxPyPzE4D_double___:P2
  // const math::XYZVector metv(met.x(), met.y(), 0.0);

  // Calculation
  // double EtSq = pow(sqrt(tau1.P()) + sqrt(tau2.P()) + sqrt(met.Mag2()), 2);
  // double PtSq = (tau1.Vect() + tau2.Vect() + metv).Dot(tau1.Vect() + tau2.Vect() + metv);
  double EtSq   = pow(tau1.Et() + tau2.Et() + met.R(), 2);
  double PtSqX  = ( tau1.px() + tau2.px() + met.x() ) * ( tau1.px() + tau2.px() + met.x() );
  double PtSqY  = ( tau1.py() + tau2.py() + met.y() ) * ( tau1.py() + tau2.py() + met.y() );
  double PtSqZ  = ( tau1.pz() + tau2.pz() ) * ( tau1.pz() + tau2.pz() );

  double PtSq   = (PtSqX + PtSqY + PtSqZ);
  double mEff   = -999.9;
  double mEffSq = EtSq - PtSq;
  if (mEffSq >= 0.0) mEff = std::sqrt(mEffSq);
  return mEff;
}

double KinematicsBkg::GetCollinearMass(const math::XYZTLorentzVector tau1, const math::XYZTLorentzVector tau2, const math::XYVector& met) {
  //
  // Must check that the two taus are not back-to-back in the lab frame:
  // | dPhi(tau1vis, tau2vis)| < 2.9 
  // for example, otherwise the approximation breaks down.
  // Mass resolution limited by missing transverse energy resolution.
  //

  // Calculation
  double mVis  = ( tau1 + tau2 ).M();
  double metEt = met.R();
  double x1    = tau1.Et()/(tau1.Et() + metEt );
  double x2    = tau2.Et()/(tau2.Et() + metEt );
  double mColl = mVis/sqrt(x1 * x2);
  return mColl;

}
