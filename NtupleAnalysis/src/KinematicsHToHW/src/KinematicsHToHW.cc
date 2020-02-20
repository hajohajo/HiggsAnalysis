// -*- c++ -*-
#include "Framework/interface/BaseSelector.h"
#include "Framework/interface/makeTH.h"

// User
#include "Auxiliary/interface/Table.h"
#include "Auxiliary/interface/Tools.h"
#include "Tools/interface/MCTools.h"
#include "Tools/interface/DirectionalCut.h"
#include "EventSelection/interface/CommonPlots.h"
#include "EventSelection/interface/EventSelections.h"

// ROOT
#include "TDirectory.h"
#include "Math/VectorUtil.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

struct PtComparator
{
  bool operator() (const genParticle p1, const genParticle p2) const { return ( p1.pt() > p2.pt() ); }
  bool operator() (const math::XYZTLorentzVector p1, const math::XYZTLorentzVector p2) const { return ( p1.pt() > p2.pt() ); }
};


class KinematicsHToHW: public BaseSelector {
public:
  explicit KinematicsHToHW(const ParameterSet& config, const TH1* skimCounters);
  virtual ~KinematicsHToHW() {}

  /// Books histograms
  virtual void book(TDirectory *dir) override;
  /// Sets up branches for reading the TTree
  virtual void setupBranches(BranchManager& branchManager) override;
  /// Called for each event
  virtual void process(Long64_t entry) override;
  virtual vector<genParticle> GetAllPreviousCopies(const vector<genParticle> genParticles, genParticle genP);
  virtual vector<genParticle> GetGenParticles(const vector<genParticle> genParticles, double ptCut, double etaCut, const int pdgId, const bool isLastCopy=true, const bool hasNoDaughters=false);
  virtual vector<GenJet> GetGenJets(const GenJetCollection& genJets, std::vector<float> ptCut, std::vector<float> etaCut);
  virtual vector<GenJet> GetGenJets(const vector<GenJet> genJets, std::vector<float> ptCut, std::vector<float> etaCut, vector<genParticle> genParticlesToMatch);
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
  const double cfg_TauPtCut;
  const double cfg_TauEtaCut;
  const ParameterSet PSet_JetSelection;
  const std::vector<float> cfg_JetPtCuts;
  const std::vector<float> cfg_JetEtaCuts;
  const DirectionalCut<int> cfg_JetNumberCut;
  const DirectionalCut<float> cfg_HtCut;
  const ParameterSet PSet_BJetSelection;
  const std::vector<float> cfg_BJetPtCuts;
  const std::vector<float> cfg_BJetEtaCuts;
  const DirectionalCut<int> cfg_BJetNumberCut;
  const ParameterSet PSet_METSelection;
  const DirectionalCut<float> cfg_METCut;
  // METSelection PSet_METSelection;
  // TopologySelection PSet_TopologySelection;
  // const DirectionalCut<double> cfg_SphericityCut;
  // const DirectionalCut<double> cfg_AplanarityCut;
  // const DirectionalCut<double> cfg_PlanarityCut;
  // const DirectionalCut<double> cfg_CircularityCut;
  // const DirectionalCut<double> cfg_Y23Cut;
  // const DirectionalCut<double> cfg_CparameterCut;
  // const DirectionalCut<double> cfg_DparameterCut;
  // const DirectionalCut<double> cfg_FoxWolframMomentCut;
  // const DirectionalCut<double> cfg_AlphaTCut;
  // const DirectionalCut<double> cfg_CentralityCut;
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
  Count cTauVeto;
  Count cJetSelection;
  Count cBJetSelection;
  Count cMETSelection;
  // Count cTopologySelection;
  // Count cTopSelection;
  Count cSelected;

  // BR Counters
  Count cInclusive;
  Count cbHt_Hpm;
  Count cbHt_HToHW_HBoson; 
  Count cbHt_HToHW_HBoson_TauTau; 
  Count cbHt_HToHW_HBoson_TauLeptonic; 
  Count cbHt_HToHW_HBoson_TauHadronic;
  Count cbHt_HToHW_HBoson_1Lep1Had;
  Count cbHt_HToHW_WBoson;
  Count cbHt_HToHW_WBoson_Quark;
  Count cbHt_HToHW_WBoson_Antiquark;
  Count cbHt_HToHW_WBoson_Leptons;
  Count cbHt_HToHW_WBoson_Neutrinos;
  Count cbHt_tWb_BQuark;
  Count cbHt_tWb_WBoson;
  Count cbHt_tWb_WBoson_Quark;
  Count cbHt_tWb_WBoson_Antiquark;
  Count cbHt_tWb_WBoson_Leptons;
  Count cbHt_tWb_WBoson_Neutrinos;

  // Event Variables
  WrappedTH1 *h_genMET_Et;
  WrappedTH1 *h_genMET_Phi;
  WrappedTH1 *h_genHT_GenParticles;
  WrappedTH1 *h_genHT_GenJets;  
  
  // GenParticles
  WrappedTH1 *h_bHt_Hpm_Pt;
  WrappedTH1 *h_bHt_Hpm_Eta;
  WrappedTH1 *h_bHt_Hpm_Rap;
  WrappedTH1 *h_bHt_TQuark_Pt;
  WrappedTH1 *h_bHt_TQuark_Eta;
  WrappedTH1 *h_bHt_TQuark_Rap;
  WrappedTH1 *h_bHt_tWb_WBoson_Pt;
  WrappedTH1 *h_bHt_tWb_WBoson_Eta;
  WrappedTH1 *h_bHt_tWb_WBoson_Rap;
  WrappedTH1 *h_bHt_tWb_BQuark_Pt;
  WrappedTH1 *h_bHt_tWb_BQuark_Eta;
  WrappedTH1 *h_bHt_tWb_BQuark_Rap;
  WrappedTH1 *h_bHt_tWb_WBoson_Quark_Pt;
  WrappedTH1 *h_bHt_tWb_WBoson_Quark_Eta;
  WrappedTH1 *h_bHt_tWb_WBoson_Quark_Rap;
  WrappedTH1 *h_bHt_tWb_WBoson_Antiquark_Pt;
  WrappedTH1 *h_bHt_tWb_WBoson_Antiquark_Eta;  
  WrappedTH1 *h_bHt_tWb_WBoson_Antiquark_Rap;  
  WrappedTH1 *h_bHt_tWb_WBoson_Leptons_Pt;
  WrappedTH1 *h_bHt_tWb_WBoson_Leptons_Eta;  
  WrappedTH1 *h_bHt_tWb_WBoson_Leptons_Rap;  
  WrappedTH1 *h_bHt_tWb_WBoson_Neutrinos_Pt;
  WrappedTH1 *h_bHt_tWb_WBoson_Neutrinos_Eta;  
  WrappedTH1 *h_bHt_tWb_WBoson_Neutrinos_Rap;  
  WrappedTH1 *h_bHt_HToHW_HBoson_Pt;
  WrappedTH1 *h_bHt_HToHW_HBoson_Eta;
  WrappedTH1 *h_bHt_HToHW_HBoson_Rap;
  WrappedTH1 *h_bHt_HToHW_WBoson_Pt;
  WrappedTH1 *h_bHt_HToHW_WBoson_Eta;
  WrappedTH1 *h_bHt_HToHW_WBoson_Rap;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau_Pt;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau_Eta;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau_Rap;
  WrappedTH1 *h_bHt_HToHW_HBoson_Antitau_Pt;
  WrappedTH1 *h_bHt_HToHW_HBoson_Antitau_Eta;
  WrappedTH1 *h_bHt_HToHW_HBoson_Antitau_Rap;
  WrappedTH1 *h_bHt_HToHW_WBoson_Quark_Pt;
  WrappedTH1 *h_bHt_HToHW_WBoson_Quark_Eta;
  WrappedTH1 *h_bHt_HToHW_WBoson_Quark_Rap;
  WrappedTH1 *h_bHt_HToHW_WBoson_Antiquark_Pt;
  WrappedTH1 *h_bHt_HToHW_WBoson_Antiquark_Eta;  
  WrappedTH1 *h_bHt_HToHW_WBoson_Antiquark_Rap;  
  WrappedTH1 *h_bHt_HToHW_WBoson_Leptons_Pt;
  WrappedTH1 *h_bHt_HToHW_WBoson_Leptons_Eta;
  WrappedTH1 *h_bHt_HToHW_WBoson_Leptons_Rap;
  WrappedTH1 *h_bHt_HToHW_WBoson_Neutrinos_Pt;
  WrappedTH1 *h_bHt_HToHW_WBoson_Neutrinos_Eta;
  WrappedTH1 *h_bHt_HToHW_WBoson_Neutrinos_Rap;

  WrappedTH1 *h_bHt_HToHW_HBoson_Taus_dR;
  WrappedTH1 *h_bHt_HToHW_HBoson_Taus_dEta;
  WrappedTH1 *h_bHt_HToHW_HBoson_Taus_dPhi;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau1_MET_dPhi;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau2_MET_dPhi;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau1_bHt_tWb_BQuark_dR;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau1_bHt_tWb_BQuark_dEta;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau1_bHt_tWb_BQuark_dPhi;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau2_bHt_tWb_BQuark_dR;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau2_bHt_tWb_BQuark_dEta;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau2_bHt_tWb_BQuark_dPhi;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Quark_dR;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Quark_dEta;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Quark_dPhi;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Quark_dR;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Quark_dEta;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Quark_dPhi;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Antiquark_dR;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Antiquark_dEta;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Antiquark_dPhi;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Antiquark_dR;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Antiquark_dEta;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Antiquark_dPhi;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau1_bHt_HToHW_WBoson_Leptons_dR;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau1_bHt_HToHW_WBoson_Leptons_dEta;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau1_bHt_HToHW_WBoson_Leptons_dPhi;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau2_bHt_HToHW_WBoson_Leptons_dR;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau2_bHt_HToHW_WBoson_Leptons_dEta;
  WrappedTH1 *h_bHt_HToHW_HBoson_Tau2_bHt_HToHW_WBoson_Leptons_dPhi;

  WrappedTH1 *h_bHt_Hpm_bHt_tWb_WBoson_dR;
  WrappedTH1 *h_bHt_Hpm_bHt_tWb_WBoson_dEta;
  WrappedTH1 *h_bHt_Hpm_bHt_tWb_WBoson_dPhi;
  WrappedTH1 *h_bHt_Hpm_bHt_tWb_WBoson_dRap;
  WrappedTH1 *h_bHt_Hpm_bHt_TQuark_dR;
  WrappedTH1 *h_bHt_Hpm_bHt_TQuark_dEta;
  WrappedTH1 *h_bHt_Hpm_bHt_TQuark_dPhi;
  WrappedTH1 *h_bHt_Hpm_bHt_TQuark_dRap;
  WrappedTH1 *h_bHt_HToHW_WBoson_bHt_tWb_WBoson_dR;
  WrappedTH1 *h_bHt_HToHW_WBoson_bHt_tWb_WBoson_dEta;
  WrappedTH1 *h_bHt_HToHW_WBoson_bHt_tWb_WBoson_dPhi;
  WrappedTH1 *h_bHt_HToHW_WBoson_bHt_tWb_WBoson_dRap;
  WrappedTH1 *h_bHt_TQuark_bHt_tWb_BQuark_dR;
  WrappedTH1 *h_bHt_TQuark_bHt_tWb_BQuark_dEta;
  WrappedTH1 *h_bHt_TQuark_bHt_tWb_BQuark_dPhi;
  WrappedTH1 *h_bHt_TQuark_bHt_tWb_BQuark_dRap;
  WrappedTH1 *h_bHt_tWb_WBoson_bHt_tWb_BQuark_dPhi;
  WrappedTH1 *h_bHt_tWb_WBoson_bHt_tWb_BQuark_dR;
  WrappedTH1 *h_bHt_tWb_WBoson_bHt_tWb_BQuark_dEta;
  WrappedTH1 *h_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark_dR;
  WrappedTH1 *h_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark_dEta;
  WrappedTH1 *h_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark_dPhi;
  WrappedTH1 *h_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark_dRap;
  WrappedTH1 *h_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark_dR;
  WrappedTH1 *h_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark_dEta;
  WrappedTH1 *h_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark_dPhi;
  WrappedTH1 *h_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark_dRap;
  WrappedTH1 *h_bHt_tWb_bqq_Pt;
  WrappedTH1 *h_bHt_tWb_bqq_Rap;
  WrappedTH1 *h_bHt_tWb_bqq_Mass;
  WrappedTH1 *h_bHt_tWb_bqq_dRMax_dRap;
  WrappedTH1 *h_bHt_tWb_bqq_dRMax_dPhi;
  WrappedTH1 *h_bHt_tWb_bqq_dRMax_dR;


  // GenJets
  WrappedTH1 *h_GenJets_N;  
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

  // TH2
  WrappedTH2 *h_bHt_tWb_bqq_dRMax_dRap_Vs_dPhi;
  WrappedTH2 *h_bHt_HToHW_HBoson_Taus_Pt_Vs_Pt;
  WrappedTH2 *h_bHt_HToHW_HBoson_Taus_Eta_Vs_Eta;
  WrappedTH2 *h_bHt_HToHW_HBoson_Taus_Phi_Vs_Phi;
  WrappedTH2 *h_bHt_HToHW_HBoson_Taus_dEta_Vs_dPhi;
  WrappedTH2 *h_bHt_HToHW_HBoson_Taus_Pt1_Vs_dR;
  WrappedTH2 *h_bHt_HToHW_HBoson_Taus_Pt2_Vs_dR;
  WrappedTH2 *h_MaxDiJetMass_dEta_Vs_dPhi;
  WrappedTH2 *h_MaxDiJetMass_dRap_Vs_dPhi;  
  WrappedTH2 *h_Bquark_Antiquark_dEta_Vs_dPhi;
  WrappedTH2 *h_Bquark_Quark_dEta_Vs_dPhi;
  WrappedTH2 *h_Quark_Antiquark_dEta_Vs_dPhi;
  WrappedTH2 *h_Jet1Jet2_dEta_Vs_Jet3Jet4_dEta;
  WrappedTH2 *h_Jet1Jet2_dPhi_Vs_Jet3Jet4_dPhi;  
  WrappedTH2 *h_Jet1Jet2_dEta_Vs_Jet1Jet2_Mass;
  WrappedTH2 *h_Jet3Jet4_dEta_Vs_Jet3Jet4_Mass;

};

#include "Framework/interface/SelectorFactory.h"
REGISTER_SELECTOR(KinematicsHToHW);

KinematicsHToHW::KinematicsHToHW(const ParameterSet& config, const TH1* skimCounters)
  : BaseSelector(config, skimCounters),
    cfg_Verbose(config.getParameter<bool>("verbose")),
    PSet_ElectronSelection(config.getParameter<ParameterSet>("ElectronSelection")),
    cfg_ElectronPtCut(config.getParameter<float>("ElectronSelection.electronPtCut")),  
    cfg_ElectronEtaCut(config.getParameter<float>("ElectronSelection.electronEtaCut")),
    PSet_MuonSelection(config.getParameter<ParameterSet>("MuonSelection")),
    cfg_MuonPtCut(config.getParameter<float>("MuonSelection.muonPtCut")),
    cfg_MuonEtaCut(config.getParameter<float>("MuonSelection.muonEtaCut")),
    PSet_TauSelection(config.getParameter<ParameterSet>("TauSelection")),
    cfg_TauPtCut(config.getParameter<float>("TauSelection.tauPtCut")),
    cfg_TauEtaCut(config.getParameter<float>("TauSelection.tauEtaCut")),
    PSet_JetSelection(config.getParameter<ParameterSet>("JetSelection")),
    cfg_JetPtCuts(config.getParameter<std::vector<float>>("JetSelection.jetPtCuts")),
    cfg_JetEtaCuts(config.getParameter<std::vector<float>>("JetSelection.jetEtaCuts")),
    cfg_JetNumberCut(config, "JetSelection.numberOfJetsCut"),
    cfg_HtCut(config, "JetSelection.HTCut"),
    PSet_BJetSelection(config.getParameter<ParameterSet>("BJetSelection")),
    cfg_BJetPtCuts(config.getParameter<std::vector<float>>("BJetSelection.jetPtCuts")),
    cfg_BJetEtaCuts(config.getParameter<std::vector<float>>("BJetSelection.jetEtaCuts")),
    cfg_BJetNumberCut(config, "BJetSelection.numberOfBJetsCut"),
    PSet_METSelection(config.getParameter<ParameterSet>("METSelection")),
    cfg_METCut(config, "METSelection.METCut"),
    // PSet_TopologySelection(config.getParameter<ParameterSet>("TopologySelection")),
    // cfg_SphericityCut(config, "TopologySelection.SphericityCut"),
    // cfg_AplanarityCut(config, "TopologySelection.AplanarityCut"),
    // cfg_PlanarityCut(config, "TopologySelection.PlanarityCut"),
    // cfg_CircularityCut(config, "TopologySelection.CircularityCut"),
    // cfg_Y23Cut(config, "TopologySelection.Y23Cut"),
    // cfg_CparameterCut(config, "TopologySelection.CparameterCut"),
    // cfg_DparameterCut(config, "TopologySelection.DparameterCut"),
    // cfg_FoxWolframMomentCut(config, "TopologySelection.FoxWolframMomentCut"),
    // cfg_AlphaTCut(config, "TopologySelection.AlphaTCut"),
    // cfg_CentralityCut(config, "TopologySelection.CentralityCut"),
    // PSet_TopSelection(config.getParameter<ParameterSet>("TopSelection")),
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
    cTauVeto(fEventCounter.addCounter("#tau-veto")),
    cJetSelection(fEventCounter.addCounter("Jets + H_{T}")),
    cBJetSelection(fEventCounter.addCounter("b-jets")),
    cMETSelection(fEventCounter.addCounter("MET")),
    // cTopologySelection(fEventCounter.addCounter("Topology")),
    // cTopSelection(fEventCounter.addCounter("Top")),
    cSelected(fEventCounter.addCounter("All Selections")),
    cInclusive(fEventCounter.addSubCounter("Branching", "All events")),
    cbHt_Hpm(fEventCounter.addSubCounter("Branching", "H+->HW")),
    cbHt_HToHW_HBoson(fEventCounter.addSubCounter("Branching", "H+->HW, H")),
    cbHt_HToHW_HBoson_TauTau(fEventCounter.addSubCounter("Branching", "H+->HW, H->tautau")),
    cbHt_HToHW_HBoson_TauLeptonic(fEventCounter.addSubCounter("Branching", "H+->HW, H->TauLeptonic")),
    cbHt_HToHW_HBoson_TauHadronic(fEventCounter.addSubCounter("Branching", "H+->HW, H->TauHadronic")),
    cbHt_HToHW_HBoson_1Lep1Had(fEventCounter.addSubCounter("Branching", "H+->HW, H->1Lep1Had")),
    cbHt_HToHW_WBoson(fEventCounter.addSubCounter("Branching", "H+->HW, W")),
    cbHt_HToHW_WBoson_Quark(fEventCounter.addSubCounter("Branching", "H+->HW, W->qq, q")),
    cbHt_HToHW_WBoson_Antiquark(fEventCounter.addSubCounter("Branching", "H+->HW, W->qq, qbar")),
    cbHt_HToHW_WBoson_Leptons(fEventCounter.addSubCounter("Branching", "H+->HW, W->lnu, l")),
    cbHt_HToHW_WBoson_Neutrinos(fEventCounter.addSubCounter("Branching", "H+->HW, W->lnu, nu")),
    cbHt_tWb_BQuark(fEventCounter.addSubCounter("Branching", "t->Wb, b")),
    cbHt_tWb_WBoson(fEventCounter.addSubCounter("Branching", "t->Wb, W")),
    cbHt_tWb_WBoson_Quark(fEventCounter.addSubCounter("Branching", "t->Wb, W->qq, q")),
    cbHt_tWb_WBoson_Antiquark(fEventCounter.addSubCounter("Branching", "t->Wb, W->qq, qbar")),
    cbHt_tWb_WBoson_Leptons(fEventCounter.addSubCounter("Branching", "t->Wb, W->lnu, l")),
    cbHt_tWb_WBoson_Neutrinos(fEventCounter.addSubCounter("Branching", "t->Wb, W->lnu, nu"))
{ }

void KinematicsHToHW::book(TDirectory *dir) {

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
  const double mindPhi = cfg_DeltaPhiBinSetting.min();
  const double maxdPhi = cfg_DeltaPhiBinSetting.max();

  const int nBinsdR  = 2*cfg_DeltaRBinSetting.bins();
  const double mindR = cfg_DeltaRBinSetting.min();
  const double maxdR = cfg_DeltaRBinSetting.max();

  TDirectory* th1 = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "TH1");
  TDirectory* th2 = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "TH2");
    
  // Event Variables
  h_genMET_Et         =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "genMET_Et"         , ";E_{T}^{miss} (GeV)", 500,  0.0,   +500.0);
  h_genMET_Phi        =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kInformative, th1, "genMET_Phi"  , ";#phi (rads)"       , nBinsPhi, minPhi, maxPhi);
  h_genHT_GenParticles=  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "genHT_GenParticles", ";GenP H_{T} (GeV)"  , nBinsM, minM, maxM   );
  h_genHT_GenJets     =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "genHT_GenJets"     , ";GenJ H_{T} (GeV)"  , nBinsM, minM, maxM   );

  // GenParticles  
  h_bHt_Hpm_Pt                    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_Hpm_Pt"                   , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_TQuark_Pt                 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_TQuark_Pt"                , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_tWb_WBoson_Pt             =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_Pt"            , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_tWb_BQuark_Pt             =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_BQuark_Pt"            , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_tWb_WBoson_Quark_Pt       =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_Quark_Pt"      , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_tWb_WBoson_Antiquark_Pt   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_Antiquark_Pt"  , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_tWb_WBoson_Leptons_Pt     =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_Leptons_Pt"    , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_tWb_WBoson_Neutrinos_Pt   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_Neutrinos_Pt"  , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_HToHW_HBoson_Pt           =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Pt"          , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_HToHW_WBoson_Pt           =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_Pt"          , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_HToHW_HBoson_Tau_Pt       =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau_Pt"      , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_HToHW_HBoson_Antitau_Pt   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Antitau_Pt"  , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_HToHW_WBoson_Leptons_Pt   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_Leptons_Pt"  , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_HToHW_WBoson_Neutrinos_Pt =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_Neutrinos_Pt", ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_HToHW_WBoson_Quark_Pt     =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_Quark_Pt"    , ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_bHt_HToHW_WBoson_Antiquark_Pt =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_Antiquark_Pt", ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);

  h_bHt_Hpm_Eta                    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_Hpm_Eta"                   , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_TQuark_Eta                 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_TQuark_Eta"                , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_tWb_WBoson_Eta             =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_Eta"            , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_tWb_BQuark_Eta             =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_BQuark_Eta"            , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_tWb_WBoson_Quark_Eta       =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_Quark_Eta"      , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_tWb_WBoson_Antiquark_Eta   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_Antiquark_Eta"  , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_tWb_WBoson_Leptons_Eta     =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_Leptons_Eta"    , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_tWb_WBoson_Neutrinos_Eta   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_Neutrinos_Eta"  , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_HToHW_HBoson_Eta           =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Eta"          , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_HToHW_WBoson_Eta           =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_Eta"          , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_HToHW_HBoson_Tau_Eta       =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau_Eta"      , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_HToHW_HBoson_Antitau_Eta   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Antitau_Eta"  , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_HToHW_WBoson_Leptons_Eta   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_Leptons_Eta"  , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_HToHW_WBoson_Neutrinos_Eta =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_Neutrinos_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_HToHW_WBoson_Quark_Eta     =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_Quark_Eta"    , ";#eta", nBinsEta, minEta, maxEta);
  h_bHt_HToHW_WBoson_Antiquark_Eta =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_Antiquark_Eta", ";#eta", nBinsEta, minEta, maxEta);

  h_bHt_Hpm_Rap                    =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_Hpm_Rap"                   , ";#omega", nBinsEta, minRap, maxRap);
  h_bHt_TQuark_Rap                 =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_TQuark_Rap"                , ";#omega", nBinsRap, minRap, maxRap);
  h_bHt_tWb_WBoson_Rap             =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_Rap"            , ";#omega", nBinsEta, minRap, maxRap);
  h_bHt_tWb_BQuark_Rap             =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_BQuark_Rap"            , ";#omega", nBinsEta, minRap, maxRap);
  h_bHt_tWb_WBoson_Quark_Rap       =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_Quark_Rap"      , ";#omega", nBinsEta, minRap, maxRap);
  h_bHt_tWb_WBoson_Antiquark_Rap   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_Antiquark_Rap"  , ";#omega", nBinsEta, minRap, maxRap);
  h_bHt_tWb_WBoson_Leptons_Rap     =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_Leptons_Rap"    , ";#omega", nBinsEta, minRap, maxRap);
  h_bHt_tWb_WBoson_Neutrinos_Rap   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_Neutrinos_Rap"  , ";#omega", nBinsEta, minRap, maxRap);
  h_bHt_HToHW_HBoson_Rap           =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Rap"          , ";#omega", nBinsEta, minRap, maxRap);
  h_bHt_HToHW_WBoson_Rap           =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_Rap"          , ";#omega", nBinsEta, minRap, maxRap);
  h_bHt_HToHW_HBoson_Tau_Rap       =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau_Rap"      , ";#omega", nBinsEta, minRap, maxRap);
  h_bHt_HToHW_HBoson_Antitau_Rap   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Antitau_Rap"  , ";#omega", nBinsEta, minRap, maxRap);
  h_bHt_HToHW_WBoson_Leptons_Rap   =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_Leptons_Rap"  , ";#omega", nBinsEta, minRap, maxRap);
  h_bHt_HToHW_WBoson_Neutrinos_Rap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_Neutrinos_Rap", ";#omega", nBinsEta, minRap, maxRap);
  h_bHt_HToHW_WBoson_Quark_Rap     =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_Quark_Rap"    , ";#omega", nBinsEta, minRap, maxRap);
  h_bHt_HToHW_WBoson_Antiquark_Rap =  fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_Antiquark_Rap", ";#omega", nBinsEta, minRap, maxRap);

  h_bHt_HToHW_HBoson_Taus_dR                   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Taus_dR"                  , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_bHt_Hpm_bHt_tWb_WBoson_dR                  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_Hpm_bHt_tWb_WBoson_dR"                 , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_bHt_Hpm_bHt_TQuark_dR                      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_Hpm_bHt_TQuark_dR"                     , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_bHt_HToHW_WBoson_bHt_tWb_WBoson_dR         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_bHt_tWb_WBoson_dR"        , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_bHt_tWb_WBoson_bHt_tWb_BQuark_dR           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_bHt_tWb_BQuark_dR"          , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_bHt_TQuark_bHt_tWb_BQuark_dR               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_TQuark_bHt_tWb_BQuark_dR"              , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark_dR     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_BQuark_bHt_tWb_WBoson_Quark_dR"    , ";#DeltaR", nBinsdR, mindR, maxdR);
  h_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark_dR = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark_dR", ";#DeltaR", nBinsdR, mindR, maxdR);
  
  h_bHt_HToHW_HBoson_Taus_dEta                   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Taus_dEta"                  , ";#Delta#eta", nBinsdEta, mindEta, maxdEta);
  h_bHt_Hpm_bHt_tWb_WBoson_dEta                  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_Hpm_bHt_tWb_WBoson_dEta"                 , ";#Delta#eta", nBinsdEta, mindEta, maxdEta);
  h_bHt_Hpm_bHt_TQuark_dEta                      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_Hpm_bHt_TQuark_dEta"                     , ";#Delta#eta", nBinsdEta, mindEta, maxdEta);
  h_bHt_HToHW_WBoson_bHt_tWb_WBoson_dEta         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_bHt_tWb_WBoson_dEta"        , ";#Delta#eta", nBinsdEta, mindEta, maxdEta);
  h_bHt_TQuark_bHt_tWb_BQuark_dEta               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_TQuark_bHt_tWb_BQuark_dEta"              , ";#Delta#eta", nBinsdEta, mindEta, maxdEta);
  h_bHt_tWb_WBoson_bHt_tWb_BQuark_dEta           = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_bHt_tWb_BQuark_dEta"          , ";#Delta#eta", nBinsdEta, mindEta, maxdEta);
  h_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark_dEta     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_BQuark_bHt_tWb_WBoson_Quark_dEta"    , ";#Delta#eta", nBinsdEta, mindEta, maxdEta);
  h_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark_dEta", ";#Delta#eta", nBinsdEta, mindEta, maxdEta);

  h_bHt_HToHW_HBoson_Taus_dPhi                  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Taus_dPhi"                  , ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_Hpm_bHt_tWb_WBoson_dPhi                 = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_Hpm_bHt_tWb_WBoson_dPhi"                 , ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_Hpm_bHt_TQuark_dPhi                     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_Hpm_bHt_TQuark_dPhi"                     , ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_TQuark_bHt_tWb_BQuark_dPhi              = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_TQuark_bHt_tWb_BQuark_dPhi"              , ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_tWb_WBoson_bHt_tWb_BQuark_dPhi          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_WBoson_bHt_tWb_BQuark_dPhi"          , ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi); 
  h_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark_dPhi    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_BQuark_bHt_tWb_WBoson_Quark_dPhi"    , ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark_dPhi= fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark_dPhi", ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_HToHW_WBoson_bHt_tWb_WBoson_dPhi        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_bHt_tWb_WBoson_dPhi"        , ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_HToHW_HBoson_Tau1_MET_dPhi        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau1_MET_dPhi"        , ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_HToHW_HBoson_Tau2_MET_dPhi        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau2_MET_dPhi"        , ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_BQuark_dR  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau1_bHt_tWb_BQuark_dR"  , ";#DeltaR"   , nBinsdR, mindR, maxdR);
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_BQuark_dEta= fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau1_bHt_tWb_BQuark_dEta", ";#Delta#eta", nBinsdEta, mindEta, maxdEta);
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_BQuark_dPhi= fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau1_bHt_tWb_BQuark_dPhi", ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_BQuark_dR  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau2_bHt_tWb_BQuark_dR"  , ";#DeltaR"   , nBinsdR, mindR, maxdR);
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_BQuark_dEta= fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau2_bHt_tWb_BQuark_dEta", ";#Delta#eta", nBinsdEta, mindEta, maxdEta); 
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_BQuark_dPhi= fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau2_bHt_tWb_BQuark_dPhi", ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Quark_dR   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Quark_dR"  , ";#DeltaR"   , nBinsdR, mindR, maxdR);
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Quark_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Quark_dEta", ";#Delta#eta", nBinsdEta, mindEta, maxdEta); 
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Quark_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Quark_dPhi", ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Quark_dR   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Quark_dR"  , ";#DeltaR"   , nBinsdR, mindR, maxdR);
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Quark_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Quark_dEta", ";#Delta#eta", nBinsdEta, mindEta, maxdEta); 
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Quark_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Quark_dPhi", ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Antiquark_dR   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Antiquark_dR"  , ";#DeltaR"   , nBinsdR, mindR, maxdR);
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Antiquark_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Antiquark_dEta", ";#Delta#eta", nBinsdEta, mindEta, maxdEta); 
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Antiquark_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Antiquark_dPhi", ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Antiquark_dR   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Antiquark_dR"  , ";#DeltaR"   , nBinsdR, mindR, maxdR);
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Antiquark_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Antiquark_dEta", ";#Delta#eta", nBinsdEta, mindEta, maxdEta); 
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Antiquark_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Antiquark_dPhi", ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_HToHW_HBoson_Tau1_bHt_HToHW_WBoson_Leptons_dR   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau1_bHt_HToHW_WBoson_Leptons_dR"  , ";#DeltaR"   , nBinsdR, mindR, maxdR); 
  h_bHt_HToHW_HBoson_Tau1_bHt_HToHW_WBoson_Leptons_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau1_bHt_HToHW_WBoson_Leptons_dEta", ";#Delta#eta", nBinsdEta, mindEta, maxdEta); 
  h_bHt_HToHW_HBoson_Tau1_bHt_HToHW_WBoson_Leptons_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau1_bHt_HToHW_WBoson_Leptons_dPhi", ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_HToHW_HBoson_Tau2_bHt_HToHW_WBoson_Leptons_dR   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau2_bHt_HToHW_WBoson_Leptons_dR"  , ";#DeltaR"   , nBinsdR, mindR, maxdR);
  h_bHt_HToHW_HBoson_Tau2_bHt_HToHW_WBoson_Leptons_dEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau2_bHt_HToHW_WBoson_Leptons_dEta", ";#Delta#eta", nBinsdEta, mindEta, maxdEta); 
  h_bHt_HToHW_HBoson_Tau2_bHt_HToHW_WBoson_Leptons_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_HBoson_Tau2_bHt_HToHW_WBoson_Leptons_dPhi", ";#Delta#phi (rads)", nBinsdPhi, mindPhi, maxdPhi);

  h_bHt_Hpm_bHt_tWb_WBoson_dRap                  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_Hpm_bHt_tWb_WBoson_dRap"                 , ";#Delta#omega", nBinsdRap, mindRap, maxdRap);
  h_bHt_Hpm_bHt_TQuark_dRap                      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_Hpm_bHt_TQuark_dRap"                     , ";#Delta#omega", nBinsdRap, mindRap, maxdRap);
  h_bHt_HToHW_WBoson_bHt_tWb_WBoson_dRap         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_HToHW_WBoson_bHt_tWb_WBoson_dRap"        , ";#Delta#omega", nBinsdRap, mindRap, maxdRap);
  h_bHt_TQuark_bHt_tWb_BQuark_dRap               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_TQuark_bHt_tWb_BQuark_dRap"              , ";#Delta#omega", nBinsdRap, mindRap, maxdRap);
  h_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark_dRap     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_BQuark_bHt_tWb_WBoson_Quark_dRap"    , ";#Delta#omega", nBinsdRap, mindRap, maxdRap);
  h_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark_dRap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark_dRap", ";#Delta#omega", nBinsdRap, mindRap, maxdRap);
  
  h_bHt_tWb_bqq_Pt         = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_bqq_Pt"        , ";p_{T} (GeV/c^{2})",  nBinsPt  , minPt  , maxPt);
  h_bHt_tWb_bqq_Rap        = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_bqq_Rap"       , ";#omega"           ,  nBinsRap , minRap , maxRap);
  h_bHt_tWb_bqq_Mass       = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_bqq_Mass"      , ";M (GeV/c^{2})"    ,  nBinsM   , minM   , maxM);
  h_bHt_tWb_bqq_dRMax_dPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_bqq_dRMax_dPhi", ";#Delta#phi (rads)",  nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_tWb_bqq_dRMax_dRap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_bqq_dRMax_dRap", ";#Delta#omega"     ,  nBinsdRap, mindRap, maxdRap);
  h_bHt_tWb_bqq_dRMax_dR   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "bHt_tWb_bqq_dRMax_dR"  , ";#DeltaR"          ,  nBinsdR  , mindR  , maxdR);  
  
  // GenJets
  h_GenJets_N   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJets_N" , ";genJet multiplicity", 30, 0.0, 30.0);
  h_GenJet1_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet1_Pt", ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenJet2_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet2_Pt", ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenJet3_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet3_Pt", ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenJet4_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet4_Pt", ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenJet5_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet5_Pt", ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenJet6_Pt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet6_Pt", ";p_{T} (GeV/c)", nBinsPt, minPt, maxPt);
  h_GenJet1_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet1_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_GenJet2_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet2_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_GenJet3_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet3_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_GenJet4_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet4_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_GenJet5_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet5_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_GenJet6_Eta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "GenJet6_Eta", ";#eta", nBinsEta, minEta, maxEta);
  h_MaxDiJetMass_Pt    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_Pt"   , ";p_{T} (GeV/c)"    , nBinsPt, minPt, maxPt);
  h_MaxDiJetMass_Eta   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_Eta"  , ";#eta"             , nBinsEta, minEta, maxEta);
  h_MaxDiJetMass_Rap   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_Rap"  , ";#omega"           , nBinsRap, minRap, maxRap);
  h_MaxDiJetMass_Mass  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_Mass" , ";M (GeV/c^{2})"    , 100,  0.0, +2000.0);  
  h_MaxDiJetMass_dEta  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_dEta" , ";#Delta#eta"       , nBinsdEta, mindEta, maxdEta);
  h_MaxDiJetMass_dRap  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_dRap" , ";#Delta#omega"     , nBinsdRap, mindRap, maxdRap);
  h_MaxDiJetMass_dPhi  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_dPhi" , ";#Delta#phi"       , nBinsdPhi, mindPhi, maxdPhi);  
  h_MaxDiJetMass_dR    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_dR"   , ";#DeltaR"          , nBinsdR, mindR, maxdR);  
  h_MaxDiJetMass_dRrap = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, th1, "MaxDiJetMass_dRrap", ";#DeltaR_{#omega}" , nBinsdR, mindR, maxdR);  

  // TH2 
  h_bHt_HToHW_HBoson_Taus_Pt_Vs_Pt     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "bHt_HToHW_HBoson_Taus_Pt_Vs_Pt", ";p_{T}^{#tau_{h,1}} (GeV/c);p_{T}^{#tau_{h,2}} (GeV/c)", nBinsPt, minPt, maxPt, nBinsPt, minPt, maxPt);
  h_bHt_HToHW_HBoson_Taus_Eta_Vs_Eta   = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "bHt_HToHW_HBoson_Taus_Eta_Vs_Eta", ";#eta^{#tau_{h,1}};#eta^{#tau_{h,2}}", nBinsEta, minEta, maxEta, nBinsEta, minEta, maxEta);
  h_bHt_HToHW_HBoson_Taus_Phi_Vs_Phi   = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "bHt_HToHW_HBoson_Taus_Phi_Vs_Phi", ";#phi^{#tau_{h,1}} (rads);#phi^{#tau_{h,2}} (rads)", nBinsPhi, minPhi, maxPhi, nBinsPhi, minPhi, maxPhi);
  h_bHt_HToHW_HBoson_Taus_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "bHt_HToHW_HBoson_Taus_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_bHt_HToHW_HBoson_Taus_Pt1_Vs_dR    = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "bHt_HToHW_HBoson_Taus_Pt1_Vs_dR", ";p_{T}^{#tau_{h,1}} (GeV/c);#DeltaR_{#tau_{h}}", nBinsPt, minPt, maxPt, nBinsdR, mindR, maxdR);  
  h_bHt_HToHW_HBoson_Taus_Pt2_Vs_dR    = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "bHt_HToHW_HBoson_Taus_Pt2_Vs_dR", ";p_{T}^{#tau_{h,2}} (GeV/c);#DeltaR_{#tau_{h}}", nBinsPt, minPt, maxPt, nBinsdR, mindR, maxdR);  

  h_bHt_tWb_bqq_dRMax_dRap_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "bHt_tWb_bqq_dRMax_dRap_Vs_dPhi", ";#Delta#omega;#Delta#phi (rads)", nBinsdRap, mindRap, maxdRap, nBinsdPhi, mindPhi, maxdPhi);
  
  h_Jet1Jet2_dEta_Vs_Jet3Jet4_dEta = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Jet1Jet2_dEta_Vs_Jet3Jet4_dEta", 
								";#Delta#eta(j_{1},j_{2});#Delta#eta(j_{3},j_{4})", nBinsdEta, mindEta, maxdEta, nBinsdEta, mindEta, maxdEta);
  h_Jet1Jet2_dPhi_Vs_Jet3Jet4_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Jet1Jet2_dPhi_Vs_Jet3Jet4_dPhi", 
								";#Delta#phi(j_{1},j_{2}) (rads);#Delta#phi(j_{3},j_{4}) (rads)", nBinsdPhi, mindPhi, maxdPhi, nBinsdPhi, mindPhi, maxdPhi);
  h_Jet1Jet2_dEta_Vs_Jet1Jet2_Mass = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Jet1Jet2_dEta_Vs_Jet1Jet2_Mass", 
								";#Delta#eta(j_{1},j_{2});M(j_{1},j_{2}) (GeV/c^{2})", nBinsdEta, mindEta, maxdEta, nBinsM, minM, maxM);
  h_Jet3Jet4_dEta_Vs_Jet3Jet4_Mass = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Jet3Jet4_dEta_Vs_Jet3Jet4_Mass", 
								";#Delta#eta(j_{3},j_{4});M(j_{4},j_{4}) (GeV/c^{2})", nBinsdEta, mindEta, maxdEta, nBinsM, minM, maxM);

  h_MaxDiJetMass_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "MaxDiJetMass_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)"  , nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_MaxDiJetMass_dRap_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "MaxDiJetMass_dRap_Vs_dPhi", ";#Delta#omega;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);

  h_Bquark_Antiquark_dEta_Vs_dPhi = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Bquark_Antiquark_dEta_Vs_dPhi", ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_Bquark_Quark_dEta_Vs_dPhi     = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Bquark_Quark_dEta_Vs_dPhi"    , ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);
  h_Quark_Antiquark_dEta_Vs_dPhi  = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, th2, "Quark_Antiquark_dEta_Vs_dPhi" , ";#Delta#eta;#Delta#phi (rads)", nBinsdEta, mindEta, maxdEta, nBinsdPhi, mindPhi, maxdPhi);

  return;
}

void KinematicsHToHW::setupBranches(BranchManager& branchManager) {
  fEvent.setupBranches(branchManager);
}


void KinematicsHToHW::process(Long64_t entry) {

  if ( !fEvent.isMC() ) return;
  
  // Create MCools object
  MCTools mcTools(fEvent);
  
  // Increment Counter
  cAllEvents.increment();

  // if (entry != 383) return;
  // if (1) std::cout << "\n=== Event " << entry << std::endl;

  //================================================================================================
  // 1) Apply trigger
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== Trigger" << std::endl;
  if ( !(fEvent.passTriggerDecision()) ) return;
  cTrigger.increment();


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
  if (cfg_Verbose)
    {
      std::cout << "\nnElectrons = " << selectedElectrons.size() << std::endl;
      for (auto& p: selectedElectrons) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads), status = " << p.status() << std::endl;
    }
  if ( selectedElectrons.size() > 0 ) return;
  cElectronVeto.increment();


  //================================================================================================
  // 5) Muon selection
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== Muon selection" << std::endl;
  vector<genParticle> selectedMuons = GetGenParticles(fEvent.genparticles().getGenParticles(), cfg_MuonPtCut, cfg_MuonEtaCut, 13, true, false);
  if (cfg_Verbose)
    {
      std::cout << "nMuons = " << selectedMuons.size() << std::endl;
      for (auto& p: selectedMuons) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads), status = " << p.status() << std::endl;
    } 
  if ( selectedMuons.size() < 1 ) return;
  cMuonVeto.increment();

  //================================================================================================
  // 6) Tau selection
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== Tau selection" << std::endl;
  vector<genParticle> selectedTaus = GetGenParticles(fEvent.genparticles().getGenParticles(), cfg_TauPtCut, cfg_TauEtaCut, 15, true, false);
  if (cfg_Verbose)
    {
      std::cout << "nTaus = " << selectedTaus.size() << std::endl;
      for (auto& p: selectedTaus) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads), status = " << p.status() << std::endl;
    }
  if ( selectedTaus.size() < 2 ) return;
  cTauVeto.increment();


  //================================================================================================
  // 7) Jet Selection
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== Jet Selection" << std::endl;
  vector<GenJet> selectedJets = GetGenJets(fEvent.genjets(), cfg_JetPtCuts, cfg_JetEtaCuts);
  if (cfg_Verbose)
    {
      std::cout << "nJets = " << selectedJets.size() << std::endl;
      for (auto& p: selectedJets) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads)" << std::endl;
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
      //if ( abs( m.pdgId() != 6) ) continue; // fixme ?
      selectedBQuarks.push_back(p);
    }
  std::sort( selectedBQuarks.begin(), selectedBQuarks.end(), PtComparator() );  

  if (cfg_Verbose)
    {
      std::cout << "\nnBQuarks = " << selectedBQuarks.size() << std::endl;
      for (auto& p: selectedBQuarks)
	{
	  mcTools.PrintGenParticle(p);
	}
      // for (auto& p: selectedBQuarks) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads), status = " << p.status() << std::endl;
    }

  // Match b-quarks with GenJets
  vector<GenJet> selectedBJets = GetGenJets(selectedJets, cfg_BJetPtCuts, cfg_BJetEtaCuts, selectedBQuarks);
  if (cfg_Verbose)
    {
      std::cout << "nBJets = " << selectedBJets.size() << std::endl;
      for (auto& p: selectedBJets) std::cout << "\tpT = " << p.pt() << " (GeV/c), eta = " << p.eta() << ", phi = " << p.phi() << " (rads)" << std::endl;
      std::cout << "" << std::endl;
    }

  // Get selected jets excluding the matched bjets
  bool isBJet = false;
  std::vector<math::XYZTLorentzVector> selJets_NoBJets_p4;
  // For-loop: Selected jets
  for (auto& jet: selectedJets) 
    {
      isBJet = false;
	    
      // For-loop: Selected bjets
      for (auto& bjet: selectedBJets) 
	{
	  double dR = ROOT::Math::VectorUtil::DeltaR(jet.p4(), bjet.p4());
	  if (dR < 0.01) isBJet = true;
	}
      if (isBJet) continue;
      jet_p4 = jet.p4();
      selJets_NoBJets_p4.push_back(jet_p4);
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
  if (cfg_Verbose) std::cout << "=== MET = " << fEvent.genMET().et() << std::endl;      
  cMETSelection.increment();


  //================================================================================================
  // 11) Topology selection 
  //================================================================================================
  /*
  float C, D, H2;
  float Circularity;
  float y23, Sphericity, SphericityT, Aplanarity, Planarity, Y; // functions to return values when properly implemented
  float HT, JT, MHT, Centrality;
  vector<float> a = GetMomentumTensorEigenValues(selJets_p4, C, D, H2);
  vector<float> b = GetMomentumTensorEigenValues2D(selJets_p4, Circularity);
  vector<float> c = GetSphericityTensorEigenValues(selJets_p4, y23, Sphericity, SphericityT, Aplanarity, Planarity, Y);
  double alphaT   = GetAlphaT(selJets_p4, HT, JT, MHT, Centrality);
  
  // Apply cuts
  if ( !cfg_CparameterCut.passedCut(C) ) return;
  if ( !cfg_DparameterCut.passedCut(D) ) return;
  if ( !cfg_FoxWolframMomentCut.passedCut(H2) ) return;
  if ( !cfg_CircularityCut.passedCut(Circularity) ) return;
  if ( !cfg_Y23Cut.passedCut(y23) ) return;
  if ( !cfg_SphericityCut.passedCut(Sphericity) ) return;
  if ( !cfg_AplanarityCut.passedCut(Aplanarity) ) return;
  if ( !cfg_PlanarityCut.passedCut(Planarity) ) return;
  if ( !cfg_CentralityCut.passedCut(Centrality) ) return;
  if ( !cfg_AlphaTCut.passedCut(alphaT) ) return;
  cTopologySelection.increment();
  */

  //================================================================================================
  // 12) Top selection 
  //================================================================================================
  // if (cfg_Verbose) std::cout << "=== Top Selection" << std::endl;
  //cTopSelection.increment();


  //================================================================================================
  // All Selections
  //================================================================================================
  if (cfg_Verbose) std::cout << "=== All Selections" << std::endl;
  cSelected.increment();

  
  if (cfg_Verbose) 
    {
      std::cout << "" << std::endl;
      std::cout << "nElectrons = " << selectedElectrons.size() << std::endl;
      std::cout << "nMuons     = " << selectedMuons.size() << std::endl;
      std::cout << "nTaus      = " << selectedTaus.size() << std::endl;
      std::cout << "nJets      = " << selectedJets.size() << std::endl;
      std::cout << "nBJets     = " << selectedBJets.size() << " (nBQuarks = " << selectedBQuarks.size() << ")" << std::endl;
      std::cout << "MET        = " << fEvent.genMET().et() << std::endl;
    }


  ///////////////////////////////////////////////////////////////////////////
  // GenParticles
  ///////////////////////////////////////////////////////////////////////////
  cInclusive.increment();
  unsigned int nGenP_Hpm              = 0;
  unsigned int nGenP_TQuarks          = 0;
  unsigned int nGenP_WBosons          = 0;
  unsigned int nGenP_BQuarks          = 0;
  unsigned int nGenP_HToHW_Quarks     = 0;
  unsigned int nGenP_HToHW_HBoson_Taus = 0;
  unsigned int nGenP_HToHW_HBoson_TauLeptonic_plus = 0;
  unsigned int nGenP_HToHW_HBoson_TauHadronic_plus = 0;
  unsigned int nGenP_HToHW_HBoson_TauLeptonic_minus = 0;
  unsigned int nGenP_HToHW_HBoson_TauHadronic_minus = 0;
  
  unsigned int nGenP_HToHW_WLeptons   = 0;
  unsigned int nGenP_HToHW_WNeutrinos = 0;
  unsigned int nGenP_HBosons          = 0;
  unsigned int nGenP_Taus             = 0;
  unsigned int nGenP_Taus_Higgs       = 0;
  unsigned int nGenP_tWb_Quarks       = 0;
  unsigned int nGenP_tWb_WLeptons     = 0;
  unsigned int nGenP_tWb_WNeutrinos   = 0;

  // 4-momenta
  math::XYZTLorentzVector bHt_Hpm_p4;
  math::XYZTLorentzVector bHt_TQuark_p4;
  math::XYZTLorentzVector bHt_HToHW_HBoson_p4;
  math::XYZTLorentzVector bHt_HToHW_WBoson_p4;
  math::XYZTLorentzVector bHt_HToHW_WBoson_Quark_p4;
  math::XYZTLorentzVector bHt_HToHW_WBoson_Antiquark_p4;
  math::XYZTLorentzVector bHt_HToHW_WBoson_Lepton_p4;
  math::XYZTLorentzVector bHt_HToHW_WBoson_Neutrino_p4;
  math::XYZTLorentzVector bHt_HToHW_HBoson_Tau_p4;
  math::XYZTLorentzVector bHt_HToHW_HBoson_Antitau_p4;
  math::XYZTLorentzVector bHt_tWb_WBoson_p4;
  math::XYZTLorentzVector bHt_tWb_WBoson_Quark_p4;
  math::XYZTLorentzVector bHt_tWb_WBoson_Antiquark_p4;
  math::XYZTLorentzVector bHt_tWb_WBoson_Lepton_p4;
  math::XYZTLorentzVector bHt_tWb_WBoson_Neutrino_p4;
  math::XYZTLorentzVector bHt_tWb_BQuark_p4;

  // Define the table
  Table table("Evt | Index | PdgId | Status | Charge | Pt | Eta | Phi | E | Vertex (mm) | D0 (mm) | Lxy (mm) | Mom | Daughters", "Text"); //LaTeX or Text

  int row = 0;
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

    // Associated genParticles
    std::vector<genParticle> genP_daughters;
    for (unsigned int i=0; i < p.daughters().size(); i++) genP_daughters.push_back(fEvent.genparticles().getGenParticles()[p.daughters().at(i)]);
    std::vector<genParticle> genP_mothers;
    for (unsigned int i=0; i < p.mothers().size(); i++) genP_mothers.push_back(fEvent.genparticles().getGenParticles()[p.mothers().at(i)]);
    std::vector<genParticle> genP_grandMothers;
    std::vector<unsigned int> genGMoms_index;
    std::vector<unsigned int> genGMoms_pdgId;
    std::vector<unsigned int> genMoms_index;
    std::vector<unsigned int> genMoms_pdgId;
    std::vector<unsigned int> genDaus_index;
    std::vector<unsigned int> genDaus_pdgId;
  
    // Removed real vertex from Tree (to save size)
    ROOT::Math::XYZPoint vtxIdeal;
    vtxIdeal.SetXYZ(0, 0, 0);
    double genP_vtxX = vtxIdeal.X(); // p.vtxX()*10; // in mm
    double genP_vtxY = vtxIdeal.Y(); // p.vtxY()*10; // in mm
    double genP_vtxZ = vtxIdeal.Z(); // p.vtxZ()*10; // in mm
  
    // Daughter, Mom and Grand-mom properties    
    genParticle m;
    genParticle g;
    genParticle d;
    
    if (genP_daughters.size() > 0) d = genP_daughters.at(0);
    if (p.mothers().size() > 0)
      {
	m = genP_mothers.at(0); // fixme
	for (unsigned int i=0; i < m.mothers().size(); i++) genP_grandMothers.push_back(fEvent.genparticles().getGenParticles()[m.mothers().at(i)]);
	if (m.mothers().size() > 0) g = genP_grandMothers.at(0); // fixme
      } 

    // For convenience, save the pdgIds in vectors
    for (unsigned int i=0; i < genP_grandMothers.size(); i++) 
      {
	if (genP_grandMothers.at(i).index() < 0) continue;
	genGMoms_index.push_back(genP_grandMothers.at(i).index());
	genGMoms_pdgId.push_back(genP_grandMothers.at(i).pdgId());
      }
    for (unsigned int i=0; i < genP_mothers.size(); i++) 
      {
	if (genP_mothers.at(i).index() < 0) continue;
	genMoms_index.push_back(genP_mothers.at(i).index());
	genMoms_pdgId.push_back(genP_mothers.at(i).pdgId());
      }
    for (unsigned int i=0; i < genP_daughters.size(); i++) 
      {
	if (genP_daughters.at(i).index() < 0) continue;
	genDaus_index.push_back(genP_daughters.at(i).index());
	genDaus_pdgId.push_back(genP_daughters.at(i).pdgId());
      }

    // Properties that need to be calculated
    bool bIsLastCopy = std::find(genDaus_pdgId.begin(), genDaus_pdgId.end(), genP_pdgId) == genDaus_pdgId.end();
    double genP_Lxy  = 0.0; 
    double genP_d0   = 0.0;
    if (genP_daughters.size() > 0 && genP_mothers.size() > 0)
      {
	genP_d0  = mcTools.GetD0 (p, m, d, vtxIdeal); // in mm
	genP_Lxy = mcTools.GetLxy(p, m, d, vtxIdeal); // in mm
      }
   
    // Print genParticle properties or decay tree ?
    if (0)
      {
    	mcTools.PrintGenParticle(p);
     	mcTools.PrintGenDaughters(p);
       }

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
    table.AddRowColumn(row, "(" + auxTools.ToString(genP_vtxX, 3) + ", " + auxTools.ToString(genP_vtxY, 3)  + ", " + auxTools.ToString(genP_vtxZ, 3) + ")" );
    table.AddRowColumn(row, auxTools.ToString(genP_d0 , 3)     );
    table.AddRowColumn(row, auxTools.ToString(genP_Lxy, 3)     );	
    if (genMoms_index.size() < 6)
      {
	table.AddRowColumn(row, auxTools.ConvertIntVectorToString(genMoms_index) );
      }
    else table.AddRowColumn(row, ".. Too many .." );
    if (genDaus_index.size() < 6)
      {
	table.AddRowColumn(row, auxTools.ConvertIntVectorToString(genDaus_index) );
	// table.AddRowColumn(row, auxTools.ConvertIntVectorToString(genDaus_pdgId) );
      }
    else table.AddRowColumn(row, ".. Too many .." );
    row++;
	       

    if (cfg_Verbose) std::cout << "=== taus" << std::endl;

    if ( abs(genP_pdgId) == 15)
      {
	if (0) cout << genP_pdgId <<" Tau found " << endl;
	// Get last copy (after all radiative corrections)
	if (!bIsLastCopy) continue;
	nGenP_Taus++;

	// N.B: Use with caution
	//if (genP_mothers.size() > 0)
	//	  {
	genParticle firstCopy;
	genParticle firstMom;
	std::vector<genParticle> genP_allCopies = GetAllPreviousCopies(fEvent.genparticles().getGenParticles(), p);
	firstCopy = genP_allCopies.at(0);
	if (firstCopy.mothers().size() > 0) firstMom = fEvent.genparticles().getGenParticles()[firstCopy.mothers().at(0)];
	//	  }

	// H0 (or h0) -> tau tau
	if ( abs(firstMom.pdgId()) == 35 ||  abs(firstMom.pdgId()) == 25 )
	  {
	    if (0) cout << genP_pdgId <<" Tau found from Higgs" << endl;
	    nGenP_Taus_Higgs++;
	    
	    if (genP_pdgId > 0) 
	      {
		
		for (auto& d: genP_daughters)
		  {
		
		    if(0) cout << "+ve Tau dauughters pt, pdgId:=   "<<  d.pt() <<" , " << d.pdgId() << endl;
		    if (abs(d.pdgId()) == 11 || abs(d.pdgId()) == 13){
			//mcTools.PrintGenParticle(p);
			if(nGenP_HToHW_HBoson_TauLeptonic_plus == 0 && genP_daughters.size() >= 1){
			  nGenP_HToHW_HBoson_TauLeptonic_plus++;
			} 
			bHt_HToHW_HBoson_Tau_p4 += d.p4();
		    }
		  }

		for (auto& d: genP_daughters)
		  {
		    if(nGenP_HToHW_HBoson_TauLeptonic_plus == 1) nGenP_HToHW_HBoson_TauHadronic_plus = 0;
		    else{
		      nGenP_HToHW_HBoson_TauHadronic_plus = 1;
		      bHt_HToHW_HBoson_Tau_p4 += d.p4();
		    }
		  }
	      }
	    else
	      {
		// For-loop: All daughters
		for (auto& d: genP_daughters)
		  {
		    if(0) cout << "-ve Tau dauughters pt, pdgId:=   "<<  d.pt() <<" , " << d.pdgId() << endl;
		    if (abs(d.pdgId()) == 11 || abs(d.pdgId()) == 13){
			//mcTools.PrintGenParticle(p);
			if(nGenP_HToHW_HBoson_TauLeptonic_minus == 0 && genP_daughters.size() >= 1){
			  nGenP_HToHW_HBoson_TauLeptonic_minus++;
			} 
			bHt_HToHW_HBoson_Antitau_p4 += d.p4();
		    }
		  }

		for (auto& d: genP_daughters)
		  {
		    if(nGenP_HToHW_HBoson_TauLeptonic_minus == 1) { 
		      nGenP_HToHW_HBoson_TauHadronic_minus = 0;
		    }
		    else{
		      nGenP_HToHW_HBoson_TauHadronic_minus = 1;
		      bHt_HToHW_HBoson_Antitau_p4 += d.p4();
		    }
		  }

	      }

	  } //from h0 or H0
      
      } // taus
  
    if (cfg_Verbose) std::cout << "=== t quarks" << std::endl;
    if ( abs(genP_pdgId) == 6)
      {

	// Get last copy (after all radiative corrections)
	if (!bIsLastCopy) continue;
	nGenP_TQuarks++;
      	bHt_TQuark_p4 = p.p4();
      } //tops
  

    if (cfg_Verbose) std::cout << "=== H+ or H- bosons-" << std::endl;
    if ( abs(genP_pdgId) == 37)
      {
	// Get last copy (after all radiative corrections)
	if (!bIsLastCopy) continue;  // if (!p.isLastCopy()) continue;
	nGenP_Hpm++;
	bHt_Hpm_p4 = p.p4();
	cbHt_Hpm.increment();
	
      } // Hpm
	    

    if (cfg_Verbose) std::cout << "=== H0 or h0 bosons" << std::endl;
    if ( (abs(genP_pdgId) == 25) || (abs(genP_pdgId) == 35) )// h0 = H0_{1}, H0 = H0_{2}
      {
	// Get last copy (after all radiative corrections)
	if (!bIsLastCopy) continue;
	nGenP_HBosons++;
	bHt_HToHW_HBoson_p4 = p.p4();
	cbHt_HToHW_HBoson.increment();
	
      } // H0 or h0
    

    if (cfg_Verbose) std::cout << "=== W bosons" << std::endl;
    if ( abs(genP_pdgId) == 24)
      {
	// Get last copy (after all radiative corrections)
	if (!bIsLastCopy) continue;

	// Use with caution
	genParticle firstCopy;
	genParticle firstMom;
	std::vector<genParticle> genP_allCopies = GetAllPreviousCopies(fEvent.genparticles().getGenParticles(), p);
	firstCopy = genP_allCopies.at(0);
	if (firstCopy.mothers().size() > 0) firstMom = fEvent.genparticles().getGenParticles()[firstCopy.mothers().at(0)];

	if ( abs(firstMom.pdgId()) == 37)
	  {
	    bHt_HToHW_WBoson_p4 = p.p4();
	    nGenP_WBosons++;
	    cbHt_HToHW_WBoson.increment();

	    // For-loop: All daughters
	    for (auto& d: genP_daughters) 
	      {
		
		// Quarks
		if ( mcTools.IsQuark(d.pdgId()) )
		  {	
		    if (d.pdgId() > 0)
		      {
			cbHt_HToHW_WBoson_Quark.increment();
			bHt_HToHW_WBoson_Quark_p4 = d.p4(); // not last copy but ok. no biggie
			nGenP_HToHW_Quarks++;
		      }
		    else
		      {
			cbHt_HToHW_WBoson_Antiquark.increment();
			bHt_HToHW_WBoson_Antiquark_p4 = d.p4();
			nGenP_HToHW_Quarks++;
		      }
		  }
		// Leptons
		else if ( mcTools.IsLepton(d.pdgId() ) )
		  {
		    cbHt_HToHW_WBoson_Leptons.increment();
		    
		    if ( abs(d.pdgId()) == 11)
		      {		      
			bHt_HToHW_WBoson_Lepton_p4 = d.p4(); // not last copy but ok. no biggie
			nGenP_HToHW_WLeptons++;
		      }
		    else if ( abs(d.pdgId()) == 13)
		      {
			bHt_HToHW_WBoson_Lepton_p4 = d.p4(); // not last copy but ok. no biggie
			nGenP_HToHW_WLeptons++;
		      }		  
		    else if ( abs(d.pdgId()) == 15)
		      {
			bHt_HToHW_WBoson_Lepton_p4 = d.p4(); // not last copy but ok. no biggie
			nGenP_HToHW_WLeptons++;
		      }		  
		    else
		      {
			nGenP_HToHW_WNeutrinos++;
			bHt_HToHW_WBoson_Neutrino_p4 = d.p4(); // not last copy but ok. no biggie
		      }
		  } // daugher IsLepton
		else
		  {
		    throw hplus::Exception("Logic") << "KinematicsHToHW::process() W daughters whose origins are not accounted for. Need to rethink this.";
		  }
	      }// for-loop: daughters
	  }
	else if ( abs(firstMom.pdgId()) == 6)
	  {
	    bHt_tWb_WBoson_p4 = p.p4();
	    nGenP_WBosons++;
	    cbHt_tWb_WBoson.increment();

	    // For-loop: All daughters
	    for (auto& d: genP_daughters) 
	      {
		
		// Quarks
		if ( mcTools.IsQuark(d.pdgId()) )
		  {	
		    if (d.pdgId() > 0)
		      {
			cbHt_tWb_WBoson_Quark.increment();
			bHt_tWb_WBoson_Quark_p4 = d.p4(); // not last copy but ok. no biggie
			nGenP_tWb_Quarks++;
		      }
		    else
		      {
			cbHt_tWb_WBoson_Antiquark.increment();
			bHt_tWb_WBoson_Antiquark_p4 = d.p4();
			nGenP_tWb_Quarks++;
		      }
		  }
		// Leptons
		else if ( mcTools.IsLepton(d.pdgId() ) )
		  {
		    cbHt_tWb_WBoson_Leptons.increment();
		    
		    if ( abs(d.pdgId()) == 11)
		      {		      
			bHt_tWb_WBoson_Lepton_p4 = d.p4(); // not last copy but ok. no biggie
			nGenP_tWb_WLeptons++;
		      }
		    else if ( abs(d.pdgId()) == 13)
		      {
			bHt_tWb_WBoson_Lepton_p4 = d.p4(); // not last copy but ok. no biggie
			nGenP_tWb_WLeptons++;
		      }		  
		    else if ( abs(d.pdgId()) == 15)
		      {
			bHt_tWb_WBoson_Lepton_p4 = d.p4(); // not last copy but ok. no biggie
			nGenP_tWb_WLeptons++;
		      }		  
		    else
		      {
			nGenP_tWb_WNeutrinos++;
			bHt_tWb_WBoson_Neutrino_p4 = d.p4(); // not last copy but ok. no biggie
		      }
		  } // daugher IsLepton
		else
		  {
		    throw hplus::Exception("Logic") << "KinematicsHToHW::process() W daughters whose origins are not accounted for. Need to rethink this.";
		  }
	      }// for-loop: daughters
	  } // Tops
	    else
	      {
		mcTools.PrintGenParticle(firstCopy);
		mcTools.PrintGenParticle(firstMom);
		mcTools.PrintGenParticle(p);
	        table.Print();
		throw hplus::Exception("Logic") << "KinematicsHToHW::process() Could not determine W boson original mother.";
	      }

      }// WBoson

    if (cfg_Verbose) std::cout << "=== b quarks" << std::endl;
    if ( abs(genP_pdgId) == 5)
      {
	// Get last copy (after all radiative corrections)
	if (!bIsLastCopy) continue;

	// Use with caution
	vector<genParticle> genP_allCopies = GetAllPreviousCopies(fEvent.genparticles().getGenParticles(), p);
	genParticle firstCopy = genP_allCopies.at(0);
	vector<int> firstMoms;
	if (firstCopy.mothers().size() > 0) 
	  {
	    for (unsigned int i=0; i < firstCopy.mothers().size(); i++) firstMoms.push_back(fEvent.genparticles().getGenParticles()[firstCopy.mothers().at(i)].pdgId());
	  }

	// b->tH+, t->bW, b
	// if ( abs(firstMom.pdgId()) == 6) 
	if ( find(firstMoms.begin(), firstMoms.end(), 6) != firstMoms.end() || find(firstMoms.begin(), firstMoms.end(), -6) != firstMoms.end() )
	  // || find(genMoms_pdgId.begin(), genMoms_pdgId.end(), 6) != genMoms_pdgId.end() || find(genMoms_pdgId.begin(), genMoms_pdgId.end(), -6) != genMoms_pdgId.end() )
	  // if ( abs(firstMom.pdgId()) != 21 && abs(firstMom.pdgId()) != 2212) 
	  {
	    bHt_tWb_BQuark_p4 = p.p4();
	    nGenP_BQuarks++;
	    cbHt_tWb_BQuark.increment();
	  }
	else
	  {
	    // mcTools.PrintGenParticle(p);
	    // for (unsigned int i=0; i < firstMoms.size(); i++) cout << i << ") index = " << genP_index << ", genP_pdgId = " << genP_pdgId << ", mom = " << firstMoms.at(i) << endl;
	    // mcTools.PrintGenParticle(p);
	    // mcTools.PrintGenParticle(firstCopy);
	    // mcTools.PrintGenParticle(firstMom);
	    // table.Print();
	  }
	
      }// b quarks

  }// for-loop: genParticles



  //cout << "nGenP_HToHW_HBoson_Taus = " << nGenP_HToHW_HBoson_Taus << endl; 

  if(0) cout << "nGenP_Taus_Higgs:= " << nGenP_Taus_Higgs << endl;

  if (nGenP_Taus_Higgs == 2){
    cbHt_HToHW_HBoson_TauTau.increment();
    if ((nGenP_HToHW_HBoson_TauLeptonic_plus + nGenP_HToHW_HBoson_TauLeptonic_minus) >= 2) cbHt_HToHW_HBoson_TauLeptonic.increment(); 
    if ((nGenP_HToHW_HBoson_TauHadronic_plus + nGenP_HToHW_HBoson_TauHadronic_minus) >= 2) cbHt_HToHW_HBoson_TauHadronic.increment(); 
    if ((nGenP_HToHW_HBoson_TauLeptonic_plus >= 1 || nGenP_HToHW_HBoson_TauLeptonic_minus >=1) && (nGenP_HToHW_HBoson_TauHadronic_plus >= 1 || nGenP_HToHW_HBoson_TauHadronic_minus >= 1)) cbHt_HToHW_HBoson_1Lep1Had.increment(); 
    if(0) cout << "nGenP_HToHW_HBoson_TauLeptonic_plus: " << nGenP_HToHW_HBoson_TauLeptonic_plus << endl;
    if(0) cout << "nGenP_HToHW_HBoson_TauLeptonic_minus: " << nGenP_HToHW_HBoson_TauLeptonic_minus << endl;
    if(0) cout << "nGenP_HToHW_HBoson_TauHadronic_plus: " << nGenP_HToHW_HBoson_TauHadronic_plus << endl;
    if(0) cout << "nGenP_HToHW_HBoson_TauHadronic_minus: " << nGenP_HToHW_HBoson_TauHadronic_minus << endl;
  }


  // Require 2 taus
  
  if (nGenP_Taus_Higgs != 2) return; // fixme - 1tau_h + 1 lepton. use as if in filling 2 tau_h histograms
  else
    {
      //cout << "nGenP_HToHW_HBoson_Taus = " << nGenP_HToHW_HBoson_Taus << endl;
    }
  
  ///////////////////////////////////////////////////////////////////////////
  // GenJets (selected)
  ///////////////////////////////////////////////////////////////////////////
  if (cfg_Verbose) std::cout << "=== GenJets" << std::endl;
  //  if (nGenP_HToHW_HBoson_Taus == 2) 
 
  double genP_HT = (bHt_HToHW_HBoson_Tau_p4.pt() + bHt_HToHW_HBoson_Antitau_p4.pt() 
		    + bHt_HToHW_WBoson_Quark_p4.pt() + bHt_HToHW_WBoson_Antiquark_p4.pt() + bHt_HToHW_WBoson_Lepton_p4.pt() 
		    + bHt_tWb_WBoson_Quark_p4.pt() + bHt_tWb_WBoson_Antiquark_p4.pt() + bHt_tWb_WBoson_Lepton_p4.pt() 
		    + bHt_tWb_BQuark_p4.pt() );

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

  // Trijet system (associated top)
  math::XYZTLorentzVector  bHt_tWb_bqq_p4 = bHt_tWb_WBoson_Quark_p4 + bHt_tWb_WBoson_Antiquark_p4 + bHt_tWb_BQuark_p4;

  // Max separation (again)
  std::vector<math::XYZTLorentzVector> bqq_p4;
  bqq_p4.push_back(bHt_tWb_BQuark_p4);
  bqq_p4.push_back(bHt_tWb_WBoson_Quark_p4);
  bqq_p4.push_back(bHt_tWb_WBoson_Antiquark_p4);
  double deltaRMax = -1.0;
  double bqq_dRap  = -1.0;
  double bqq_dPhi  = -1.0;
  int deltaRMax_i  = -1;
  int deltaRMax_j  = -1;

  // For-loop: All p4 of bqq system
  for (size_t i = 0; i < bqq_p4.size(); i++)    
    {
      for (size_t j = i+1; j < bqq_p4.size(); j++)
	{
	  double deltaR = ROOT::Math::VectorUtil::DeltaR(bqq_p4.at(i), bqq_p4.at(j));
	  if (deltaR > deltaRMax)
	    {
	      deltaRMax   = deltaR;
	      deltaRMax_i = i;
	      deltaRMax_j = j;
	    }
	}
    } // For-loop: All p4 of bqq system

  // bqq_dEta = abs(bqq_p4.at(deltaRMax_i).eta() - bqq_p4.at(deltaRMax_j).eta());
  bqq_dRap = abs(mcTools.GetRapidity(bqq_p4.at(deltaRMax_i) ) - mcTools.GetRapidity(bqq_p4.at(deltaRMax_j) ) );
  bqq_dPhi = abs(ROOT::Math::VectorUtil::DeltaPhi(bqq_p4.at(deltaRMax_i), bqq_p4.at(deltaRMax_j)));
    
  // Print the table with genP info
  if (cfg_Verbose) table.Print();
    
  // Sanity check! Sometimes pruned genParticles are all messed up. And/Or I must improve my code a bit
  //bool bNotOk = (nGenP_Hpm != 1) || (nGenP_HBosons != 1) || (nGenP_TQuarks != 1) || (nGenP_WBosons != 2) || (nGenP_BQuarks != 1) || (nGenP_HToHW_HBoson_Taus != 2);
  bool bNotOk = (nGenP_Hpm != 1) || (nGenP_HBosons != 1) || (nGenP_TQuarks != 1) || (nGenP_WBosons != 2) || (nGenP_BQuarks != 1) || (nGenP_Taus_Higgs != 2);

  if (bNotOk)
    {      
      if (cfg_Verbose)
	{
	  cout << "\n=== Event " << entry << endl;
	  cout << "nGenP_Hpm              = " << nGenP_Hpm               << endl;
	  cout << "nGenP_HBosons          = " << nGenP_HBosons           << endl;
	  cout << "nGenP_TQuarks          = " << nGenP_TQuarks           << endl;
	  cout << "nGenP_WBosons          = " << nGenP_WBosons           << endl;
	  cout << "nGenP_BQuarks          = " << nGenP_BQuarks           << endl;
	  cout << "nGenP_HToHW_HBoson_Taus= " << nGenP_HToHW_HBoson_Taus << endl; // old counter

	  // cout << "nGenP_HToHW_HBoson_TauLeptonic:=   " << nGenP_HToHW_HBoson_TauLeptonic << endl;
	  // cout << "nGenP_HToHW_HBoson_TauHadronic:=   " << nGenP_HToHW_HBoson_TauHadronic << endl;
	  // cout << "nGenP_Taus             = " << nGenP_Taus              << endl;
	  // cout << "nGenP_HToHW_tQuarks    = " << nGenP_HToHW_Quarks      << endl;
	  // cout << "nGenP_HToHW_WLeptons   = " << nGenP_HToHW_WLeptons    << endl;
	  // cout << "nGenP_HToHW_WNeutrinos = " << nGenP_HToHW_WNeutrinos  << endl;
	  // cout << "nGenP_tWb_Quarks       = " << nGenP_tWb_Quarks        << endl;
	  // cout << "nGenP_tWb_WLeptons     = " << nGenP_tWb_WLeptons      << endl;
	  // cout << "nGenP_tWb_WNeutrinos   = " << nGenP_tWb_WNeutrinos    << endl;
	  // table.Print();
	}
      return;
    }

  // Calculations
  double dR_bHt_Hpm_bHt_tWb_WBoson                  = ROOT::Math::VectorUtil::DeltaR( bHt_Hpm_p4, bHt_tWb_WBoson_p4 );
  double dR_bHt_Hpm_bHt_TQuark                      = ROOT::Math::VectorUtil::DeltaR( bHt_Hpm_p4, bHt_TQuark_p4 );
  double dR_bHt_HToHW_WBoson_bHt_tWb_WBoson         = ROOT::Math::VectorUtil::DeltaR( bHt_HToHW_WBoson_p4, bHt_tWb_WBoson_p4 );
  double dR_bHt_TQuark_bHt_tWb_BQuark               = ROOT::Math::VectorUtil::DeltaR( bHt_TQuark_p4, bHt_tWb_BQuark_p4 ); 
  double dR_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark     = ROOT::Math::VectorUtil::DeltaR( bHt_tWb_BQuark_p4, bHt_tWb_WBoson_Quark_p4 );
  double dR_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark = ROOT::Math::VectorUtil::DeltaR( bHt_tWb_BQuark_p4, bHt_tWb_WBoson_Antiquark_p4 );

  double dEta_bHt_Hpm_bHt_tWb_WBoson                  = abs( bHt_Hpm_p4.eta() - bHt_tWb_WBoson_p4.eta() );
  double dEta_bHt_Hpm_bHt_TQuark                      = abs( bHt_Hpm_p4.eta() - bHt_TQuark_p4.eta() );
  double dEta_bHt_HToHW_WBoson_bHt_tWb_WBoson         = abs( bHt_HToHW_WBoson_p4.eta() - bHt_tWb_WBoson_p4.eta() );
  double dEta_bHt_tWb_WBoson_bHt_tWb_BQuark           = abs( bHt_tWb_WBoson_p4.eta() - bHt_tWb_BQuark_p4.eta() );
  double dEta_bHt_TQuark_bHt_tWb_BQuark               = abs( bHt_TQuark_p4.eta()     - bHt_tWb_BQuark_p4.eta() ); 
  double dEta_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark     = abs( bHt_tWb_BQuark_p4.eta() - bHt_tWb_WBoson_Quark_p4.eta() );
  double dEta_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark = abs( bHt_tWb_BQuark_p4.eta() - bHt_tWb_WBoson_Antiquark_p4.eta() );

  double dPhi_bHt_Hpm_bHt_tWb_WBoson                  = abs(ROOT::Math::VectorUtil::DeltaPhi( bHt_Hpm_p4, bHt_tWb_WBoson_p4 ) );
  double dPhi_bHt_Hpm_bHt_TQuark                      = abs(ROOT::Math::VectorUtil::DeltaPhi( bHt_Hpm_p4, bHt_TQuark_p4 ) );
  double dPhi_bHt_HToHW_WBoson_bHt_tWb_WBoson         = abs(ROOT::Math::VectorUtil::DeltaPhi( bHt_HToHW_WBoson_p4, bHt_tWb_WBoson_p4 ) );
  double dPhi_bHt_tWb_WBoson_bHt_tWb_BQuark           = abs(ROOT::Math::VectorUtil::DeltaPhi( bHt_tWb_WBoson_p4, bHt_tWb_BQuark_p4) );
  double dPhi_bHt_TQuark_bHt_tWb_BQuark               = abs(ROOT::Math::VectorUtil::DeltaPhi( bHt_TQuark_p4    , bHt_tWb_BQuark_p4 ) ); 
  double dPhi_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark     = abs(ROOT::Math::VectorUtil::DeltaPhi( bHt_tWb_BQuark_p4, bHt_tWb_WBoson_Quark_p4 ) );
  double dPhi_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark = abs(ROOT::Math::VectorUtil::DeltaPhi( bHt_tWb_BQuark_p4, bHt_tWb_WBoson_Antiquark_p4 ) );

  double dRap_bHt_Hpm_bHt_tWb_WBoson                  = abs( mcTools.GetRapidity(bHt_Hpm_p4) - mcTools.GetRapidity(bHt_tWb_WBoson_p4) );
  double dRap_bHt_Hpm_bHt_TQuark                      = abs( mcTools.GetRapidity(bHt_Hpm_p4) - mcTools.GetRapidity(bHt_TQuark_p4) );
  double dRap_bHt_HToHW_WBoson_bHt_tWb_WBoson         = abs( mcTools.GetRapidity(bHt_HToHW_WBoson_p4) - mcTools.GetRapidity(bHt_tWb_WBoson_p4) );
  double dRap_bHt_TQuark_bHt_tWb_BQuark               = abs( mcTools.GetRapidity(bHt_TQuark_p4) - mcTools.GetRapidity(bHt_tWb_BQuark_p4) ); 
  double dRap_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark     = abs( mcTools.GetRapidity(bHt_tWb_BQuark_p4) - mcTools.GetRapidity(bHt_tWb_WBoson_Quark_p4) );
  double dRap_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark = abs( mcTools.GetRapidity(bHt_tWb_BQuark_p4) - mcTools.GetRapidity(bHt_tWb_WBoson_Antiquark_p4) );

  double dEta_bqbar = abs(bHt_tWb_BQuark_p4.eta() - bHt_tWb_WBoson_Antiquark_p4.eta());
  double dEta_bq    = abs(bHt_tWb_BQuark_p4.eta() - bHt_tWb_WBoson_Quark_p4.eta());
  double dEta_qq    = abs(bHt_tWb_WBoson_Antiquark_p4.eta() - bHt_tWb_WBoson_Quark_p4.eta());
  double dPhi_bqbar = abs( ROOT::Math::VectorUtil::DeltaPhi( bHt_tWb_BQuark_p4, bHt_tWb_WBoson_Antiquark_p4) );
  double dPhi_bq    = abs( ROOT::Math::VectorUtil::DeltaPhi( bHt_tWb_BQuark_p4, bHt_tWb_WBoson_Quark_p4) );
  double dPhi_qq    = abs( ROOT::Math::VectorUtil::DeltaPhi( bHt_tWb_WBoson_Antiquark_p4, bHt_tWb_WBoson_Quark_p4) );


  math::XYZTLorentzVector Tau1_p4;
  math::XYZTLorentzVector Tau2_p4;
  if (bHt_HToHW_HBoson_Tau_p4.pt() > bHt_HToHW_HBoson_Antitau_p4.pt())
    {
      Tau1_p4 = bHt_HToHW_HBoson_Tau_p4;
      Tau2_p4 = bHt_HToHW_HBoson_Antitau_p4;
    }
  else
    {
      Tau1_p4 = bHt_HToHW_HBoson_Antitau_p4;
      Tau2_p4 = bHt_HToHW_HBoson_Tau_p4;
    }
  // cout << "Tau1_p4.pt() = " << Tau1_p4.pt() << ", Tau2_p4.pt() = " << Tau2_p4.pt() << endl;
  // cout << "bHt_HToHW_HBoson_Tau_p4.pt() = " << bHt_HToHW_HBoson_Tau_p4.pt() << ", bHt_HToHW_HBoson_Antitau_p4.pt() = " << bHt_HToHW_HBoson_Antitau_p4.pt() << endl;

  double Taus_dEta     = abs( Tau1_p4.eta() - Tau2_p4.eta());
  double Taus_dPhi     = abs(ROOT::Math::VectorUtil::DeltaPhi( Tau1_p4, Tau2_p4 ) );
  double Taus_dR       = ROOT::Math::VectorUtil::DeltaR( Tau1_p4, Tau2_p4 );
  double Tau1_MET_dPhi = abs( auxTools.DeltaPhi( Tau1_p4.phi(), fEvent.genMET().Phi()) );
  double Tau2_MET_dPhi = abs( auxTools.DeltaPhi( Tau2_p4.phi(), fEvent.genMET().Phi()) );
  double Tau1_bHt_tWb_BQuark_dR   = ROOT::Math::VectorUtil::DeltaR( Tau1_p4, bHt_tWb_BQuark_p4);
  double Tau1_bHt_tWb_BQuark_dEta = abs(Tau1_p4.eta() - bHt_tWb_BQuark_p4.eta());
  double Tau1_bHt_tWb_BQuark_dPhi = abs(ROOT::Math::VectorUtil::DeltaPhi( Tau1_p4, bHt_tWb_BQuark_p4));
  double Tau2_bHt_tWb_BQuark_dR   = ROOT::Math::VectorUtil::DeltaR( Tau2_p4, bHt_tWb_BQuark_p4);
  double Tau2_bHt_tWb_BQuark_dEta = abs(Tau2_p4.eta() - bHt_tWb_BQuark_p4.eta());
  double Tau2_bHt_tWb_BQuark_dPhi = abs(ROOT::Math::VectorUtil::DeltaPhi( Tau2_p4, bHt_tWb_BQuark_p4));
  double Tau1_bHt_tWb_WBoson_Quark_dR       = ROOT::Math::VectorUtil::DeltaR( Tau1_p4, bHt_tWb_WBoson_Quark_p4);
  double Tau1_bHt_tWb_WBoson_Quark_dEta     = abs(Tau1_p4.eta() - bHt_tWb_WBoson_Quark_p4.eta());
  double Tau1_bHt_tWb_WBoson_Quark_dPhi     = abs(ROOT::Math::VectorUtil::DeltaPhi( Tau1_p4, bHt_tWb_WBoson_Quark_p4) );
  double Tau1_bHt_tWb_WBoson_Antiquark_dR   = ROOT::Math::VectorUtil::DeltaR( Tau1_p4, bHt_tWb_WBoson_Antiquark_p4);
  double Tau1_bHt_tWb_WBoson_Antiquark_dEta = abs(Tau1_p4.eta() - bHt_tWb_WBoson_Antiquark_p4.eta() );
  double Tau1_bHt_tWb_WBoson_Antiquark_dPhi = abs(ROOT::Math::VectorUtil::DeltaPhi( Tau1_p4, bHt_tWb_WBoson_Antiquark_p4) );
  double Tau2_bHt_tWb_WBoson_Quark_dR       = ROOT::Math::VectorUtil::DeltaR( Tau2_p4, bHt_tWb_WBoson_Quark_p4);
  double Tau2_bHt_tWb_WBoson_Quark_dEta     = abs(Tau2_p4.eta() - bHt_tWb_WBoson_Quark_p4.eta());
  double Tau2_bHt_tWb_WBoson_Quark_dPhi     = abs(ROOT::Math::VectorUtil::DeltaPhi( Tau2_p4, bHt_tWb_WBoson_Quark_p4) );
  double Tau2_bHt_tWb_WBoson_Antiquark_dR   = ROOT::Math::VectorUtil::DeltaR( Tau2_p4, bHt_tWb_WBoson_Antiquark_p4);
  double Tau2_bHt_tWb_WBoson_Antiquark_dEta = abs(Tau2_p4.eta() - bHt_tWb_WBoson_Antiquark_p4.eta());
  double Tau2_bHt_tWb_WBoson_Antiquark_dPhi = abs(ROOT::Math::VectorUtil::DeltaPhi( Tau2_p4, bHt_tWb_WBoson_Antiquark_p4) );
  double Tau1_bHt_tWb_WBoson_Leptons_dR     = ROOT::Math::VectorUtil::DeltaR( Tau1_p4, bHt_tWb_WBoson_Lepton_p4);
  double Tau1_bHt_tWb_WBoson_Leptons_dEta   = abs(Tau1_p4.eta() - bHt_tWb_WBoson_Lepton_p4.eta());
  double Tau1_bHt_tWb_WBoson_Leptons_dPhi   = abs(ROOT::Math::VectorUtil::DeltaPhi( Tau1_p4, bHt_tWb_WBoson_Lepton_p4) );
  double Tau2_bHt_tWb_WBoson_Leptons_dR     = ROOT::Math::VectorUtil::DeltaR( Tau2_p4, bHt_tWb_WBoson_Lepton_p4);
  double Tau2_bHt_tWb_WBoson_Leptons_dEta   = abs(Tau2_p4.eta() - bHt_tWb_WBoson_Lepton_p4.eta());
  double Tau2_bHt_tWb_WBoson_Leptons_dPhi   = abs(ROOT::Math::VectorUtil::DeltaPhi( Tau2_p4, bHt_tWb_WBoson_Lepton_p4) );

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Fill TH1
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (cfg_Verbose) std::cout << "=== Fill TH1" << std::endl;
  if (selJets_p4.size() > 3) 
    {
      double jet1_Eta = selJets_p4.at(0).eta();
      double jet2_Eta = selJets_p4.at(1).eta();
      double jet3_Eta = selJets_p4.at(2).eta();
      double jet4_Eta = selJets_p4.at(3).eta();

      double jet1_Phi = selJets_p4.at(0).phi();
      double jet2_Phi = selJets_p4.at(1).phi();
      double jet3_Phi = selJets_p4.at(2).phi();
      double jet4_Phi = selJets_p4.at(3).phi();
      
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
  h_genHT_GenParticles->Fill(genP_HT);
  h_genHT_GenJets     ->Fill(genJ_HT);
  h_genMET_Et         ->Fill(fEvent.genMET().et()); 
  h_genMET_Phi        ->Fill(fEvent.genMET().Phi());
  h_GenJets_N         ->Fill(selectedJets.size());
  h_MaxDiJetMass_Mass ->Fill( maxDijetMass_mass     );
  h_MaxDiJetMass_Pt   ->Fill( maxDijetMass_p4.pt()  );
  h_MaxDiJetMass_Eta  ->Fill( maxDijetMass_p4.eta() );
  h_MaxDiJetMass_Rap  ->Fill( maxDijetMass_rapidity );
  h_MaxDiJetMass_dR   ->Fill( maxDijetMass_dR       );
  h_MaxDiJetMass_dRrap->Fill( maxDijetMass_dRrap    );
  h_MaxDiJetMass_dEta ->Fill( maxDijetMass_dEta     );
  h_MaxDiJetMass_dPhi ->Fill( maxDijetMass_dPhi     );
  h_MaxDiJetMass_dRap ->Fill( maxDijetMass_dRap     );

  h_bHt_HToHW_HBoson_Taus_dR  ->Fill(Taus_dR);
  h_bHt_HToHW_HBoson_Taus_dEta->Fill(Taus_dEta);
  h_bHt_HToHW_HBoson_Taus_dPhi->Fill(Taus_dPhi);
  h_bHt_HToHW_HBoson_Tau1_MET_dPhi->Fill(Tau1_MET_dPhi);
  h_bHt_HToHW_HBoson_Tau2_MET_dPhi->Fill(Tau2_MET_dPhi);										
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_BQuark_dR->Fill(Tau1_bHt_tWb_BQuark_dR);
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_BQuark_dEta->Fill(Tau1_bHt_tWb_BQuark_dEta);
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_BQuark_dPhi->Fill(Tau1_bHt_tWb_BQuark_dPhi);
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_BQuark_dR->Fill(Tau2_bHt_tWb_BQuark_dR);
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_BQuark_dEta->Fill(Tau2_bHt_tWb_BQuark_dEta);
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_BQuark_dPhi->Fill(Tau2_bHt_tWb_BQuark_dPhi);
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Quark_dR->Fill(Tau1_bHt_tWb_WBoson_Quark_dR);
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Quark_dEta->Fill(Tau1_bHt_tWb_WBoson_Quark_dEta);
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Quark_dPhi->Fill(Tau1_bHt_tWb_WBoson_Quark_dPhi);
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Quark_dR->Fill(Tau2_bHt_tWb_WBoson_Quark_dR);
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Quark_dEta->Fill(Tau2_bHt_tWb_WBoson_Quark_dEta);
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Quark_dPhi->Fill(Tau2_bHt_tWb_WBoson_Quark_dPhi);
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Antiquark_dR->Fill(Tau1_bHt_tWb_WBoson_Antiquark_dR);
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Antiquark_dEta->Fill(Tau1_bHt_tWb_WBoson_Antiquark_dEta);
  h_bHt_HToHW_HBoson_Tau1_bHt_tWb_WBoson_Antiquark_dPhi->Fill(Tau1_bHt_tWb_WBoson_Antiquark_dPhi);
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Antiquark_dR->Fill(Tau2_bHt_tWb_WBoson_Antiquark_dR);
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Antiquark_dEta->Fill(Tau2_bHt_tWb_WBoson_Antiquark_dEta);
  h_bHt_HToHW_HBoson_Tau2_bHt_tWb_WBoson_Antiquark_dPhi->Fill(Tau2_bHt_tWb_WBoson_Antiquark_dPhi);
  h_bHt_HToHW_HBoson_Tau1_bHt_HToHW_WBoson_Leptons_dR->Fill(Tau1_bHt_tWb_WBoson_Leptons_dR);
  h_bHt_HToHW_HBoson_Tau1_bHt_HToHW_WBoson_Leptons_dEta->Fill(Tau1_bHt_tWb_WBoson_Leptons_dEta);
  h_bHt_HToHW_HBoson_Tau1_bHt_HToHW_WBoson_Leptons_dPhi->Fill(Tau1_bHt_tWb_WBoson_Leptons_dPhi);
  h_bHt_HToHW_HBoson_Tau2_bHt_HToHW_WBoson_Leptons_dR->Fill(Tau2_bHt_tWb_WBoson_Leptons_dR);
  h_bHt_HToHW_HBoson_Tau2_bHt_HToHW_WBoson_Leptons_dEta->Fill(Tau2_bHt_tWb_WBoson_Leptons_dEta);
  h_bHt_HToHW_HBoson_Tau2_bHt_HToHW_WBoson_Leptons_dPhi->Fill(Tau2_bHt_tWb_WBoson_Leptons_dPhi);

  h_bHt_Hpm_Pt         ->Fill( bHt_Hpm_p4.pt() );
  h_bHt_Hpm_Eta        ->Fill( bHt_Hpm_p4.eta() );
  h_bHt_Hpm_Rap        ->Fill( mcTools.GetRapidity(bHt_Hpm_p4 ) );
  h_bHt_TQuark_Pt      ->Fill( bHt_TQuark_p4.pt()  );
  h_bHt_TQuark_Eta     ->Fill( bHt_TQuark_p4.eta() );
  h_bHt_TQuark_Rap     ->Fill( mcTools.GetRapidity(bHt_TQuark_p4) );
  h_bHt_tWb_WBoson_Pt  ->Fill( bHt_tWb_WBoson_p4.pt()  );
  h_bHt_tWb_WBoson_Eta ->Fill( bHt_tWb_WBoson_p4.eta() );
  h_bHt_tWb_WBoson_Rap ->Fill( mcTools.GetRapidity(bHt_tWb_WBoson_p4 ) );
  h_bHt_tWb_BQuark_Pt  ->Fill( bHt_tWb_BQuark_p4.pt()  ); 
  h_bHt_tWb_BQuark_Eta ->Fill( bHt_tWb_BQuark_p4.eta() );
  h_bHt_tWb_BQuark_Rap ->Fill( mcTools.GetRapidity(bHt_tWb_BQuark_p4) );

  if ( (nGenP_tWb_WNeutrinos > 0) && (nGenP_tWb_WLeptons > 0) )
    {
      h_bHt_tWb_WBoson_Leptons_Pt   ->Fill( bHt_tWb_WBoson_Lepton_p4.pt()  );
      h_bHt_tWb_WBoson_Leptons_Eta  ->Fill( bHt_tWb_WBoson_Lepton_p4.eta() );
      h_bHt_tWb_WBoson_Leptons_Rap  ->Fill( mcTools.GetRapidity(bHt_tWb_WBoson_Lepton_p4) );
      h_bHt_tWb_WBoson_Neutrinos_Pt ->Fill( bHt_tWb_WBoson_Neutrino_p4.pt()  );
      h_bHt_tWb_WBoson_Neutrinos_Eta->Fill( bHt_tWb_WBoson_Neutrino_p4.eta() );
      h_bHt_tWb_WBoson_Neutrinos_Rap->Fill( mcTools.GetRapidity(bHt_tWb_WBoson_Neutrino_p4) );
    }
  else
    {
      h_bHt_tWb_WBoson_Quark_Pt     ->Fill( bHt_tWb_WBoson_Quark_p4.pt()  );
      h_bHt_tWb_WBoson_Quark_Eta    ->Fill( bHt_tWb_WBoson_Quark_p4.eta() );
      h_bHt_tWb_WBoson_Quark_Rap    ->Fill( mcTools.GetRapidity(bHt_tWb_WBoson_Quark_p4) );
      h_bHt_tWb_WBoson_Antiquark_Pt ->Fill( bHt_tWb_WBoson_Antiquark_p4.pt()  );
      h_bHt_tWb_WBoson_Antiquark_Eta->Fill( bHt_tWb_WBoson_Antiquark_p4.eta() );
      h_bHt_tWb_WBoson_Antiquark_Rap->Fill( mcTools.GetRapidity(bHt_tWb_WBoson_Antiquark_p4) );

      h_bHt_tWb_bqq_Pt        ->Fill( bHt_tWb_bqq_p4.pt() );
      h_bHt_tWb_bqq_Rap       ->Fill( mcTools.GetRapidity(bHt_tWb_bqq_p4) );
      h_bHt_tWb_bqq_Mass      ->Fill( bHt_tWb_bqq_p4.mass() );
      h_bHt_tWb_bqq_dRMax_dR  ->Fill( deltaRMax );
      h_bHt_tWb_bqq_dRMax_dRap->Fill( bqq_dRap );
      h_bHt_tWb_bqq_dRMax_dPhi->Fill( bqq_dPhi );      
    }

  h_bHt_HToHW_HBoson_Pt          ->Fill( bHt_HToHW_HBoson_p4.pt()  );
  h_bHt_HToHW_HBoson_Eta         ->Fill( bHt_HToHW_HBoson_p4.eta() );
  h_bHt_HToHW_HBoson_Rap         ->Fill( mcTools.GetRapidity(bHt_HToHW_HBoson_p4) );
  // Tau1_p4.phi()  // xenios
  h_bHt_HToHW_HBoson_Tau_Pt      ->Fill( bHt_HToHW_HBoson_Tau_p4.pt()  );
  h_bHt_HToHW_HBoson_Tau_Eta     ->Fill( bHt_HToHW_HBoson_Tau_p4.eta() );
  h_bHt_HToHW_HBoson_Tau_Rap     ->Fill( mcTools.GetRapidity(bHt_HToHW_HBoson_Tau_p4) );
  h_bHt_HToHW_HBoson_Antitau_Pt  ->Fill( bHt_HToHW_HBoson_Antitau_p4.pt()  );
  h_bHt_HToHW_HBoson_Antitau_Eta ->Fill( bHt_HToHW_HBoson_Antitau_p4.eta() );
  h_bHt_HToHW_HBoson_Antitau_Rap ->Fill( mcTools.GetRapidity(bHt_HToHW_HBoson_Antitau_p4) );
  h_bHt_HToHW_WBoson_Pt          ->Fill( bHt_HToHW_WBoson_p4.pt()  );
  h_bHt_HToHW_WBoson_Eta         ->Fill( bHt_HToHW_WBoson_p4.eta() );
  h_bHt_HToHW_WBoson_Rap         ->Fill( mcTools.GetRapidity(bHt_HToHW_WBoson_p4) );
  if ( (nGenP_HToHW_WNeutrinos > 0) && (nGenP_HToHW_WLeptons > 0) )
    {
      h_bHt_HToHW_WBoson_Leptons_Pt   ->Fill( bHt_HToHW_WBoson_Lepton_p4.pt() );
      h_bHt_HToHW_WBoson_Leptons_Eta  ->Fill( bHt_HToHW_WBoson_Lepton_p4.eta() );
      h_bHt_HToHW_WBoson_Leptons_Rap  ->Fill( mcTools.GetRapidity(bHt_HToHW_WBoson_Lepton_p4) );
      h_bHt_HToHW_WBoson_Neutrinos_Pt ->Fill( bHt_HToHW_WBoson_Neutrino_p4.pt() );
      h_bHt_HToHW_WBoson_Neutrinos_Eta->Fill( bHt_HToHW_WBoson_Neutrino_p4.eta() );
      h_bHt_HToHW_WBoson_Neutrinos_Rap->Fill( mcTools.GetRapidity(bHt_HToHW_WBoson_Neutrino_p4) );
    }
  else
    {
      h_bHt_HToHW_WBoson_Quark_Pt     ->Fill( bHt_HToHW_WBoson_Quark_p4.pt()  );
      h_bHt_HToHW_WBoson_Quark_Eta    ->Fill( bHt_HToHW_WBoson_Quark_p4.eta() );
      h_bHt_HToHW_WBoson_Quark_Rap    ->Fill( mcTools.GetRapidity(bHt_HToHW_WBoson_Quark_p4) );
      h_bHt_HToHW_WBoson_Antiquark_Pt ->Fill( bHt_HToHW_WBoson_Antiquark_p4.pt()  );
      h_bHt_HToHW_WBoson_Antiquark_Eta->Fill( bHt_HToHW_WBoson_Antiquark_p4.eta() );
      h_bHt_HToHW_WBoson_Antiquark_Rap->Fill( mcTools.GetRapidity(bHt_HToHW_WBoson_Antiquark_p4) );
    }
  h_bHt_Hpm_bHt_tWb_WBoson_dR         ->Fill(dR_bHt_Hpm_bHt_tWb_WBoson);
  h_bHt_Hpm_bHt_TQuark_dR             ->Fill(dR_bHt_Hpm_bHt_TQuark);
  h_bHt_HToHW_WBoson_bHt_tWb_WBoson_dR->Fill(dR_bHt_HToHW_WBoson_bHt_tWb_WBoson); 
  h_bHt_tWb_WBoson_bHt_tWb_BQuark_dR  ->Fill(dEta_bHt_tWb_WBoson_bHt_tWb_BQuark);
  h_bHt_TQuark_bHt_tWb_BQuark_dR      ->Fill(dR_bHt_TQuark_bHt_tWb_BQuark);

  if ( (nGenP_tWb_WNeutrinos > 0) && (nGenP_tWb_WLeptons > 0) ) {}
  else
    {
      h_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark_dEta     ->Fill(dEta_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark);
      h_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark_dEta ->Fill(dEta_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark);
      h_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark_dPhi     ->Fill(dPhi_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark);
      h_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark_dPhi ->Fill(dPhi_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark);
      h_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark_dR       ->Fill(dR_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark);
      h_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark_dR   ->Fill(dR_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark);
      h_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark_dRap     ->Fill(dRap_bHt_tWb_BQuark_bHt_tWb_WBoson_Quark);
      h_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark_dRap ->Fill(dRap_bHt_tWb_BQuark_bHt_tWb_WBoson_Antiquark);
    }

  h_bHt_Hpm_bHt_tWb_WBoson_dEta         ->Fill(dEta_bHt_Hpm_bHt_tWb_WBoson);
  h_bHt_Hpm_bHt_TQuark_dEta             ->Fill(dEta_bHt_Hpm_bHt_TQuark);
  h_bHt_HToHW_WBoson_bHt_tWb_WBoson_dEta->Fill(dEta_bHt_HToHW_WBoson_bHt_tWb_WBoson); 
  h_bHt_tWb_WBoson_bHt_tWb_BQuark_dEta  ->Fill(dEta_bHt_tWb_WBoson_bHt_tWb_BQuark);
  h_bHt_TQuark_bHt_tWb_BQuark_dEta      ->Fill(dEta_bHt_TQuark_bHt_tWb_BQuark);
  h_bHt_Hpm_bHt_tWb_WBoson_dPhi         ->Fill(dPhi_bHt_Hpm_bHt_tWb_WBoson);
  h_bHt_Hpm_bHt_TQuark_dPhi             ->Fill(dPhi_bHt_Hpm_bHt_TQuark);
  h_bHt_TQuark_bHt_tWb_BQuark_dPhi      ->Fill(dPhi_bHt_TQuark_bHt_tWb_BQuark);
  h_bHt_HToHW_WBoson_bHt_tWb_WBoson_dPhi->Fill(dPhi_bHt_HToHW_WBoson_bHt_tWb_WBoson); 
  h_bHt_tWb_WBoson_bHt_tWb_BQuark_dPhi  ->Fill(dPhi_bHt_tWb_WBoson_bHt_tWb_BQuark);
  h_bHt_Hpm_bHt_tWb_WBoson_dRap         ->Fill(dRap_bHt_Hpm_bHt_tWb_WBoson);
  h_bHt_Hpm_bHt_TQuark_dRap             ->Fill(dRap_bHt_Hpm_bHt_TQuark);
  h_bHt_HToHW_WBoson_bHt_tWb_WBoson_dRap->Fill(dRap_bHt_HToHW_WBoson_bHt_tWb_WBoson); 
  h_bHt_TQuark_bHt_tWb_BQuark_dRap      ->Fill(dRap_bHt_TQuark_bHt_tWb_BQuark);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Fill TH2
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (cfg_Verbose) std::cout << "=== Fill TH2" << std::endl;
  h_bHt_tWb_bqq_dRMax_dRap_Vs_dPhi ->Fill( bqq_dRap, bqq_dPhi );

  h_bHt_HToHW_HBoson_Taus_Pt_Vs_Pt     ->Fill( Tau1_p4.pt(), Tau2_p4.pt());
  h_bHt_HToHW_HBoson_Taus_Eta_Vs_Eta   ->Fill( Tau1_p4.eta(), Tau2_p4.eta());
  h_bHt_HToHW_HBoson_Taus_Phi_Vs_Phi   ->Fill( Tau1_p4.phi(), Tau2_p4.phi());
  h_bHt_HToHW_HBoson_Taus_dEta_Vs_dPhi ->Fill( Taus_dEta, Taus_dPhi);
  h_bHt_HToHW_HBoson_Taus_Pt1_Vs_dR    ->Fill( Tau1_p4.pt(), Taus_dR);
  h_bHt_HToHW_HBoson_Taus_Pt2_Vs_dR    ->Fill( Tau2_p4.pt(), Taus_dR);

  h_MaxDiJetMass_dEta_Vs_dPhi->Fill( maxDijetMass_dEta, maxDijetMass_dPhi );
  h_MaxDiJetMass_dRap_Vs_dPhi->Fill( maxDijetMass_dRap, maxDijetMass_dPhi );
 
  h_Bquark_Antiquark_dEta_Vs_dPhi->Fill( dEta_bqbar, dPhi_bqbar );
  h_Bquark_Quark_dEta_Vs_dPhi    ->Fill( dEta_bq   , dPhi_bq    );
  h_Quark_Antiquark_dEta_Vs_dPhi ->Fill( dEta_qq   , dPhi_qq    );

  return;
}



vector<float> KinematicsHToHW::GetMomentumTensorEigenValues(std::vector<math::XYZTLorentzVector> jets,
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


vector<float> KinematicsHToHW::GetMomentumTensorEigenValues2D(std::vector<math::XYZTLorentzVector> jets,
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



vector<float> KinematicsHToHW::GetSphericityTensorEigenValues(std::vector<math::XYZTLorentzVector> jets,
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

TMatrixDSym KinematicsHToHW::ComputeMomentumTensor(std::vector<math::XYZTLorentzVector> jets, double r)
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



TMatrixDSym KinematicsHToHW::ComputeMomentumTensor2D(std::vector<math::XYZTLorentzVector> jets)
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


double KinematicsHToHW::GetAlphaT(std::vector<math::XYZTLorentzVector> jets,
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

vector<GenJet> KinematicsHToHW::GetGenJets(const vector<GenJet> genJets, std::vector<float> ptCuts, std::vector<float> etaCuts, vector<genParticle> genParticlesToMatch)
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
	  
      // For-loop: All Generated Jets
      for (auto jet: genJets) 
	{
	  
	  // Jet index (for pT and eta cuts)
	  jet_index++;
	  
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
	  if (0) std::cout << "dR = " << dR << ": dPt = " << dPt << ", dEta = " << dEta << ", dPhi = " << dPhi << std::endl;

	  // Increment cut index only. Cannot be bigger than the size of the cut list provided
	  if (ptCut_index  < ptCuts.size()-1  ) ptCut_index++;
	  if (etaCut_index < etaCuts.size()-1  ) etaCut_index++;
	  break;
	}
    }
  if (0) std::cout << "bjets.size() = " << jets.size() << std::endl;
  return jets;
}

vector<GenJet> KinematicsHToHW::GetGenJets(const GenJetCollection& genJets, std::vector<float> ptCuts, std::vector<float> etaCuts)
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


vector<genParticle> KinematicsHToHW::GetAllPreviousCopies(const vector<genParticle> genParticles, genParticle genP)
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
	      continue;// iro - fixme
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


vector<genParticle> KinematicsHToHW::GetGenParticles(const vector<genParticle> genParticles, double ptCut, double etaCut, const int pdgId, const bool isLastCopy, const bool hasNoDaughters)
  {
  
  std::vector<genParticle> particles;

  // For-loop: All genParticles
  for (auto& p: genParticles) 
    {

      // std::cout << "p.pdgId() = " << p.pdgId() << std::endl;

      // Find last copy of a given particle
      // if (isLastCopy) if (!p.isLastCopy()) continue; // crashes
      if (isLastCopy)
	{
	  // fixme
	  if (abs(p.pdgId()) == 5){ if (p.status() != 23) continue;} // Consider only status=23 (outgoing) particles
	  else{ if (p.status() != 1 and p.status() != 2) continue;}
	}
      
      // Commonly enables for parton-based jet flavour definition
      if (hasNoDaughters) if (p.daughters().size() > 0) continue;

      // Consider only particles
      if (std::abs(p.pdgId()) != pdgId) continue;

      // Apply cuts
      if ( p.pt() < ptCut ) continue;
      if (std::abs(p.eta()) > etaCut) continue;
      
      // Save this particle
      particles.push_back(p);
    }

  // std::cout << "Returning " << particles.size() << " particles." << std::endl;
  return particles;
  }
