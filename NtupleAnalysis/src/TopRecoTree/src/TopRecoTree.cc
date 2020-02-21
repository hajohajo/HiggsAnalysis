// -*- c++ -*-
#include "Framework/interface/BaseSelector.h"
#include "Framework/interface/makeTH.h"

#include "EventSelection/interface/CommonPlots.h"
#include "Tools/interface/DirectionalCut.h"
#include "EventSelection/interface/EventSelections.h"
#include "EventSelection/interface/BJetSelection.h"
#include "Tools/interface/MCTools.h"
#include "Auxiliary/interface/Tools.h"
#include "Auxiliary/interface/Table.h"
#include "Tools/interface/DirectionalCut.h"

#include "TDirectory.h"
#include "Math/VectorUtil.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorT.h"
#include "TVectorD.h"

#include "TTree.h"
#include "TBranch.h"
#include <TLorentzVector.h>

#include "../work/TopRecoTree.h"

struct PtComparator
{
  bool operator() (const genParticle p1, const genParticle p2) const { return ( p1.pt() > p2.pt() ); }
  bool operator() (const math::XYZTLorentzVector p1, const math::XYZTLorentzVector p2) const { return ( p1.pt() > p2.pt() ); }
};

struct TrijetSelections{
  std::vector<Jet> Jet1;
  std::vector<Jet> Jet2;
  std::vector<Jet> BJet;
  std::vector <double> MVA;
  std::vector<math::XYZTLorentzVector> TrijetP4;
  std::vector<math::XYZTLorentzVector> DijetP4;
};


class TopRecoTree: public BaseSelector {
public:
  explicit TopRecoTree(const ParameterSet& config, const TH1* skimCounters);
  virtual ~TopRecoTree() {}

  Jet getLeadingSubleadingJet(const Jet& jet0, const Jet& jet1, string selectedJet);
  std::vector<int> SortInPt(std::vector<int> Vector);
  //Vector sorting according to the pt. - Descending order
  std::vector<math::XYZTLorentzVector> SortInPt(std::vector<math::XYZTLorentzVector> Vector);
  //returns the last copy of a gen particle
  genParticle findLastCopy(int index);
  //is Bjet
  bool isBJet(const Jet& jet, const std::vector<Jet>& bjets);
  bool isMatchedJet(const Jet& jet, const std::vector<Jet>& jets);
  
  bool isWsubjet(const Jet& jet, const std::vector<Jet>& jets1, const std::vector<Jet>& jets2);
  bool isLepton(const Jet& jet, const std::vector<Electron>& selectedElectrons, const std::vector<Muon>& selectedMuons);
  ///Are same Jets
  bool areSameJets(const Jet& jet1, const Jet& jet2);
  /// Books histograms
  virtual void book(TDirectory *dir ) override;
  /// Sets up branches for reading the TTree
  virtual void setupBranches(BranchManager& branchManager) override;
  /// Called for each event
  virtual void process(Long64_t entry) override;
  void getTopDecayProducts(const Event& fEvent, genParticle top, vector<genParticle> &quarks, vector<genParticle> &bquarks);
  void getTopDecayProducts(const Event& fEvent, genParticle top, genParticle &quark1, genParticle &quark2, genParticle &bquark); 
  const genParticle GetLastCopy(const vector<genParticle> genParticles, const genParticle &p);
  vector<genParticle> GetGenParticles(const vector<genParticle> genParticles, const int pdgId);


private:
  // Input parameters
  const std::string cfg_SelectionsType;
  const DirectionalCut<double> cfg_MiniIsoCut;
  const DirectionalCut<double> cfg_METCut;
  const DirectionalCut<double> cfg_LepBJetDRCut;
  const float cfg_DeltaPtOverPtCut;
  const float cfg_DeltaRCut;

  const float cfg_MuonPtCut;
  const float cfg_MuonEtaCut;
  const float cfg_ElePtCut;
  const float cfg_EleEtaCut;

  const HistogramSettings cfg_PtBinSetting;
  const HistogramSettings cfg_EtaBinSetting;
  const HistogramSettings cfg_PhiBinSetting;
  const HistogramSettings cfg_MassBinSetting;
  const HistogramSettings cfg_DeltaEtaBinSetting;
  const HistogramSettings cfg_DeltaPhiBinSetting;
  const HistogramSettings cfg_DeltaRBinSetting;
  

  Tools auxTools;


  // Common plots
  CommonPlots fCommonPlots;
  // Event selection classes and event counters (in same order like they are applied)
  Count              cAllEvents;
  Count              cTrigger;
  METFilterSelection fMETFilterSelection;
  Count              cVertexSelection;
  ElectronSelection  fElectronSelection;
  MuonSelection      fMuonSelection;
  TauSelection       fTauSelection;
  JetSelection       fJetSelection;
  BJetSelection      fBJetSelection;
  Count              cBTaggingSFCounter;
  METSelection       fMETSelection;
  Count              cSelected;

  // Histograms  
// Trijet Distributions
  WrappedTH1Triplet *hTrijetPt;
  WrappedTH1Triplet *hTrijetEta;
  WrappedTH1Triplet *hTrijetPhi;
  WrappedTH1Triplet *hTrijetMass;
  WrappedTH1Triplet *hTrijetPtDr;
  //WrappedTH1Triplet *hTrijetPFCharge;
  WrappedTH1Triplet *hTrijetCombinedCvsL;
  //WrappedTH1Triplet *hTrijetDeepCvsL;
  WrappedTH1Triplet *hTrijetPtD;
  //WrappedTH1Triplet *hTrijetAxis1;
  WrappedTH1Triplet *hTrijetAxis2;
  WrappedTH1Triplet *hTrijetMult;
  WrappedTH1Triplet *hTrijetQGLikelihood;
  WrappedTH1Triplet *hTrijetQGLikelihood_avg;
  WrappedTH1Triplet *hTrijetChiSquared;
  
  // Dijet distributions
  WrappedTH1Triplet *hDijetPt;
  WrappedTH1Triplet *hDijetEta;
  WrappedTH1Triplet *hDijetPhi;
  WrappedTH1Triplet *hDijetMass;
  WrappedTH1Triplet *hDijetPtDr;
  //WrappedTH1Triplet *hDijetPFCharge;
  WrappedTH1Triplet *hDijetCombinedCvsL;
  //WrappedTH1Triplet *hDijetDeepCvsL;
  WrappedTH1Triplet *hDijetPtD;
  //WrappedTH1Triplet *hDijetAxis1;
  WrappedTH1Triplet *hDijetAxis2;
  WrappedTH1Triplet *hDijetMult;
  WrappedTH1Triplet *hDijetQGLikelihood;
  WrappedTH1Triplet *hDijetQGLikelihood_avg;
  WrappedTH1Triplet *hDijetMassOverTrijetMass;
  WrappedTH1Triplet *hDijetChiSquared;
  
  // Leading jet from dijet distributions
  WrappedTH1Triplet *hLdgJetPt;
  WrappedTH1Triplet *hLdgJetEta;
  WrappedTH1Triplet *hLdgJetPhi;
  WrappedTH1Triplet *hLdgJetMass;
  //WrappedTH1Triplet *hLdgJetPFCharge;
  WrappedTH1Triplet *hLdgJetCombinedCvsL;
  //WrappedTH1Triplet *hLdgJetDeepCvsL;
  WrappedTH1Triplet *hLdgJetBdisc;
  WrappedTH1Triplet *hLdgJetPtD;
  //WrappedTH1Triplet *hLdgJetAxis1;
  WrappedTH1Triplet *hLdgJetAxis2;
  WrappedTH1Triplet *hLdgJetMult;
  WrappedTH1Triplet *hLdgJetQGLikelihood;
  //WrappedTH1Triplet *hLdgJetPullMagnitude;
  
  // Subleading jet from dijet distributions
  WrappedTH1Triplet *hSubldgJetPt;
  WrappedTH1Triplet *hSubldgJetEta;
  WrappedTH1Triplet *hSubldgJetPhi;
  WrappedTH1Triplet *hSubldgJetMass;
  //WrappedTH1Triplet *hSubldgJetPFCharge;
  WrappedTH1Triplet *hSubldgJetCombinedCvsL;
  //WrappedTH1Triplet *hSubldgJetDeepCvsL;
  WrappedTH1Triplet *hSubldgJetBdisc;
  WrappedTH1Triplet *hSubldgJetPtD;
  //WrappedTH1Triplet *hSubldgJetAxis1;
  WrappedTH1Triplet *hSubldgJetAxis2;
  WrappedTH1Triplet *hSubldgJetMult;
  WrappedTH1Triplet *hSubldgJetQGLikelihood;
  //WrappedTH1Triplet *hSubldgJetPullMagnitude;
  
  // b-jet distributions
  WrappedTH1Triplet *hBJetPt;
  WrappedTH1Triplet *hBJetEta;
  WrappedTH1Triplet *hBJetPhi;
  WrappedTH1Triplet *hBJetMass;
  //WrappedTH1Triplet *hBJetPFCharge;
  WrappedTH1Triplet *hBJetCombinedCvsL;
  //WrappedTH1Triplet *hBJetDeepCvsL;
  WrappedTH1Triplet *hBJetBdisc;
  WrappedTH1Triplet *hBJetPtD;
  //WrappedTH1Triplet *hBJetAxis1;
  WrappedTH1Triplet *hBJetAxis2;
  WrappedTH1Triplet *hBJetMult;
  WrappedTH1Triplet *hBJetQGLikelihood;
  
  // b-jet + (sub)leading jet from dijet system
  WrappedTH1Triplet *hBJetLdgJet_Mass;
  //WrappedTH1Triplet *hBJetLdgJet_PFCharge;
  WrappedTH1Triplet *hBJetSubldgJet_Mass;
  //WrappedTH1Triplet *hBJetSubldgJet_PFCharge;
  
  WrappedTH1Triplet *hSoftDrop_n2;
  //WrappedTH1Triplet *hPullAngleJ1J2;
  //WrappedTH1Triplet *hPullAngleJ2J1;
  
  // Distances
  //DEta
  WrappedTH1Triplet *hDEtaJ1withJ2;
  WrappedTH1Triplet *hDEtaJ1withBJet;
  WrappedTH1Triplet *hDEtaJ2withBJet;
  WrappedTH1Triplet *hDEtaDijetwithBJet;
  WrappedTH1Triplet *hDEtaJ1BJetwithJ2;
  WrappedTH1Triplet *hDEtaJ2BJetwithJ1;
  //DPhi
  WrappedTH1Triplet *hDPhiJ1withJ2;
  WrappedTH1Triplet *hDPhiJ1withBJet;
  WrappedTH1Triplet *hDPhiJ2withBJet;
  WrappedTH1Triplet *hDPhiDijetwithBJet;
  WrappedTH1Triplet *hDPhiJ1BJetwithJ2;
  WrappedTH1Triplet *hDPhiJ2BJetwithJ1;
  //DR
  WrappedTH1Triplet *hDRJ1withJ2;
  WrappedTH1Triplet *hDRJ1withBJet;
  WrappedTH1Triplet *hDRJ2withBJet;
  WrappedTH1Triplet *hDRDijetwithBJet;
  WrappedTH1Triplet *hDRJ1BJetwithJ2;
  WrappedTH1Triplet *hDRJ2BJetwithJ1;
  //
  WrappedTH1Triplet *hDRJ1withElectron;
  WrappedTH1Triplet *hDRJ2withElectron;
  WrappedTH1Triplet *hDRBwithElectron;
  WrappedTH1Triplet *hDRJ1withElectron_ptD0p80;
  WrappedTH1Triplet *hDRJ2withElectron_ptD0p80;
  WrappedTH1Triplet *hDRBwithElectron_ptD0p80;

  WrappedTH1Triplet *hDRJ1withMuon;
  WrappedTH1Triplet *hDRJ2withMuon;
  WrappedTH1Triplet *hDRBwithMuon;
  WrappedTH1Triplet *hDRJ1withMuon_ptD0p80;
  WrappedTH1Triplet *hDRJ2withMuon_ptD0p80;
  WrappedTH1Triplet *hDRBwithMuon_ptD0p80;

  //Ctrl plots for semi-leptonic selections
  WrappedTH1 *hCtrl_JetMult;
  WrappedTH1 *hCtrl_JetPt;
  WrappedTH1 *hCtrl_BJetMult;
  WrappedTH1 *hCtrl_muonMult;
  WrappedTH1 *hCtrl_muonIso;
  WrappedTH1 *hCtrl_eleMult;
  WrappedTH1 *hCtrl_eleIso;
  WrappedTH1 *hCtrl_BJetPt;
  WrappedTH1 *hCtrl_MET;
  WrappedTH1 *hCtrl_HT;
  WrappedTH1 *hCtrl_NhadronicTops;
  WrappedTH1 *hCtrl_hadrGenTopPt;
  WrappedTH1 *hCtrl_lepGenTopPt;

  // All jets 
  WrappedTH1 *hAllJetCvsL;
  WrappedTH1 *hAllJetPtD;
  WrappedTH1 *hAllJetAxis2;
  WrappedTH1 *hAllJetMult;
  WrappedTH1 *hAllJetBdisc;
  WrappedTH1 *hAllJetQGLikelihood;
  // CJets
  WrappedTH1 *hCJetCvsL;
  WrappedTH1 *hCJetPtD;
  WrappedTH1 *hCJetAxis2;
  WrappedTH1 *hCJetMult;
  WrappedTH1 *hCJetBdisc;

  //Matching check
  WrappedTH1        *hNmatchedTop;
  WrappedTH1        *hNmatchedTrijets;
  WrappedTH1Triplet *hGenTop_Pt;
  WrappedTH1Triplet *hGenQuark_Pt;
  WrappedTH1        *hTrijetDrMin;
  WrappedTH1        *hTrijetDPtOverGenPt;
  WrappedTH1        *hTrijetDEtaOverGenEta;
  WrappedTH1        *hTrijetDPhiOverGenPhi;
  WrappedTH1        *hTrijetDPt_matched;
  WrappedTH1        *hTrijetDEta_matched;
  WrappedTH1        *hTrijetDPhi_matched;
  WrappedTH1        *hJetsDeltaRmin;
  
  // Dpt/pt (quark pt bins)
  WrappedTH1        *hQuarkJetMinDr03_DeltaPtOverPt;
  WrappedTH1        *hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt0To40GeV;
  WrappedTH1        *hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt40To60GeV;
  WrappedTH1        *hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt60To80GeV;
  WrappedTH1        *hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt80To100GeV;
  WrappedTH1        *hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt100To120GeV;
  WrappedTH1        *hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt120To140GeV;
  WrappedTH1        *hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt140To160GeV;
  WrappedTH1        *hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt160To180GeV;
  WrappedTH1        *hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt180To200GeV;
  WrappedTH1        *hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt200ToInfGeV;

  //Correlations
  WrappedTH2 *hQuarkJetMinDr03_DeltaPtOverPt_vs_QuarkPt;
  WrappedTH2 *hQuarkJetMinDr03_DeltaR_vs_QuarkPt;
  WrappedTH2 *hQuarkJetMinDr03_DeltaPtOverPt_vs_DeltaRmin;
  
  //next WrappedTH1Triplet

  // TTree - TBranches       
  TTree *treeS;
  TTree *treeB;
  //TBranch

  // === SIGNAL
  //Variables from Analysis Note AN-16-437 (branch names will be renamed - kept here temporary)
  TBranch *TrijetPtDR_S;
  TBranch *TrijetDijetPtDR_S;
  TBranch *TrijetBjetMass_S;
  TBranch *TrijetLdgJetBDisc_S;
  TBranch *TrijetSubldgJetBDisc_S;
  TBranch *TrijetBJetLdgJetMass_S;
  TBranch *TrijetBJetSubldgJetMass_S;
  TBranch *TrijetDijetMass_S;
  TBranch *TrijetBJetBDisc_S;
  TBranch *TrijetMass_S;
  TBranch *TrijetSoftDrop_n2_S;
  TBranch *TrijetLdgJetCvsL_S;
  TBranch *TrijetSubldgJetCvsL_S;
  TBranch *TrijetLdgJetPtD_S;
  TBranch *TrijetSubldgJetPtD_S;
  TBranch *TrijetLdgJetAxis2_S;
  TBranch *TrijetSubldgJetAxis2_S;
  TBranch *TrijetLdgJetMult_S;
  TBranch *TrijetSubldgJetMult_S;

  // === BACKGROUND
  //Variables from Analysis Note AN-16-437 (branch names will be renamed - kept here temporary)
  TBranch *TrijetPtDR_B;
  TBranch *TrijetDijetPtDR_B;
  TBranch *TrijetBjetMass_B;
  TBranch *TrijetLdgJetBDisc_B;
  TBranch *TrijetSubldgJetBDisc_B;
  TBranch *TrijetBJetLdgJetMass_B;
  TBranch *TrijetBJetSubldgJetMass_B;
  TBranch *TrijetDijetMass_B;
  TBranch *TrijetBJetBDisc_B;
  TBranch *TrijetMass_B;
  TBranch *TrijetSoftDrop_n2_B;
  TBranch *TrijetLdgJetCvsL_B;
  TBranch *TrijetSubldgJetCvsL_B;
  TBranch *TrijetLdgJetPtD_B;
  TBranch *TrijetSubldgJetPtD_B;
  TBranch *TrijetLdgJetAxis2_B;
  TBranch *TrijetSubldgJetAxis2_B;
  TBranch *TrijetLdgJetMult_B;
  TBranch *TrijetSubldgJetMult_B;
  
  //New TTree branches
  TBranch *branch_eventWeight_S;
  TBranch *branch_trijetPt_S;
  TBranch *branch_trijetEta_S;
  TBranch *branch_trijetPhi_S;
  TBranch *branch_trijetMass_S;
  TBranch *branch_trijetPtDR_S;
  //TBranch *branch_trijetPFCharge_S;
  TBranch *branch_trijetCombinedCvsL_S;
  //TBranch *branch_trijetDeepCvsL_S;
  TBranch *branch_trijetPtD_S;
  //TBranch *branch_trijetAxis1_S;
  TBranch *branch_trijetAxis2_S;
  TBranch *branch_trijetMult_S;
  TBranch *branch_trijetQGLikelihood_S;
  TBranch *branch_trijetQGLikelihood_avg_S;
  TBranch *branch_trijetChiSquared_S;
  TBranch *branch_dijetPt_S;
  TBranch *branch_dijetEta_S;
  TBranch *branch_dijetPhi_S;
  TBranch *branch_dijetMass_S;
  TBranch *branch_dijetPtDR_S;
  //TBranch *branch_dijetPFCharge_S;
  TBranch *branch_dijetCombinedCvsL_S;
  //TBranch *branch_dijetDeepCvsL_S;
  TBranch *branch_dijetPtD_S;
  //TBranch *branch_dijetAxis1_S;
  TBranch *branch_dijetAxis2_S;
  TBranch *branch_dijetMult_S;
  TBranch *branch_dijetQGLikelihood_S;
  TBranch *branch_dijetQGLikelihood_avg_S;
  TBranch *branch_dijetChiSquared_S;
  TBranch *branch_LdgJetPt_S;
  TBranch *branch_LdgJetEta_S;
  TBranch *branch_LdgJetPhi_S;
  TBranch *branch_LdgJetMass_S;
  //TBranch *branch_LdgJetPFCharge_S;
  TBranch *branch_LdgJetBdisc_S;
  TBranch *branch_LdgJetCombinedCvsL_S;
  //TBranch *branch_LdgJetDeepCvsL_S;
  TBranch *branch_LdgJetPtD_S;
  TBranch *branch_LdgJetAxis2_S;
  //TBranch *branch_LdgJetAxis1_S;
  TBranch *branch_LdgJetMult_S;
  TBranch *branch_LdgJetQGLikelihood_S;
  //TBranch *branch_LdgJetPullMagnitude_S;
  TBranch *branch_SubldgJetPt_S;
  TBranch *branch_SubldgJetEta_S;
  TBranch *branch_SubldgJetPhi_S;
  TBranch *branch_SubldgJetMass_S;
  //TBranch *branch_SubldgJetPFCharge_S;
  TBranch *branch_SubldgJetBdisc_S;
  TBranch *branch_SubldgJetCombinedCvsL_S;
  //TBranch *branch_SubldgJetDeepCvsL_S;
  TBranch *branch_SubldgJetPtD_S;
  TBranch *branch_SubldgJetAxis2_S;
  //TBranch *branch_SubldgJetAxis1_S;
  TBranch *branch_SubldgJetMult_S;
  TBranch *branch_SubldgJetQGLikelihood_S;
  //TBranch *branch_SubldgJetPullMagnitude_S;
  TBranch *branch_bjetPt_S;
  TBranch *branch_bjetEta_S;
  TBranch *branch_bjetPhi_S;
  TBranch *branch_bjetBdisc_S;
  TBranch *branch_bjetMass_S;
  TBranch *branch_bjetQGLikelihood_S;
  TBranch *branch_bjetCombinedCvsL_S;
  //TBranch *branch_bjetDeepCvsL_S;
  TBranch *branch_bjetPtD_S;
  TBranch *branch_bjetAxis2_S;
  //TBranch *branch_bjetAxis1_S;
  TBranch *branch_bjetMult_S;
  //TBranch *branch_bjetPFCharge_S;
  TBranch *branch_bjetLdgJetMass_S;
  TBranch *branch_bjetSubldgJetMass_S;
  TBranch *branch_SoftDrop_n2_S;
  //TBranch *branch_PullAngleJ1J2_S;
  //TBranch *branch_PullAngleJ2J1_S;
  TBranch *branch_DEtaJ1withJ2_S;
  TBranch *branch_DEtaJ1withBJet_S;
  TBranch *branch_DEtaJ2withBJet_S;
  TBranch *branch_DEtaDijetwithBJet_S;
  TBranch *branch_DEtaJ1BJetwithJ2_S;
  TBranch *branch_DEtaJ2BJetwithJ1_S;
  TBranch *branch_DPhiJ1withJ2_S;
  TBranch *branch_DPhiJ1withBJet_S;
  TBranch *branch_DPhiJ2withBJet_S;
  TBranch *branch_DPhiDijetwithBJet_S;
  TBranch *branch_DPhiJ1BJetwithJ2_S;
  TBranch *branch_DPhiJ2BJetwithJ1_S;
  TBranch *branch_DRJ1withJ2_S;
  TBranch *branch_DRJ1withBJet_S;
  TBranch *branch_DRJ2withBJet_S;
  TBranch *branch_DRDijetwithBJet_S;
  TBranch *branch_DRJ1BJetwithJ2_S;
  TBranch *branch_DRJ2BJetwithJ1_S;
  TBranch *branch_dijetMassOverTrijetMass_S;
  
  TBranch *branch_eventWeight_B;
  TBranch *branch_trijetPt_B;
  TBranch *branch_trijetEta_B;
  TBranch *branch_trijetPhi_B;
  TBranch *branch_trijetMass_B;
  TBranch *branch_trijetPtDR_B;
  //TBranch *branch_trijetPFCharge_B;
  TBranch *branch_trijetCombinedCvsL_B;
  //TBranch *branch_trijetDeepCvsL_B;
  TBranch *branch_trijetPtD_B;
  //TBranch *branch_trijetAxis1_B;
  TBranch *branch_trijetAxis2_B;
  TBranch *branch_trijetMult_B;
  TBranch *branch_trijetQGLikelihood_B;
  TBranch *branch_trijetQGLikelihood_avg_B;
  TBranch *branch_trijetChiSquared_B;
  TBranch *branch_dijetPt_B;
  TBranch *branch_dijetEta_B;
  TBranch *branch_dijetPhi_B;
  TBranch *branch_dijetMass_B;
  TBranch *branch_dijetPtDR_B;
  //TBranch *branch_dijetPFCharge_B;
  TBranch *branch_dijetCombinedCvsL_B;
  //TBranch *branch_dijetDeepCvsL_B;
  TBranch *branch_dijetPtD_B;
  //TBranch *branch_dijetAxis1_B;
  TBranch *branch_dijetAxis2_B;
  TBranch *branch_dijetMult_B;
  TBranch *branch_dijetQGLikelihood_B;
  TBranch *branch_dijetQGLikelihood_avg_B;
  TBranch *branch_dijetChiSquared_B;
  TBranch *branch_LdgJetPt_B;
  TBranch *branch_LdgJetEta_B;
  TBranch *branch_LdgJetPhi_B;
  TBranch *branch_LdgJetMass_B;
  //TBranch *branch_LdgJetPFCharge_B;
  TBranch *branch_LdgJetBdisc_B;
  TBranch *branch_LdgJetCombinedCvsL_B;
  //TBranch *branch_LdgJetDeepCvsL_B;
  TBranch *branch_LdgJetPtD_B;
  TBranch *branch_LdgJetAxis2_B;
  //TBranch *branch_LdgJetAxis1_B;
  TBranch *branch_LdgJetMult_B;
  TBranch *branch_LdgJetQGLikelihood_B;
  //TBranch *branch_LdgJetPullMagnitude_B;
  TBranch *branch_SubldgJetPt_B;
  TBranch *branch_SubldgJetEta_B;
  TBranch *branch_SubldgJetPhi_B;
  TBranch *branch_SubldgJetMass_B;
  //TBranch *branch_SubldgJetPFCharge_B;
  TBranch *branch_SubldgJetBdisc_B;
  TBranch *branch_SubldgJetCombinedCvsL_B;
  //TBranch *branch_SubldgJetDeepCvsL_B;
  TBranch *branch_SubldgJetPtD_B;
  TBranch *branch_SubldgJetAxis2_B;
  //TBranch *branch_SubldgJetAxis1_B;
  TBranch *branch_SubldgJetMult_B;
  TBranch *branch_SubldgJetQGLikelihood_B;
  //TBranch *branch_SubldgJetPullMagnitude_B;
  TBranch *branch_bjetPt_B;
  TBranch *branch_bjetEta_B;
  TBranch *branch_bjetPhi_B;
  TBranch *branch_bjetBdisc_B;
  TBranch *branch_bjetMass_B;
  TBranch *branch_bjetQGLikelihood_B;
  TBranch *branch_bjetCombinedCvsL_B;
  //TBranch *branch_bjetDeepCvsL_B;
  TBranch *branch_bjetPtD_B;
  TBranch *branch_bjetAxis2_B;
  //TBranch *branch_bjetAxis1_B;
  TBranch *branch_bjetMult_B;
  //TBranch *branch_bjetPFCharge_B;
  TBranch *branch_bjetLdgJetMass_B;
  TBranch *branch_bjetSubldgJetMass_B;
  TBranch *branch_SoftDrop_n2_B;
  //TBranch *branch_PullAngleJ1J2_B;
  //TBranch *branch_PullAngleJ2J1_B;
  TBranch *branch_DEtaJ1withJ2_B;
  TBranch *branch_DEtaJ1withBJet_B;
  TBranch *branch_DEtaJ2withBJet_B;
  TBranch *branch_DEtaDijetwithBJet_B;
  TBranch *branch_DEtaJ1BJetwithJ2_B;
  TBranch *branch_DEtaJ2BJetwithJ1_B;
  TBranch *branch_DPhiJ1withJ2_B;
  TBranch *branch_DPhiJ1withBJet_B;
  TBranch *branch_DPhiJ2withBJet_B;
  TBranch *branch_DPhiDijetwithBJet_B;
  TBranch *branch_DPhiJ1BJetwithJ2_B;
  TBranch *branch_DPhiJ2BJetwithJ1_B;
  TBranch *branch_DRJ1withJ2_B;
  TBranch *branch_DRJ1withBJet_B;
  TBranch *branch_DRJ2withBJet_B;
  TBranch *branch_DRDijetwithBJet_B;
  TBranch *branch_DRJ1BJetwithJ2_B;
  TBranch *branch_DRJ2BJetwithJ1_B;
  TBranch *branch_dijetMassOverTrijetMass_B;

  
};

#include "Framework/interface/SelectorFactory.h"
REGISTER_SELECTOR(TopRecoTree);

TopRecoTree::TopRecoTree(const ParameterSet& config, const TH1* skimCounters)
  : BaseSelector(config, skimCounters),

    cfg_SelectionsType(config.getParameter<std::string>("RecoTopMVASelection.SelectionsType")),
    cfg_MiniIsoCut(config,  "RecoTopMVASelection.MiniIsoCut"),
    cfg_METCut(config,      "RecoTopMVASelection.METCut"),
    cfg_LepBJetDRCut(config, "RecoTopMVASelection.LepBJetDRCut"),
    cfg_DeltaPtOverPtCut(config.getParameter<float>("RecoTopMVASelection.DeltaPtOverPtCutValue")),
    cfg_DeltaRCut(config.getParameter<float>("RecoTopMVASelection.DeltaRCutValue")),

    // Muon Selection Cuts
    cfg_MuonPtCut(config.getParameter<float>("MuonSelection.muonPtCut")),
    cfg_MuonEtaCut(config.getParameter<float>("MuonSelection.muonEtaCut")),
    cfg_ElePtCut(config.getParameter<float>("ElectronSelection.electronPtCut")),
    cfg_EleEtaCut(config.getParameter<float>("ElectronSelection.electronEtaCut")),

    cfg_PtBinSetting(config.getParameter<ParameterSet>("CommonPlots.ptBins")),
    cfg_EtaBinSetting(config.getParameter<ParameterSet>("CommonPlots.etaBins")),
    cfg_PhiBinSetting(config.getParameter<ParameterSet>("CommonPlots.phiBins")),

    cfg_MassBinSetting(config.getParameter<ParameterSet>("CommonPlots.invMassBins")),
    cfg_DeltaEtaBinSetting(config.getParameter<ParameterSet>("CommonPlots.deltaEtaBins")),
    cfg_DeltaPhiBinSetting(config.getParameter<ParameterSet>("CommonPlots.deltaPhiBins")),
    cfg_DeltaRBinSetting(config.getParameter<ParameterSet>("CommonPlots.deltaRBins")),
    fCommonPlots(config.getParameter<ParameterSet>("CommonPlots"), CommonPlots::kTopReco, fHistoWrapper),
    cAllEvents(fEventCounter.addCounter("all events")),
    cTrigger(fEventCounter.addCounter("passed trigger")),
    fMETFilterSelection(config.getParameter<ParameterSet>("METFilter"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cVertexSelection(fEventCounter.addCounter("passed PV")),
    fElectronSelection(config.getParameter<ParameterSet>("ElectronSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    fMuonSelection(config.getParameter<ParameterSet>("MuonSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    fTauSelection(config.getParameter<ParameterSet>("TauSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, "Veto"),
    fJetSelection(config.getParameter<ParameterSet>("JetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    fBJetSelection(config.getParameter<ParameterSet>("BJetSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    cBTaggingSFCounter(fEventCounter.addCounter("b tag SF")),
    fMETSelection(config.getParameter<ParameterSet>("METSelection"), fEventCounter, fHistoWrapper, &fCommonPlots, ""),
    
    cSelected(fEventCounter.addCounter("Selected Events"))

{ }


void TopRecoTree::book(TDirectory *dir) {
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
  
  // Fixed-binning
  const int nBinsDEta  = cfg_DeltaEtaBinSetting.bins();
  const double minDEta = cfg_DeltaEtaBinSetting.min();
  const double maxDEta = cfg_DeltaEtaBinSetting.max();
  
  const int nBinsDPhi  = cfg_DeltaPhiBinSetting.bins();
  const double minDPhi = cfg_DeltaPhiBinSetting.min();
  const double maxDPhi = cfg_DeltaPhiBinSetting.max();
  
  const int nBinsDR  = cfg_DeltaRBinSetting.bins();
  const double minDR = cfg_DeltaRBinSetting.min();
  const double maxDR = cfg_DeltaRBinSetting.max();
  
  const int nBinsPt    = cfg_PtBinSetting.bins();
  const double minPt   = cfg_PtBinSetting.min();
  const double maxPt   = cfg_PtBinSetting.max();
  
  const int nBinsEta   = cfg_EtaBinSetting.bins();
  const double minEta  = cfg_EtaBinSetting.min();
  const double maxEta  = cfg_EtaBinSetting.max();
 
  const int nBinsPhi   = cfg_PhiBinSetting.bins();
  const double minPhi  = cfg_PhiBinSetting.min();
  const double maxPhi  = cfg_PhiBinSetting.max();
  
  // const int nBinsAxis1 = 100; 
  // const double minAxis1 = 0.0;
  // const double maxAxis1 = 0.40;
  
  const int nBinsAxis2 = 50; 
  const double minAxis2 = 0.0;
  const double maxAxis2 = 0.20; 
  
  const int nBinsMult = 50;
  const double minMult = 0.0;
  const double maxMult = 50.0;
  
  // const int nBinsCharge  = 100;
  // const double minCharge = -2.0;
  // const double maxCharge = +2.0;
  
  const int nBinsQGL = 100;
  const double minQGL = 0.0;
  const double maxQGL = 1.0;
  
  // const int nBinsPullMag = 800;
  // const double minPullMag = 0.0;
  // const double maxPullMag = 0.4;
  
  // const int nBinsPullAngle = 200;
  // const double minPullAngle = -4.0;
  // const double maxPullAngle = +4.0;
  
  const int nNBins     = fCommonPlots.getNjetsBinSettings().bins();
  const float fNMin    = fCommonPlots.getNjetsBinSettings().min();
  const float fNMax    = fCommonPlots.getNjetsBinSettings().max();

  const int nMetBins  = fCommonPlots.getMetBinSettings().bins();
  const float fMetMin = fCommonPlots.getMetBinSettings().min();
  const float fMetMax = 2*fCommonPlots.getMetBinSettings().max();

  const int nHtBins    = fCommonPlots.getHtBinSettings().bins();
  const float fHtMin   = fCommonPlots.getHtBinSettings().min();
  const float fHtMax   = fCommonPlots.getHtBinSettings().max();

  // Create directories for normalization                                                                                                                                                
  std::string myInclusiveLabel  = "TrijetCandidate";
  std::string myFakeLabel       = myInclusiveLabel+"Fake";
  std::string myGenuineLabel    = myInclusiveLabel+"Genuine";
  TDirectory* myInclusiveDir    = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myInclusiveLabel);
  TDirectory* myFakeDir    = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myFakeLabel);
  TDirectory* myGenuineDir = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myGenuineLabel);
  std::vector<TDirectory*> myDirs = {myInclusiveDir, myFakeDir, myGenuineDir};
  
  std::string myInclusiveLabel_2d  = "Scatterplots_";
  std::string myFakeLabel_2d       = myInclusiveLabel_2d+"Fake";
  std::string myGenuineLabel_2d    = myInclusiveLabel_2d+"Genuine";
  TDirectory* myInclusiveDir_2d    = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myInclusiveLabel_2d);
  TDirectory* myFakeDir_2d    = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myFakeLabel_2d);
  TDirectory* myGenuineDir_2d = fHistoWrapper.mkdir(HistoLevel::kSystematics, dir, myGenuineLabel_2d);
  std::vector<TDirectory*> myDirs_2d = {myInclusiveDir_2d, myFakeDir_2d, myGenuineDir_2d};


  // Trijet Distributions
  hTrijetPt           = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "TrijetPt", ";p_{T} (GeV/c)", 2*nBinsPt, minPt , 2*maxPt);
  hTrijetEta          = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "TrijetEta", ";|#eta|", nBinsEta/2, minEta, maxEta);
  hTrijetPhi          = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "TrijetPhi", ";#phi (rads)", nBinsPhi , minPhi , maxPhi );
  hTrijetMass         = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "TrijetMass", ";m_{jjb} (GeV/c^{2})",150,0.0,1500);
  hTrijetPtDr         = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "TrijetPtDr", ";p_{T}#Delta R",150,0.0,1500);
  //hTrijetPFCharge     = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "TrijetPFCharge", "PF-charge", nBinsCharge, minCharge, maxCharge);
  hTrijetCombinedCvsL = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "TrijetAvgCombinedCvsL", ";avg CombinedCvsL discr", 200,-1,1);
  //hTrijetDeepCvsL     = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "TrijetAvgDeepCvsL", ";avg DeepCvsL discr", 200,-1,1);
  hTrijetPtD          = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "TrijetAvgPtD", ";avg p_{T}D",100,0.0,1.0);
  //hTrijetAxis1        = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "TrijetAvgAxis1", ";avg axis1", nBinsAxis1, minAxis1, maxAxis1);
  hTrijetAxis2        = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "TrijetAvgAxis2", ";avg axis2", nBinsAxis2, minAxis2, maxAxis2);
  hTrijetMult         = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "TrijetAvgMult", ";avg mult", nBinsMult, minMult, maxMult);
  hTrijetQGLikelihood = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "TrijetQGLikelihood",";Quark-Gluon Likelihood", nBinsQGL, minQGL, maxQGL);
  hTrijetQGLikelihood_avg = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "TrijetQGLikelihood_avg",";Quark-Gluon Likelihood", nBinsQGL, minQGL, maxQGL);
  hTrijetChiSquared       = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "TrijetChiSquared", ";#chi^{2}", 1000,  0.0, 1000.0);
  
  // Dijet Distributions
  hDijetPt            = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DijetPt", ";p_{T} (GeV/c)", 2*nBinsPt, minPt , 2*maxPt);
  hDijetEta           = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DijetEta", ";|#eta|", nBinsEta/2, minEta, maxEta);
  hDijetPhi           = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DijetPhi", ";#phi (rads)", nBinsPhi , minPhi , maxPhi );
  hDijetMass          = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DijetMass", ";m_{W} (GeV/c^{2})",100,0.0,1000);
  hDijetPtDr          = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DijetPtDr", ";p_{T}#Delta R",150,0.0,1500);
  //hDijetPFCharge      = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DijetPFCharge", "PF-charge", nBinsCharge, minCharge, maxCharge);
  hDijetCombinedCvsL  = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DijetAvgCombinedCvsL", ";avg CombinedCvsL discr", 200,-1,1);
  //hDijetDeepCvsL      = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DijetAvgDeepCvsL", ";avg DeepCvsL discr", 200,-1,1);
  hDijetPtD           = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DijetAvgPtD",   ";avg p_{T}D",100,0.0,1.0);
  //hDijetAxis1         = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DijetAvgAxis1",  ";avg axis1", nBinsAxis1, minAxis1, maxAxis1);
  hDijetAxis2         = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DijetAvgAxis2",  ";avg axis2", nBinsAxis2, minAxis2, maxAxis2);
  hDijetMult          = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DijetAvgMult",   ";avg mult",50,0,50);
  hDijetQGLikelihood  = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DijetQGLikelihood",";Quark-Gluon Likelihood", nBinsQGL, minQGL, maxQGL);
  hDijetQGLikelihood_avg   = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DijetQGLikelihood_avg",";Quark-Gluon Likelihood", nBinsQGL, minQGL, maxQGL);
  hDijetMassOverTrijetMass = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DijetMassOverTrijetMass", ";m_{W}/m_{top}", 100,0.0,1.0);
  hDijetChiSquared         = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DijetChiSquared", ";#chi^{2}", 1000,  0.0, 1000.0);
  
  // Leading Jet Distributions
  hLdgJetPt           = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "LdgJetPt", ";p_{T} (GeV/c)", 2*nBinsPt, minPt , 2*maxPt);
  hLdgJetEta          = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "LdgJetEta", ";|#eta|", nBinsEta/2, minEta, maxEta);
  hLdgJetPhi          = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "LdgJetPhi", ";#phi (rads)", nBinsPhi , minPhi , maxPhi);
  hLdgJetMass         = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "LdgJetMass", ";m_{j} (GeV)", 100,0.0,1000);
  //hLdgJetPFCharge     = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "LdgJetPFCharge", "PF-charge", nBinsCharge, minCharge, maxCharge);
  hLdgJetCombinedCvsL = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "LdgJetCombinedCvsL",";CombinedCvsL discr", 200,-1,1);
  //hLdgJetDeepCvsL     = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "LdgJetDeepCvsL",";DeepCvsL discr", 200,-1,1);
  hLdgJetBdisc        = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "LdgJetBdisc", ";b-tag discr",100,0.0,1.0);
  hLdgJetPtD          = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "LdgJetPtD",";p_{T}D",100,0.0,1.0);
  //hLdgJetAxis1        = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "LdgJetAxis1",";axis1", nBinsAxis1, minAxis1, maxAxis1);
  hLdgJetAxis2        = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "LdgJetAxis2",";axis2", nBinsAxis2, minAxis2, maxAxis2);
  hLdgJetMult         = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "LdgJetMult",";mult", nBinsMult, minMult, maxMult);
  hLdgJetQGLikelihood = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "LdgJetQGLikelihood",";Quark-Gluon Likelihood", nBinsQGL, minQGL, maxQGL);
  //hLdgJetPullMagnitude= fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "LdgJetPullMagnitude", ";pull magnitude", nBinsPullMag, minPullMag, maxPullMag);
  
  // Subleading Jet Distributions
  hSubldgJetPt        = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "SubldgJetPt", ";p_{T} (GeV/c)", 2*nBinsPt, minPt , 2*maxPt);
  hSubldgJetEta       = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "SubldgJetEta",  ";|#eta|", nBinsEta/2, minEta, maxEta);
  hSubldgJetPhi       = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "SubldgJetPhi", ";#phi (rads)", nBinsPhi , minPhi , maxPhi);
  hSubldgJetMass      = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "SubldgJetMass", ";m_{j} (GeV)", 100,0.0,1000);
  //hSubldgJetPFCharge  = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "SubldgJetPFCharge", "PF-charge", nBinsCharge, minCharge, maxCharge);
  hSubldgJetBdisc     = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "SubldgJetBdisc", ";b-tag discr",100,0.0,1.0);
  hSubldgJetCombinedCvsL = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "SubldgJetCombinedCvsL",";CombinedCvsL discr", 200,-1,1);
  //hSubldgJetDeepCvsL  = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "SubldgJetDeepCvsL",";DeepCvsL discr", 200,-1,1);
  hSubldgJetPtD       = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "SubldgJetPtD",";p_{T}D",100,0.0,1.0);
  //hSubldgJetAxis1     = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "SubldgJetAxis1",";axis1", nBinsAxis1, minAxis1, maxAxis1);
  hSubldgJetAxis2     = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "SubldgJetAxis2",";axis2", nBinsAxis2, minAxis2, maxAxis2);
  hSubldgJetMult      = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "SubldgJetMult",";mult", nBinsMult, minMult, maxMult);
  hSubldgJetQGLikelihood  = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "SubldgJetQGLikelihood",";Quark-Gluon Likelihood", nBinsQGL, minQGL, maxQGL);
  // hSubldgJetPullMagnitude = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "SubldgJetPullMagnitude", ";pull magnitude", nBinsPullMag, minPullMag, maxPullMag);
  
  // B-Jet Distributions
  hBJetPt             = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "BJetPt", ";p_{T} (GeV/c)", 2*nBinsPt, minPt , 2*maxPt);
  hBJetEta            = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "BJetEta", ";|#eta|", nBinsEta/2, minEta, maxEta);
  hBJetPhi            = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "BJetPhi", ";#phi (rads)", nBinsPhi , minPhi , maxPhi);
  hBJetBdisc          = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "BJetBdisc", ";b-tag discr",100,0.0,1.0);
  hBJetMass           = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "BJetMass", ";m_{b} (GeV/c^{2})", 120,0,120);
  hBJetLdgJet_Mass    = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "BJetLdgJet_Mass", ";M (GeV/c^{2})",100,0.0,1000);
  //hBJetLdgJet_PFCharge= fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "BJetLdgJet_PFCharge", "PF-charge", nBinsCharge, minCharge, maxCharge);
  hBJetSubldgJet_Mass = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "BJetSubldgJet_Mass", ";M (GeV/c^{2})",100,0.0,1000);
  //hBJetSubldgJet_PFCharge = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "BJetSubldgJet_PFCharge", "PF-charge", nBinsCharge, minCharge, maxCharge);
  hBJetCombinedCvsL   = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "BJetCombinedCvsL",";CombinedCvsL discr", 200,-1,1);
  //hBJetDeepCvsL       = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "BJetDeepCvsL",";DeepCvsL discr", 200,-1,1);
  hBJetPtD            = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "BJetPtD",";p_{T}D",100,0.0,1.0);
  //hBJetAxis1          = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "BJetAxis1",";axis1", nBinsAxis1, minAxis1, maxAxis1);
  hBJetAxis2          = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "BJetAxis2",";axis2", nBinsAxis2, minAxis2, maxAxis2);
  //hBJetPFCharge       = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "BJetPFCharge", ";PF-charge", nBinsCharge, minCharge, maxCharge);
  hBJetMult           = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "BJetMult",";mult", nBinsMult, minMult, maxMult);
  hBJetQGLikelihood   = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "BJetQGLikelihood",";Quark-Gluon Likelihood", nBinsQGL, minQGL, maxQGL);
  
  hSoftDrop_n2        = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "SoftDrop_n2", ";SoftDrop_n2", 50, 0, 2);
  //hPullAngleJ1J2      = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "PullAngleJ1J2", ";#theta(j1, j2)", nBinsPullAngle, minPullAngle, maxPullAngle);
  //hPullAngleJ2J1      = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "PullAngleJ2J1", ";#theta(j2, j1)", nBinsPullAngle, minPullAngle, maxPullAngle);
  
  hDEtaJ1withJ2      = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DEtaJ1withJ2"  , ";#Delta#eta", nBinsDEta, minDEta, maxDEta);
  hDEtaJ1withBJet    = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DEtaJ1withBJet", ";#Delta#eta", nBinsDEta, minDEta, maxDEta);
  hDEtaJ2withBJet    = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DEtaJ2withBJet", ";#Delta#eta", nBinsDEta, minDEta, maxDEta);
  hDEtaDijetwithBJet = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DEtaDijetwithBJet", ";#Delta#eta", nBinsDEta, minDEta, maxDEta);
  hDEtaJ1BJetwithJ2  = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DEtaJ1BJetwithJ2" , ";#Delta#eta", nBinsDEta, minDEta, maxDEta);
  hDEtaJ2BJetwithJ1  = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DEtaJ2BJetwithJ1" , ";#Delta#eta", nBinsDEta, minDEta, maxDEta);
  
  hDPhiJ1withJ2      = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DPhiJ1withJ2"  , ";#Delta#phi", nBinsDPhi, minDPhi, maxDPhi);
  hDPhiJ1withBJet    = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DPhiJ1withBJet", ";#Delta#phi", nBinsDPhi, minDPhi, maxDPhi);
  hDPhiJ2withBJet    = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DPhiJ2withBJet", ";#Delta#phi", nBinsDPhi, minDPhi, maxDPhi);
  hDPhiDijetwithBJet = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DPhiDijetwithBJet", ";#Delta#phi", nBinsDPhi, minDPhi, maxDPhi);
  hDPhiJ1BJetwithJ2  = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DPhiJ1BJetwithJ2" , ";#Delta#phi", nBinsDPhi, minDPhi, maxDPhi);
  hDPhiJ2BJetwithJ1  = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DPhiJ2BJetwithJ1" , ";#Delta#phi", nBinsDPhi, minDPhi, maxDPhi);
  
  hDRJ1withJ2      = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRJ1withJ2"  , ";#Delta R", nBinsDR, minDR, maxDR);
  hDRJ1withBJet    = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRJ1withBJet", ";#Delta R", nBinsDR, minDR, maxDR);
  hDRJ2withBJet    = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRJ2withBJet", ";#Delta R", nBinsDR, minDR, maxDR);
  hDRDijetwithBJet = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRDijetwithBJet", ";#Delta R", nBinsDR, minDR, maxDR);
  hDRJ1BJetwithJ2  = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRJ1BJetwithJ2" , ";#Delta R", nBinsDR, minDR, maxDR);
  hDRJ2BJetwithJ1  = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRJ2BJetwithJ1" , ";#Delta R", nBinsDR, minDR, maxDR);

  hDRJ1withElectron = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRJ1withElectron", ";#Delta R", nBinsDR, minDR, maxDR);
  hDRJ2withElectron = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRJ2withElectron", ";#Delta R", nBinsDR, minDR, maxDR);
  hDRBwithElectron  = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRBwithElectron", ";#Delta R", nBinsDR, minDR, maxDR);
  hDRJ1withElectron_ptD0p80 = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRJ1withElectron_ptD0p80", ";#Delta R", nBinsDR, minDR, maxDR);
  hDRJ2withElectron_ptD0p80 = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRJ2withElectron_ptD0p80", ";#Delta R", nBinsDR, minDR, maxDR);
  hDRBwithElectron_ptD0p80  = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRBwithElectron_ptD0p80", ";#Delta R", nBinsDR, minDR, maxDR);

  hDRJ1withMuon = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRJ1withMuon", ";#Delta R", nBinsDR, minDR, maxDR);
  hDRJ2withMuon = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRJ2withMuon", ";#Delta R", nBinsDR, minDR, maxDR);
  hDRBwithMuon  = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRBwithMuon", ";#Delta R", nBinsDR, minDR, maxDR);
  hDRJ1withMuon_ptD0p80 = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRJ1withMuon_ptD0p80", ";#Delta R", nBinsDR, minDR, maxDR);
  hDRJ2withMuon_ptD0p80 = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRJ2withMuon_ptD0p80", ";#Delta R", nBinsDR, minDR, maxDR);
  hDRBwithMuon_ptD0p80  = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs, "DRBwithMuon_ptD0p80", ";#Delta R", nBinsDR, minDR, maxDR);

  // Ctrl plots for semi-leptonic selections
  hCtrl_JetMult  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug,myInclusiveDir,"Ctrl_JetMult", ";Jet multiplicity", nNBins, fNMin, fNMax);
  hCtrl_JetPt    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug,myInclusiveDir,"Ctrl_JetPt", ";Jet p_{T}",  nBinsPt, minPt , maxPt);
  hCtrl_BJetMult = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug,myInclusiveDir,"Ctrl_BJetMult", ";BJet multiplicity", nNBins, fNMin, fNMax);
  hCtrl_muonMult = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug,myInclusiveDir,"Ctrl_muonMult", ";Muon multiplicity", nNBins, fNMin, fNMax);
  hCtrl_muonIso  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug,myInclusiveDir,"Ctrl_muonIso", ";Muon Isolation", 20, 0, 1);
  hCtrl_eleMult = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug,myInclusiveDir,"Ctrl_eleMult", ";Electron multiplicity", nNBins, fNMin, fNMax);
  hCtrl_eleIso  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug,myInclusiveDir,"Ctrl_eleIso", ";Electron Isolation", 20, 0, 1);
  hCtrl_BJetPt   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug,myInclusiveDir,"Ctrl_BJetPt", ";BJet p_{T}",  nBinsPt, minPt , maxPt);
  hCtrl_MET              = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug,myInclusiveDir,"Ctrl_MET", "Ctrl_MET", nMetBins, fMetMin, fMetMax);
  hCtrl_HT               = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug,myInclusiveDir,"Ctrl_HT", "Ctrl_HT", nHtBins, fHtMin, fHtMax);
  hCtrl_NhadronicTops    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug,myInclusiveDir,"Ctrl_NhadronicTops", "Ctrl_NhadronicTops", 4, -0.5, 3.5);
  hCtrl_hadrGenTopPt     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, myInclusiveDir, "Ctrl_hadrGenTop_Pt", ";p_{T} (GeV/c)", 2*nBinsPt, minPt , 2*maxPt);
  hCtrl_lepGenTopPt     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, myInclusiveDir, "Ctrl_lepGenTop_Pt", ";p_{T} (GeV/c)", 2*nBinsPt, minPt , 2*maxPt);

  // All jets 
  hAllJetCvsL     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, myInclusiveDir, "AllJetCvsL",";CvsL discr", 200,-1,1);
  hAllJetPtD      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, myInclusiveDir, "AllJetPtD",";p_{T}D",100,0.0,1.0);
  hAllJetAxis2    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, myInclusiveDir, "AllJetAxis2",";axis2",50,0,0.2);
  hAllJetMult     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, myInclusiveDir, "AllJetMult",";mult",50,0,50);
  hAllJetBdisc    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, myInclusiveDir,"AllJetBdisc",";b-tag discr",100,0.0,1.0);
  hAllJetQGLikelihood  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, myInclusiveDir, "AllJetQGLikelihood",";Quark-Gluon Likelihood",100,0.0,1.0);
  // CJets
  hCJetCvsL     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, myInclusiveDir, "CJetCvsL",";CvsL discr", 200,-1,1);
  hCJetPtD      = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, myInclusiveDir, "CJetPtD",";p_{T}D",100,0.0,1.0);
  hCJetAxis2    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, myInclusiveDir, "CJetAxis2",";axis2",50,0,0.2);
  hCJetMult     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, myInclusiveDir, "CJetMult",";mult",50,0,50);
  hCJetBdisc    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, myInclusiveDir, "CJetBdisc",";b-tag discr",100,0.0,1.0);

  //Matching check
  hNmatchedTop     = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug,myInclusiveDir,"NmatchedTrijets",";NTrijet_{matched}",4,-0.5,3.5);
  hNmatchedTrijets = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug,myInclusiveDir,"NmatchedTrijetCand",";NTrijet_{matched}",4,-0.5,3.5);
  hGenTop_Pt       = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs,"GenTop_Pt", ";p_{T} (GeV/c)", 2*nBinsPt, minPt , 2*maxPt);
  hGenQuark_Pt     = fHistoWrapper.makeTHTriplet<TH1F>(true, HistoLevel::kVital, myDirs,"GenQuark_Pt", ";p_{T} (GeV/c)", 2*nBinsPt, minPt , 2*maxPt);

  hTrijetDrMin          = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital,myInclusiveDir,"TrijetDrMin",";#Delta R", 50,0.0,0.5);
  hTrijetDPtOverGenPt   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital,myInclusiveDir,"TrijetDPtOverGenPt",";#Delta P_{T}/P_{T,qqb}", 500,-2.5,2.5);
  hTrijetDEtaOverGenEta = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug,myInclusiveDir,"TrijetDEtaOverGenEta",";#Delta #eta/#eta_{qqb}", 200,-1.0,1.0);
  hTrijetDPhiOverGenPhi = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug,myInclusiveDir,"TrijetDPhiOverGenPhi",";#Delta #phi/#phi_{qqb}", 200,-1.0,1.0);
  hTrijetDPt_matched    = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug,myInclusiveDir,"TrijetDPt_matched",";#Delta P_{T}", 2*nBinsPt, -maxPt , maxPt);
  hTrijetDEta_matched   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital,myInclusiveDir,"TrijetDEta_matched",";#Delta #eta", 80,-0.4,0.4);
  hTrijetDPhi_matched   = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital,myInclusiveDir,"TrijetDPhi_matched",";#Delta #phi", 80,-0.4,0.4);

  hJetsDeltaRmin =   fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital,myInclusiveDir,"JetsDeltaRmin",";#Delta R", 50,0.0,0.5);

  hQuarkJetMinDr03_DeltaPtOverPt  = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital,myInclusiveDir,"QuarkJetMinDr03_DeltaPtOverPt","#;Delta P_{T}(jet-quark)/P_{T,q}", 500,-2.5,2.5);

  hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt0To40GeV=fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital,myInclusiveDir,
									    "QuarkJetMinDr03_DeltaPtOverPt_0To40GeV","#;Delta P_{T}(jet-quark)/P_{T,q}",500,-2.5,2.5);
  hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt40To60GeV=fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital,myInclusiveDir,
									     "QuarkJetMinDr03_DeltaPtOverPt_40To60GeV","#;Delta P_{T}(jet-quark)/P_{T,q}",500,-2.5,2.5);
  hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt60To80GeV=fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital,myInclusiveDir,
									     "QuarkJetMinDr03_DeltaPtOverPt_60To80GeV","#;Delta P_{T}(jet-quark)/P_{T,q}",500,-2.5,2.5);
  hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt80To100GeV=fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital,myInclusiveDir,
									      "QuarkJetMinDr03_DeltaPtOverPt_80To100GeV","#;Delta P_{T}(jet-quark)/P_{T,q}",500,-2.5,2.5);
  hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt100To120GeV=fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital,myInclusiveDir,
									       "QuarkJetMinDr03_DeltaPtOverPt_100To120GeV","#;Delta P_{T}(jet-quark)/P_{T,q}",500,-2.5,2.5);
  hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt120To140GeV=fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital,myInclusiveDir,
									       "QuarkJetMinDr03_DeltaPtOverPt_120To140GeV","#;Delta P_{T}(jet-quark)/P_{T,q}",500,-2.5,2.5);
  hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt140To160GeV=fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital,myInclusiveDir,
									       "QuarkJetMinDr03_DeltaPtOverPt_140To160GeV","#;Delta P_{T}(jet-quark)/P_{T,q}",500,-2.5,2.5);
  hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt160To180GeV=fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital,myInclusiveDir,
									       "QuarkJetMinDr03_DeltaPtOverPt_160To180GeV","#;Delta P_{T}(jet-quark)/P_{T,q}",500,-2.5,2.5);
  hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt180To200GeV=fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital,myInclusiveDir,
									       "QuarkJetMinDr03_DeltaPtOverPt_180To200GeV","#;Delta P_{T}(jet-quark)/P_{T,q}",500,-2.5,2.5);
  hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt200ToInfGeV=fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital,myInclusiveDir,
									       "QuarkJetMinDr03_DeltaPtOverPt_2000ToInfGeV","#;Delta P_{T}(jet-quark)/P_{T,q}",500,-2.5,2.5);
  
  //Correlation plots
  hQuarkJetMinDr03_DeltaPtOverPt_vs_QuarkPt = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, myInclusiveDir_2d,
									 "QuarkJetMinDr03_DeltaPtOverPt_vs_QuarkPt", ";#Delta P_{T}/P_{T};P_{T} (GeV/c)", 
									 500,-2.5,2.5, 2*nBinsPt, minPt , 2*maxPt);
  hQuarkJetMinDr03_DeltaR_vs_QuarkPt        = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, myInclusiveDir_2d,
									 "QuarkJetMinDr03_DeltaR_vs_QuarkPt", ";#Delta R;P_{T} (GeV/c)",
									 100,0.0,1.0, 2*nBinsPt, minPt , 2*maxPt);
  hQuarkJetMinDr03_DeltaPtOverPt_vs_DeltaRmin   = fHistoWrapper.makeTH<TH2F>(HistoLevel::kVital, myInclusiveDir_2d,
									     "QuarkJetMinDr03_DeltaPtOverPt_vs_DeltaRmin", ";#Delta P_{T}/P_{T};#Delta R",
									     500,-2.5,2.5, 100,0.0,1.0);


  // TTree
  treeS = new TTree("treeS", "TTree");
  treeB = new TTree("treeB", "TTree");

  //OLD branches 
  TrijetPtDR_S              = treeS -> Branch ("TrijetPtDR",              &TrijetPtDR_S,             "TrijetPtDR_S/F"              );
  TrijetDijetPtDR_S         = treeS -> Branch ("TrijetDijetPtDR",         &TrijetDijetPtDR_S,        "TrijetDijetPtDR_S/F"         );
  TrijetBjetMass_S          = treeS -> Branch ("TrijetBjetMass",          &TrijetBjetMass_S,         "TrijetBjetMass_S/F"          );
  TrijetLdgJetBDisc_S       = treeS -> Branch ("TrijetLdgJetBDisc",       &TrijetLdgJetBDisc_S,      "TrijetLdgJetBDisc_S/F"       );
  TrijetSubldgJetBDisc_S    = treeS -> Branch ("TrijetSubldgJetBDisc",    &TrijetSubldgJetBDisc_S,   "TrijetSubldgJetBDisc_S/F"    );
  TrijetBJetLdgJetMass_S    = treeS -> Branch ("TrijetBJetLdgJetMass",    &TrijetBJetLdgJetMass_S,   "TrijetBJetLdgJetMass_S/F"    );
  TrijetBJetSubldgJetMass_S = treeS -> Branch ("TrijetBJetSubldgJetMass", &TrijetBJetSubldgJetMass_S,"TrijetBJetSubldgJetMass_S/F" );
  TrijetDijetMass_S         = treeS -> Branch ("TrijetDijetMass",         &TrijetDijetMass_S,        "TrijetDijetMass_S/F"         );
  TrijetBJetBDisc_S         = treeS -> Branch ("TrijetBJetBDisc",         &TrijetBJetBDisc_S,        "TrijetBJetBDisc_S/F"         );
  TrijetMass_S              = treeS -> Branch ("TrijetMass",              &TrijetMass_S,             "TrijetMass_S/F"              );
  TrijetSoftDrop_n2_S       = treeS -> Branch ("TrijetSoftDrop_n2",       &TrijetSoftDrop_n2_S,      "TrijetSoftDrop_n2_S/F"       );
  TrijetLdgJetCvsL_S       = treeS -> Branch ("TrijetLdgJetCvsL",         &TrijetLdgJetCvsL_S,       "TrijetLdgJetCvsL_S/F"        );
  TrijetSubldgJetCvsL_S    = treeS -> Branch ("TrijetSubldgJetCvsL",      &TrijetSubldgJetCvsL_S,     "TrijetSubldgJetCvsL_S/F"    );
  TrijetLdgJetPtD_S        = treeS -> Branch ("TrijetLdgJetPtD",          &TrijetLdgJetPtD_S,        "TrijetLdgJetPtD_S/F"         );
  TrijetSubldgJetPtD_S     = treeS -> Branch ("TrijetSubldgJetPtD",       &TrijetSubldgJetPtD_S,     "TrijetSubldgJetPtD_S/F"      );
  TrijetLdgJetAxis2_S      = treeS -> Branch ("TrijetLdgJetAxis2",        &TrijetLdgJetAxis2_S,      "TrijetLdgJetAxis2_S/F"       );
  TrijetSubldgJetAxis2_S   = treeS -> Branch ("TrijetSubldgJetAxis2",     &TrijetSubldgJetAxis2_S,   "TrijetSubldgJetAxis2_S/F"    );
  TrijetLdgJetMult_S       = treeS -> Branch ("TrijetLdgJetMult",         &TrijetLdgJetMult_S,       "TrijetLdgJetMult_S/I"        );
  TrijetSubldgJetMult_S    = treeS -> Branch ("TrijetSubldgJetMult",      &TrijetSubldgJetMult_S,    "TrijetSubldgJetMult_S/I"     );
  //
  TrijetPtDR_B              = treeB -> Branch ("TrijetPtDR",              &TrijetPtDR_B,             "TrijetPtDR_B/F"              );
  TrijetDijetPtDR_B         = treeB -> Branch ("TrijetDijetPtDR",         &TrijetDijetPtDR_B,        "TrijetDijetPtDR_B/F"         );
  TrijetBjetMass_B          = treeB -> Branch ("TrijetBjetMass",          &TrijetBjetMass_B,         "TrijetBjetMass_B/F"          );
  TrijetLdgJetBDisc_B       = treeB -> Branch ("TrijetLdgJetBDisc",       &TrijetLdgJetBDisc_B,      "TrijetLdgJetBDisc_B/F"       );
  TrijetSubldgJetBDisc_B    = treeB -> Branch ("TrijetSubldgJetBDisc",    &TrijetSubldgJetBDisc_B,   "TrijetSubldgJetBDisc_B/F"    );
  TrijetBJetLdgJetMass_B    = treeB -> Branch ("TrijetBJetLdgJetMass",    &TrijetBJetLdgJetMass_B,   "TrijetBJetLdgJetMass_B/F"    );
  TrijetBJetSubldgJetMass_B = treeB -> Branch ("TrijetBJetSubldgJetMass", &TrijetBJetSubldgJetMass_B,"TrijetBJetSubldgJetMass_B/F" );
  TrijetDijetMass_B         = treeB -> Branch ("TrijetDijetMass",         &TrijetDijetMass_B,        "TrijetDijetMass_B/F"         );
  TrijetBJetBDisc_B         = treeB -> Branch ("TrijetBJetBDisc",         &TrijetBJetBDisc_B,        "TrijetBJetBDisc_B/F"         );
  TrijetMass_B              = treeB -> Branch ("TrijetMass",              &TrijetMass_B,             "TrijetMass_B/F"              );
  TrijetSoftDrop_n2_B       = treeB -> Branch ("TrijetSoftDrop_n2",       &TrijetSoftDrop_n2_B,      "TrijetSoftDrop_n2_B/F"       );
  TrijetLdgJetCvsL_B       = treeB -> Branch ("TrijetLdgJetCvsL",         &TrijetLdgJetCvsL_B,       "TrijetLdgJetCvsL_B/F"        );
  TrijetSubldgJetCvsL_B    = treeB -> Branch ("TrijetSubldgJetCvsL",      &TrijetSubldgJetCvsL_B,     "TrijetSubldgJetCvsL_B/F"    );
  TrijetLdgJetPtD_B        = treeB -> Branch ("TrijetLdgJetPtD",          &TrijetLdgJetPtD_B,        "TrijetLdgJetPtD_B/F"         );
  TrijetSubldgJetPtD_B     = treeB -> Branch ("TrijetSubldgJetPtD",       &TrijetSubldgJetPtD_B,     "TrijetSubldgJetPtD_B/F"      );
  TrijetLdgJetAxis2_B      = treeB -> Branch ("TrijetLdgJetAxis2",        &TrijetLdgJetAxis2_B,      "TrijetLdgJetAxis2_B/F"       );
  TrijetSubldgJetAxis2_B   = treeB -> Branch ("TrijetSubldgJetAxis2",     &TrijetSubldgJetAxis2_B,   "TrijetSubldgJetAxis2_B/F"    );
  TrijetLdgJetMult_B       = treeB -> Branch ("TrijetLdgJetMult",         &TrijetLdgJetMult_B,       "TrijetLdgJetMult_B/I"        );
  TrijetSubldgJetMult_B    = treeB -> Branch ("TrijetSubldgJetMult",      &TrijetSubldgJetMult_B,    "TrijetSubldgJetMult_B/I"     );

  //new branch
  branch_eventWeight_S		   = treeS -> Branch("eventWeight", &branch_eventWeight_S, "eventWeight_S/F");
  branch_trijetPt_S		   = treeS -> Branch("trijetPt", &branch_trijetPt_S, "trijetPt_S/F");
  branch_trijetEta_S		   = treeS -> Branch("trijetEta", &branch_trijetEta_S, "trijetEta_S/F");
  branch_trijetPhi_S		   = treeS -> Branch("trijetPhi", &branch_trijetPhi_S, "trijetPhi_S/F");
  branch_trijetMass_S		   = treeS -> Branch("trijetMass", &branch_trijetMass_S, "trijetMass_S/F");
  branch_trijetPtDR_S		   = treeS -> Branch("trijetPtDR", &branch_trijetPtDR_S, "trijetPtDR_S/F");
  //branch_trijetPFCharge_S	   = treeS -> Branch("trijetPFCharge", &branch_trijetPFCharge_S, "trijetPFCharge_S/F");
  branch_trijetCombinedCvsL_S	   = treeS -> Branch("trijetCombinedCvsL", &branch_trijetCombinedCvsL_S, "trijetCombinedCvsL_S/F");
  //branch_trijetDeepCvsL_S	   = treeS -> Branch("trijetDeepCvsL", &branch_trijetDeepCvsL_S, "trijetDeepCvsL_S/F");
  branch_trijetPtD_S		   = treeS -> Branch("trijetPtD", &branch_trijetPtD_S, "trijetPtD_S/F");
  //branch_trijetAxis1_S		   = treeS -> Branch("trijetAxis1", &branch_trijetAxis1_S, "trijetAxis1_S/F");
  branch_trijetAxis2_S		   = treeS -> Branch("trijetAxis2", &branch_trijetAxis2_S, "trijetAxis2_S/F");
  branch_trijetMult_S		   = treeS -> Branch("trijetMult", &branch_trijetMult_S, "trijetMult_S/F");
  branch_trijetQGLikelihood_S	   = treeS -> Branch("trijetQGLikelihood", &branch_trijetQGLikelihood_S, "trijetQGLikelihood_S/F");
  branch_trijetQGLikelihood_avg_S  = treeS -> Branch("trijetQGLikelihood_avg", &branch_trijetQGLikelihood_avg_S, "trijetQGLikelihood_avg_S/F");
  branch_trijetChiSquared_S	   = treeS -> Branch("trijetChiSquared", &branch_trijetChiSquared_S, "trijetChiSquared_S/F");
  branch_dijetPt_S		   = treeS -> Branch("dijetPt", &branch_dijetPt_S, "dijetPt_S/F");
  branch_dijetEta_S		   = treeS -> Branch("dijetEta", &branch_dijetEta_S, "dijetEta_S/F");
  branch_dijetPhi_S		   = treeS -> Branch("dijetPhi", &branch_dijetPhi_S, "dijetPhi_S/F");
  branch_dijetMass_S		   = treeS -> Branch("dijetMass", &branch_dijetMass_S, "dijetMass_S/F");
  branch_dijetPtDR_S		   = treeS -> Branch("dijetPtDR", &branch_dijetPtDR_S, "dijetPtDR_S/F");
  //branch_dijetPFCharge_S	   = treeS -> Branch("dijetPFCharge", &branch_dijetPFCharge_S, "dijetPFCharge_S/F");
  branch_dijetCombinedCvsL_S	   = treeS -> Branch("dijetCombinedCvsL", &branch_dijetCombinedCvsL_S, "dijetCombinedCvsL_S/F");
  //branch_dijetDeepCvsL_S	   = treeS -> Branch("dijetDeepCvsL", &branch_dijetDeepCvsL_S, "dijetDeepCvsL_S/F");
  branch_dijetPtD_S		   = treeS -> Branch("dijetPtD", &branch_dijetPtD_S, "dijetPtD_S/F");
  //branch_dijetAxis1_S		   = treeS -> Branch("dijetAxis1", &branch_dijetAxis1_S, "dijetAxis1_S/F");
  branch_dijetAxis2_S		   = treeS -> Branch("dijetAxis2", &branch_dijetAxis2_S, "dijetAxis2_S/F");
  branch_dijetMult_S		   = treeS -> Branch("dijetMult", &branch_dijetMult_S, "dijetMult_S/F");
  branch_dijetQGLikelihood_S	   = treeS -> Branch("dijetQGLikelihood", &branch_dijetQGLikelihood_S, "dijetQGLikelihood_S/F");
  branch_dijetQGLikelihood_avg_S   = treeS -> Branch("dijetQGLikelihood_avg", &branch_dijetQGLikelihood_avg_S, "dijetQGLikelihood_avg_S/F");
  branch_dijetChiSquared_S	   = treeS -> Branch("dijetChiSquared", &branch_dijetChiSquared_S, "dijetChiSquared_S/F");
  branch_LdgJetPt_S		   = treeS -> Branch("LdgJetPt", &branch_LdgJetPt_S, "LdgJetPt_S/F");
  branch_LdgJetEta_S		   = treeS -> Branch("LdgJetEta", &branch_LdgJetEta_S, "LdgJetEta_S/F");
  branch_LdgJetPhi_S		   = treeS -> Branch("LdgJetPhi", &branch_LdgJetPhi_S, "LdgJetPhi_S/F");
  branch_LdgJetMass_S		   = treeS -> Branch("LdgJetMass", &branch_LdgJetMass_S, "LdgJetMass_S/F");
  //branch_LdgJetPFCharge_S	   = treeS -> Branch("LdgJetPFCharge", &branch_LdgJetPFCharge_S, "LdgJetPFCharge_S/F");
  branch_LdgJetBdisc_S		   = treeS -> Branch("LdgJetBdisc", &branch_LdgJetBdisc_S, "LdgJetBdisc_S/F");
  branch_LdgJetCombinedCvsL_S	   = treeS -> Branch("LdgJetCombinedCvsL", &branch_LdgJetCombinedCvsL_S, "LdgJetCombinedCvsL_S/F");
  //branch_LdgJetDeepCvsL_S	   = treeS -> Branch("LdgJetDeepCvsL", &branch_LdgJetDeepCvsL_S, "LdgJetDeepCvsL_S/F");
  branch_LdgJetPtD_S		   = treeS -> Branch("LdgJetPtD", &branch_LdgJetPtD_S, "LdgJetPtD_S/F");
  branch_LdgJetAxis2_S		   = treeS -> Branch("LdgJetAxis2", &branch_LdgJetAxis2_S, "LdgJetAxis2_S/F");
  //branch_LdgJetAxis1_S		   = treeS -> Branch("LdgJetAxis1", &branch_LdgJetAxis1_S, "LdgJetAxis1_S/F");
  branch_LdgJetMult_S		   = treeS -> Branch("LdgJetMult", &branch_LdgJetMult_S, "LdgJetMult_S/I");
  branch_LdgJetQGLikelihood_S	   = treeS -> Branch("LdgJetQGLikelihood", &branch_LdgJetQGLikelihood_S, "LdgJetQGLikelihood_S/F");
  //branch_LdgJetPullMagnitude_S	   = treeS -> Branch("LdgJetPullMagnitude", &branch_LdgJetPullMagnitude_S, "LdgJetPullMagnitude_S/F");
  branch_SubldgJetPt_S		   = treeS -> Branch("SubldgJetPt", &branch_SubldgJetPt_S, "SubldgJetPt_S/F");
  branch_SubldgJetEta_S		   = treeS -> Branch("SubldgJetEta", &branch_SubldgJetEta_S, "SubldgJetEta_S/F");
  branch_SubldgJetPhi_S		   = treeS -> Branch("SubldgJetPhi", &branch_SubldgJetPhi_S, "SubldgJetPhi_S/F");
  branch_SubldgJetMass_S	   = treeS -> Branch("SubldgJetMass", &branch_SubldgJetMass_S, "SubldgJetMass_S/F");
  //branch_SubldgJetPFCharge_S	   = treeS -> Branch("SubldgJetPFCharge", &branch_SubldgJetPFCharge_S, "SubldgJetPFCharge_S/F");
  branch_SubldgJetBdisc_S	   = treeS -> Branch("SubldgJetBdisc", &branch_SubldgJetBdisc_S, "SubldgJetBdisc_S/F");
  branch_SubldgJetCombinedCvsL_S   = treeS -> Branch("SubldgJetCombinedCvsL", &branch_SubldgJetCombinedCvsL_S, "SubldgJetCombinedCvsL_S/F");
  //branch_SubldgJetDeepCvsL_S	   = treeS -> Branch("SubldgJetDeepCvsL", &branch_SubldgJetDeepCvsL_S, "SubldgJetDeepCvsL_S/F");
  branch_SubldgJetPtD_S		   = treeS -> Branch("SubldgJetPtD", &branch_SubldgJetPtD_S, "SubldgJetPtD_S/F");
  branch_SubldgJetAxis2_S	   = treeS -> Branch("SubldgJetAxis2", &branch_SubldgJetAxis2_S, "SubldgJetAxis2_S/F");
  //branch_SubldgJetAxis1_S	   = treeS -> Branch("SubldgJetAxis1", &branch_SubldgJetAxis1_S, "SubldgJetAxis1_S/F");
  branch_SubldgJetMult_S	   = treeS -> Branch("SubldgJetMult", &branch_SubldgJetMult_S, "SubldgJetMult_S/I");
  branch_SubldgJetQGLikelihood_S   = treeS -> Branch("SubldgJetQGLikelihood", &branch_SubldgJetQGLikelihood_S, "SubldgJetQGLikelihood_S/F");
  //branch_SubldgJetPullMagnitude_S  = treeS -> Branch("SubldgJetPullMagnitude", &branch_SubldgJetPullMagnitude_S, "SubldgJetPullMagnitude_S/F");
  branch_bjetPt_S		   = treeS -> Branch("bjetPt", &branch_bjetPt_S, "bjetPt_S/F");
  branch_bjetEta_S		   = treeS -> Branch("bjetEta", &branch_bjetEta_S, "bjetEta_S/F");
  branch_bjetPhi_S		   = treeS -> Branch("bjetPhi", &branch_bjetPhi_S, "bjetPhi_S/F");
  branch_bjetBdisc_S		   = treeS -> Branch("bjetBdisc", &branch_bjetBdisc_S, "bjetBdisc_S/F");
  branch_bjetMass_S		   = treeS -> Branch("bjetMass", &branch_bjetMass_S, "bjetMass_S/F");
  branch_bjetQGLikelihood_S	   = treeS -> Branch("bjetQGLikelihood", &branch_bjetQGLikelihood_S, "bjetQGLikelihood_S/F");
  branch_bjetCombinedCvsL_S	   = treeS -> Branch("bjetCombinedCvsL", &branch_bjetCombinedCvsL_S, "bjetCombinedCvsL_S/F");
  //branch_bjetDeepCvsL_S		   = treeS -> Branch("bjetDeepCvsL", &branch_bjetDeepCvsL_S, "bjetDeepCvsL_S/F");
  branch_bjetPtD_S		   = treeS -> Branch("bjetPtD", &branch_bjetPtD_S, "bjetPtD_S/F");
  branch_bjetAxis2_S		   = treeS -> Branch("bjetAxis2", &branch_bjetAxis2_S, "bjetAxis2_S/F");
  //branch_bjetAxis1_S		   = treeS -> Branch("bjetAxis1", &branch_bjetAxis1_S, "bjetAxis1_S/F");
  branch_bjetMult_S		   = treeS -> Branch("bjetMult", &branch_bjetMult_S, "bjetMult_S/I");
  //branch_bjetPFCharge_S		   = treeS -> Branch("bjetPFCharge", &branch_bjetPFCharge_S, "bjetPFCharge_S/F");
  branch_bjetLdgJetMass_S          = treeS -> Branch("bjetLdgJetMass", &branch_bjetLdgJetMass_S, "bjetLdgJetMass/F");          
  branch_bjetSubldgJetMass_S       = treeS -> Branch("bjetSubldgJetMass", &branch_bjetSubldgJetMass_S, "bjetSubldgJetMass_S/F");
  branch_SoftDrop_n2_S		   = treeS -> Branch("SoftDrop_n2", &branch_SoftDrop_n2_S, "SoftDrop_n2_S/F");
  //branch_PullAngleJ1J2_S	   = treeS -> Branch("PullAngleJ1J2", &branch_PullAngleJ1J2_S, "PullAngleJ1J2_S/F");
  //branch_PullAngleJ2J1_S	   = treeS -> Branch("PullAngleJ2J1", &branch_PullAngleJ2J1_S, "PullAngleJ2J1_S/F");
  branch_DEtaJ1withJ2_S		   = treeS -> Branch("DEtaJ1withJ2", &branch_DEtaJ1withJ2_S, "DEtaJ1withJ2_S/F");
  branch_DEtaJ1withBJet_S	   = treeS -> Branch("DEtaJ1withBJet", &branch_DEtaJ1withBJet_S, "DEtaJ1withBJet_S/F");
  branch_DEtaJ2withBJet_S	   = treeS -> Branch("DEtaJ2withBJet", &branch_DEtaJ2withBJet_S, "DEtaJ2withBJet_S/F");
  branch_DEtaDijetwithBJet_S	   = treeS -> Branch("DEtaDijetwithBJet", &branch_DEtaDijetwithBJet_S, "DEtaDijetwithBJet_S/F"); 
  branch_DEtaJ1BJetwithJ2_S	   = treeS -> Branch("DEtaJ1BJetwithJ2", &branch_DEtaJ1BJetwithJ2_S, "DEtaJ1BJetwithJ2_S/F");
  branch_DEtaJ2BJetwithJ1_S	   = treeS -> Branch("DEtaJ2BJetwithJ1", &branch_DEtaJ2BJetwithJ1_S, "DEtaJ2BJetwithJ1_S/F");
  branch_DPhiJ1withJ2_S		   = treeS -> Branch("DPhiJ1withJ2", &branch_DPhiJ1withJ2_S, "DPhiJ1withJ2_S/F");
  branch_DPhiJ1withBJet_S	   = treeS -> Branch("DPhiJ1withBJet", &branch_DPhiJ1withBJet_S, "DPhiJ1withBJet_S/F");
  branch_DPhiJ2withBJet_S	   = treeS -> Branch("DPhiJ2withBJet", &branch_DPhiJ2withBJet_S, "DPhiJ2withBJet_S/F");
  branch_DPhiDijetwithBJet_S	   = treeS -> Branch("DPhiDijetwithBJet", &branch_DPhiDijetwithBJet_S, "DPhiDijetwithBJet_S/F");
  branch_DPhiJ1BJetwithJ2_S	   = treeS -> Branch("DPhiJ1BJetwithJ2", &branch_DPhiJ1BJetwithJ2_S, "DPhiJ1BJetwithJ2_S/F");
  branch_DPhiJ2BJetwithJ1_S	   = treeS -> Branch("DPhiJ2BJetwithJ1", &branch_DPhiJ2BJetwithJ1_S, "DPhiJ2BJetwithJ1_S/F");
  branch_DRJ1withJ2_S		   = treeS -> Branch("DRJ1withJ2", &branch_DRJ1withJ2_S, "DRJ1withJ2_S/F");
  branch_DRJ1withBJet_S		   = treeS -> Branch("DRJ1withBJet", &branch_DRJ1withBJet_S, "DRJ1withBJet_S/F");
  branch_DRJ2withBJet_S		   = treeS -> Branch("DRJ2withBJet", &branch_DRJ2withBJet_S, "DRJ2withBJet_S/F");
  branch_DRDijetwithBJet_S	   = treeS -> Branch("DRDijetwithBJet", &branch_DRDijetwithBJet_S, "DRDijetwithBJet_S/F");
  branch_DRJ1BJetwithJ2_S	   = treeS -> Branch("DRJ1BJetwithJ2", &branch_DRJ1BJetwithJ2_S, "DRJ1BJetwithJ2_S/F");
  branch_DRJ2BJetwithJ1_S	   = treeS -> Branch("DRJ2BJetwithJ1", &branch_DRJ2BJetwithJ1_S, "DRJ2BJetwithJ1_S/F");
  branch_dijetMassOverTrijetMass_S = treeS -> Branch("dijetMassOverTrijetMass", &branch_dijetMassOverTrijetMass_S, "dijetMassOverTrijetMass_S/F");
  
  // Background branches
  branch_eventWeight_B		   = treeB -> Branch("eventWeight", &branch_eventWeight_B, "eventWeight_B/F");
  branch_trijetPt_B		   = treeB -> Branch("trijetPt", &branch_trijetPt_B, "trijetPt_B/F");
  branch_trijetEta_B		   = treeB -> Branch("trijetEta", &branch_trijetEta_B, "trijetEta_B/F");
  branch_trijetPhi_B		   = treeB -> Branch("trijetPhi", &branch_trijetPhi_B, "trijetPhi_B/F");
  branch_trijetMass_B		   = treeB -> Branch("trijetMass", &branch_trijetMass_B, "trijetMass_B/F");
  branch_trijetPtDR_B		   = treeB -> Branch("trijetPtDR", &branch_trijetPtDR_B, "trijetPtDR_B/F");
  //branch_trijetPFCharge_B	   = treeB -> Branch("trijetPFCharge", &branch_trijetPFCharge_B, "trijetPFCharge_B/F");
  branch_trijetCombinedCvsL_B	   = treeB -> Branch("trijetCombinedCvsL", &branch_trijetCombinedCvsL_B, "trijetCombinedCvsL_B/F");
  //branch_trijetDeepCvsL_B	   = treeB -> Branch("trijetDeepCvsL", &branch_trijetDeepCvsL_B, "trijetDeepCvsL_B/F");
  branch_trijetPtD_B		   = treeB -> Branch("trijetPtD", &branch_trijetPtD_B, "trijetPtD_B/F");
  //branch_trijetAxis1_B		   = treeB -> Branch("trijetAxis1", &branch_trijetAxis1_B, "trijetAxis1_B/F");
  branch_trijetAxis2_B		   = treeB -> Branch("trijetAxis2", &branch_trijetAxis2_B, "trijetAxis2_B/F");
  branch_trijetMult_B		   = treeB -> Branch("trijetMult", &branch_trijetMult_B, "trijetMult_B/F");
  branch_trijetQGLikelihood_B	   = treeB -> Branch("trijetQGLikelihood", &branch_trijetQGLikelihood_B, "trijetQGLikelihood_B/F");
  branch_trijetQGLikelihood_avg_B  = treeB -> Branch("trijetQGLikelihood_avg", &branch_trijetQGLikelihood_avg_B, "trijetQGLikelihood_avg_B/F");
  branch_trijetChiSquared_B	   = treeB -> Branch("trijetChiSquared", &branch_trijetChiSquared_B, "trijetChiSquared_B/F");
  branch_dijetPt_B		   = treeB -> Branch("dijetPt", &branch_dijetPt_B, "dijetPt_B/F");
  branch_dijetEta_B		   = treeB -> Branch("dijetEta", &branch_dijetEta_B, "dijetEta_B/F");
  branch_dijetPhi_B		   = treeB -> Branch("dijetPhi", &branch_dijetPhi_B, "dijetPhi_B/F");
  branch_dijetMass_B		   = treeB -> Branch("dijetMass", &branch_dijetMass_B, "dijetMass_B/F");
  branch_dijetPtDR_B		   = treeB -> Branch("dijetPtDR", &branch_dijetPtDR_B, "dijetPtDR_B/F");
  //branch_dijetPFCharge_B	   = treeB -> Branch("dijetPFCharge", &branch_dijetPFCharge_B, "dijetPFCharge_B/F");
  branch_dijetCombinedCvsL_B	   = treeB -> Branch("dijetCombinedCvsL", &branch_dijetCombinedCvsL_B, "dijetCombinedCvsL_B/F");
  //branch_dijetDeepCvsL_B	   = treeB -> Branch("dijetDeepCvsL", &branch_dijetDeepCvsL_B, "dijetDeepCvsL_B/F");
  branch_dijetPtD_B		   = treeB -> Branch("dijetPtD", &branch_dijetPtD_B, "dijetPtD_B/F");
  //branch_dijetAxis1_B		   = treeB -> Branch("dijetAxis1", &branch_dijetAxis1_B, "dijetAxis1_B/F");
  branch_dijetAxis2_B		   = treeB -> Branch("dijetAxis2", &branch_dijetAxis2_B, "dijetAxis2_B/F");
  branch_dijetMult_B		   = treeB -> Branch("dijetMult", &branch_dijetMult_B, "dijetMult_B/F");
  branch_dijetQGLikelihood_B	   = treeB -> Branch("dijetQGLikelihood", &branch_dijetQGLikelihood_B, "dijetQGLikelihood_B/F");
  branch_dijetQGLikelihood_avg_B   = treeB -> Branch("dijetQGLikelihood_avg", &branch_dijetQGLikelihood_avg_B, "dijetQGLikelihood_avg_B/F");
  branch_dijetChiSquared_B	   = treeB -> Branch("dijetChiSquared", &branch_dijetChiSquared_B, "dijetChiSquared_B/F");
  branch_LdgJetPt_B		   = treeB -> Branch("LdgJetPt", &branch_LdgJetPt_B, "LdgJetPt_B/F");
  branch_LdgJetEta_B		   = treeB -> Branch("LdgJetEta", &branch_LdgJetEta_B, "LdgJetEta_B/F");
  branch_LdgJetPhi_B		   = treeB -> Branch("LdgJetPhi", &branch_LdgJetPhi_B, "LdgJetPhi_B/F");
  branch_LdgJetMass_B		   = treeB -> Branch("LdgJetMass", &branch_LdgJetMass_B, "LdgJetMass_B/F");
  //branch_LdgJetPFCharge_B	   = treeB -> Branch("LdgJetPFCharge", &branch_LdgJetPFCharge_B, "LdgJetPFCharge_B/F");
  branch_LdgJetBdisc_B		   = treeB -> Branch("LdgJetBdisc", &branch_LdgJetBdisc_B, "LdgJetBdisc_B/F");
  branch_LdgJetCombinedCvsL_B	   = treeB -> Branch("LdgJetCombinedCvsL", &branch_LdgJetCombinedCvsL_B, "LdgJetCombinedCvsL_B/F");
  //branch_LdgJetDeepCvsL_B	   = treeB -> Branch("LdgJetDeepCvsL", &branch_LdgJetDeepCvsL_B, "LdgJetDeepCvsL_B/F");
  branch_LdgJetPtD_B		   = treeB -> Branch("LdgJetPtD", &branch_LdgJetPtD_B, "LdgJetPtD_B/F");
  branch_LdgJetAxis2_B		   = treeB -> Branch("LdgJetAxis2", &branch_LdgJetAxis2_B, "LdgJetAxis2_B/F");
  //branch_LdgJetAxis1_B		   = treeB -> Branch("LdgJetAxis1", &branch_LdgJetAxis1_B, "LdgJetAxis1_B/F");
  branch_LdgJetMult_B		   = treeB -> Branch("LdgJetMult", &branch_LdgJetMult_B, "LdgJetMult_B/I");
  branch_LdgJetQGLikelihood_B	   = treeB -> Branch("LdgJetQGLikelihood", &branch_LdgJetQGLikelihood_B, "LdgJetQGLikelihood_B/F");
  //branch_LdgJetPullMagnitude_B	   = treeB -> Branch("LdgJetPullMagnitude", &branch_LdgJetPullMagnitude_B, "LdgJetPullMagnitude_B/F");
  branch_SubldgJetPt_B		   = treeB -> Branch("SubldgJetPt", &branch_SubldgJetPt_B, "SubldgJetPt_B/F");
  branch_SubldgJetEta_B		   = treeB -> Branch("SubldgJetEta", &branch_SubldgJetEta_B, "SubldgJetEta_B/F");
  branch_SubldgJetPhi_B		   = treeB -> Branch("SubldgJetPhi", &branch_SubldgJetPhi_B, "SubldgJetPhi_B/F");
  branch_SubldgJetMass_B	   = treeB -> Branch("SubldgJetMass", &branch_SubldgJetMass_B, "SubldgJetMass_B/F");
  //branch_SubldgJetPFCharge_B	   = treeB -> Branch("SubldgJetPFCharge", &branch_SubldgJetPFCharge_B, "SubldgJetPFCharge_B/F");
  branch_SubldgJetBdisc_B	   = treeB -> Branch("SubldgJetBdisc", &branch_SubldgJetBdisc_B, "SubldgJetBdisc_B/F");
  branch_SubldgJetCombinedCvsL_B   = treeB -> Branch("SubldgJetCombinedCvsL", &branch_SubldgJetCombinedCvsL_B, "SubldgJetCombinedCvsL_B/F");
  //branch_SubldgJetDeepCvsL_B	   = treeB -> Branch("SubldgJetDeepCvsL", &branch_SubldgJetDeepCvsL_B, "SubldgJetDeepCvsL_B/F");
  branch_SubldgJetPtD_B		   = treeB -> Branch("SubldgJetPtD", &branch_SubldgJetPtD_B, "SubldgJetPtD_B/F");
  branch_SubldgJetAxis2_B	   = treeB -> Branch("SubldgJetAxis2", &branch_SubldgJetAxis2_B, "SubldgJetAxis2_B/F");
  //branch_SubldgJetAxis1_B	   = treeB -> Branch("SubldgJetAxis1", &branch_SubldgJetAxis1_B, "SubldgJetAxis1_B/F");
  branch_SubldgJetMult_B	   = treeB -> Branch("SubldgJetMult", &branch_SubldgJetMult_B, "SubldgJetMult_B/I");
  branch_SubldgJetQGLikelihood_B   = treeB -> Branch("SubldgJetQGLikelihood", &branch_SubldgJetQGLikelihood_B, "SubldgJetQGLikelihood_B/F");
  //branch_SubldgJetPullMagnitude_B  = treeB -> Branch("SubldgJetPullMagnitude", &branch_SubldgJetPullMagnitude_B, "SubldgJetPullMagnitude_B/F");
  branch_bjetPt_B		   = treeB -> Branch("bjetPt", &branch_bjetPt_B, "bjetPt_B/F");
  branch_bjetEta_B		   = treeB -> Branch("bjetEta", &branch_bjetEta_B, "bjetEta_B/F");
  branch_bjetPhi_B		   = treeB -> Branch("bjetPhi", &branch_bjetPhi_B, "bjetPhi_B/F");
  branch_bjetBdisc_B		   = treeB -> Branch("bjetBdisc", &branch_bjetBdisc_B, "bjetBdisc_B/F");
  branch_bjetMass_B		   = treeB -> Branch("bjetMass", &branch_bjetMass_B, "bjetMass_B/F");
  branch_bjetQGLikelihood_B	   = treeB -> Branch("bjetQGLikelihood", &branch_bjetQGLikelihood_B, "bjetQGLikelihood_B/F");
  branch_bjetCombinedCvsL_B	   = treeB -> Branch("bjetCombinedCvsL", &branch_bjetCombinedCvsL_B, "bjetCombinedCvsL_B/F");
  //branch_bjetDeepCvsL_B		   = treeB -> Branch("bjetDeepCvsL", &branch_bjetDeepCvsL_B, "bjetDeepCvsL_B/F");
  branch_bjetPtD_B		   = treeB -> Branch("bjetPtD", &branch_bjetPtD_B, "bjetPtD_B/F");
  branch_bjetAxis2_B		   = treeB -> Branch("bjetAxis2", &branch_bjetAxis2_B, "bjetAxis2_B/F");
  //branch_bjetAxis1_B		   = treeB -> Branch("bjetAxis1", &branch_bjetAxis1_B, "bjetAxis1_B/F");
  branch_bjetMult_B		   = treeB -> Branch("bjetMult", &branch_bjetMult_B, "bjetMult_B/I");
  //branch_bjetPFCharge_B		   = treeB -> Branch("bjetPFCharge", &branch_bjetPFCharge_B, "bjetPFCharge_B/F");
  branch_bjetLdgJetMass_B          = treeB -> Branch("bjetLdgJetMass", &branch_bjetLdgJetMass_B, "bjetLdgJetMass_B/F");
  branch_bjetSubldgJetMass_B       = treeB -> Branch("bjetSubldgJetMass", &branch_bjetSubldgJetMass_B, "bjetSubldgJetMass_B/F");
  branch_SoftDrop_n2_B		   = treeB -> Branch("SoftDrop_n2", &branch_SoftDrop_n2_B, "SoftDrop_n2_B/F");
  //branch_PullAngleJ1J2_B	   = treeB -> Branch("PullAngleJ1J2", &branch_PullAngleJ1J2_B, "PullAngleJ1J2_B/F");
  //branch_PullAngleJ2J1_B	   = treeB -> Branch("PullAngleJ2J1", &branch_PullAngleJ2J1_B, "PullAngleJ2J1_B/F");
  branch_DEtaJ1withJ2_B		   = treeB -> Branch("DEtaJ1withJ2", &branch_DEtaJ1withJ2_B, "DEtaJ1withJ2_B/F");
  branch_DEtaJ1withBJet_B	   = treeB -> Branch("DEtaJ1withBJet", &branch_DEtaJ1withBJet_B, "DEtaJ1withBJet_B/F");
  branch_DEtaJ2withBJet_B	   = treeB -> Branch("DEtaJ2withBJet", &branch_DEtaJ2withBJet_B, "DEtaJ2withBJet_B/F");
  branch_DEtaDijetwithBJet_B	   = treeB -> Branch("DEtaDijetwithBJet", &branch_DEtaDijetwithBJet_B, "DEtaDijetwithBJet_B/F"); 
  branch_DEtaJ1BJetwithJ2_B	   = treeB -> Branch("DEtaJ1BJetwithJ2", &branch_DEtaJ1BJetwithJ2_B, "DEtaJ1BJetwithJ2_B/F");
  branch_DEtaJ2BJetwithJ1_B	   = treeB -> Branch("DEtaJ2BJetwithJ1", &branch_DEtaJ2BJetwithJ1_B, "DEtaJ2BJetwithJ1_B/F");
  branch_DPhiJ1withJ2_B		   = treeB -> Branch("DPhiJ1withJ2", &branch_DPhiJ1withJ2_B, "DPhiJ1withJ2_B/F");
  branch_DPhiJ1withBJet_B	   = treeB -> Branch("DPhiJ1withBJet", &branch_DPhiJ1withBJet_B, "DPhiJ1withBJet_B/F");
  branch_DPhiJ2withBJet_B	   = treeB -> Branch("DPhiJ2withBJet", &branch_DPhiJ2withBJet_B, "DPhiJ2withBJet_B/F");
  branch_DPhiDijetwithBJet_B	   = treeB -> Branch("DPhiDijetwithBJet", &branch_DPhiDijetwithBJet_B, "DPhiDijetwithBJet_B/F");
  branch_DPhiJ1BJetwithJ2_B	   = treeB -> Branch("DPhiJ1BJetwithJ2", &branch_DPhiJ1BJetwithJ2_B, "DPhiJ1BJetwithJ2_B/F");
  branch_DPhiJ2BJetwithJ1_B	   = treeB -> Branch("DPhiJ2BJetwithJ1", &branch_DPhiJ2BJetwithJ1_B, "DPhiJ2BJetwithJ1_B/F");
  branch_DRJ1withJ2_B		   = treeB -> Branch("DRJ1withJ2", &branch_DRJ1withJ2_B, "DRJ1withJ2_B/F");
  branch_DRJ1withBJet_B		   = treeB -> Branch("DRJ1withBJet", &branch_DRJ1withBJet_B, "DRJ1withBJet_B/F");
  branch_DRJ2withBJet_B		   = treeB -> Branch("DRJ2withBJet", &branch_DRJ2withBJet_B, "DRJ2withBJet_B/F");
  branch_DRDijetwithBJet_B	   = treeB -> Branch("DRDijetwithBJet", &branch_DRDijetwithBJet_B, "DRDijetwithBJet_B/F");
  branch_DRJ1BJetwithJ2_B	   = treeB -> Branch("DRJ1BJetwithJ2", &branch_DRJ1BJetwithJ2_B, "DRJ1BJetwithJ2_B/F");
  branch_DRJ2BJetwithJ1_B	   = treeB -> Branch("DRJ2BJetwithJ1", &branch_DRJ2BJetwithJ1_B, "DRJ2BJetwithJ1_B/F");
  branch_dijetMassOverTrijetMass_B = treeB -> Branch("dijetMassOverTrijetMass", &branch_dijetMassOverTrijetMass_B, "dijetMassOverTrijetMass_B/F");

  return;
}


//  Get all gen particles by pdgId                                                                                                                
vector<genParticle> TopRecoTree::GetGenParticles(const vector<genParticle> genParticles, const int pdgId)
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


void TopRecoTree::getTopDecayProducts(const Event& fEvent, genParticle top, vector<genParticle> &quarks, vector<genParticle> &bquarks){
  
  for (size_t i=0; i<top.daughters().size(); i++){
    int dau_index = top.daughters().at(i);
    genParticle dau = fEvent.genparticles().getGenParticles()[dau_index]; //dau: first copy

    // B-Quark                                                                                                                                                                       
    if (std::abs(dau.pdgId()) ==  5) bquarks.push_back(dau);
    
    // W-Boson                                                                                                                                                                       
    if (std::abs(dau.pdgId()) == 24){
      // Get the last copy                                                                                                        
      genParticle W = GetLastCopy(fEvent.genparticles().getGenParticles(), dau);      

      // Find the decay products of W                                                                               
      for (size_t idau=0; idau<W.daughters().size(); idau++){	
	int Wdau_index = W.daughters().at(idau);
	genParticle Wdau = fEvent.genparticles().getGenParticles()[Wdau_index];	
	//std::cout<<"pdgId: "<<Wdau.pdgId()<<std::endl;
	// Consider only quarks as decaying products                                                                    
	if (std::abs(Wdau.pdgId()) > 5) continue;
	quarks.push_back(Wdau);
      }//for (size_t idau=0; idau<W.daughters().size(); idau++)
    }//if (std::abs(dau.pdgId()) == 24)
  }//for (size_t i=0; i<top.daughters().size(); i++)
  
  //sort quarks in pt
  std::sort( quarks.begin(),  quarks.end(),  PtComparator() );
  std::sort( bquarks.begin(), bquarks.end(), PtComparator() );
  //Debug
  if (quarks.size() == 2)
    if (quarks.at(0).pt() < quarks.at(1).pt()) std::cout<<"not sorted in pt?"<<std::endl;
}


void TopRecoTree::getTopDecayProducts(const Event& fEvent, genParticle top, genParticle &quark1, genParticle &quark2, genParticle &bquark){
  vector<genParticle> quarks, bquarks;
  getTopDecayProducts(fEvent, top, quarks, bquarks);
  
  if (quarks.size() > 0)  quark1 = quarks.at(0);
  if (quarks.size() > 1)  quark2 = quarks.at(1);
  if (bquarks.size() > 0) bquark = bquarks.at(0);
}

//  Get the last copy of a particle.
const genParticle TopRecoTree::GetLastCopy(const vector<genParticle> genParticles, const genParticle &p){

  int gen_pdgId = p.pdgId();

  for (size_t i=0; i<p.daughters().size(); i++){

    const genParticle genDau = genParticles[p.daughters().at(i)];
    int genDau_pdgId   = genDau.pdgId();

    if (gen_pdgId == genDau_pdgId)  return GetLastCopy(genParticles, genDau);
  }
  return p;
}

Jet TopRecoTree::getLeadingSubleadingJet(const Jet& jet0, const Jet& jet1, string selectedJet){
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

std::vector<int> TopRecoTree::SortInPt(std::vector<int> Vector)
{
  int size = Vector.size();
  for (int i=0; i<size-1; i++){
    genParticle genPart1 = fEvent.genparticles().getGenParticles()[Vector.at(i)];
    for (int j=i+1;  j<size; j++){
      genParticle genPart2 = fEvent.genparticles().getGenParticles()[Vector.at(j)];
      if (genPart1.pt() > genPart2.pt()) continue;
      int temp = Vector.at(i);
      Vector.at(i) = Vector.at(j);
      Vector.at(j) = temp;
    }
  }
  return Vector;
}
std::vector<math::XYZTLorentzVector> TopRecoTree::SortInPt(std::vector<math::XYZTLorentzVector> Vector)
{
  int size = Vector.size();
  for (int i=0; i<size-1; i++){
    math::XYZTLorentzVector p4_i = Vector.at(i);
    for (int j=i+1;  j<size; j++){
      math::XYZTLorentzVector p4_j = Vector.at(j);
      if (p4_i.pt() > p4_j.pt()) continue;
      Vector.at(i) = p4_j;
      Vector.at(j) = p4_i;
    }
  }
  return Vector;
}


genParticle TopRecoTree::findLastCopy(int index){
  genParticle gen_particle = fEvent.genparticles().getGenParticles()[index];
  int gen_pdgId = gen_particle.pdgId();
  for (size_t i=0; i<gen_particle.daughters().size(); i++){
    
    genParticle genDau = fEvent.genparticles().getGenParticles()[gen_particle.daughters().at(i)];
    int genDau_index   = genDau.index();
    int genDau_pdgId   = genDau.pdgId();
    if (gen_pdgId == genDau_pdgId) return findLastCopy(genDau_index);
  }
  return gen_particle;
}


bool TopRecoTree::isWsubjet(const Jet& jet , const std::vector<Jet>& jets1 , const std::vector<Jet>& jets2){
  return  (isMatchedJet(jet,jets1)||isMatchedJet(jet,jets2));
}

bool TopRecoTree::isLepton(const Jet& jet, const std::vector<Electron>& selectedElectrons, const std::vector<Muon>& selectedMuons){
  bool passEle = (selectedElectrons.size() > 0);
  bool passMu = (selectedMuons.size() > 0);
  float dR = 999.999;
  float dR_match = 0.3;
  
  if (passEle) dR = ROOT::Math::VectorUtil::DeltaR(selectedElectrons.at(0).p4(), jet.p4());
  else if (passMu)  dR = ROOT::Math::VectorUtil::DeltaR(selectedMuons.at(0).p4(), jet.p4());  
  return (dR <= dR_match);
}

  
bool TopRecoTree::isBJet(const Jet& jet, const std::vector<Jet>& bjets) {
  for (auto bjet: bjets)
    {
      if (areSameJets(jet, bjet)) return true;
    }
  return false;
}

bool TopRecoTree::isMatchedJet(const Jet& jet, const std::vector<Jet>& jets) {
  for (auto Jet: jets)
    {
      if (areSameJets(jet, Jet)) return true;
    }
  return false;
}

bool TopRecoTree::areSameJets(const Jet& jet1, const Jet& jet2) {
  float dR = ROOT::Math::VectorUtil::DeltaR(jet1.p4(), jet2.p4());
  float dR_match = 0.1;
  if (dR <= dR_match) return true;
  else return false;
}

void TopRecoTree::setupBranches(BranchManager& branchManager) {
  fEvent.setupBranches(branchManager);
  return;
}


void TopRecoTree::process(Long64_t entry) {

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
  // 4) Trigger SF
  //================================================================================================   
  // if(0) std::cout << "=== MET Trigger SF" << std::endl;
  // const METSelection::Data silentMETData = fMETSelection.silentAnalyze(fEvent, nVertices);
  // if (fEvent.isMC()) {
  //   fEventWeight.multiplyWeight(silentMETData.getMETTriggerSF());
  // }
  // cMetTriggerSFCounter.increment();
  // fCommonPlots.fillControlPlotsAfterMETTriggerScaleFactor(fEvent);
  

  //================================================================================================   
  // 5) Electron Selection
  //================================================================================================   
  if (0) std::cout << "=== Electron veto" << std::endl;
  const ElectronSelection::Data eData = fElectronSelection.analyze(fEvent);
  std::vector<Electron> selectedElectrons;

  if (cfg_SelectionsType == "hadronic"){
    // Electron veto
    if (eData.hasIdentifiedElectrons()) return;
  }
  else{ //semiLeptonic
    
    // For-loop: All electrons
    for (Electron ele: fEvent.electrons()){
      // Apply cut on pt
      if (ele.pt() < cfg_ElePtCut) continue;
      // Apply cut on abs(eta)
      if (std::fabs(ele.eta()) > cfg_EleEtaCut) continue;
      //=== Apply cut on electron ID (MVA)
      //bool passedIDCut   = getEleIdMVADecision(ele);

      //== fixme!
      bool passedIDCut = false;
      double AbsEta = std::abs(ele.eta());
      if (AbsEta<=0.8 && ele.MVA()>=-0.041) passedIDCut = true;
      if (AbsEta>0.8 && AbsEta<1.479 && ele.MVA()>=0.383) passedIDCut =true;
      if (AbsEta>=1.479 && ele.MVA()>=-0.515) passedIDCut =true;
      //=== fixme!
      
      if (!passedIDCut) continue;
      hCtrl_eleIso -> Fill(ele.effAreaMiniIso());
      //=== Apply cut on electron isolation
      bool ElePass_Iso = cfg_MiniIsoCut.passedCut(ele.effAreaMiniIso());
      if (! ElePass_Iso) continue;
      selectedElectrons.push_back(ele);
    }
  } //semiLeptonic
  //================================================================================================
  // 6) Muon Selection
  //================================================================================================
  if (0) std::cout << "=== Muon veto" << std::endl;
  const MuonSelection::Data muData = fMuonSelection.analyze(fEvent);
  std::vector<Muon> selectedMuons;

  if (cfg_SelectionsType == "hadronic"){
    // Muon veto
    if (muData.hasIdentifiedMuons()) return;
  }
  else{ //semiLeptonic

    // For-loop: All muons
    for(Muon muon: fEvent.muons()) {
      // Apply cut on pt
      if (muon.pt() < cfg_MuonPtCut) continue;
      // Apply cut on abs(eta)
      if (std::fabs(muon.eta()) > cfg_MuonEtaCut) continue;
      // Apply cut on muon ID
      if (!muon.muonIDDiscriminator()) continue;
      hCtrl_muonIso -> Fill(muon.effAreaMiniIso());
      //Apply cut on mini isolation
      bool MuPass_Iso       = cfg_MiniIsoCut.passedCut(muon.effAreaMiniIso());
      if (!MuPass_Iso) continue;
      selectedMuons.push_back(muon);
    }
  }

  bool passElectron = (selectedElectrons.size() == 1);
  bool passMuon     = (selectedMuons.size() == 1);
  bool passLepton   = (passElectron && !passMuon) || (!passElectron && passMuon);
  
  //if (passElectron) std::cout<<"Electron: "<<passElectron<<" "<<selectedElectrons.size()<<" "<<selectedElectrons.at(0).effAreaMiniIso()<<std::endl;
  //if (passMuon) std::cout<<"Muons: "<<passMuon<<" "<<selectedMuons.size()<<" "<<selectedMuons.at(0).effAreaMiniIso()<<std::endl;

  //SemiLeptonic: Select exactly one lepton (electron or muon)
  if (cfg_SelectionsType != "hadronic"){
    if (!passLepton) return;
  }

  //================================================================================================   
  // 7) Tau Veto (HToTauNu Orthogonality)
  //================================================================================================   
  if (0) std::cout << "=== Tau-Veto" << std::endl;
  const TauSelection::Data tauData = fTauSelection.analyze(fEvent);
  if (tauData.hasIdentifiedTaus() ) return;

  //================================================================================================
  // 8) Jet selection
  //================================================================================================
  if (0) std::cout << "=== Jet selection" << std::endl;
  const JetSelection::Data jetData = fJetSelection.analyzeWithoutTau(fEvent);
  if (!jetData.passedSelection()) return;

  
  //================================================================================================
  // Standard Selections
  //================================================================================================
  if (0) std::cout << "=== Standard selection" << std::endl;
  //fCommonPlots.fillControlPlotsAfterTopologicalSelections(fEvent, true);  

  //================================================================================================  
  // 9) BJet selection
  //================================================================================================
  if (0) std::cout << "=== BJet selection" << std::endl;
  const BJetSelection::Data bjetData = fBJetSelection.analyze(fEvent, jetData);
  if (!bjetData.passedSelection()) return;
  
  
  //================================================================================================  
  // 10) BJet SF  
  //================================================================================================
  if (0) std::cout << "=== BJet SF" << std::endl;
  if (fEvent.isMC()) 
    {
      fEventWeight.multiplyWeight(bjetData.getBTaggingScaleFactorEventWeight());
    }
  cBTaggingSFCounter.increment();


  //================================================================================================
  // 11) MET selection
  //================================================================================================
  if (0) std::cout << "=== MET selection" << std::endl;
  const METSelection::Data METData = fMETSelection.analyze(fEvent, nVertices);
  //if (!METData.passedSelection()) return;
  
  if (cfg_SelectionsType != "hadronic"){
    bool MetPassCut = cfg_METCut.passedCut(METData.getMET().R());
    if (!MetPassCut) return;
  }

  //================================================================================================
  // All cuts passed
  //================================================================================================
  if (0) std::cout << "=== All cuts passed" << std::endl;
  cSelected.increment();

  //Ctrl plots for semi-leptonic selections
  hCtrl_MET      -> Fill(METData.getMET().R());
  hCtrl_JetMult  -> Fill(jetData.getSelectedJets().size());
  hCtrl_BJetMult -> Fill(bjetData.getSelectedBJets().size());
  if (passMuon)     hCtrl_muonMult -> Fill(selectedMuons.size());
  if (passElectron) hCtrl_eleMult -> Fill(selectedElectrons.size());
  for (auto& jet: jetData.getSelectedJets()){
    hCtrl_JetPt  -> Fill(jet.pt());
  }
  for (auto& bjet: bjetData.getSelectedBJets()){
    hCtrl_BJetPt -> Fill(bjet.pt());
  }

  // TTree Variables  
  //OLD
  TrijetPtDR_S                 -> SetAddress(&trijetPtDR_S);
  TrijetDijetPtDR_S            -> SetAddress(&dijetPtDR_S);
  TrijetBjetMass_S             -> SetAddress(&bjetMass_S);
  TrijetLdgJetBDisc_S          -> SetAddress(&LdgJetBdisc_S);
  TrijetSubldgJetBDisc_S       -> SetAddress(&SubldgJetBdisc_S);
  TrijetBJetLdgJetMass_S       -> SetAddress(&bjetLdgJetMass_S);
  TrijetBJetSubldgJetMass_S    -> SetAddress(&bjetSubldgJetMass_S);
  TrijetDijetMass_S            -> SetAddress(&dijetMass_S);
  TrijetBJetBDisc_S            -> SetAddress(&bjetBdisc_S);
  TrijetMass_S                 -> SetAddress(&trijetMass_S);
  TrijetSoftDrop_n2_S          -> SetAddress(&SoftDrop_n2_S);
  TrijetLdgJetCvsL_S           -> SetAddress(&LdgJetCombinedCvsL_S);
  TrijetSubldgJetCvsL_S        -> SetAddress(&SubldgJetCombinedCvsL_S);
  TrijetLdgJetPtD_S            -> SetAddress(&LdgJetPtD_S);
  TrijetSubldgJetPtD_S         -> SetAddress(&SubldgJetPtD_S);
  TrijetLdgJetAxis2_S          -> SetAddress(&LdgJetAxis2_S);
  TrijetSubldgJetAxis2_S       -> SetAddress(&SubldgJetAxis2_S);
  TrijetLdgJetMult_S           -> SetAddress(&LdgJetMult_S);
  TrijetSubldgJetMult_S        -> SetAddress(&SubldgJetMult_S);
  //
  TrijetPtDR_B                 -> SetAddress(&trijetPtDR_B);
  TrijetDijetPtDR_B            -> SetAddress(&dijetPtDR_B);
  TrijetBjetMass_B             -> SetAddress(&bjetMass_B);
  TrijetLdgJetBDisc_B          -> SetAddress(&LdgJetBdisc_B);
  TrijetSubldgJetBDisc_B       -> SetAddress(&SubldgJetBdisc_B);
  TrijetBJetLdgJetMass_B       -> SetAddress(&bjetLdgJetMass_B);
  TrijetBJetSubldgJetMass_B    -> SetAddress(&bjetSubldgJetMass_B);
  TrijetDijetMass_B            -> SetAddress(&dijetMass_B);
  TrijetBJetBDisc_B            -> SetAddress(&bjetBdisc_B);
  TrijetMass_B                 -> SetAddress(&trijetMass_B);
  TrijetSoftDrop_n2_B          -> SetAddress(&SoftDrop_n2_B);
  TrijetLdgJetCvsL_B           -> SetAddress(&LdgJetCombinedCvsL_B);
  TrijetSubldgJetCvsL_B        -> SetAddress(&SubldgJetCombinedCvsL_B);
  TrijetLdgJetPtD_B            -> SetAddress(&LdgJetPtD_B);
  TrijetSubldgJetPtD_B         -> SetAddress(&SubldgJetPtD_B);
  TrijetLdgJetAxis2_B          -> SetAddress(&LdgJetAxis2_B);
  TrijetSubldgJetAxis2_B       -> SetAddress(&SubldgJetAxis2_B);
  TrijetLdgJetMult_B           -> SetAddress(&LdgJetMult_B);
  TrijetSubldgJetMult_B        -> SetAddress(&SubldgJetMult_B);
  
  // TTree Branches
  branch_eventWeight_S -> SetAddress(&eventWeight_S);
  branch_trijetPt_S    -> SetAddress(&trijetPt_S);
  branch_trijetEta_S   -> SetAddress(&trijetEta_S);
  branch_trijetPhi_S   -> SetAddress(&trijetPhi_S);
  branch_trijetMass_S  -> SetAddress(&trijetMass_S);
  branch_trijetPtDR_S  -> SetAddress(&trijetPtDR_S);
  //branch_trijetPFCharge_S -> SetAddress(&trijetPFCharge_S);
  branch_trijetCombinedCvsL_S -> SetAddress(&trijetCombinedCvsL_S);
  //branch_trijetDeepCvsL_S -> SetAddress(&trijetDeepCvsL_S);
  branch_trijetPtD_S -> SetAddress(&trijetPtD_S);
  //branch_trijetAxis1_S -> SetAddress(&trijetAxis1_S);
  branch_trijetAxis2_S -> SetAddress(&trijetAxis2_S);
  branch_trijetMult_S -> SetAddress(&trijetMult_S);
  branch_trijetQGLikelihood_S -> SetAddress(&trijetQGLikelihood_S);
  branch_trijetQGLikelihood_avg_S -> SetAddress(&trijetQGLikelihood_avg_S);
  branch_trijetChiSquared_S -> SetAddress(&trijetChiSquared_S);
  branch_dijetPt_S -> SetAddress(&dijetPt_S);
  branch_dijetEta_S -> SetAddress(&dijetEta_S);
  branch_dijetPhi_S -> SetAddress(&dijetPhi_S);
  branch_dijetMass_S -> SetAddress(&dijetMass_S);
  branch_dijetPtDR_S -> SetAddress(&dijetPtDR_S);
  //branch_dijetPFCharge_S -> SetAddress(&dijetPFCharge_S);
  branch_dijetCombinedCvsL_S -> SetAddress(&dijetCombinedCvsL_S);
  //branch_dijetDeepCvsL_S -> SetAddress(&dijetDeepCvsL_S);
  branch_dijetPtD_S -> SetAddress(&dijetPtD_S);
  //branch_dijetAxis1_S -> SetAddress(&dijetAxis1_S);
  branch_dijetAxis2_S -> SetAddress(&dijetAxis2_S);
  branch_dijetMult_S -> SetAddress(&dijetMult_S);
  branch_dijetQGLikelihood_S -> SetAddress(&dijetQGLikelihood_S);
  branch_dijetQGLikelihood_avg_S -> SetAddress(&dijetQGLikelihood_avg_S);
  branch_dijetChiSquared_S -> SetAddress(&dijetChiSquared_S);
  branch_LdgJetPt_S -> SetAddress(&LdgJetPt_S);
  branch_LdgJetEta_S -> SetAddress(&LdgJetEta_S);
  branch_LdgJetPhi_S -> SetAddress(&LdgJetPhi_S);
  branch_LdgJetMass_S -> SetAddress(&LdgJetMass_S);
  //branch_LdgJetPFCharge_S -> SetAddress(&LdgJetPFCharge_S);
  branch_LdgJetBdisc_S -> SetAddress(&LdgJetBdisc_S);
  branch_LdgJetCombinedCvsL_S -> SetAddress(&LdgJetCombinedCvsL_S);
  //branch_LdgJetDeepCvsL_S -> SetAddress(&LdgJetDeepCvsL_S);
  branch_LdgJetPtD_S -> SetAddress(&LdgJetPtD_S);
  branch_LdgJetAxis2_S -> SetAddress(&LdgJetAxis2_S);
  //branch_LdgJetAxis1_S -> SetAddress(&LdgJetAxis1_S);
  branch_LdgJetMult_S -> SetAddress(&LdgJetMult_S);
  branch_LdgJetQGLikelihood_S -> SetAddress(&LdgJetQGLikelihood_S);
  //branch_LdgJetPullMagnitude_S -> SetAddress(&LdgJetPullMagnitude_S);
  branch_SubldgJetPt_S -> SetAddress(&SubldgJetPt_S);
  branch_SubldgJetEta_S -> SetAddress(&SubldgJetEta_S);
  branch_SubldgJetPhi_S -> SetAddress(&SubldgJetPhi_S);
  branch_SubldgJetMass_S -> SetAddress(&SubldgJetMass_S);
  //branch_SubldgJetPFCharge_S -> SetAddress(&SubldgJetPFCharge_S);
  branch_SubldgJetBdisc_S -> SetAddress(&SubldgJetBdisc_S);
  branch_SubldgJetCombinedCvsL_S -> SetAddress(&SubldgJetCombinedCvsL_S);
  //branch_SubldgJetDeepCvsL_S -> SetAddress(&SubldgJetDeepCvsL_S);
  branch_SubldgJetPtD_S -> SetAddress(&SubldgJetPtD_S);
  branch_SubldgJetAxis2_S -> SetAddress(&SubldgJetAxis2_S);
  //branch_SubldgJetAxis1_S -> SetAddress(&SubldgJetAxis1_S);
  branch_SubldgJetMult_S -> SetAddress(&SubldgJetMult_S);
  branch_SubldgJetQGLikelihood_S -> SetAddress(&SubldgJetQGLikelihood_S);
  //branch_SubldgJetPullMagnitude_S -> SetAddress(&SubldgJetPullMagnitude_S);
  branch_bjetPt_S -> SetAddress(&bjetPt_S);
  branch_bjetEta_S -> SetAddress(&bjetEta_S);
  branch_bjetPhi_S -> SetAddress(&bjetPhi_S);
  branch_bjetBdisc_S -> SetAddress(&bjetBdisc_S);
  branch_bjetMass_S -> SetAddress(&bjetMass_S);
  branch_bjetQGLikelihood_S -> SetAddress(&bjetQGLikelihood_S);
  branch_bjetCombinedCvsL_S -> SetAddress(&bjetCombinedCvsL_S);
  //branch_bjetDeepCvsL_S -> SetAddress(&bjetDeepCvsL_S);
  branch_bjetPtD_S -> SetAddress(&bjetPtD_S);
  branch_bjetAxis2_S -> SetAddress(&bjetAxis2_S);
  //branch_bjetAxis1_S -> SetAddress(&bjetAxis1_S);
  branch_bjetMult_S -> SetAddress(&bjetMult_S);
  //branch_bjetPFCharge_S -> SetAddress(&bjetPFCharge_S);
  branch_bjetLdgJetMass_S -> SetAddress(&bjetLdgJetMass_S);
  branch_bjetSubldgJetMass_S -> SetAddress(&bjetSubldgJetMass_S);
  branch_SoftDrop_n2_S -> SetAddress(&SoftDrop_n2_S);
  //branch_PullAngleJ1J2_S -> SetAddress(&PullAngleJ1J2_S);
  //branch_PullAngleJ2J1_S -> SetAddress(&PullAngleJ2J1_S);
  branch_DEtaJ1withJ2_S -> SetAddress(&DEtaJ1withJ2_S);
  branch_DEtaJ1withBJet_S -> SetAddress(&DEtaJ1withBJet_S);
  branch_DEtaJ2withBJet_S -> SetAddress(&DEtaJ2withBJet_S);
  branch_DEtaDijetwithBJet_S -> SetAddress(&DEtaDijetwithBJet_S);
  branch_DEtaJ1BJetwithJ2_S -> SetAddress(&DEtaJ1BJetwithJ2_S);
  branch_DEtaJ2BJetwithJ1_S -> SetAddress(&DEtaJ2BJetwithJ1_S);
  branch_DPhiJ1withJ2_S -> SetAddress(&DPhiJ1withJ2_S);
  branch_DPhiJ1withBJet_S -> SetAddress(&DPhiJ1withBJet_S);
  branch_DPhiJ2withBJet_S -> SetAddress(&DPhiJ2withBJet_S);
  branch_DPhiDijetwithBJet_S -> SetAddress(&DPhiDijetwithBJet_S);
  branch_DPhiJ1BJetwithJ2_S -> SetAddress(&DPhiJ1BJetwithJ2_S);
  branch_DPhiJ2BJetwithJ1_S -> SetAddress(&DPhiJ2BJetwithJ1_S);
  branch_DRJ1withJ2_S -> SetAddress(&DRJ1withJ2_S);
  branch_DRJ1withBJet_S -> SetAddress(&DRJ1withBJet_S);
  branch_DRJ2withBJet_S -> SetAddress(&DRJ2withBJet_S);
  branch_DRDijetwithBJet_S -> SetAddress(&DRDijetwithBJet_S);
  branch_DRJ1BJetwithJ2_S -> SetAddress(&DRJ1BJetwithJ2_S);
  branch_DRJ2BJetwithJ1_S -> SetAddress(&DRJ2BJetwithJ1_S);
  branch_dijetMassOverTrijetMass_S -> SetAddress(&dijetMassOverTrijetMass_S);
  
  // Background branches
  branch_eventWeight_B -> SetAddress(& eventWeight_B);
  branch_trijetPt_B -> SetAddress(&trijetPt_B);
  branch_trijetEta_B -> SetAddress(&trijetEta_B);
  branch_trijetPhi_B -> SetAddress(&trijetPhi_B);
  branch_trijetMass_B -> SetAddress(&trijetMass_B);
  branch_trijetPtDR_B -> SetAddress(&trijetPtDR_B);
  //branch_trijetPFCharge_B -> SetAddress(&trijetPFCharge_B);
  branch_trijetCombinedCvsL_B -> SetAddress(&trijetCombinedCvsL_B);
  //branch_trijetDeepCvsL_B -> SetAddress(&trijetDeepCvsL_B);
  branch_trijetPtD_B -> SetAddress(&trijetPtD_B);
  //branch_trijetAxis1_B -> SetAddress(&trijetAxis1_B);
  branch_trijetAxis2_B -> SetAddress(&trijetAxis2_B);
  branch_trijetMult_B -> SetAddress(&trijetMult_B);
  branch_trijetQGLikelihood_B -> SetAddress(&trijetQGLikelihood_B);
  branch_trijetQGLikelihood_avg_B -> SetAddress(&trijetQGLikelihood_avg_B);
  branch_trijetChiSquared_B -> SetAddress(&trijetChiSquared_B);
  branch_dijetPt_B -> SetAddress(&dijetPt_B);
  branch_dijetEta_B -> SetAddress(&dijetEta_B);
  branch_dijetPhi_B -> SetAddress(&dijetPhi_B);
  branch_dijetMass_B -> SetAddress(&dijetMass_B);
  branch_dijetPtDR_B -> SetAddress(&dijetPtDR_B);
  //branch_dijetPFCharge_B -> SetAddress(&dijetPFCharge_B);
  branch_dijetCombinedCvsL_B -> SetAddress(&dijetCombinedCvsL_B);
  //branch_dijetDeepCvsL_B -> SetAddress(&dijetDeepCvsL_B);
  branch_dijetPtD_B -> SetAddress(&dijetPtD_B);
  //branch_dijetAxis1_B -> SetAddress(&dijetAxis1_B);
  branch_dijetAxis2_B -> SetAddress(&dijetAxis2_B);
  branch_dijetMult_B -> SetAddress(&dijetMult_B);
  branch_dijetQGLikelihood_B -> SetAddress(&dijetQGLikelihood_B);
  branch_dijetQGLikelihood_avg_B -> SetAddress(&dijetQGLikelihood_avg_B);
  branch_dijetChiSquared_B -> SetAddress(&dijetChiSquared_B);
  branch_LdgJetPt_B -> SetAddress(&LdgJetPt_B);
  branch_LdgJetEta_B -> SetAddress(&LdgJetEta_B);
  branch_LdgJetPhi_B -> SetAddress(&LdgJetPhi_B);
  branch_LdgJetMass_B -> SetAddress(&LdgJetMass_B);
  //branch_LdgJetPFCharge_B -> SetAddress(&LdgJetPFCharge_B);
  branch_LdgJetBdisc_B -> SetAddress(&LdgJetBdisc_B);
  branch_LdgJetCombinedCvsL_B -> SetAddress(&LdgJetCombinedCvsL_B);
  //branch_LdgJetDeepCvsL_B -> SetAddress(&LdgJetDeepCvsL_B);
  branch_LdgJetPtD_B -> SetAddress(&LdgJetPtD_B);
  branch_LdgJetAxis2_B -> SetAddress(&LdgJetAxis2_B);
  //branch_LdgJetAxis1_B -> SetAddress(&LdgJetAxis1_B);
  branch_LdgJetMult_B -> SetAddress(&LdgJetMult_B);
  branch_LdgJetQGLikelihood_B -> SetAddress(&LdgJetQGLikelihood_B);
  //branch_LdgJetPullMagnitude_B -> SetAddress(&LdgJetPullMagnitude_B);
  branch_SubldgJetPt_B -> SetAddress(&SubldgJetPt_B);
  branch_SubldgJetEta_B -> SetAddress(&SubldgJetEta_B);
  branch_SubldgJetPhi_B -> SetAddress(&SubldgJetPhi_B);
  branch_SubldgJetMass_B -> SetAddress(&SubldgJetMass_B);
  //branch_SubldgJetPFCharge_B -> SetAddress(&SubldgJetPFCharge_B);
  branch_SubldgJetBdisc_B -> SetAddress(&SubldgJetBdisc_B);
  branch_SubldgJetCombinedCvsL_B -> SetAddress(&SubldgJetCombinedCvsL_B);
  //branch_SubldgJetDeepCvsL_B -> SetAddress(&SubldgJetDeepCvsL_B);
  branch_SubldgJetPtD_B -> SetAddress(&SubldgJetPtD_B);
  branch_SubldgJetAxis2_B -> SetAddress(&SubldgJetAxis2_B);
  //branch_SubldgJetAxis1_B -> SetAddress(&SubldgJetAxis1_B);
  branch_SubldgJetMult_B -> SetAddress(&SubldgJetMult_B);
  branch_SubldgJetQGLikelihood_B -> SetAddress(&SubldgJetQGLikelihood_B);
  //branch_SubldgJetPullMagnitude_B -> SetAddress(&SubldgJetPullMagnitude_B);
  branch_bjetPt_B -> SetAddress(&bjetPt_B);
  branch_bjetEta_B -> SetAddress(&bjetEta_B);
  branch_bjetPhi_B -> SetAddress(&bjetPhi_B);
  branch_bjetBdisc_B -> SetAddress(&bjetBdisc_B);
  branch_bjetMass_B -> SetAddress(&bjetMass_B);
  branch_bjetQGLikelihood_B -> SetAddress(&bjetQGLikelihood_B);
  branch_bjetCombinedCvsL_B -> SetAddress(&bjetCombinedCvsL_B);
  //branch_bjetDeepCvsL_B -> SetAddress(&bjetDeepCvsL_B);
  branch_bjetPtD_B -> SetAddress(&bjetPtD_B);
  branch_bjetAxis2_B -> SetAddress(&bjetAxis2_B);
  //branch_bjetAxis1_B -> SetAddress(&bjetAxis1_B);
  branch_bjetMult_B -> SetAddress(&bjetMult_B);
  //branch_bjetPFCharge_B -> SetAddress(&bjetPFCharge_B);
  branch_bjetLdgJetMass_B -> SetAddress(&bjetLdgJetMass_B);
  branch_bjetSubldgJetMass_B -> SetAddress(&bjetSubldgJetMass_B);
  branch_SoftDrop_n2_B -> SetAddress(&SoftDrop_n2_B);
  //branch_PullAngleJ1J2_B -> SetAddress(&PullAngleJ1J2_B);
  //branch_PullAngleJ2J1_B -> SetAddress(&PullAngleJ2J1_B);
  branch_DEtaJ1withJ2_B -> SetAddress(&DEtaJ1withJ2_B);
  branch_DEtaJ1withBJet_B -> SetAddress(&DEtaJ1withBJet_B);
  branch_DEtaJ2withBJet_B -> SetAddress(&DEtaJ2withBJet_B);
  branch_DEtaDijetwithBJet_B -> SetAddress(&DEtaDijetwithBJet_B);
  branch_DEtaJ1BJetwithJ2_B -> SetAddress(&DEtaJ1BJetwithJ2_B);
  branch_DEtaJ2BJetwithJ1_B -> SetAddress(&DEtaJ2BJetwithJ1_B);
  branch_DPhiJ1withJ2_B -> SetAddress(&DPhiJ1withJ2_B);
  branch_DPhiJ1withBJet_B -> SetAddress(&DPhiJ1withBJet_B);
  branch_DPhiJ2withBJet_B -> SetAddress(&DPhiJ2withBJet_B);
  branch_DPhiDijetwithBJet_B -> SetAddress(&DPhiDijetwithBJet_B);
  branch_DPhiJ1BJetwithJ2_B -> SetAddress(&DPhiJ1BJetwithJ2_B);
  branch_DPhiJ2BJetwithJ1_B -> SetAddress(&DPhiJ2BJetwithJ1_B);
  branch_DRJ1withJ2_B -> SetAddress(&DRJ1withJ2_B);
  branch_DRJ1withBJet_B -> SetAddress(&DRJ1withBJet_B);
  branch_DRJ2withBJet_B -> SetAddress(&DRJ2withBJet_B);
  branch_DRDijetwithBJet_B -> SetAddress(&DRDijetwithBJet_B);
  branch_DRJ1BJetwithJ2_B -> SetAddress(&DRJ1BJetwithJ2_B);
  branch_DRJ2BJetwithJ1_B -> SetAddress(&DRJ2BJetwithJ1_B);
  branch_dijetMassOverTrijetMass_B -> SetAddress(&dijetMassOverTrijetMass_B);


  //================================================================================================//
  //                            Gen Trijet subjets selection                                        //
  //================================================================================================//

  //Soti
  if (!fEvent.isMC()) return;
  //Get top quarks
  vector<genParticle> GenTops = GetGenParticles(fEvent.genparticles().getGenParticles(), 6);
  vector<genParticle> GenTops_BQuark, GenTops_LdgQuark, GenTops_SubldgQuark, GenTops_Quarks;

  //For loop: top quarks
  for (auto& top: GenTops){
    vector<genParticle> quarks, bquarks;

    // Get top decay products (Returns direct b quarks, quarks from W boson)
    getTopDecayProducts(fEvent, top, quarks, bquarks);
    bool foundB = (bquarks.size() == 1);
    
    //Fill control plots (gen top pt)
    if (quarks.size() == 2) hCtrl_hadrGenTopPt -> Fill(top.pt());
    else hCtrl_lepGenTopPt -> Fill(top.pt());
    
    if (cfg_SelectionsType == "hadronic"){
      // Skip if at least one top does not decay hadronically
      if (!(quarks.size() == 2 && foundB)) return;
    }
    else{
      // If top does not decay into W+b return
      if (!foundB) return;
    }
    
    if (quarks.size() == 2){
      // Fill vectors for b-quarks, leading and subleading quarks coming from tops                                                                                                   
      GenTops_Quarks.push_back(bquarks.at(0));
      GenTops_Quarks.push_back(quarks.at(0));
      GenTops_Quarks.push_back(quarks.at(1));
      
      GenTops_LdgQuark.push_back(quarks.at(0));
      GenTops_SubldgQuark.push_back(quarks.at(1));
      GenTops_BQuark.push_back(bquarks.at(0));
    }
  }//for (auto& top: GenTops)                                                                                                                                                      
  
  //Definitions
  int imatched =0;
  vector <Jet> MCtrue_LdgJet, MCtrue_SubldgJet, MCtrue_Bjet;   // Keep matched jets
  vector <genParticle> MGen_LdgJet, MGen_SubldgJet, MGen_Bjet; // Keep get particles
  vector <Jet> BJetCand, MC_Jets;
  vector <double> dRminB;
  Jet firstBjet;

  //========================================================================================================
  //======= B jet matching (Loop over all Jets)
  //========================================================================================================  
  for (size_t i=0; i<GenTops_BQuark.size(); i++){
    genParticle BQuark      = GenTops_BQuark.at(i);
    Jet mcMatched_BJet;
    double dRmin  = 99999.9;

    for (auto& bjet: jetData.getSelectedJets()){
      double dR        = ROOT::Math::VectorUtil::DeltaR( bjet.p4(), BQuark.p4());
      double dPtOverPt = std::abs((bjet.pt() - BQuark.pt())/BQuark.pt());
      if (dR > cfg_DeltaRCut)   continue;
      if (dR > dRmin)           continue;
      if (dPtOverPt > cfg_DeltaPtOverPtCut) continue;

      dRmin  = dR;
      mcMatched_BJet = bjet;
    }
    dRminB.push_back(dRmin);
    BJetCand.push_back(mcMatched_BJet);
  } //for (size_t i=0; i<GenTops.size(); i++){

  if (0) std::cout<<"N hadronic tops "<<GenTops_LdgQuark.size()<<std::endl;
  hCtrl_NhadronicTops -> Fill(GenTops_LdgQuark.size());
    
  //========================================================================================================
  //======= Dijet matching (Loop over all Jets)
  //========================================================================================================

  for (size_t i=0; i<GenTops_BQuark.size(); i++){
    genParticle LdgQuark    = GenTops_LdgQuark.at(i);
    genParticle SubldgQuark = GenTops_SubldgQuark.at(i);
      
    Jet mcMatched_LdgJet;
    Jet mcMatched_SubldgJet;
      
    double dR1min, dR2min, dPtOverPt1min, dPtOverPt2min;
    dR1min = dR2min = dPtOverPt1min = dPtOverPt2min = 99999.9;

    for (auto& jet: jetData.getSelectedJets()){
      bool same = false;
      for (size_t k=0; k<GenTops_BQuark.size(); k++){
	if (dRminB.at(k) <= cfg_DeltaRCut){
	  if( areSameJets(jet,BJetCand.at(k))) same = true; //Skip the jets that are matched with bquarks
	}
      }

      if (same) continue;

      double dR1 = ROOT::Math::VectorUtil::DeltaR(jet.p4(), LdgQuark.p4());
      double dR2 = ROOT::Math::VectorUtil::DeltaR(jet.p4(), SubldgQuark.p4());
	
      if (std::min(dR1, dR2) > cfg_DeltaRCut) continue;
	
      double dPtOverPt1 = std::abs((jet.pt() - LdgQuark.pt())/LdgQuark.pt());
      double dPtOverPt2 = std::abs( (jet.pt() - SubldgQuark.pt())/SubldgQuark.pt());

      if (dR1 < dR2) // jet is closer to the leading quark
	{
	  if (dR1 < dR1min)
	    {
	      if(dPtOverPt1 < cfg_DeltaPtOverPtCut)
		{
		  dR1min = dR1;
		  dPtOverPt1min= dPtOverPt1;
		  mcMatched_LdgJet = jet;
		}
	    }	
	  else if (dR2 <= cfg_DeltaRCut && dR2 < dR2min)
	    {
	      if (dPtOverPt2 < cfg_DeltaPtOverPtCut)
		{
		  dR2min  = dR2;
		  dPtOverPt2min = dPtOverPt2;
		  mcMatched_SubldgJet = jet;
		}
	    }
	} //if (dR1 < dR2)

      else // jet is closer to the subleading quark
	{
	  if (dR2 < dR2min)
	    {
	      if(dPtOverPt2 < cfg_DeltaPtOverPtCut)
		{
		  dR2min  = dR2;
		  dPtOverPt2min = dPtOverPt2;
		  mcMatched_SubldgJet = jet;
		}
	    }
	  else if (dR1 <= cfg_DeltaRCut && dR1 < dR1min)
	    {
	      if (dPtOverPt1 < cfg_DeltaPtOverPtCut)
		{
		  dR1min  = dR1;
		  dPtOverPt1min = dPtOverPt1;
		  mcMatched_LdgJet = jet;
		}
	    }
	}// else   
    } //for (auto& jet: jetData.getSelectedJets()){
  
    //Gen Quarks: Pt
    hGenQuark_Pt -> Fill((dR1min<= cfg_DeltaRCut), LdgQuark.pt());
    hGenQuark_Pt -> Fill((dR2min <= cfg_DeltaRCut), SubldgQuark.pt());
    hGenQuark_Pt -> Fill((dRminB.at(i) <= cfg_DeltaRCut), GenTops_BQuark.at(i).pt());
    //Genuine if all the top quarks are matched
    genParticle top = GenTops.at(i);
    bool genuine = (dR1min<= cfg_DeltaRCut && dR2min <= cfg_DeltaRCut && dRminB.at(i) <= cfg_DeltaRCut);
      
    if (genuine){
      imatched ++;
      hGenTop_Pt->Fill(true, top.pt());
      MCtrue_LdgJet.push_back(mcMatched_LdgJet);
      MCtrue_SubldgJet.push_back(mcMatched_SubldgJet);
      MCtrue_Bjet.push_back(BJetCand.at(i));
	
      MC_Jets.push_back(mcMatched_LdgJet);
      MC_Jets.push_back(mcMatched_SubldgJet);
      MC_Jets.push_back(BJetCand.at(i));

      MGen_LdgJet.push_back(GenTops_LdgQuark.at(i));
      MGen_SubldgJet.push_back(GenTops_SubldgQuark.at(i));
      MGen_Bjet.push_back(GenTops_BQuark.at(i));
    }
    else {
      hGenTop_Pt->Fill(false, top.pt());
    }
  } // for (size_t i=0; i<GenTops_LdgQuark.size(); i++){
  hNmatchedTop ->Fill(imatched);
  
  if (0) std::cout<<"nmatched top candidates: "<<imatched<<std::endl;

  // USED in 2016 analysis
  // Matching criterion: Select events with DeltaR(q,q') > 0.8
  // Otherwise: Event not used for the training    
  // size_t nTops = GenTops_BQuark.size();
  // for (size_t i=0; i<nTops; i++)
  //   {
  //  	double dR12 = ROOT::Math::VectorUtil::DeltaR(GenTops_LdgQuark.at(i).p4(),    GenTops_SubldgQuark.at(i).p4());
  //  	double dR1b = ROOT::Math::VectorUtil::DeltaR(GenTops_LdgQuark.at(i).p4(),    GenTops_BQuark.at(i).p4());
  //  	double dR2b = ROOT::Math::VectorUtil::DeltaR(GenTops_SubldgQuark.at(i).p4(), GenTops_BQuark.at(i).p4());
  //  	double dRmin = min(min(dR12, dR1b), dR2b);
  //  	if (dRmin < 0.8) return;
  //   }

  //================================================================================================//                      
  //                                    Top Candidates                                              //
  //================================================================================================//
  bool applyBjetPtCut = true;
  TrijetSelections TopCandidates;
  int index0 = -1;
  for (auto& jet0: jetData.getSelectedJets()){
    index0++;
    if (isLepton(jet0, selectedElectrons, selectedMuons)) continue;
    int index1 = -1;
    for (auto& jet1: jetData.getSelectedJets()){
      index1++;
      if (isLepton(jet1, selectedElectrons, selectedMuons)) continue;
      if (index1 < index0) continue;
      if (areSameJets(jet1, jet0)) continue;
      int index2 = -1;
      for (auto& jet2: jetData.getSelectedJets()){
	index2++;
	if (isLepton(jet2, selectedElectrons, selectedMuons)) continue;
	if (index2 < index1) continue;
	if (areSameJets(jet2,  jet1)) continue;
	if (areSameJets(jet2,  jet0)) continue;

	  
	//********************************** Bjet Matched OR dijet matched*********************************//
	if ( isBJet(jet0, MCtrue_Bjet)    ||    (isWsubjet(jet1,MCtrue_LdgJet, MCtrue_SubldgJet)  &&  isWsubjet(jet2,MCtrue_LdgJet, MCtrue_SubldgJet))){
	  //getTrijet(jet0, jet1, jet2, applyBjetPtCut, 40);
	  if (applyBjetPtCut && (jet0.pt() < 40)) continue;
	  TopCandidates.BJet.push_back(jet0);
	  TopCandidates.Jet1.push_back(getLeadingSubleadingJet(jet1,jet2,"leading"));
	  TopCandidates.Jet2.push_back(getLeadingSubleadingJet(jet1,jet2,"subleading"));
	}
	  
	else if ( isBJet(jet1, MCtrue_Bjet)    ||  (isWsubjet(jet0,MCtrue_LdgJet, MCtrue_SubldgJet)  &&  isWsubjet(jet2,MCtrue_LdgJet, MCtrue_SubldgJet))){
	  if (applyBjetPtCut && (jet1.pt() < 40)) continue;
	  TopCandidates.BJet.push_back(jet1);
	  TopCandidates.Jet1.push_back(getLeadingSubleadingJet(jet0,jet2,"leading"));
	  TopCandidates.Jet2.push_back(getLeadingSubleadingJet(jet0,jet2,"subleading"));

	}
	else if ( isBJet(jet2, MCtrue_Bjet)    ||  (isWsubjet(jet0,MCtrue_LdgJet, MCtrue_SubldgJet)  &&  isWsubjet(jet1,MCtrue_LdgJet, MCtrue_SubldgJet))){
	  if (applyBjetPtCut && (jet2.pt() < 40)) continue;
	  TopCandidates.BJet.push_back(jet2);
	  TopCandidates.Jet1.push_back(getLeadingSubleadingJet(jet0,jet1,"leading"));
	  TopCandidates.Jet2.push_back(getLeadingSubleadingJet(jet0,jet1,"subleading"));
	}
	  
	//********************************** One of the dijet subjets matched*********************************//
	else if (isWsubjet(jet0,MCtrue_LdgJet, MCtrue_SubldgJet)){

	  if (jet1.bjetDiscriminator() > jet2.bjetDiscriminator()){
	    if (applyBjetPtCut && (jet1.pt() < 40)) continue;
	    TopCandidates.BJet.push_back(jet1);
	    TopCandidates.Jet1.push_back(getLeadingSubleadingJet(jet0,jet2,"leading"));
	    TopCandidates.Jet2.push_back(getLeadingSubleadingJet(jet0,jet2,"subleading"));
	  }
	  else{
	    if (applyBjetPtCut && (jet2.pt() < 40)) continue;
	    TopCandidates.BJet.push_back(jet2);
	    TopCandidates.Jet1.push_back(getLeadingSubleadingJet(jet0,jet1,"leading"));
	    TopCandidates.Jet2.push_back(getLeadingSubleadingJet(jet0,jet1,"subleading"));
	  }
	}
	else if (isWsubjet(jet1,MCtrue_LdgJet, MCtrue_SubldgJet)){

	  if (jet0.bjetDiscriminator() > jet2.bjetDiscriminator()){
	    if (applyBjetPtCut && (jet0.pt() < 40)) continue;
	    TopCandidates.BJet.push_back(jet0);
	    TopCandidates.Jet1.push_back(getLeadingSubleadingJet(jet1,jet2,"leading"));
	    TopCandidates.Jet2.push_back(getLeadingSubleadingJet(jet1,jet2,"subleading"));
	  }
	  else{
	    if (applyBjetPtCut && (jet2.pt() < 40)) continue;
	    TopCandidates.BJet.push_back(jet2);
	    TopCandidates.Jet1.push_back(getLeadingSubleadingJet(jet0,jet1,"leading"));
	    TopCandidates.Jet2.push_back(getLeadingSubleadingJet(jet0,jet1,"subleading"));
	  }
	}
	else if (isWsubjet(jet2,MCtrue_LdgJet, MCtrue_SubldgJet)){

	  if (jet0.bjetDiscriminator() > jet1.bjetDiscriminator()){
	    if (applyBjetPtCut && (jet0.pt() < 40)) continue;
	    TopCandidates.BJet.push_back(jet0);
	    TopCandidates.Jet1.push_back(getLeadingSubleadingJet(jet1,jet2,"leading"));
	    TopCandidates.Jet2.push_back(getLeadingSubleadingJet(jet1,jet2,"subleading"));
	  }
	  else{
	    if (applyBjetPtCut && (jet1.pt() < 40)) continue;
	    TopCandidates.BJet.push_back(jet1);
	    TopCandidates.Jet1.push_back(getLeadingSubleadingJet(jet0,jet2,"leading"));
	    TopCandidates.Jet2.push_back(getLeadingSubleadingJet(jet0,jet2,"subleading"));
	  }
	}	  	 	  
	//********************** Non of the three subjets is matched************************//
	  
	else if (jet0.bjetDiscriminator() > jet1.bjetDiscriminator() && jet0.bjetDiscriminator() > jet2.bjetDiscriminator()){            
	  if (applyBjetPtCut && (jet0.pt() < 40)) continue;
	  TopCandidates.BJet.push_back(jet0);
	  TopCandidates.Jet1.push_back(getLeadingSubleadingJet(jet1,jet2,"leading"));
	  TopCandidates.Jet2.push_back(getLeadingSubleadingJet(jet1,jet2,"subleading"));	    
	}
	else if (jet1.bjetDiscriminator() > jet0.bjetDiscriminator() && jet1.bjetDiscriminator() > jet2.bjetDiscriminator()){
	  if (applyBjetPtCut && (jet1.pt() < 40)) continue;
	  TopCandidates.BJet.push_back(jet1);	    
	  TopCandidates.Jet1.push_back(getLeadingSubleadingJet(jet0,jet2,"leading"));
	  TopCandidates.Jet2.push_back(getLeadingSubleadingJet(jet0,jet2,"subleading"));
	}
	else if (jet2.bjetDiscriminator() > jet0.bjetDiscriminator() && jet2.bjetDiscriminator() > jet1.bjetDiscriminator()){
	  if (applyBjetPtCut && (jet2.pt() < 40)) continue;
	  TopCandidates.BJet.push_back(jet2);
	  TopCandidates.Jet1.push_back(getLeadingSubleadingJet(jet0,jet1,"leading"));
	  TopCandidates.Jet2.push_back(getLeadingSubleadingJet(jet0,jet1,"subleading"));
	}
      } //for (auto& jet2: jetData.getSelectedJets()){
    } //for (auto& jet1: jetData.getSelectedJets()){
  } //for (auto& jet0: jetData.getSelectedJets()){

    //========================================================================================================
    //                       Identification of fake, genuine TopCandidates
    //========================================================================================================
  vector <bool> GenuineTop; 
  int jmatched = 0;
  //    std::cout<<"==="<<std::endl;
  for (size_t i=0; i<TopCandidates.BJet.size(); i++){
    bool genuine = false;
    for (size_t j=0; j<MCtrue_Bjet.size(); j++){
      Jet jet1 = TopCandidates.Jet1.at(i);
      Jet jet2 = TopCandidates.Jet2.at(i);
      Jet bjet = TopCandidates.BJet.at(i);
      
      bool same1 = areSameJets(jet1, MCtrue_LdgJet.at(j))    && areSameJets(jet2, MCtrue_SubldgJet.at(j)) && areSameJets(bjet,  MCtrue_Bjet.at(j));
      bool same2 = areSameJets(jet1, MCtrue_SubldgJet.at(j)) && areSameJets(jet2, MCtrue_LdgJet.at(j))    && areSameJets(bjet,  MCtrue_Bjet.at(j));
      if (same1 || same2){
	genuine = true;
      }
    }
    if (genuine) jmatched++;
    GenuineTop.push_back(genuine);
  }
  hNmatchedTrijets ->Fill(jmatched);

  //========================================================================================================
  //                                 Fill trees
  //========================================================================================================
  // Definitions
  Bool_t isGenuineTop = false;
  // Chi Squared
  const double mTopMass = 172.34;  // arXiv:1812.10534
  const double mWMass   = 80.385;    // TopSelection (Chi-Squared Method)
  const double sigmaTrijet = 27.5; // from Mikela's fitting (Chi-Squared Method)
  const double sigmaDijet  = 10.59; // from Mikela's fitting (Chi-Squared Method)
    
  for (size_t i=0; i<TopCandidates.BJet.size(); i++){
    Jet jet1 = TopCandidates.Jet1.at(i);
    Jet jet2 = TopCandidates.Jet2.at(i);
    Jet bjet = TopCandidates.BJet.at(i);
      
    math::XYZTLorentzVector TrijetP4, DijetP4;
    TrijetP4 = bjet.p4()+jet1.p4()+jet2.p4();
    DijetP4  = jet1.p4()+jet2.p4();
      
    // Keep only tops in the low/high pT range
    //if (!cfg_TopCandPtCut.passedCut(TrijetP4.pt())) continue;
      
    //Debug:
    if (applyBjetPtCut && TopCandidates.BJet.at(i).pt() < 40) std::cout<<"never reach here"<<std::endl;
    isGenuineTop = GenuineTop.at(i);
      
    // Compute distances
    double dRJ12 = ROOT::Math::VectorUtil::DeltaR(jet1.p4(), jet2.p4());      
    // Compute Soft drop value N2
    double softDrop_n2 = min(jet2.pt(), jet1.pt())/( (jet2.pt() + jet1.pt() )*dRJ12*dRJ12);
      
    hSoftDrop_n2   -> Fill(isGenuineTop, softDrop_n2);
    //hPullAngleJ1J2 -> Fill(isGenuineTop, getJetPullAngle(jet1, jet2));  
    //hPullAngleJ2J1 -> Fill(isGenuineTop, getJetPullAngle(jet2, jet1));
    hDEtaJ1withJ2      -> Fill(isGenuineTop, std::abs(jet1.eta() - jet2.eta()));
    hDEtaJ1withBJet    -> Fill(isGenuineTop, std::abs(jet1.eta() - bjet.eta()));
    hDEtaJ2withBJet    -> Fill(isGenuineTop, std::abs(jet2.eta() - bjet.eta()));
    hDEtaDijetwithBJet -> Fill(isGenuineTop, std::abs(DijetP4.Eta() - bjet.eta()));
    hDEtaJ1BJetwithJ2  -> Fill(isGenuineTop, std::abs((jet1.p4() + bjet.p4()).Eta() - jet2.eta()));
    hDEtaJ2BJetwithJ1  -> Fill(isGenuineTop, std::abs((jet2.p4() + bjet.p4()).Eta() - jet1.eta()));
    
    hDPhiJ1withJ2      -> Fill(isGenuineTop, std::abs(ROOT::Math::VectorUtil::DeltaPhi(jet1.p4(), jet2.p4()))); 
    hDPhiJ1withBJet    -> Fill(isGenuineTop, std::abs(ROOT::Math::VectorUtil::DeltaPhi(jet1.p4(), bjet.p4())));
    hDPhiJ2withBJet    -> Fill(isGenuineTop, std::abs(ROOT::Math::VectorUtil::DeltaPhi(jet2.p4(), bjet.p4())));
    hDPhiDijetwithBJet -> Fill(isGenuineTop, std::abs(ROOT::Math::VectorUtil::DeltaPhi(DijetP4, bjet.p4())));
    hDPhiJ1BJetwithJ2  -> Fill(isGenuineTop, std::abs(ROOT::Math::VectorUtil::DeltaPhi((jet1.p4()+bjet.p4()), jet2.p4())));
    hDPhiJ2BJetwithJ1  -> Fill(isGenuineTop, std::abs(ROOT::Math::VectorUtil::DeltaPhi((jet2.p4()+bjet.p4()), jet1.p4())));
      
    hDRJ1withJ2       -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR(jet1.p4(), jet2.p4()));
    hDRJ1withBJet     -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR(jet1.p4(), bjet.p4()));
    hDRJ2withBJet     -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR(jet2.p4(), bjet.p4()));
    hDRDijetwithBJet  -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR(DijetP4, bjet.p4()));
    hDRJ1BJetwithJ2   -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR((jet1.p4()+bjet.p4()), jet2.p4()));
    hDRJ2BJetwithJ1   -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR((jet2.p4()+bjet.p4()), jet1.p4()));
      
    if (passElectron){
      Electron ele = selectedElectrons.at(0);
      hDRJ1withElectron -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR(ele.p4(), jet1.p4()));
      hDRJ2withElectron -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR(ele.p4(), jet2.p4()));
      hDRBwithElectron  -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR(ele.p4(), bjet.p4()));      
      if (jet1.QGTaggerAK4PFCHSptD() > 0.80) hDRJ1withElectron_ptD0p80 -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR(ele.p4(), jet1.p4()));
      if (jet2.QGTaggerAK4PFCHSptD() > 0.80) hDRJ2withElectron_ptD0p80 -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR(ele.p4(), jet2.p4()));
      if (bjet.QGTaggerAK4PFCHSptD() > 0.80) hDRBwithElectron_ptD0p80  -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR(ele.p4(), bjet.p4()));
    }
    if (passMuon){
      Muon mu = selectedMuons.at(0);
      hDRJ1withMuon -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR(mu.p4(), jet1.p4()));
      hDRJ2withMuon -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR(mu.p4(), jet2.p4()));
      hDRBwithMuon  -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR(mu.p4(), bjet.p4()));      
      if (jet1.QGTaggerAK4PFCHSptD() > 0.80) hDRJ1withMuon_ptD0p80 -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR(mu.p4(), jet1.p4()));
      if (jet2.QGTaggerAK4PFCHSptD() > 0.80) hDRJ2withMuon_ptD0p80 -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR(mu.p4(), jet2.p4()));
      if (bjet.QGTaggerAK4PFCHSptD() > 0.80) hDRBwithMuon_ptD0p80  -> Fill(isGenuineTop, ROOT::Math::VectorUtil::DeltaR(mu.p4(), bjet.p4()));
    }
	
    // Trijet system
    hTrijetPt           -> Fill(isGenuineTop, TrijetP4.Pt());
    hTrijetEta          -> Fill(isGenuineTop, TrijetP4.Eta());
    hTrijetPhi          -> Fill(isGenuineTop, TrijetP4.Phi());
    hTrijetMass         -> Fill(isGenuineTop, TrijetP4.M());
    hTrijetPtDr         -> Fill(isGenuineTop, TrijetP4.Pt() * ROOT::Math::VectorUtil::DeltaR(DijetP4, bjet.p4()));
    //hTrijetPFCharge     -> Fill(isGenuineTop, jet1.pfcharge()+jet2.pfcharge()+bjet.pfcharge());
    hTrijetCombinedCvsL -> Fill(isGenuineTop, (jet1.pfCombinedCvsLJetTags() + jet2.pfCombinedCvsLJetTags() + bjet.pfCombinedCvsLJetTags())/3.);
    //hTrijetDeepCvsL     -> Fill(isGenuineTop, (jet1.pfDeepCSVCvsLJetTags()  + jet2.pfDeepCSVCvsLJetTags()  + bjet.pfDeepCSVCvsLJetTags())/3.);
    hTrijetPtD              -> Fill(isGenuineTop, (jet1.QGTaggerAK4PFCHSptD() + jet2.QGTaggerAK4PFCHSptD() + bjet.QGTaggerAK4PFCHSptD())/3.);
    hTrijetAxis2            -> Fill(isGenuineTop, (jet1.QGTaggerAK4PFCHSaxis2() + jet2.QGTaggerAK4PFCHSaxis2() + bjet.QGTaggerAK4PFCHSaxis2())/3.);
    //hTrijetAxis1            -> Fill(isGenuineTop, (jet1.axis1() + jet2.axis1() + bjet.axis1())/3.);
    hTrijetMult             -> Fill(isGenuineTop, (jet1.QGTaggerAK4PFCHSmult() + jet2.QGTaggerAK4PFCHSmult()  + bjet.QGTaggerAK4PFCHSmult())/3.);
    hTrijetQGLikelihood     -> Fill(isGenuineTop, jet1.QGTaggerAK4PFCHSqgLikelihood()*jet2.QGTaggerAK4PFCHSqgLikelihood()*bjet.QGTaggerAK4PFCHSqgLikelihood()); //.QGTaggerAK4PFCHSqgLikelihood()
    hTrijetQGLikelihood_avg -> Fill(isGenuineTop, (jet1.QGTaggerAK4PFCHSqgLikelihood()+jet2.QGTaggerAK4PFCHSqgLikelihood()+bjet.QGTaggerAK4PFCHSqgLikelihood())/3.);
    hTrijetChiSquared -> Fill(isGenuineTop, (TrijetP4.M() - mTopMass)/sigmaTrijet); 
      
    // Dijet system
    hDijetPt                 -> Fill(isGenuineTop, DijetP4.Pt());
    hDijetEta                -> Fill(isGenuineTop, DijetP4.Eta());
    hDijetPhi                -> Fill(isGenuineTop, DijetP4.Phi());
    hDijetMass               -> Fill(isGenuineTop, DijetP4.M());
    //hDijetPFCharge           -> Fill(isGenuineTop, jet1.pfcharge()+jet2.pfcharge());
    hDijetPtDr               -> Fill(isGenuineTop, DijetP4.Pt() * dRJ12);
    hDijetCombinedCvsL       -> Fill(isGenuineTop, (jet1.pfCombinedCvsLJetTags() + jet2.pfCombinedCvsLJetTags())/2.);
    //hDijetDeepCvsL           -> Fill(isGenuineTop, (jet1.pfDeepCSVCvsLJetTags() + jet2.pfDeepCSVCvsLJetTags())/2.);
    hDijetPtD                -> Fill(isGenuineTop, (jet1.QGTaggerAK4PFCHSptD() + jet2.QGTaggerAK4PFCHSptD())/2.);
    hDijetAxis2              -> Fill(isGenuineTop, (jet1.QGTaggerAK4PFCHSaxis2() + jet2.QGTaggerAK4PFCHSaxis2())/2.);
    //hDijetAxis1              -> Fill(isGenuineTop, (jet1.axis1() + jet2.axis1())/2.);
    hDijetMult               -> Fill(isGenuineTop, (jet1.QGTaggerAK4PFCHSmult() + jet2.QGTaggerAK4PFCHSmult())/2.);
    hDijetQGLikelihood       -> Fill(isGenuineTop, jet1.QGTaggerAK4PFCHSqgLikelihood()*jet2.QGTaggerAK4PFCHSqgLikelihood());
    hDijetQGLikelihood_avg   -> Fill(isGenuineTop, (jet1.QGTaggerAK4PFCHSqgLikelihood()+jet2.QGTaggerAK4PFCHSqgLikelihood())/2.);
    hDijetMassOverTrijetMass -> Fill(isGenuineTop, DijetP4.M()/TrijetP4.M());
    hDijetChiSquared         -> Fill(isGenuineTop, (DijetP4.M() - mWMass)/sigmaDijet);
      
      
    // Leading jet from dijet system
    hLdgJetPt                -> Fill(isGenuineTop, jet1.pt());
    hLdgJetEta               -> Fill(isGenuineTop, jet1.eta());
    hLdgJetPhi               -> Fill(isGenuineTop, jet1.phi());
    hLdgJetMass              -> Fill(isGenuineTop, jet1.p4().M());
    //hLdgJetPFCharge          -> Fill(isGenuineTop, jet1.pfcharge());
    hLdgJetBdisc             -> Fill(isGenuineTop, jet1.bjetDiscriminator());
    hLdgJetCombinedCvsL      -> Fill(isGenuineTop, jet1.pfCombinedCvsLJetTags());
    //hLdgJetDeepCvsL          -> Fill(isGenuineTop, jet1.pfDeepCSVCvsLJetTags());
    hLdgJetPtD               -> Fill(isGenuineTop, jet1.QGTaggerAK4PFCHSptD());
    hLdgJetAxis2             -> Fill(isGenuineTop, jet1.QGTaggerAK4PFCHSaxis2());
    //hLdgJetAxis1             -> Fill(isGenuineTop, jet1.axis1());
    hLdgJetMult              -> Fill(isGenuineTop, jet1.QGTaggerAK4PFCHSmult());
    hLdgJetQGLikelihood      -> Fill(isGenuineTop, jet1.QGTaggerAK4PFCHSqgLikelihood());
    //hLdgJetPullMagnitude     -> Fill(isGenuineTop, getJetPullMagnitude(jet1));      
      
    // Subleading jet from dijet system
    hSubldgJetPt             -> Fill(isGenuineTop, jet2.pt());
    hSubldgJetEta            -> Fill(isGenuineTop, jet2.eta());
    hSubldgJetPhi            -> Fill(isGenuineTop, jet2.phi());
    hSubldgJetBdisc          -> Fill(isGenuineTop, jet2.bjetDiscriminator());
      
    hSubldgJetCombinedCvsL   -> Fill(isGenuineTop, jet2.pfCombinedCvsLJetTags());
    //hSubldgJetDeepCvsL       -> Fill(isGenuineTop, jet2.pfDeepCSVCvsLJetTags());
    hSubldgJetPtD            -> Fill(isGenuineTop, jet2.QGTaggerAK4PFCHSptD());
    hSubldgJetAxis2          -> Fill(isGenuineTop, jet2.QGTaggerAK4PFCHSaxis2());
    //hSubldgJetAxis1          -> Fill(isGenuineTop, jet2.axis1());
    hSubldgJetMult           -> Fill(isGenuineTop, jet2.QGTaggerAK4PFCHSmult());
    hSubldgJetQGLikelihood   -> Fill(isGenuineTop, jet2.QGTaggerAK4PFCHSqgLikelihood());
    //hSubldgJetPFCharge       -> Fill(isGenuineTop, jet2.pfcharge());
    //hSubldgJetPullMagnitude  -> Fill(isGenuineTop, getJetPullMagnitude(jet2));
      
    hBJetPt                  -> Fill(isGenuineTop, bjet.pt());
    hBJetEta                 -> Fill(isGenuineTop, bjet.eta());
    hBJetPhi                 -> Fill(isGenuineTop, bjet.phi());
    hBJetBdisc               -> Fill(isGenuineTop, bjet.bjetDiscriminator());
    hBJetMass                -> Fill(isGenuineTop, bjet.p4().M());
    hBJetQGLikelihood        -> Fill(isGenuineTop, bjet.QGTaggerAK4PFCHSqgLikelihood());
    hBJetCombinedCvsL        -> Fill(isGenuineTop, bjet.pfCombinedCvsLJetTags());
    //hBJetDeepCvsL            -> Fill(isGenuineTop, bjet.pfDeepCSVCvsLJetTags());
    hBJetPtD                 -> Fill(isGenuineTop, bjet.QGTaggerAK4PFCHSptD());
    hBJetAxis2               -> Fill(isGenuineTop, bjet.QGTaggerAK4PFCHSaxis2());
    //hBJetAxis1               -> Fill(isGenuineTop, bjet.axis1());
    hBJetMult                -> Fill(isGenuineTop, bjet.QGTaggerAK4PFCHSmult());
    //hBJetPFCharge            -> Fill(isGenuineTop, bjet.pfcharge());
      
    hBJetLdgJet_Mass         -> Fill(isGenuineTop, (bjet.p4()+jet1.p4()).M());
    //hBJetLdgJet_PFCharge     -> Fill(isGenuineTop, bjet.pfcharge() + jet1.pfcharge());
    hBJetSubldgJet_Mass      -> Fill(isGenuineTop, (bjet.p4()+jet2.p4()).M());
    //hBJetSubldgJet_PFCharge  -> Fill(isGenuineTop, bjet.pfcharge()+jet2.pfcharge());

    if (isGenuineTop){

      // Fill signal background
      eventWeight_S = fEventWeight.getWeight();
	  
      // Trijet Variables
      trijetPt_S		   = TrijetP4.Pt();
      trijetEta_S		   = TrijetP4.Eta();
      trijetPhi_S		   = TrijetP4.Phi();
      trijetMass_S		   = TrijetP4.M();
      trijetPtDR_S		   = TrijetP4.Pt() * ROOT::Math::VectorUtil::DeltaR(DijetP4, bjet.p4());
      //trijetPFCharge_S	   = jet1.pfcharge() + jet2.pfcharge() + bjet.pfcharge();
      trijetCombinedCvsL_S	   = (jet1.pfCombinedCvsLJetTags() + jet2.pfCombinedCvsLJetTags() + bjet.pfCombinedCvsLJetTags())/3.;
      //trijetDeepCvsL_S	   = (jet1.pfDeepCSVCvsLJetTags() + jet2.pfDeepCSVCvsLJetTags() + bjet.pfDeepCSVCvsLJetTags())/3.;
      trijetPtD_S		   = (jet1.QGTaggerAK4PFCHSptD()   + jet2.QGTaggerAK4PFCHSptD()   + bjet.QGTaggerAK4PFCHSptD())/3.;
      //trijetAxis1_S		   = (jet1.axis1() + jet2.axis1() + bjet.axis1())/3.;
      trijetAxis2_S		   = (jet1.QGTaggerAK4PFCHSaxis2() + jet2.QGTaggerAK4PFCHSaxis2() + bjet.QGTaggerAK4PFCHSaxis2())/3.;
      trijetMult_S		   = (jet1.QGTaggerAK4PFCHSmult()  + jet2.QGTaggerAK4PFCHSmult()  + bjet.QGTaggerAK4PFCHSmult())/3.;
      trijetQGLikelihood_S	   = jet1.QGTaggerAK4PFCHSqgLikelihood()*jet2.QGTaggerAK4PFCHSqgLikelihood()*bjet.QGTaggerAK4PFCHSqgLikelihood();
      trijetQGLikelihood_avg_S = (jet1.QGTaggerAK4PFCHSqgLikelihood()+jet2.QGTaggerAK4PFCHSqgLikelihood()+bjet.QGTaggerAK4PFCHSqgLikelihood())/3.;
      trijetChiSquared_S	   = (TrijetP4.M() - mTopMass)/sigmaTrijet;
	  
      // Dijet Variables
      dijetPt_S		   = DijetP4.Pt();
      dijetEta_S		   = DijetP4.Eta();
      dijetPhi_S		   = DijetP4.Phi();
      dijetMass_S		   = DijetP4.M();
      dijetPtDR_S		   = DijetP4.Pt() * ROOT::Math::VectorUtil::DeltaR(jet1.p4(), jet2.p4());
      //dijetPFCharge_S	   = jet1.pfcharge()+jet2.pfcharge();
      dijetCombinedCvsL_S	   = (jet1.pfCombinedCvsLJetTags() + jet2.pfCombinedCvsLJetTags())/2.0;
      //dijetDeepCvsL_S	   = (jet1.pfDeepCSVCvsLJetTags() + jet2.pfDeepCSVCvsLJetTags())/2.0;
      dijetPtD_S		   = (jet1.QGTaggerAK4PFCHSptD() + jet2.QGTaggerAK4PFCHSptD())/2.0;
      //dijetAxis1_S		   = (jet1.axis1() + jet2.axis1())/2.;
      dijetAxis2_S		   = (jet1.QGTaggerAK4PFCHSaxis2() + jet2.QGTaggerAK4PFCHSaxis2())/2.;
      dijetMult_S		   = (jet1.QGTaggerAK4PFCHSmult()  + jet2.QGTaggerAK4PFCHSmult())/2.;
      dijetQGLikelihood_S	   = jet1.QGTaggerAK4PFCHSqgLikelihood()*jet2.QGTaggerAK4PFCHSqgLikelihood();
      dijetQGLikelihood_avg_S  = (jet1.QGTaggerAK4PFCHSqgLikelihood()+jet2.QGTaggerAK4PFCHSqgLikelihood())/2.;
      dijetChiSquared_S	   = (DijetP4.M() - mWMass)/sigmaDijet;
	  
      // Leading jet from dijet system
      LdgJetPt_S		    = jet1.pt();
      LdgJetEta_S		    = jet1.eta();
      LdgJetPhi_S		    = jet1.phi();
      LdgJetMass_S		    = jet1.p4().M();
      //LdgJetPFCharge_S	    = jet1.pfcharge();
      if (jet1.bjetDiscriminator() < 0.0)
	{
	  LdgJetBdisc_S	    = -1.0;
	}
      else
	{
	  LdgJetBdisc_S	    = jet1.bjetDiscriminator();
	}
      LdgJetCombinedCvsL_S	    = jet1.pfCombinedCvsLJetTags();
      //LdgJetDeepCvsL_S	    = jet1.pfDeepCSVCvsLJetTags();
      LdgJetPtD_S		    = jet1.QGTaggerAK4PFCHSptD();
      LdgJetAxis2_S		    = jet1.QGTaggerAK4PFCHSaxis2();
      //LdgJetAxis1_S		    = jet1.axis1();
      LdgJetMult_S		    = jet1.QGTaggerAK4PFCHSmult();
      LdgJetQGLikelihood_S	    = jet1.QGTaggerAK4PFCHSqgLikelihood();
      //LdgJetPullMagnitude_S	    = getJetPullMagnitude(jet1);
	  
      // Subleading jet from dijet system
      SubldgJetPt_S		    = jet2.pt();
      SubldgJetEta_S	    = jet2.eta();
      SubldgJetPhi_S	    = jet2.phi();
      SubldgJetMass_S	    = jet2.p4().M();
      //SubldgJetPFCharge_S	    = jet2.pfcharge();
      if (jet2.bjetDiscriminator() < 0.0)
	{
	  SubldgJetBdisc_S	    = -1.0;
	}
      else
	{
	  SubldgJetBdisc_S	    = jet2.bjetDiscriminator();
	}
      SubldgJetCombinedCvsL_S   = jet2.pfCombinedCvsLJetTags();
      //SubldgJetDeepCvsL_S	    = jet2.pfDeepCSVCvsLJetTags();
      SubldgJetPtD_S	    = jet2.QGTaggerAK4PFCHSptD();
      SubldgJetAxis2_S	    = jet2.QGTaggerAK4PFCHSaxis2();
      //SubldgJetAxis1_S	    = jet2.axis1();
      SubldgJetMult_S	    = jet2.QGTaggerAK4PFCHSmult();
      SubldgJetQGLikelihood_S   = jet2.QGTaggerAK4PFCHSqgLikelihood();
      //SubldgJetPullMagnitude_S  = getJetPullMagnitude(jet2);
	  
      // b-jet 
      bjetPt_S		    = bjet.pt();
      bjetEta_S		    = bjet.eta();
      bjetPhi_S		    = bjet.phi();
      bjetBdisc_S		    = bjet.bjetDiscriminator();
      bjetMass_S		    = bjet.p4().M();
      bjetQGLikelihood_S	    = bjet.QGTaggerAK4PFCHSqgLikelihood();
      bjetCombinedCvsL_S	    = bjet.pfCombinedCvsLJetTags();
      //bjetDeepCvsL_S	    = bjet.pfDeepCSVCvsLJetTags();
      bjetPtD_S		    = bjet.QGTaggerAK4PFCHSptD();
      bjetAxis2_S		    = bjet.QGTaggerAK4PFCHSaxis2();
      //bjetAxis1_S		    = bjet.axis1();
      bjetMult_S		    = bjet.QGTaggerAK4PFCHSmult();
      //bjetPFCharge_S	    = bjet.pfcharge();
      bjetLdgJetMass_S          = (bjet.p4() + jet1.p4()).M();
      bjetSubldgJetMass_S       = (bjet.p4() + jet2.p4()).M();
	  
      // Others (Soft-drop, distances, pull variables, etc)
      SoftDrop_n2_S		    = softDrop_n2;
      //PullAngleJ1J2_S	    = getJetPullAngle(jet1, jet2);  
      //PullAngleJ2J1_S	    = getJetPullAngle(jet2, jet1);
	  
      DEtaJ1withJ2_S	    = std::abs(jet1.eta() - jet2.eta());
      DEtaJ1withBJet_S	    = std::abs(jet1.eta() - bjet.eta());
      DEtaJ2withBJet_S	    = std::abs(jet2.eta() - bjet.eta());
      DEtaDijetwithBJet_S	    = std::abs(DijetP4.Eta() - bjet.eta());
      DEtaJ1BJetwithJ2_S	    = std::abs((jet1.p4() + bjet.p4()).Eta() - jet2.eta());
      DEtaJ2BJetwithJ1_S	    = std::abs((jet2.p4() + bjet.p4()).Eta() - jet1.eta());
      DPhiJ1withJ2_S	    = std::abs(ROOT::Math::VectorUtil::DeltaPhi(jet1.p4(), jet2.p4())); 
      DPhiJ1withBJet_S	    = std::abs(ROOT::Math::VectorUtil::DeltaPhi(jet1.p4(), bjet.p4()));
      DPhiJ2withBJet_S	    = std::abs(ROOT::Math::VectorUtil::DeltaPhi(jet2.p4(), bjet.p4()));
      DPhiDijetwithBJet_S	    = std::abs(ROOT::Math::VectorUtil::DeltaPhi(DijetP4, bjet.p4()));
      DPhiJ1BJetwithJ2_S	    = std::abs(ROOT::Math::VectorUtil::DeltaPhi((jet1.p4()+bjet.p4()), jet2.p4()));
      DPhiJ2BJetwithJ1_S	    = std::abs(ROOT::Math::VectorUtil::DeltaPhi((jet2.p4()+bjet.p4()), jet1.p4()));
      DRJ1withJ2_S		    = ROOT::Math::VectorUtil::DeltaR(jet1.p4(), jet2.p4());
      DRJ1withBJet_S	    = ROOT::Math::VectorUtil::DeltaR(jet1.p4(), bjet.p4());
      DRJ2withBJet_S	    = ROOT::Math::VectorUtil::DeltaR(jet2.p4(), bjet.p4());
      DRDijetwithBJet_S	    = ROOT::Math::VectorUtil::DeltaR(DijetP4, bjet.p4());
      DRJ1BJetwithJ2_S	    = ROOT::Math::VectorUtil::DeltaR((jet1.p4()+bjet.p4()), jet2.p4());
      DRJ2BJetwithJ1_S	    = ROOT::Math::VectorUtil::DeltaR((jet2.p4()+bjet.p4()), jet1.p4());
	  
      dijetMassOverTrijetMass_S = DijetP4.M()/TrijetP4.M();	  

      treeS -> Fill();
    }
    else{
	
      // Fill signal background
      eventWeight_B = fEventWeight.getWeight();
	  
      // Trijet Variables
      trijetPt_B		   = TrijetP4.Pt();
      trijetEta_B		   = TrijetP4.Eta();
      trijetPhi_B		   = TrijetP4.Phi();
      trijetMass_B		   = TrijetP4.M();
      trijetPtDR_B		   = TrijetP4.Pt() * ROOT::Math::VectorUtil::DeltaR(DijetP4, bjet.p4());
      //trijetPFCharge_B	   = jet1.pfcharge() + jet2.pfcharge() + bjet.pfcharge();
      trijetCombinedCvsL_B	   = (jet1.pfCombinedCvsLJetTags() + jet2.pfCombinedCvsLJetTags() + bjet.pfCombinedCvsLJetTags())/3.;
      //trijetDeepCvsL_B	   = (jet1.pfDeepCSVCvsLJetTags() + jet2.pfDeepCSVCvsLJetTags() + bjet.pfDeepCSVCvsLJetTags())/3.;
      trijetPtD_B		   = (jet1.QGTaggerAK4PFCHSptD()   + jet2.QGTaggerAK4PFCHSptD()   + bjet.QGTaggerAK4PFCHSptD())/3.;
      //trijetAxis1_B		   = (jet1.axis1() + jet2.axis1() + bjet.axis1())/3.;
      trijetAxis2_B		   = (jet1.QGTaggerAK4PFCHSaxis2() + jet2.QGTaggerAK4PFCHSaxis2() + bjet.QGTaggerAK4PFCHSaxis2())/3.;
      trijetMult_B		   = (jet1.QGTaggerAK4PFCHSmult()  + jet2.QGTaggerAK4PFCHSmult()  + bjet.QGTaggerAK4PFCHSmult())/3.;
      trijetQGLikelihood_B	   = jet1.QGTaggerAK4PFCHSqgLikelihood()*jet2.QGTaggerAK4PFCHSqgLikelihood()*bjet.QGTaggerAK4PFCHSqgLikelihood();
      trijetQGLikelihood_avg_B = (jet1.QGTaggerAK4PFCHSqgLikelihood()+jet2.QGTaggerAK4PFCHSqgLikelihood()+bjet.QGTaggerAK4PFCHSqgLikelihood())/3.;
      trijetChiSquared_B	   = (TrijetP4.M() - mTopMass)/sigmaTrijet;
	  
      // Dijet Variables
      dijetPt_B		   = DijetP4.Pt();
      dijetEta_B		   = DijetP4.Eta();
      dijetPhi_B		   = DijetP4.Phi();
      dijetMass_B		   = DijetP4.M();
      dijetPtDR_B		   = DijetP4.Pt() * ROOT::Math::VectorUtil::DeltaR(jet1.p4(), jet2.p4());
      //dijetPFCharge_B	   = jet1.pfcharge()+jet2.pfcharge();
      dijetCombinedCvsL_B	   = (jet1.pfCombinedCvsLJetTags() + jet2.pfCombinedCvsLJetTags())/2.0;
      //dijetDeepCvsL_B	   = (jet1.pfDeepCSVCvsLJetTags() + jet2.pfDeepCSVCvsLJetTags())/2.0;
      dijetPtD_B		   = (jet1.QGTaggerAK4PFCHSptD() + jet2.QGTaggerAK4PFCHSptD())/2.0;
      //dijetAxis1_B		   = (jet1.axis1() + jet2.axis1())/2.;
      dijetAxis2_B		   = (jet1.QGTaggerAK4PFCHSaxis2() + jet2.QGTaggerAK4PFCHSaxis2())/2.;
      dijetMult_B		   = (jet1.QGTaggerAK4PFCHSmult()  + jet2.QGTaggerAK4PFCHSmult())/2.;
      dijetQGLikelihood_B	   = jet1.QGTaggerAK4PFCHSqgLikelihood()*jet2.QGTaggerAK4PFCHSqgLikelihood();
      dijetQGLikelihood_avg_B  = (jet1.QGTaggerAK4PFCHSqgLikelihood()+jet2.QGTaggerAK4PFCHSqgLikelihood())/2.;
      dijetChiSquared_B	   = (DijetP4.M() - mWMass)/sigmaDijet;
	  
      // Leading jet from dijet system
      LdgJetPt_B		    = jet1.pt();
      LdgJetEta_B		    = jet1.eta();
      LdgJetPhi_B		    = jet1.phi();
      LdgJetMass_B		    = jet1.p4().M();
      //LdgJetPFCharge_B	    = jet1.pfcharge();
      if (jet1.bjetDiscriminator() < 0.0)
	{
	  LdgJetBdisc_B	    = -1.0;
	}
      else
	{
	  LdgJetBdisc_B	    = jet1.bjetDiscriminator();
	}
      LdgJetCombinedCvsL_B	    = jet1.pfCombinedCvsLJetTags();
      //LdgJetDeepCvsL_B	    = jet1.pfDeepCSVCvsLJetTags();
      LdgJetPtD_B		    = jet1.QGTaggerAK4PFCHSptD();
      LdgJetAxis2_B		    = jet1.QGTaggerAK4PFCHSaxis2();
      //LdgJetAxis1_B		    = jet1.axis1();
      LdgJetMult_B		    = jet1.QGTaggerAK4PFCHSmult();
      LdgJetQGLikelihood_B	    = jet1.QGTaggerAK4PFCHSqgLikelihood();
      //LdgJetPullMagnitude_B	    = getJetPullMagnitude(jet1);
	  
      // Subleading jet from dijet system
      SubldgJetPt_B		    = jet2.pt();
      SubldgJetEta_B	    = jet2.eta();
      SubldgJetPhi_B	    = jet2.phi();
      SubldgJetMass_B	    = jet2.p4().M();
      //SubldgJetPFCharge_B	    = jet2.pfcharge();
      if (jet2.bjetDiscriminator() < 0.0)
	{
	  SubldgJetBdisc_B	    = -1.0;
	}
      else
	{
	  SubldgJetBdisc_B	    = jet2.bjetDiscriminator();
	}
      SubldgJetCombinedCvsL_B   = jet2.pfCombinedCvsLJetTags();
      //SubldgJetDeepCvsL_B	    = jet2.pfDeepCSVCvsLJetTags();
      SubldgJetPtD_B	    = jet2.QGTaggerAK4PFCHSptD();
      SubldgJetAxis2_B	    = jet2.QGTaggerAK4PFCHSaxis2();
      //SubldgJetAxis1_B	    = jet2.axis1();
      SubldgJetMult_B	    = jet2.QGTaggerAK4PFCHSmult();
      SubldgJetQGLikelihood_B   = jet2.QGTaggerAK4PFCHSqgLikelihood();
      //SubldgJetPullMagnitude_B  = getJetPullMagnitude(jet2);
	  
      // b-jet 
      bjetPt_B		    = bjet.pt();
      bjetEta_B		    = bjet.eta();
      bjetPhi_B		    = bjet.phi();
      bjetBdisc_B		    = bjet.bjetDiscriminator();
      bjetMass_B		    = bjet.p4().M();
      bjetQGLikelihood_B	    = bjet.QGTaggerAK4PFCHSqgLikelihood();
      bjetCombinedCvsL_B	    = bjet.pfCombinedCvsLJetTags();
      //bjetDeepCvsL_B	    = bjet.pfDeepCSVCvsLJetTags();
      bjetPtD_B		    = bjet.QGTaggerAK4PFCHSptD();
      bjetAxis2_B		    = bjet.QGTaggerAK4PFCHSaxis2();
      //bjetAxis1_B		    = bjet.axis1();
      bjetMult_B		    = bjet.QGTaggerAK4PFCHSmult();
      //bjetPFCharge_B	    = bjet.pfcharge();
      bjetLdgJetMass_B          = (bjet.p4() + jet1.p4()).M();
      bjetSubldgJetMass_B       = (bjet.p4() + jet2.p4()).M();
	  
      // Others (Soft-drop, distances, pull variables, etc)
      SoftDrop_n2_B		    = softDrop_n2;
      //PullAngleJ1J2_B	    = getJetPullAngle(jet1, jet2);  
      //PullAngleJ2J1_B	    = getJetPullAngle(jet2, jet1);
	  
      DEtaJ1withJ2_B	    = std::abs(jet1.eta() - jet2.eta());
      DEtaJ1withBJet_B	    = std::abs(jet1.eta() - bjet.eta());
      DEtaJ2withBJet_B	    = std::abs(jet2.eta() - bjet.eta());
      DEtaDijetwithBJet_B	    = std::abs(DijetP4.Eta() - bjet.eta());
      DEtaJ1BJetwithJ2_B	    = std::abs((jet1.p4() + bjet.p4()).Eta() - jet2.eta());
      DEtaJ2BJetwithJ1_B	    = std::abs((jet2.p4() + bjet.p4()).Eta() - jet1.eta());
      DPhiJ1withJ2_B	    = std::abs(ROOT::Math::VectorUtil::DeltaPhi(jet1.p4(), jet2.p4())); 
      DPhiJ1withBJet_B	    = std::abs(ROOT::Math::VectorUtil::DeltaPhi(jet1.p4(), bjet.p4()));
      DPhiJ2withBJet_B	    = std::abs(ROOT::Math::VectorUtil::DeltaPhi(jet2.p4(), bjet.p4()));
      DPhiDijetwithBJet_B	    = std::abs(ROOT::Math::VectorUtil::DeltaPhi(DijetP4, bjet.p4()));
      DPhiJ1BJetwithJ2_B	    = std::abs(ROOT::Math::VectorUtil::DeltaPhi((jet1.p4()+bjet.p4()), jet2.p4()));
      DPhiJ2BJetwithJ1_B	    = std::abs(ROOT::Math::VectorUtil::DeltaPhi((jet2.p4()+bjet.p4()), jet1.p4()));
      DRJ1withJ2_B		    = ROOT::Math::VectorUtil::DeltaR(jet1.p4(), jet2.p4());
      DRJ1withBJet_B	    = ROOT::Math::VectorUtil::DeltaR(jet1.p4(), bjet.p4());
      DRJ2withBJet_B	    = ROOT::Math::VectorUtil::DeltaR(jet2.p4(), bjet.p4());
      DRDijetwithBJet_B	    = ROOT::Math::VectorUtil::DeltaR(DijetP4, bjet.p4());
      DRJ1BJetwithJ2_B	    = ROOT::Math::VectorUtil::DeltaR((jet1.p4()+bjet.p4()), jet2.p4());
      DRJ2BJetwithJ1_B	    = ROOT::Math::VectorUtil::DeltaR((jet2.p4()+bjet.p4()), jet1.p4());
	  
      dijetMassOverTrijetMass_B = DijetP4.M()/TrijetP4.M();	  

      treeB -> Fill();
	
    }            
  }

  //========================================================================================================
  // Pseudo-matching: Used to calculate the dPt/Pt cut
  //========================================================================================================  
  for (size_t i=0; i<GenTops_Quarks.size(); i++)
    {
      genParticle Quark      = GenTops_Quarks.at(i);
      Jet JetClosest;
      double dRmin  = 99999.9;
      for (auto& jet: jetData.getSelectedJets())
	{
	  double dR  = ROOT::Math::VectorUtil::DeltaR( jet.p4(), Quark.p4());
	  if (dR > cfg_DeltaRCut)   continue;
	  if (dR > dRmin)   continue;
	  dRmin  = dR;
	  JetClosest = jet;
	}
      if (dRmin > cfg_DeltaRCut) continue;
      double dPtOverPt = (JetClosest.pt() - Quark.pt())/Quark.pt();
      hQuarkJetMinDr03_DeltaPtOverPt            -> Fill (dPtOverPt);   //2sigma = 2*0.16 = 0.32
      hQuarkJetMinDr03_DeltaPtOverPt_vs_QuarkPt -> Fill(dPtOverPt, Quark.pt());
      hQuarkJetMinDr03_DeltaPtOverPt_vs_DeltaRmin -> Fill(dPtOverPt, dRmin);

      if (Quark.pt() < 40)       hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt0To40GeV    -> Fill(dPtOverPt);
      else if (Quark.pt() < 60)  hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt40To60GeV   -> Fill(dPtOverPt);
      else if (Quark.pt() < 80)  hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt60To80GeV   -> Fill(dPtOverPt);
      else if (Quark.pt() < 100) hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt80To100GeV  -> Fill(dPtOverPt);
      else if (Quark.pt() < 120) hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt100To120GeV -> Fill(dPtOverPt);
      else if (Quark.pt() < 140) hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt120To140GeV -> Fill(dPtOverPt);
      else if (Quark.pt() < 160) hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt140To160GeV -> Fill(dPtOverPt);
      else if (Quark.pt() < 180) hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt160To180GeV -> Fill(dPtOverPt);
      else if (Quark.pt() < 200) hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt180To200GeV -> Fill(dPtOverPt);
      else hQuarkJetMinDr03_DeltaPtOverPt_QuarkPt200ToInfGeV -> Fill(dPtOverPt);
      hQuarkJetMinDr03_DeltaR_vs_QuarkPt        -> Fill(dRmin, Quark.pt());
    }

  //========================================================================================================
  // Matching:  Minimum DeltaR and DeltaPt check
  //========================================================================================================
  
  //Sanity chack
  double DeltaR_min = 999.999;
  for (auto& jet1: jetData.getSelectedJets())
    {
      for (auto& jet2: jetData.getSelectedJets())
	{
	  if (areSameJets(jet1, jet2)) continue;
	  double DeltaR = ROOT::Math::VectorUtil::DeltaR( jet1.p4(), jet2.p4());
	  if (DeltaR > DeltaR_min) continue;
	  DeltaR_min = DeltaR;
	}
    }

  hJetsDeltaRmin -> Fill(DeltaR_min);

  //========================================================================================================
  // Check Properties of matched objects
  //========================================================================================================
  for (size_t j=0; j<MCtrue_Bjet.size(); j++){

    hTrijetDrMin -> Fill( ROOT::Math::VectorUtil::DeltaR(MCtrue_Bjet.at(j).p4(),      MGen_Bjet.at(j).p4()));
    hTrijetDrMin -> Fill( ROOT::Math::VectorUtil::DeltaR(MCtrue_LdgJet.at(j).p4(),    MGen_LdgJet.at(j).p4()));
    hTrijetDrMin -> Fill( ROOT::Math::VectorUtil::DeltaR(MCtrue_SubldgJet.at(j).p4(), MGen_SubldgJet.at(j).p4()));
      
    hTrijetDPtOverGenPt -> Fill( (MCtrue_Bjet.at(j).pt()      - MGen_Bjet.at(j).pt())/MGen_Bjet.at(j).pt());
    hTrijetDPtOverGenPt -> Fill( (MCtrue_LdgJet.at(j).pt()    - MGen_LdgJet.at(j).pt())/MGen_LdgJet.at(j).pt());
    hTrijetDPtOverGenPt -> Fill( (MCtrue_SubldgJet.at(j).pt() - MGen_SubldgJet.at(j).pt())/MGen_SubldgJet.at(j).pt());
      
    hTrijetDPt_matched -> Fill( (MCtrue_Bjet.at(j).pt()      - MGen_Bjet.at(j).pt()));
    hTrijetDPt_matched -> Fill( (MCtrue_LdgJet.at(j).pt()    - MGen_LdgJet.at(j).pt()));
    hTrijetDPt_matched -> Fill( (MCtrue_SubldgJet.at(j).pt() - MGen_SubldgJet.at(j).pt()));
      
    hTrijetDEtaOverGenEta -> Fill( (MCtrue_Bjet.at(j).eta()      - MGen_Bjet.at(j).eta())/MGen_Bjet.at(j).eta());
    hTrijetDEtaOverGenEta -> Fill( (MCtrue_LdgJet.at(j).eta()    - MGen_LdgJet.at(j).eta())/MGen_LdgJet.at(j).eta());
    hTrijetDEtaOverGenEta -> Fill( (MCtrue_SubldgJet.at(j).eta() - MGen_SubldgJet.at(j).eta())/MGen_SubldgJet.at(j).eta());
      
    hTrijetDEta_matched -> Fill( (MCtrue_Bjet.at(j).eta()      - MGen_Bjet.at(j).eta()));
    hTrijetDEta_matched -> Fill( (MCtrue_LdgJet.at(j).eta()    - MGen_LdgJet.at(j).eta()));
    hTrijetDEta_matched -> Fill( (MCtrue_SubldgJet.at(j).eta() - MGen_SubldgJet.at(j).eta()));
      
    hTrijetDPhiOverGenPhi -> Fill( (ROOT::Math::VectorUtil::DeltaPhi(MCtrue_Bjet.at(j).p4(),      MGen_Bjet.at(j).p4()))/MGen_Bjet.at(j).phi());
    hTrijetDPhiOverGenPhi -> Fill( (ROOT::Math::VectorUtil::DeltaPhi(MCtrue_LdgJet.at(j).p4(),    MGen_LdgJet.at(j).p4()))/MGen_LdgJet.at(j).phi());
    hTrijetDPhiOverGenPhi -> Fill( (ROOT::Math::VectorUtil::DeltaPhi(MCtrue_SubldgJet.at(j).p4(), MGen_SubldgJet.at(j).p4()))/MGen_SubldgJet.at(j).phi());
      
    hTrijetDPhi_matched -> Fill( (ROOT::Math::VectorUtil::DeltaPhi(MCtrue_Bjet.at(j).p4(),      MGen_Bjet.at(j).p4())));
    hTrijetDPhi_matched -> Fill( (ROOT::Math::VectorUtil::DeltaPhi(MCtrue_LdgJet.at(j).p4(),    MGen_LdgJet.at(j).p4())));
    hTrijetDPhi_matched -> Fill( (ROOT::Math::VectorUtil::DeltaPhi(MCtrue_SubldgJet.at(j).p4(), MGen_SubldgJet.at(j).p4())));
      
  }

  //========================================================================================================
  // Sanity check
  //========================================================================================================
    
  for (auto& jet: jetData.getSelectedJets()){
    hAllJetCvsL         -> Fill(jet.pfCombinedCvsLJetTags());
    hAllJetPtD          -> Fill(jet.QGTaggerAK4PFCHSptD());
    hAllJetAxis2        -> Fill(jet.QGTaggerAK4PFCHSaxis2());
    hAllJetMult         -> Fill(jet.QGTaggerAK4PFCHSmult());
    hAllJetBdisc        -> Fill(jet.bjetDiscriminator());
    hAllJetQGLikelihood -> Fill(jet.QGTaggerAK4PFCHSqgLikelihood());
  }
  
  vector<genParticle> GenCharm = GetGenParticles(fEvent.genparticles().getGenParticles(), 4);
  vector<Jet> CJets;

  for (size_t i=0; i<GenCharm.size(); i++){
    double dRmin = 10000.0;

    // double dPtOverPtmin = 10000.0;

    Jet mcMatched_CJet;
    for (auto& jet: jetData.getSelectedJets()){
      double dR  = ROOT::Math::VectorUtil::DeltaR( jet.p4(), GenCharm.at(i).p4());
      double dPtOverPt = std::abs((jet.pt() - GenCharm.at(i).pt())/ GenCharm.at(i).pt());
      if (dR > cfg_DeltaRCut) continue;
      if (dR > dRmin) continue;
      //if (dPtOverPt > dPtOverPtmin) continue;
      if (dPtOverPt > cfg_DeltaPtOverPtCut) continue;
      dRmin = dR;
      // dPtOverPtmin = dPtOverPt;
      mcMatched_CJet = jet;
    }
    if (dRmin < cfg_DeltaRCut){
      hCJetCvsL         -> Fill(mcMatched_CJet.pfCombinedCvsLJetTags());
      hCJetPtD          -> Fill(mcMatched_CJet.QGTaggerAK4PFCHSptD());
      hCJetAxis2        -> Fill(mcMatched_CJet.QGTaggerAK4PFCHSaxis2());
      hCJetMult         -> Fill(mcMatched_CJet.QGTaggerAK4PFCHSmult());
      hCJetBdisc        -> Fill(mcMatched_CJet.bjetDiscriminator());
    }
  }
  
  //================================================================================================
  // Finalize
  //================================================================================================
  fEventSaver.save();
}

  
