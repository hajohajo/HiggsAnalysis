// -*- c++ -*-
#ifndef EventSelection_TopSelectionBDT_h
#define EventSelection_TopSelectionBDT_h

#include "EventSelection/interface/BaseSelection.h"
#include "EventSelection/interface/JetSelection.h"
#include "EventSelection/interface/BJetSelection.h"
#include "EventSelection/interface/TopTagSFCalculator.h"
#include "DataFormat/interface/Jet.h"
#include "Framework/interface/EventCounter.h"
#include "Tools/interface/DirectionalCut.h"

#include "Tools/interface/MCTools.h"

#include <string>
#include <vector>

#include <TDirectory.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TROOT.h>

#include <TMVA/Factory.h>
#include <TMVA/Tools.h>
#include <TMVA/TMVAGui.h>
#include <TMVA/Reader.h>

// Forward declarations
class ParameterSet;
class CommonPlots;
class Event;
class EventCounter;
class HistoWrapper;
class WrappedTH1;
class WrappedTH2;


struct TrijetSelection{
  std::vector<Jet> Jet1;
  std::vector<Jet> Jet2;
  std::vector<Jet> BJet;
  std::vector<double> MVA;
  std::vector<math::XYZTLorentzVector> TrijetP4;
  std::vector<math::XYZTLorentzVector> DijetP4; 
  std::vector<bool> isGenuine;
  std::vector<bool> isTagged;
};

struct SelectedTrijets{
  Jet Jet1;
  Jet Jet2;
  Jet BJet;
  double MVA;
  math::XYZTLorentzVector TrijetP4;
  math::XYZTLorentzVector DijetP4;
  bool isGenuine;
  bool isTagged;
};


class TopSelectionBDT: public BaseSelection {
public:
  class Data {
  public:
    // The reason for pointer instead of reference is that const
    // reference allows temporaries, while const pointer does not.
    // Here the object pointed-to must live longer than this object.
    Data();
    ~Data();

    // Status of passing event selection
    bool passedSelection() const { return bPassedSelection; }    

    // Trijet
    const float getTopMVA() const { return fTopMVA; }
    const Jet getTopJet1() const { return fTopJet1; } 
    const Jet getTopJet2() const { return fTopJet2; } 
    const Jet getTopBJet() const { return fTopBJet; } 
    //    const math::XYZTLorentzVector getTopDijetP4() const {return fTopDijet_p4; }
    const math::XYZTLorentzVector getTopDijet() const {return fTopDijet_p4; }
    const math::XYZTLorentzVector getTop() const {return fTop_p4; }

    const std::vector<Jet>& getSelectedTopsJet1() const { return fSelectedTopsJet1; }
    const std::vector<Jet>& getSelectedTopsJet2() const { return fSelectedTopsJet2; }
    const std::vector<Jet>& getSelectedTopsBJet() const { return fSelectedTopsBJet; }
    const std::vector<float>& getSelectedTopsMVA() const { return fSelectedTopsMVA; }
    const size_t getSelectedTopsSize() const { return fSelectedTopsMVA.size(); }

    const std::vector<Jet>& getNotSelectedTopsJet1() const { return fNotSelectedTopsJet1; }
    const std::vector<Jet>& getNotSelectedTopsJet2() const { return fNotSelectedTopsJet2; }
    const std::vector<Jet>& getNotSelectedTopsBJet() const { return fNotSelectedTopsBJet; }
    const std::vector<float>& getNotSelectedTopsMVA() const { return fNotSelectedTopsMVA; }
    const size_t getNotSelectedTopsSize() const { return fNotSelectedTopsMVA.size(); }

    const std::vector<Jet>& getAllTopsJet1() const { return fAllTopsJet1; }
    const std::vector<Jet>& getAllTopsJet2() const { return fAllTopsJet2; }
    const std::vector<Jet>& getAllTopsBJet() const { return fAllTopsBJet; }
    const std::vector<float>& getAllTopsMVA() const { return fAllTopsMVA; }
    const size_t getAllTopsSize() const { return fAllTopsMVA.size(); }

    const std::vector<Jet>& getSelectedCleanedTopsJet1() const { return fSelectedCleanedTopsJet1; }
    const std::vector<Jet>& getSelectedCleanedTopsJet2() const { return fSelectedCleanedTopsJet2; }
    const std::vector<Jet>& getSelectedCleanedTopsBJet() const { return fSelectedCleanedTopsBJet; }
    const std::vector<float>& getSelectedCleanedTopsMVA() const { return fSelectedCleanedTopsMVA; }
    const size_t getSelectedCleanedTopsSize() const { return fSelectedCleanedTopsMVA.size(); }

    const std::vector<Jet>& getAllCleanedTopsJet1() const { return fAllCleanedTopsJet1; }
    const std::vector<Jet>& getAllCleanedTopsJet2() const { return fAllCleanedTopsJet2; }
    const std::vector<Jet>& getAllCleanedTopsBJet() const { return fAllCleanedTopsBJet; }
    const std::vector<float>& getAllCleanedTopsMVA() const { return fAllCleanedTopsMVA; }
    const size_t getAllCleanedTopsSize() const { return fAllCleanedTopsMVA.size(); }
    
    const double getTopMassWMassRatio() const
    { 
      double R = fTop_p4.mass()/fTopDijet_p4.mass();
      return R;
    }


    /// Obtain the b-tagging event weight 
    const double getTopTaggingScaleFactorEventWeight() const { return fTopTaggingScaleFactorEventWeight; }

    friend class TopSelectionBDT;

  private:
    /// Boolean for passing selection
    bool bPassedSelection;
    std::vector<Jet> fJetsUsedAsBJets;
    std::vector<Jet> fFailedBJetsUsedAsBJets;
    /// Trijet-1
    float fTopMVA;
    Jet fTopJet1;
    Jet fTopJet2;
    Jet fTopBJet;
    math::XYZTLorentzVector fTopDijet_p4;
    math::XYZTLorentzVector fTop_p4;

    std::vector<Jet> fSelectedTopsJet1;
    std::vector<Jet> fSelectedTopsJet2;
    std::vector<Jet> fSelectedTopsBJet;
    std::vector<float> fSelectedTopsMVA;

    std::vector<Jet> fNotSelectedTopsJet1;
    std::vector<Jet> fNotSelectedTopsJet2;
    std::vector<Jet> fNotSelectedTopsBJet;
    std::vector<float> fNotSelectedTopsMVA;

    std::vector<Jet> fAllTopsJet1;
    std::vector<Jet> fAllTopsJet2;
    std::vector<Jet> fAllTopsBJet;
    std::vector<float> fAllTopsMVA;

    std::vector<Jet> fSelectedCleanedTopsJet1;
    std::vector<Jet> fSelectedCleanedTopsJet2;
    std::vector<Jet> fSelectedCleanedTopsBJet;
    std::vector<float> fSelectedCleanedTopsMVA;

    std::vector<Jet> fAllCleanedTopsJet1;
    std::vector<Jet> fAllCleanedTopsJet2;
    std::vector<Jet> fAllCleanedTopsBJet;
    std::vector<float> fAllCleanedTopsMVA;

    // top-tagging scale factor event weight
    double fTopTaggingScaleFactorEventWeight;

  };
  
  // Main class
  /// Constructor with histogramming
  explicit TopSelectionBDT(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& postfix = "");
  /// Constructor without histogramming
  explicit TopSelectionBDT(const ParameterSet& config);
  virtual ~TopSelectionBDT();

  virtual void bookHistograms(TDirectory* dir);
  
  /// Use silentAnalyze if you do not want to fill histograms or increment counters
  Data silentAnalyze(const Event& event, const JetSelection::Data& jetData, const BJetSelection::Data& bjetData);
  /// analyze does fill histograms and incrementes counters
  Data analyze(const Event& event, const JetSelection::Data& jetData, const BJetSelection::Data& bjetData);

  TMVA::Reader *reader;
  
  Float_t TrijetPtDR;
  Float_t TrijetDijetPtDR;
  Float_t TrijetBjetMass;
  Float_t TrijetLdgJetBDisc;
  Float_t TrijetSubldgJetBDisc;
  Float_t TrijetBJetLdgJetMass;
  Float_t TrijetBJetSubldgJetMass;
  Float_t TrijetMass;
  Float_t TrijetDijetMass;
  Float_t TrijetBJetBDisc;
  Float_t TrijetSoftDrop_n2;
  Float_t TrijetLdgJetCvsL;    
  Float_t TrijetSubldgJetCvsL;
  Float_t TrijetLdgJetPtD;
  Float_t TrijetSubldgJetPtD;
  Float_t TrijetLdgJetAxis2;
  Float_t TrijetSubldgJetAxis2;
  Float_t TrijetLdgJetMult;
  Float_t TrijetSubldgJetMult;
  Float_t TrijetLdgJetQGLikelihood;
  Float_t TrijetSubldgJetQGLikelihood;

private:
  /// Initialisation called from constructor
  void initialize(const ParameterSet& config);
  /// The actual selection
  Data privateAnalyze(const Event& event, const std::vector<Jet> selectedJets, const std::vector<Jet> selectedBjets);
  /// Returns true if the two jets are the same
  bool areSameJets(const Jet& jet1, const Jet& jet2);
  /// Return true if a selected jet matches a selected bjet
  bool isBJet(const Jet& jet1, const std::vector<Jet>& bjets);
  /// Determine if top candidate is MC matched
  bool _getIsMatchedTop(bool isMC, Jet bjet, Jet jet1, Jet jet2, TrijetSelection mcTrueTrijets);
  Jet getLeadingSubleadingJet(const Jet& jet0, const Jet& jet1, string selectedJet);
  bool isMatchedJet(const Jet& jet, const TrijetSelection& myTops, const unsigned int index);
  TrijetSelection SortInMVAvalue(TrijetSelection TopCand);
  SelectedTrijets GetSelectedTopCandidate(TrijetSelection TopCand, int index);
  bool TopIsCrossCleaned(int Index, TrijetSelection TopCand, const std::vector<Jet>& bjets);
  vector<genParticle> GetGenParticles(const vector<genParticle> genParticles, const int pdgId);
  const genParticle GetLastCopy(const vector<genParticle> genParticles, const genParticle &p);
  genParticle getLeadingSubleadingParton(const genParticle& quark0, const genParticle& quark1, string selectedParton);
  vector<genParticle> GetWpartons( genParticle daughter, const Event& event);
  bool IsGenuineTop(const Jet& trijetJet1, const Jet& trijetJet2, const Jet& trijetBJet,
		    const std::vector<Jet>& MCtrue_LdgJet,  const std::vector<Jet>& MCtrue_SubldgJet, const std::vector<Jet>& MCtrue_Bjet);

 
  // Input parameters
  const DirectionalCut<int> cfg_NumberOfTopsCut;
  const DirectionalCut<double> cfg_TopMVACut;
  const DirectionalCut<double> cfg_TopMassLowCut;
  const DirectionalCut<double> cfg_TopMassUppCut;
  const DirectionalCut<double> cfg_CSV_bDiscCut;

  // Event counter for passing selection
  Count cPassedTopSelectionBDT;

  // Sub counters
  Count cSubAll;
  Count cSubPassedBjetsCut;
  Count cSubPassedMVACut;
  Count cSubPassedNTopsCut;
  //
  Count cTopsAll;
  Count cTopsPassTopMassLowCut;
  Count cTopsPassTopMassUppCut;
  Count cTopsPassBDiscCut;
  Count cTopsPassBDTCut;
  Count cTopsPassCrossCleanCut;

  // Scalefactor calculator
  TopTagSFCalculator fTopTagSFCalculator;

  // Histograms (1D)
  WrappedTH1 *hTopBDT_AllCandidates;
  WrappedTH1 *hTopMass_AllCandidates;
  WrappedTH1 *hTopPt_AllCandidates;
  WrappedTH1 *hTopMultiplicity_AllCandidates;

  WrappedTH1 *hTopBDT_SelectedCandidates;
  WrappedTH1 *hTopMass_SelectedCandidates;
  WrappedTH1 *hTopPt_SelectedCandidates;
  WrappedTH1 *hTopMultiplicity_SelectedCandidates;

  WrappedTH1 *hTopBDT_SelectedCleanedCandidates;
  WrappedTH1 *hTopMass_SelectedCleanedCandidates;
  WrappedTH1 *hTopPt_SelectedCleanedCandidates;
  WrappedTH1 *hTopMultiplicity_SelectedCleanedCandidates;

  WrappedTH1 *hTopBDT_NotSelectedCandidates;
  WrappedTH1 *hTopMass_NotSelectedCandidates;
  WrappedTH1 *hTopPt_NotSelectedCandidates;
  WrappedTH1 *hTopMultiplicity_NotSelectedCandidates;

  WrappedTH1 *hTopBDT_AllCleanedCandidates;
  WrappedTH1 *hTopMass_AllCleanedCandidates;
  WrappedTH1 *hTopPt_AllCleanedCandidates;
  WrappedTH1 *hTopMultiplicity_AllCleanedCandidates;

  WrappedTH1  *hTopPt;
  WrappedTH1  *hTopMass;
  WrappedTH1  *hTopJet1Pt;
  WrappedTH1  *hTopJet1Eta;
  WrappedTH1  *hTopJet1BDisc;
  WrappedTH1  *hTopJet2Pt;
  WrappedTH1  *hTopJet2Eta;
  WrappedTH1  *hTopJet2BDisc;
  WrappedTH1  *hTopBJetPt;
  WrappedTH1  *hTopBJetEta;
  WrappedTH1  *hTopBJetBDisc;
  WrappedTH1  *hTopDiJetPt;
  WrappedTH1  *hTopDiJetEta;
  WrappedTH1  *hTopDiJetMass;
  WrappedTH1  *hTopDijetDeltaR;
  WrappedTH1  *hTopMassWMassRatio;
  
};

#endif
