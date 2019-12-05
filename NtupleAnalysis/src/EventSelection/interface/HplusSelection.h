// -*- c++ -*-
#ifndef EventSelection_HplusSelection_h
#define EventSelection_HplusSelection_h

#include "EventSelection/interface/BaseSelection.h"
#include "EventSelection/interface/JetSelection.h"
#include "EventSelection/interface/BJetSelection.h"
#include "EventSelection/interface/TopSelectionMVA.h"
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

class HplusSelection: public BaseSelection {
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
    // Status of passing presence of any two tops (with given MVA cut) and free b-jet
    bool passedAnyTwoTopsAndFreeB() const { return bHasTwoTopsAndFreeB;}
    // Status of passing presence of at least n free bjets
    bool hasFreeBJet() const {return bPass_FreeBjet;}
    const size_t getAllCleanedTopsSize() const { return nAllCleanedTops;}
    
    // Trijet-1
    const float getMVAmax1() const { return fMVAmax1; }
    const Jet getTrijet1Jet1() const { return fTrijet1Jet1; } 
    const Jet getTrijet1Jet2() const { return fTrijet1Jet2; } 
    const Jet getTrijet1BJet() const { return fTrijet1BJet; } 
    const math::XYZTLorentzVector getTrijet1DijetP4() const {return fTrijet1Dijet_p4; }
    const math::XYZTLorentzVector getTriJet1() const {return fTrijet1_p4; }
    // Trijet-2
    const float getMVAmax2() const { return fMVAmax2; }
    const Jet getTrijet2Jet1() const { return fTrijet2Jet1; } 
    const Jet getTrijet2Jet2() const { return fTrijet2Jet2; } 
    const Jet getTrijet2BJet() const { return fTrijet2BJet; }     

    // Leading Trijet
    const math::XYZTLorentzVector getLdgTrijet() const
    {       
      if (fTrijet1_p4.pt() > fTrijet2_p4.pt()) return fTrijet1_p4;
      else return fTrijet2_p4;
    }
    const float getMVALdgInPt() const
    {
      if (fTrijet1_p4.pt() > fTrijet2_p4.pt()) return fMVAmax1;
      else return fMVAmax2;
    }
    const Jet getLdgTrijetBJet() const
    {
      if (fTrijet1_p4.pt() > fTrijet2_p4.pt()) return fTrijet1BJet;
      else return fTrijet2BJet;
    }     
    const Jet getLdgTrijetJet1() const
    {
      if (fTrijet1_p4.pt() > fTrijet2_p4.pt()) return fTrijet1Jet1;
      else return fTrijet2Jet1;
    }
    const Jet getLdgTrijetJet2() const
    {
      if (fTrijet1_p4.pt() > fTrijet2_p4.pt()) return fTrijet1Jet2;
      else return fTrijet2Jet2;
    }     
    const math::XYZTLorentzVector getLdgTrijetDijet() const
    {
      if (fTrijet1_p4.pt() > fTrijet2_p4.pt()) return fTrijet1Dijet_p4;
      else return fTrijet2Dijet_p4;
    }
    //Subleading Trijet
    const math::XYZTLorentzVector getSubldgTrijet() const
    {
      if (fTrijet1_p4.pt() > fTrijet2_p4.pt()) return fTrijet2_p4;
      else return fTrijet1_p4;
    }
    const float getMVASubldgInPt() const
    {
      if (fTrijet1_p4.pt() > fTrijet2_p4.pt()) return fMVAmax2;
      else return fMVAmax1;
    }    
    const Jet getSubldgTrijetBJet() const
    {
      if (fTrijet1_p4.pt() < fTrijet2_p4.pt()) return fTrijet1BJet;
      else return fTrijet2BJet;
    }
    const Jet getSubldgTrijetJet1() const
    {
      if (fTrijet1_p4.pt() < fTrijet2_p4.pt()) return fTrijet1Jet1;
      else return fTrijet2Jet1;
    }
    const Jet getSubldgTrijetJet2() const
    {
      if (fTrijet1_p4.pt() < fTrijet2_p4.pt()) return fTrijet1Jet2;
      else return fTrijet2Jet2;
    }
    const math::XYZTLorentzVector getSubldgTrijetDijet() const
    {
      if (fTrijet1_p4.pt() > fTrijet2_p4.pt()) return fTrijet2Dijet_p4;
      else return fTrijet1Dijet_p4;
    }

    //Tetrajet
    const math::XYZTLorentzVector getLdgTetrajet() const {return fLdgTetrajet_p4;} // uses ldg-trijet and tetrajetBjet (NOT the tetrajet with largest pt)
    const Jet getTetrajetBJet() const {return fTetrajetBJet;}

    /// Obtain the top-tagging event weight     
    const double getTopTaggingScaleFactorEventWeight() const { return fTopTaggingScaleFactorEventWeight; }

    friend class HplusSelection;

  private:
    /// Boolean for passing selection
    bool bPassedSelection;
    bool bHasTwoTopsAndFreeB;
    bool bPass_FreeBjet;
    size_t nAllCleanedTops;
    float fMVAmax1;
    float fMVAmax2;
    math::XYZTLorentzVector fTrijet1_p4;
    math::XYZTLorentzVector fTrijet2_p4;
    Jet fTrijet1BJet;
    Jet fTrijet1Jet1;
    Jet fTrijet1Jet2;
    Jet fTrijet2BJet;
    Jet fTrijet2Jet1;
    Jet fTrijet2Jet2;
    math::XYZTLorentzVector fTrijet1Dijet_p4;
    math::XYZTLorentzVector fTrijet2Dijet_p4;
    math::XYZTLorentzVector fLdgTetrajet_p4;
    Jet fTetrajetBJet;
    // top-tagging scale factor event weight
    double fTopTaggingScaleFactorEventWeight;  
  };
  
  // Main class
  /// Constructor with histogramming
  explicit HplusSelection(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& postfix = "");
  /// Constructor without histogramming
  explicit HplusSelection(const ParameterSet& config);
  virtual ~HplusSelection();

  virtual void bookHistograms(TDirectory* dir);
  
  /// Use silentAnalyze if you do not want to fill histograms or increment counters
  Data silentAnalyze(const Event& event, const JetSelection::Data& jetData, const BJetSelection::Data& bjetData, const TopSelectionMVA::Data& topData);
  /// analyze does fill histograms and incrementes counters
  Data analyze(const Event& event, const JetSelection::Data& jetData, const BJetSelection::Data& bjetData, const TopSelectionMVA::Data& topData);

private:
  /// Initialisation called from constructor
  void initialize(const ParameterSet& config);
  /// The actual selection
  Data privateAnalyze(const Event& event, const std::vector<Jet> selectedJets, const std::vector<Jet> selectedBjets, const TopSelectionMVA::Data& topData);
  /// Returns true if the two jets are the same
  bool areSameJets(const Jet& jet1, const Jet& jet2);
  /// Return true if a selected jet matches a selected bjet
  bool isBJet(const Jet& jet1, const std::vector<Jet>& bjets);

  bool isMatchedJet(const Jet& jet, std::vector<Jet> Jet1, std::vector<Jet> Jet2, std::vector<Jet> BJet, const unsigned int index);
  bool TopIsCrossCleaned( int Index, std::vector<Jet> jet1, std::vector<Jet> jet2, std::vector<Jet> bjet, std::vector<Jet> bjets);
  
  // Input parameters
  const DirectionalCut<float> cfg_TopMVACut;
  const DirectionalCut<float> cfg_AnyTopMVACut;
  const DirectionalCut<int>   cfg_FreeBjetsCut;

  // Scalefactor calculator
  TopTagSFCalculator fTopTagSFCalculator;

  // Histograms (1D)
  WrappedTH1 *hTopMVA_AllCandidates;
  WrappedTH1 *hTopMass_AllCandidates;
  WrappedTH1 *hTopPt_AllCandidates;
  WrappedTH1 *hTopMultiplicity_AllCandidates;
  WrappedTH1 *hTopMVA_SelectedCandidates;
  WrappedTH1 *hTopMass_SelectedCandidates;
  WrappedTH1 *hTopPt_SelectedCandidates;
  WrappedTH1 *hTopMultiplicity_SelectedCandidates;
  WrappedTH1 *hTopMVA_SelectedCleanedCandidates;
  WrappedTH1 *hTopMass_SelectedCleanedCandidates;
  WrappedTH1 *hTopPt_SelectedCleanedCandidates;
  WrappedTH1 *hTopMultiplicity_SelectedCleanedCandidates;
  WrappedTH1 *hTopMVA_NotSelectedCandidates;
  WrappedTH1 *hTopMass_NotSelectedCandidates;
  WrappedTH1 *hTopPt_NotSelectedCandidates;
  WrappedTH1 *hTopMultiplicity_NotSelectedCandidates;
  WrappedTH1 *hTopMVA_AllCleanedCandidates;
  WrappedTH1 *hTopMass_AllCleanedCandidates;
  WrappedTH1 *hTopPt_AllCleanedCandidates;
  WrappedTH1 *hTopMultiplicity_AllCleanedCandidates;
  WrappedTH1 *hTetrajetBJetPt;
  WrappedTH1 *hTetrajetBJetEta;
  WrappedTH1 *hTetrajetBJetBDisc;
  WrappedTH1 *hTetrajetPt;
  WrappedTH1 *hTetrajetMass;
  WrappedTH1 *hTetrajetEta;
  WrappedTH1 *hLdgTrijetPt;
  WrappedTH1 *hLdgTrijetMass;
  WrappedTH1 *hLdgTrijetJet1Pt;
  WrappedTH1 *hLdgTrijetJet1Eta;
  WrappedTH1 *hLdgTrijetJet1BDisc;
  WrappedTH1 *hLdgTrijetJet2Pt;
  WrappedTH1 *hLdgTrijetJet2Eta;
  WrappedTH1 *hLdgTrijetJet2BDisc;
  WrappedTH1 *hLdgTrijetBJetPt;
  WrappedTH1 *hLdgTrijetBJetEta;
  WrappedTH1 *hLdgTrijetBJetBDisc;
  WrappedTH1 *hLdgTrijetDiJetPt;
  WrappedTH1 *hLdgTrijetDiJetEta;
  WrappedTH1 *hLdgTrijetDiJetMass;
  WrappedTH1 *hLdgTrijetDijetDeltaR;
  WrappedTH1 *hLdgTrijetTopMassWMassRatio;
  WrappedTH1 *hLdgTrijet_DeltaR_Trijet_TetrajetBjet;
  WrappedTH1 *hLdgTrijet_DeltaEta_Trijet_TetrajetBjet;
  WrappedTH1 *hLdgTrijet_DeltaPhi_Trijet_TetrajetBjet;
  WrappedTH1 *hLdgTrijet_DeltaY_Trijet_TetrajetBjet;
  WrappedTH1 *hSubldgTrijetPt;
  WrappedTH1 *hSubldgTrijetMass;
  WrappedTH1 *hSubldgTrijetJet1Pt;
  WrappedTH1 *hSubldgTrijetJet1Eta;
  WrappedTH1 *hSubldgTrijetJet1BDisc;
  WrappedTH1 *hSubldgTrijetJet2Pt;
  WrappedTH1 *hSubldgTrijetJet2Eta;
  WrappedTH1 *hSubldgTrijetJet2BDisc;
  WrappedTH1 *hSubldgTrijetBJetPt;
  WrappedTH1 *hSubldgTrijetBJetEta;
  WrappedTH1 *hSubldgTrijetBJetBDisc;
  WrappedTH1 *hSubldgTrijetDiJetPt;
  WrappedTH1 *hSubldgTrijetDiJetEta;
  WrappedTH1 *hSubldgTrijetDiJetMass;
  WrappedTH1 *hSubldgTrijetDijetDeltaR;
  WrappedTH1 *hSubldgTrijetTopMassWMassRatio;
  WrappedTH1 *hSubldgTrijet_DeltaR_Trijet_TetrajetBjet;
  WrappedTH1 *hSubldgTrijet_DeltaEta_Trijet_TetrajetBjet;
  WrappedTH1 *hSubldgTrijet_DeltaPhi_Trijet_TetrajetBjet;
  WrappedTH1 *hSubldgTrijet_DeltaY_Trijet_TetrajetBjet;
  WrappedTH2 *hDeltaR_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet;
  WrappedTH2 *hDeltaEta_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet;
  WrappedTH2 *hDeltaPhi_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet;
  WrappedTH2 *hDeltaY_LdgTrijet_TetrajetBjet_Vs_SubldgTrijet_TetrajetBjet;
};

#endif
