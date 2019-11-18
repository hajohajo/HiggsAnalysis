// -*- c++ -*-
#ifndef EventSelection_HplusSelectionNN_h
#define EventSelection_HplusSelectionNN_h

#include "EventSelection/interface/BaseSelection.h"
#include "EventSelection/interface/JetSelection.h"
#include "EventSelection/interface/BJetSelection.h"
#include "EventSelection/interface/TopSelectionNN.h"
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

class HplusSelectionNN: public BaseSelection {
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
    const size_t getAllCleanedTopsSize() const { return nAllCleanedTops;}
    
    // Trijet-1
    const float getMVAmax1() const { return fMVAmax1; }
    // Trijet-2
    const float getMVAmax2() const { return fMVAmax2; }

    // Leading Trijet
    const math::XYZTLorentzVector getLdgTrijet() const
    {       
      if (fTrijet1_p4.pt() > fTrijet2_p4.pt()) return fTrijet1_p4;
      else return fTrijet2_p4;
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

    friend class HplusSelectionNN;

  private:
    /// Boolean for passing selection
    bool bPassedSelection;
    bool bHasTwoTopsAndFreeB;
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
  explicit HplusSelectionNN(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& postfix = "");
  /// Constructor without histogramming
  explicit HplusSelectionNN(const ParameterSet& config);
  virtual ~HplusSelectionNN();

  virtual void bookHistograms(TDirectory* dir);
  
  /// Use silentAnalyze if you do not want to fill histograms or increment counters
  Data silentAnalyze(const Event& event, const JetSelection::Data& jetData, const BJetSelection::Data& bjetData, const TopSelectionNN::Data& topData);
  /// analyze does fill histograms and incrementes counters
  Data analyze(const Event& event, const JetSelection::Data& jetData, const BJetSelection::Data& bjetData, const TopSelectionNN::Data& topData);

private:
  /// Initialisation called from constructor
  void initialize(const ParameterSet& config);
  /// The actual selection
  Data privateAnalyze(const Event& event, const std::vector<Jet> selectedJets, const std::vector<Jet> selectedBjets, const TopSelectionNN::Data& topData);
  /// Returns true if the two jets are the same
  bool areSameJets(const Jet& jet1, const Jet& jet2);
  /// Return true if a selected jet matches a selected bjet
  bool isBJet(const Jet& jet1, const std::vector<Jet>& bjets);

  bool isMatchedJet(const Jet& jet, std::vector<Jet> Jet1, std::vector<Jet> Jet2, std::vector<Jet> BJet, const unsigned int index);
  bool TopIsCrossCleaned( int Index, std::vector<Jet> jet1, std::vector<Jet> jet2, std::vector<Jet> bjet);
  
  // Input parameters
  const DirectionalCut<float> cfg_TopMVACut;
  const DirectionalCut<float> cfg_AnyTopMVACut;
  const DirectionalCut<int>   cfg_FreeBjetsCut;
  //insert here

  // Scalefactor calculator
  TopTagSFCalculator fTopTagSFCalculator;

  // Histograms (1D)
  //WrappedTH1 *hTopMVA_AllCandidates;
  // insert here  
};

#endif
