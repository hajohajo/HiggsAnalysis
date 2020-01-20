// -*- c++ -*-
#ifndef EventSelection_MuonSelection_h
#define EventSelection_MuonSelection_h

#include "EventSelection/interface/BaseSelection.h"
#include "DataFormat/interface/Muon.h"
#include "Framework/interface/EventCounter.h"
#include "Framework/interface/GenericScaleFactor.h"

#include <string>
#include <vector>

class ParameterSet;
class CommonPlots;
class CommonPlots_ttm;
class Event;
class EventCounter;
class HistoWrapper;
class WrappedTH1;
class WrappedTH2;

class MuonSelection: public BaseSelection {
public:
    /**
    * Class to encapsulate the access to the data members of
    * TauSelection. If you want to add a new accessor, add it here
    * and keep all the data of TauSelection private.
    */
  class Data {
  public:
    // The reason for pointer instead of reference is that const
    // reference allows temporaries, while const pointer does not.
    // Here the object pointed-to must live longer than this object.
    Data();
    ~Data();

    const bool hasIdentifiedMuons() const { return (fSelectedMuons.size() > 0); }
    const std::vector<Muon>& getSelectedMuons() const { return fSelectedMuons; }
    const float getHighestSelectedMuonPt() const { return fHighestSelectedMuonPt; }
    const float getHighestSelectedMuonEta() const { return fHighestSelectedMuonEta; }
    const float getHighestSelectedMuonPhi() const { return fHighestSelectedMuonPhi; }
    const float getHighestSelectedMuonPtBeforePtCut() const { return fHighestSelectedMuonPtBeforePtCut; }
    const float getMuonTriggerSF() const { return fMuonTriggerSF; }
    const float getMuonIDSF() const { return fMuonIDSF; }
    const int getHLTMuonCharge() const { return fHLTMuonCharge; }
    // FIXME: Add MC information if deemed necessary
//     const bool eventContainsMuonFromCJet() const { return fHasMuonFromCjetStatus; }
//     const bool eventContainsMuonFromBJet() const { return fHasMuonFromBjetStatus; }
//     const bool eventContainsMuonFromCorBJet() const { return eventContainsMuonFromCJet() || eventContainsMuonFromBJet(); }

    // Getters for anti-isolated muons (i.e. passed other cuts but not isolation)
    const bool isAntiIsolated() const { return !hasIdentifiedMuons(); }
    const bool hasAntiIsolatedMuons() const { return (fAntiIsolatedMuons.size() > 0); }
    const std::vector<Muon>& getAntiIsolatedMuons() const { return fAntiIsolatedMuons; }

    friend class MuonSelection;

  private:
    /// pt and eta of highest pt muon passing the selection
    float fHighestSelectedMuonPt;
    float fHighestSelectedMuonEta;
    float fHighestSelectedMuonPhi;
    float fHighestSelectedMuonPtBeforePtCut;
    /// Cache muon identification scale factor
    float fMuonIDSF;
    /// Cache for muon trigger SF
    float fMuonTriggerSF;
    int fHLTMuonCharge;
    /// MC info about non-isolated muons
    //bool fHasMuonFromCjetStatus;
    //bool fHasMuonFromBjetStatus;
    /// Muon collection after all selections
    std::vector<Muon> fSelectedMuons;
    std::vector<Muon> fAntiIsolatedMuons;
  };
  
  // Main class
  /// Constructor with histogramming
  explicit MuonSelection(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& postfix);
  explicit MuonSelection(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots_ttm* commonPlots, const std::string& postfix);
  explicit MuonSelection(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, std::nullptr_t, const std::string& postfix);
  /// Constructor without histogramming
  explicit MuonSelection(const ParameterSet& config, const std::string& postfix);
  virtual ~MuonSelection();

  virtual void bookHistograms(TDirectory* dir);
  
  /// Use silentAnalyze if you do not want to fill histograms or increment counters
  Data silentAnalyze(const Event& event);
  /// analyze does fill histograms and incrementes counters
  Data analyze(const Event& event);
  Data analyzeLoose(const Event& event);
private:
  /// Initialisation called from constructor
  void initialize(const ParameterSet& config, const std::string& postfix);
  /// The actual selection
  Data privateAnalyze(const Event& iEvent);
  Data privateAnalyzeLoose(const Event& iEvent);
  bool passTrgMatching(const Muon& muon, std::vector<math::LorentzVectorT<double>>& trgMuons) const;

  // Input parameters
  bool cfg_ApplyTriggerMatching;
  float cfg_TriggerMatchingCone;
  const double cfg_MuonPtCut;
  const double cfg_MuonEtaCut;
  float fRelIsoCut;
  float fMiniIsoCut;
  bool fVetoMode;
  bool fMiniIsol;
  
  // muon identification SF
  GenericScaleFactor fMuonIDSFReader;
  // muon trigger SF
  GenericScaleFactor fMuonTriggerSFReader;

  // Event counter for passing selection
  Count cPassedMuonSelection;
  // Sub counters
  Count cSubAll;
  Count cSubPassedIsPresent;
  Count cSubPassedTriggerMatching;
  Count cSubPassedPt;
  Count cSubPassedEta;
  Count cSubPassedID;
  Count cSubPassedIsolation;
  Count cSubPassedSelection;
  Count cSubPassedVeto;
  
  // Histograms
  WrappedTH1 *hTriggerMatchDeltaR;
  WrappedTH1 *hMuonNAll;
  WrappedTH1 *hMuonPtAll;
  WrappedTH1 *hMuonEtaAll;
  WrappedTH1 *hMuonRelIsoAll;
  WrappedTH1 *hMuonMiniIsoAll;
  
  WrappedTH1 *hMuonNPassed;
  WrappedTH1 *hMuonPtPassed;
  WrappedTH1 *hMuonEtaPassed;
  WrappedTH1 *hMuonRelIsoPassed;
  WrappedTH1 *hMuonMiniIsoPassed;
  
  WrappedTH1 *hPtResolution;
  WrappedTH1 *hEtaResolution;
  WrappedTH1 *hPhiResolution;
  
  WrappedTH1 *hIsolPtBefore;
  WrappedTH1 *hIsolEtaBefore;
  WrappedTH1 *hIsolVtxBefore;
  WrappedTH1 *hIsolRelIsoBefore;
  WrappedTH1 *hIsolMiniIsoBefore;
  
  WrappedTH1 *hIsolPtAfter;
  WrappedTH1 *hIsolEtaAfter;
  WrappedTH1 *hIsolVtxAfter;
  WrappedTH1 *hIsolRelIsoAfter;
  WrappedTH1 *hIsolMiniIsoAfter;
};

#endif
