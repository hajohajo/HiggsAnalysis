#ifndef EventSelection_MVASelection_h
#define EventSelection_MVASelection_h


#include "EventSelection/interface/BaseSelection.h"
#include "EventSelection/interface/TauSelection.h"
#include "EventSelection/interface/BJetSelection.h"


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

class ParameterSet;
class CommonPlots;
class Event;
class EventCounter;
class HistoWrapper;
class WrappedTH1;
class WrappedTH2;

class MVASelection: public BaseSelection {
public:

  class Data {
  public:
    Data();
    ~Data();
    bool passedSelection() const { return bPassedSelection;}
    float mvaValue() const { return fMVAOutputValue;}
    void setTrue() {this->bPassedSelection=true;}
    void setValue(float value) {this->fMVAOutputValue=value;}

    friend class MVASelection;

  private:
    float fMVAOutputValue;
    bool bPassedSelection;

  };

  explicit MVASelection(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& prefix, const std::string& postfix);
  explicit MVASelection(const ParameterSet& config, const std::string& postfix);
  virtual ~MVASelection();

  virtual void bookHistograms(TDirectory* dir);

  Data analyze(const Event& event, TMVA::Reader& reader, const BJetSelection::Data& bJetData);
  Data silentAnalyze(const Event& event, TMVA::Reader& reader, const BJetSelection::Data& bJetData);

  TMVA::Reader *reader;
  Float_t MET, TransMass, met_x, met_y;
  Float_t Tau_pt, Tau_eta, Tau_phi, Tau_lChTrkPt;
  Float_t Bjet_pt, Bjet_eta, Bjet_phi;
  Float_t DPhi_tau_miss, DPhi_bjet_miss, Dist_tau_bjet;
  Float_t Transmass, Upsilon;


private:
  void initialize(const ParameterSet& config, const std::string& postfix);

  Data privateAnalyze(const Event& iEvent, const BJetSelection::Data& bJetData);

  Count cPassedMVASelection;
  Count cSubAll;

  WrappedTH1 *hMVAValueAll;

};

#endif
