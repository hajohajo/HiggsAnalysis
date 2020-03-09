#ifndef EventSelection_DnnSelection_h
#define EventSelection_DnnSelection_h

#include "EventSelection/interface/BaseSelection.h"
#include "EventSelection/interface/TauSelection.h"
#include "EventSelection/interface/METSelection.h"
#include "EventSelection/interface/BJetSelection.h"
#include "EventSelection/interface/TransverseMass.h"

#include <Math/Vector2D.h>
#include <Math/VectorUtil.h>

#include <tensorflow/c/c_api.h>

class ParameterSet;
class CommonPlots;
class Event;
class EventCounter;
class HistoWrapper;
class WrappedTH1;

class DnnSelection: public BaseSelection
{
public:
    class Data {
    public:
        Data();
        ~Data();

        bool passedSelection() const {return bPassedDnnSelection;}
        float getClassifierOutput() const {return fDnnOutputValue;}

        friend class DnnSelection;

    private:
        bool bPassedDnnSelection;
        float fDnnOutputValue;
        float fMt;
    };

    explicit DnnSelection(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& postfix = "");
    explicit DnnSelection(const ParameterSet& config);

    virtual ~DnnSelection();
    virtual void bookHistograms(TDirectory* dir);

    Data silentAnalyze(const Event& event, const Tau& tau, const METSelection::Data& metData, const BJetSelection::Data& bjetData);
    Data analyze(const Event& event, const Tau& tau, const METSelection::Data& metData, const BJetSelection::Data& bjetData);

    TF_Graph* graph;
    TF_Status* status;
    TF_SessionOptions* sessionOptions;
    TF_Session* session;
    TF_Buffer* runOptions;
    TF_Tensor* inputTensor;
    TF_Tensor* outputTensor;
    TF_Tensor** inputValues;
    TF_Tensor** outputValues;
    TF_Operation* inputOperation;
    TF_Operation* outputOperation;
    TF_Output* runInputs;
    TF_Output* runOutputs;

    int ntags;
    std::vector<std::int64_t> dims;
    std::int64_t dataSize;


    float MET;
    float tauPt;
    float ldgTrkPtFrac;
    float deltaPhiTauMet;
    float deltaPhiTauBjet;
    float bjetPt;
    float deltaPhiBjetMet;
    float TransverseMass;

    float looseCut = 0.3;
    float mediumCut = 0.5;
    float tightCut = 0.7;

private:
    void initialize(const ParameterSet& config);
    Data privateAnalyze(const Event& event, const Tau& selectedTau, const math::XYVectorD& METVector, const Jet& bjet);

    Count cSubAll;
    Count cSubPassedDnnSelection;
    Count cPassedDnnSelection;

    WrappedTH1 *hDnnOutput;
    WrappedTH1 *hDnnTransverseMass;
    WrappedTH1 *hDnnTransverseMassLoose;
    WrappedTH1 *hDnnTransverseMassMedium;
    WrappedTH1 *hDnnTransverseMassTight;
};

#endif
