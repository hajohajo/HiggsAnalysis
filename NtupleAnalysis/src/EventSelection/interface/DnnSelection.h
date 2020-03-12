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

        const bool passedSelection() const {return bPassedDnnSelection;}
        const float getDnnOutput() const {return fDnnOutputValue;}
        const float getMt() const {return fMt;}
        const bool passedLooseSelection() const {return bPassedLooseCut;}
        const bool passedMediumSelection() const {return bPassedMediumCut;}
        const bool passedTightSelection() const {return bPassedTightCut;}

        friend class DnnSelection;

    private:
        float fMt;
        float fDnnOutputValue;
        bool bPassedDnnSelection;
        bool bPassedLooseCut;
        bool bPassedMediumCut;
        bool bPassedTightCut;
    };

    explicit DnnSelection(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& postfix = "");
    explicit DnnSelection(const ParameterSet& config);

    virtual ~DnnSelection();
    virtual void bookHistograms(TDirectory* dir);

    Data silentAnalyze(const Event& event, const Tau& tau, const METSelection::Data& metData, const BJetSelection::Data& bjetData);
    Data analyze(const Event& event, const Tau& tau, const METSelection::Data& metData, const BJetSelection::Data& bjetData);

    TF_Graph* graphEven;
    TF_Graph* graphOdd;
    TF_Status* status;
    TF_SessionOptions* sessionOptions;
    TF_Session* sessionEven;
    TF_Session* sessionOdd;
    TF_Buffer* runOptions;
    TF_Tensor* inputTensor;
    TF_Tensor* outputTensor;
    TF_Tensor** inputValues;
    TF_Tensor** outputValues;
    TF_Operation* inputOperationEven;
    TF_Operation* outputOperationEven;
    TF_Output* runInputsEven;
    TF_Output* runOutputsEven;
    TF_Operation* inputOperationOdd;
    TF_Operation* outputOperationOdd;
    TF_Output* runInputsOdd;
    TF_Output* runOutputsOdd;


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

private:
    void initialize(const ParameterSet& config);
    Data privateAnalyze(const Event& event, const Tau& selectedTau, const math::XYVectorD& METVector, const Jet& bjet);

    float fLooseCut;
    float fMediumCut;
    float fTightCut;

    Count cPassedDnnSelection;
    Count cSubAll;
    Count cSubPassedDnnSelection;

    WrappedTH1 *hDnnOutput;
    WrappedTH1 *hDnnTransverseMass;
    WrappedTH1 *hDnnTransverseMassLoose;
    WrappedTH1 *hDnnTransverseMassMedium;
    WrappedTH1 *hDnnTransverseMassTight;

};

#endif
