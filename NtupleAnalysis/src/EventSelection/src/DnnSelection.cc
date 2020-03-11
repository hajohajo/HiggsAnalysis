#include "EventSelection/interface/DnnSelection.h"
#include "EventSelection/interface/CommonPlots.h"

#include "Framework/interface/HistoWrapper.h"
#include "Framework/interface/Exception.h"

#include "DataFormat/interface/Event.h"

void Deallocator(void* data, size_t length, void* arg) 
{
    std::free(data);
}

DnnSelection::Data::Data()
:
    bPassedDnnSelection(false),
    fDnnOutputValue(-99.9)
{}

DnnSelection::Data::~Data() { }

DnnSelection::DnnSelection(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& postfix)
: BaseSelection(eventCounter, histoWrapper, commonPlots, postfix),
    cPassedDnnSelection(fEventCounter.addCounter("passed DNN event selection (" +postfix+")")),
    cSubAll(fEventCounter.addSubCounter("DNN event selection ("+postfix+")", "All")),
    cSubPassedDnnSelection(fEventCounter.addSubCounter("DNN event selection ("+postfix+")", "#geq DNN cut value"))

{
    initialize(config);
}

DnnSelection::DnnSelection(const ParameterSet& config)
: BaseSelection(),
    cPassedDnnSelection(fEventCounter.addCounter("passed DNN event selection")),
    cSubAll(fEventCounter.addSubCounter("DNN event selection", "All")),
    cSubPassedDnnSelection(fEventCounter.addSubCounter("DNN event selection", "#geq DNN cut value"))
{
    initialize(config);
    bookHistograms(new TDirectory());
}

DnnSelection::~DnnSelection()
{
    delete hDnnOutput;
    delete hDnnTransverseMass;
    delete hDnnTransverseMassLoose;
    delete hDnnTransverseMassMedium;
    delete hDnnTransverseMassTight;

    TF_DeleteGraph(graph);
    TF_DeleteSession(session, status);
    TF_DeleteSessionOptions(sessionOptions);
    TF_DeleteStatus(status);


}

void DnnSelection::initialize(const ParameterSet& config)
{
    graph = TF_NewGraph();
    status = TF_NewStatus();
    sessionOptions = TF_NewSessionOptions();
    runOptions = nullptr;

    const char* savedModelDir = "/work/hajohajo/MVAtoChargedHiggs/HiggsAnalysis/NtupleAnalysis/src/EventSelection/src/savedModel";
    const char* tags = "serve";
    ntags = 1;

    session = TF_LoadSessionFromSavedModel(sessionOptions, runOptions, savedModelDir, &tags, ntags, graph, nullptr, status);

    inputOperation = TF_GraphOperationByName(graph, "serving_default_inputClassifier");
    outputOperation = TF_GraphOperationByName(graph, "StatefulPartitionedCall");

    runInputs = (TF_Output*)malloc(1 * sizeof(TF_Output));
    runOutputs = (TF_Output*)malloc(1 * sizeof(TF_Output));

    runInputs[0].oper = inputOperation;
    runInputs[0].index = 0;

    runOutputs[0].oper = outputOperation;
    runOutputs[0].index = 0;

    inputValues = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*1);
    outputValues = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*1);

}

void DnnSelection::bookHistograms(TDirectory *dir)
{
    const int nMtBin = 50;
    const float fMtMin = 0.0;
    const float fMtMax = 500.0;

    const int nFracBin = 20;
    const float fFracMin = 0.0;
    const float fFracMax = 1.0;

    TDirectory* subdir = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "DnnEventSelection_" + sPostfix);

    hDnnOutput = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Dnn_Output", ";Dnn output value", nFracBin, fFracMin, fFracMax);
    hDnnTransverseMass = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Dnn_TransverseMass", ";m_{T}", nMtBin, fMtMin, fMtMax);
    hDnnTransverseMassLoose = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Dnn_TransverseMassLoose", ";m_{T}", nMtBin, fMtMin, fMtMax);
    hDnnTransverseMassMedium = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Dnn_TransverseMassMedium", ";m_{T}", nMtBin, fMtMin, fMtMax);
    hDnnTransverseMassTight = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Dnn_TransverseMassTight", ";m_{T}", nMtBin, fMtMin, fMtMax);

}

DnnSelection::Data DnnSelection::silentAnalyze(const Event& event, const Tau& tau, const METSelection::Data& metData, const BJetSelection::Data& bjetData)
{
    if(bjetData.getSelectedBJets().size() == 0) {
        DnnSelection::Data myData = DnnSelection::Data();
        return myData;
    }

    ensureAnalyzeAllowed(event.eventID());
    disableHistogramsAndCounters();
    DnnSelection::Data myData = privateAnalyze(event, tau, metData.getMET(), (const Jet&)bjetData.getSelectedBJets().at(0));
    enableHistogramsAndCounters();

    return myData;
}

DnnSelection::Data DnnSelection::analyze(const Event& event, const Tau& tau, const METSelection::Data& metData, const BJetSelection::Data& bjetData) {
    ensureAnalyzeAllowed(event.eventID());
    DnnSelection::Data myData = privateAnalyze(event, tau, metData.getMET(), (const Jet&)bjetData.getSelectedBJets().at(0));

    if (fCommonPlotsIsEnabled()){
//        fCommonPlots->fillControlPlotsAtDnnSelection(event, myData);
    }

    return myData;
}


DnnSelection::Data DnnSelection::privateAnalyze(const Event& event, const Tau& selectedTau, const math::XYVectorD& METVector, const Jet& selectedBjet)
{
    Data output;
    cSubAll.increment();

    MET = (float)std::sqrt(METVector.mag2());
    tauPt = (float)selectedTau.pt();
    ldgTrkPtFrac = (float)2.0*selectedTau.lChTrkPt()/selectedTau.pt() - 1.0;
    deltaPhiTauMet = (float)ROOT::Math::VectorUtil::DeltaPhi(selectedTau.p4(), METVector);
    deltaPhiTauBjet = (float)ROOT::Math::VectorUtil::DeltaPhi(selectedTau.p4(), selectedBjet.p4());
    bjetPt = (float)selectedBjet.pt();
    deltaPhiBjetMet = (float)ROOT::Math::VectorUtil::DeltaPhi(selectedBjet.p4(), METVector);
    TransverseMass = (float)TransverseMass::reconstruct(selectedTau, METVector);
    output.fMt = TransverseMass;

    //Ensuring the inputs are given in the exact order the network expects
    std::vector<float> inputVector = {MET, tauPt, ldgTrkPtFrac, deltaPhiTauMet, deltaPhiTauBjet, bjetPt, deltaPhiBjetMet, TransverseMass};

    dims = {1, 8};
    dataSize = std::accumulate(dims.begin(), dims.end(), sizeof(float), std::multiplies<std::int64_t>{});
    float* data = static_cast<float*>(std::malloc(dataSize));
    std::copy(inputVector.begin(), inputVector.end(), data);

    inputValues = (TF_Tensor**)malloc(1 * sizeof(TF_Tensor*));
    outputValues = (TF_Tensor**)malloc(1 * sizeof(TF_Tensor*));

    inputTensor = TF_NewTensor(TF_FLOAT, dims.data(), static_cast<int>(dims.size()), data, dataSize, &Deallocator, nullptr);
    inputValues[0] = inputTensor;

    TF_SessionRun(session, nullptr,
                runInputs, inputValues, 1,
                runOutputs, outputValues, 1,
                nullptr, 0, nullptr,
                status);

    float* outputData = (float*)TF_TensorData(outputValues[0]);
    float DnnOutput = outputData[0];

    hDnnOutput->Fill(DnnOutput);

    hDnnTransverseMass->Fill(TransverseMass);
    if(DnnOutput >= looseCut) {
        hDnnTransverseMassLoose->Fill(TransverseMass);
    }

    if(DnnOutput >= mediumCut) {
        hDnnTransverseMassMedium->Fill(TransverseMass);
    }

    if(DnnOutput >= tightCut) {
        hDnnTransverseMassTight->Fill(TransverseMass);
    }

    if(DnnOutput >= mediumCut) {
        output.bPassedDnnSelection = true;
    }

    output.fDnnOutputValue = DnnOutput;
    return output;

}
