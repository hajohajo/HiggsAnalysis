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
    fMt(-99.9),
    fDnnOutputValue(-99.9),
    bPassedDnnSelection(false),
    bPassedLooseCut(false),
    bPassedMediumCut(false),
    bPassedTightCut(false)
{}

DnnSelection::Data::~Data() { }

DnnSelection::DnnSelection(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& postfix)
: BaseSelection(eventCounter, histoWrapper, commonPlots, postfix),
    fLooseCut(config.getParameter<float>("looseCut")),
    fMediumCut(config.getParameter<float>("mediumCut")),
    fTightCut(config.getParameter<float>("tightCut")),
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

//    TF_CloseSession(sessionEven, status);
//    TF_CloseSession(sessionOdd, status);
    TF_DeleteGraph(graphEven);
    TF_DeleteSession(sessionEven, status);
    TF_DeleteGraph(graphOdd);
    TF_DeleteSession(sessionOdd, status);
    TF_DeleteSessionOptions(sessionOptions);
    TF_DeleteStatus(status);

}

void DnnSelection::initialize(const ParameterSet& config)
{
    std::srand(std::time(0)); 
    graphEven = TF_NewGraph();
    graphOdd = TF_NewGraph();
    status = TF_NewStatus();
    sessionOptions = TF_NewSessionOptions();
    runOptions = nullptr;

    const char* savedModelDirEven = "/home/joona/Documents/CodeRepos/HiggsAnalysis/NtupleAnalysis/data/DnnClassifierModels/evenModel";
    const char* savedModelDirOdd = "/home/joona/Documents/CodeRepos/HiggsAnalysis/NtupleAnalysis/data/DnnClassifierModels/oddModel";
    const char* tags = "serve";
    ntags = 1;

    sessionEven = TF_LoadSessionFromSavedModel(sessionOptions, runOptions, savedModelDirEven, &tags, ntags, graphEven, nullptr, status);
    sessionOdd = TF_LoadSessionFromSavedModel(sessionOptions, runOptions, savedModelDirOdd, &tags, ntags, graphOdd, nullptr, status);


    inputOperationEven = TF_GraphOperationByName(graphEven, "serving_default_input_layer");
    outputOperationEven = TF_GraphOperationByName(graphEven, "StatefulPartitionedCall");
    inputOperationOdd = TF_GraphOperationByName(graphOdd, "serving_default_input_layer");
    outputOperationOdd = TF_GraphOperationByName(graphOdd, "StatefulPartitionedCall");

    runInputsEven = (TF_Output*)malloc(1 * sizeof(TF_Output));
    runOutputsEven = (TF_Output*)malloc(1 * sizeof(TF_Output));

    runInputsEven[0].oper = inputOperationEven;
    runInputsEven[0].index = 0;

    runOutputsEven[0].oper = outputOperationEven;
    runOutputsEven[0].index = 0;

    runInputsOdd = (TF_Output*)malloc(1 * sizeof(TF_Output));
    runOutputsOdd = (TF_Output*)malloc(1 * sizeof(TF_Output));

    runInputsOdd[0].oper = inputOperationOdd;
    runInputsOdd[0].index = 0;

    runOutputsOdd[0].oper = outputOperationOdd;
    runOutputsOdd[0].index = 0;



    inputValues = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*1);
    outputValues = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*1);

    lookup[190.0] = 180.0;
    lookup[210.0] = 200.0;
    lookup[235.0] = 220.0;
    lookup[275.0] = 250.0;
    lookup[350.0] = 300.0;
    lookup[450.0] = 400.0;
    lookup[625.0] = 500.0;
    lookup[775.0] = 750.0;
    lookup[900.0] = 800.0;
    lookup[1250.0] = 1000.0;
    lookup[1750.0] = 1500.0;
    lookup[2500.0] = 2000.0;
    lookup[3500.0] = 3000.0;
}

void DnnSelection::bookHistograms(TDirectory *dir)
{
    const int nFracBin = 20;
    const float fFracMin = 0.0;
    const float fFracMax = 1.0;

    TDirectory* subdir = fHistoWrapper.mkdir(HistoLevel::kVital, dir, "DnnEventSelection_" + sPostfix);

    hDnnOutput = fHistoWrapper.makeTH<TH1F>(HistoLevel::kVital, subdir, "Dnn_Output", ";Dnn output value", nFracBin, fFracMin, fFracMax);
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
        fCommonPlots->fillControlPlotsAtDnnSelection(event, myData);
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
    
    int ind = std::rand() % true_masses.size(); 
    TrueMass = lookup.lower_bound(boost::algorithm::clamp(TransverseMass, 0.0, 3200.0))->second;
//    std::cout<<"MT: "<<TransverseMass<<std::endl;
//    std::cout<<"True: "<<TrueMass<<std::endl;
//    TrueMass =  true_masses[ind];

    //Ensuring the inputs are given in the exact order the network expects
    std::vector<float> inputVector = {ldgTrkPtFrac, deltaPhiBjetMet, deltaPhiTauMet, deltaPhiTauBjet, MET, tauPt, bjetPt,  TransverseMass, TrueMass};
  dims = {1, 9};

//    std::vector<float> inputVector = {ldgTrkPtFrac, deltaPhiBjetMet, deltaPhiTauMet, deltaPhiTauBjet, MET, tauPt, bjetPt,  TransverseMass};
//    dims = {1, 8};
    dataSize = std::accumulate(dims.begin(), dims.end(), sizeof(float), std::multiplies<std::int64_t>{});
    float* data = static_cast<float*>(std::malloc(dataSize));
    std::copy(inputVector.begin(), inputVector.end(), data);

    inputValues = (TF_Tensor**)malloc(1 * sizeof(TF_Tensor*));
    outputValues = (TF_Tensor**)malloc(1 * sizeof(TF_Tensor*));

    inputTensor = TF_NewTensor(TF_FLOAT, dims.data(), static_cast<int>(dims.size()), data, dataSize, &Deallocator, nullptr);
    inputValues[0] = inputTensor;

    //NOTE: "even" means the network should be used on even events when evaluating...!
    if(event.eventID().event()%2==0) {
        TF_SessionRun(sessionEven, nullptr,
                    runInputsEven, inputValues, 1,
                    runOutputsEven, outputValues, 1,
                    nullptr, 0, nullptr,
                    status);
    }else{
        TF_SessionRun(sessionOdd, nullptr,
                    runInputsOdd, inputValues, 1,
                    runOutputsOdd, outputValues, 1,
                    nullptr, 0, nullptr,
                    status);
    }

    float* outputData = (float*)TF_TensorData(outputValues[0]);
    output.fDnnOutputValue = outputData[0];
    TF_DeleteTensor(*inputValues);
    TF_DeleteTensor(*outputValues);

    hDnnOutput->Fill(output.getDnnOutput());
    //std::cout<<"Output: "<<output.getDnnOutput()<<std::endl;

    if(output.getDnnOutput() >= fLooseCut) {
        output.bPassedLooseCut = true;
    }

    if(output.getDnnOutput() >= fMediumCut) {
        output.bPassedMediumCut = true;
    }

    if(output.getDnnOutput() >= fTightCut) {
        output.bPassedTightCut = true;
    }

    if(output.bPassedMediumCut) {
        output.bPassedDnnSelection = true;
    }
//    else{
//	for (auto i = inputVector.begin(); i != inputVector.end(); ++i){
//    		std::cout << *i << ' ';
//	}
//    }

    return output;

}
