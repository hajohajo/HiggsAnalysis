#include "KerasModel.h"
#include <iostream>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "test_main.h"

// EXAMPLE:
// root -l -b -q test_main.cc\(\"../Keras_4Layers_32relu_32relu_32relu_1sigmoid_500Epochs_32BatchSize_500000Entrystop_19Inputs_07Feb2020_11h29m27s\",\"../histograms-TT_19var.root\"\)

// LAST USE:
// root -l -b -q test_main.cc\(\"../Keras_4Layers_32relu_32relu_32relu_1sigmoid_500Epochs_32BatchSize_500000Entrystop_19Inputs_07Feb2020_11h29m27s\",\"../histograms-TT_19var.root\"\)

void test_main(std::string folder, TString fileName) {
  std::string model = folder+"/model.txt";
  std::cout<<"===  test_main.cc"<<std::endl;
  std::cout<<"    Opening ROOT file "<<fileName<<std::endl;  
  TFile *file = TFile::Open(fileName, "R");
  //Read signal and background trees
  TTree *sigTree = (TTree*)file->Get("treeS");
  TTree *bkgTree = (TTree*)file->Get("treeB");

  int Nevts = 20; //sigTree->GetEntries();
  if (0) std::cout<<"Number of events "<<Nevts<<std::endl;

  //Keras Model
  KerasModel kerasModel;  
  std::cout<<"===  test_main.cc"<<std::endl;
  std::cout<<"    Loading model "<<model<<std::endl;
  //load keras model weights
  kerasModel.load_weights(model);

  //Signal
  SetBranchAddresses(sigTree);

  std::cout<<"===  test_main.cc"<<std::endl;
  std::cout<<"    Get the DNN score"<<std::endl;
  std::cout<<"    predicted signal:"<<std::endl;

  //https://stackoverflow.com/questions/28329457/cern-root-how-to-read-contents-from-a-ttree-root-file-into-array
  for (int i = 0; i < Nevts; i++) {
    sigTree->GetEntry(i);
    //std::cout<<"TrijetMass "<<TrijetMass<<std::endl;
    vector<float> inputs = GetListOfInputs(kerasModel.GetInputNames());
    float MVAoutput = kerasModel.predict(inputs).at(0);
    std::cout<<"    ["<<MVAoutput<<"]"<<std::endl;
  }
  
  //Background
  SetBranchAddresses(bkgTree);
  std::cout<<"    predicted background:"<<std::endl;
  for (int i = 0; i < Nevts; i++) {
    bkgTree->GetEntry(i);
    //std::cout<<"TrijetMass "<<TrijetMass<<std::endl;
    vector<float> inputs = GetListOfInputs(kerasModel.GetInputNames());
    float MVAoutput = kerasModel.predict(inputs).at(0);
    std::cout<<"    ["<<MVAoutput<<"]"<<std::endl;
  }  
}
