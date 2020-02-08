#import <vector>
//Variable definition
float TrijetPtDR;
float TrijetDijetPtDR;
float TrijetBjetMass;
float TrijetLdgJetBDisc;
float TrijetSubldgJetBDisc;
float TrijetBJetLdgJetMass;
float TrijetBJetSubldgJetMass;
float TrijetMass;
float TrijetDijetMass;
float TrijetBJetBDisc;
float TrijetSoftDrop_n2;
float TrijetLdgJetCvsL;
float TrijetSubldgJetCvsL;
float TrijetLdgJetPtD;
float TrijetSubldgJetPtD;
float TrijetLdgJetAxis2; 
float TrijetSubldgJetAxis2;
float TrijetLdgJetMult;
float TrijetSubldgJetMult;

void SetBranchAddresses(TTree *tree){
  tree->SetBranchAddress("TrijetPtDR", &TrijetPtDR);
  tree->SetBranchAddress("TrijetDijetPtDR", &TrijetDijetPtDR);
  tree->SetBranchAddress("TrijetBjetMass", &TrijetBjetMass);
  tree->SetBranchAddress("TrijetLdgJetBDisc", &TrijetLdgJetBDisc);
  tree->SetBranchAddress("TrijetSubldgJetBDisc", &TrijetSubldgJetBDisc);
  tree->SetBranchAddress("TrijetBJetLdgJetMass", &TrijetBJetLdgJetMass);
  tree->SetBranchAddress("TrijetBJetSubldgJetMass", &TrijetBJetSubldgJetMass);
  tree->SetBranchAddress("TrijetMass", &TrijetMass);
  tree->SetBranchAddress("TrijetDijetMass", &TrijetDijetMass);
  tree->SetBranchAddress("TrijetBJetBDisc", &TrijetBJetBDisc);
  tree->SetBranchAddress("TrijetSoftDrop_n2", &TrijetSoftDrop_n2);
  tree->SetBranchAddress("TrijetLdgJetCvsL", &TrijetLdgJetCvsL);
  tree->SetBranchAddress("TrijetSubldgJetCvsL", &TrijetSubldgJetCvsL);
  tree->SetBranchAddress("TrijetLdgJetPtD", &TrijetLdgJetPtD);
  tree->SetBranchAddress("TrijetSubldgJetPtD", &TrijetSubldgJetPtD);
  tree->SetBranchAddress("TrijetLdgJetAxis2", &TrijetLdgJetAxis2);
  tree->SetBranchAddress("TrijetSubldgJetAxis2", &TrijetSubldgJetAxis2);
  tree->SetBranchAddress("TrijetLdgJetMult", &TrijetLdgJetMult);
  tree->SetBranchAddress("TrijetSubldgJetMult", &TrijetSubldgJetMult); 
}

vector<float> GetListOfInputs(vector<std::string> inputNames){

  vector<float> inputs;
  
  // For Loop: Protects form giving the inputs in wrong order (They should be given in the same order as they appear in the top tagger!)
  for (size_t i=0; i < inputNames.size(); i++){
    std::string var = inputNames.at(i);
    if (0) std::cout<<i<<": "<<var<<std::endl;
    if (var == "TrijetPtDR")                   inputs.push_back(TrijetPtDR);
    else if (var == "TrijetDijetPtDR")         inputs.push_back(TrijetDijetPtDR);
    else if (var == "TrijetBjetMass")          inputs.push_back(TrijetBjetMass);
    else if (var == "TrijetLdgJetBDisc")       inputs.push_back(TrijetLdgJetBDisc);
    else if (var == "TrijetSubldgJetBDisc")    inputs.push_back(TrijetSubldgJetBDisc);
    else if (var == "TrijetBJetLdgJetMass")    inputs.push_back(TrijetBJetLdgJetMass);
    else if (var == "TrijetBJetSubldgJetMass") inputs.push_back(TrijetBJetSubldgJetMass);
    else if (var == "TrijetMass")              inputs.push_back(TrijetMass);
    else if (var == "TrijetDijetMass")         inputs.push_back(TrijetDijetMass);
    else if (var == "TrijetBJetBDisc")         inputs.push_back(TrijetBJetBDisc);
    else if (var == "TrijetSoftDrop_n2")       inputs.push_back(TrijetSoftDrop_n2);
    else if (var == "TrijetLdgJetCvsL")        inputs.push_back(TrijetLdgJetCvsL);
    else if (var == "TrijetSubldgJetCvsL")     inputs.push_back(TrijetSubldgJetCvsL);
    else if (var == "TrijetLdgJetPtD")         inputs.push_back(TrijetLdgJetPtD);
    else if (var == "TrijetSubldgJetPtD")      inputs.push_back(TrijetSubldgJetPtD);
    else if (var == "TrijetLdgJetAxis2")       inputs.push_back(TrijetLdgJetAxis2);
    else if (var == "TrijetSubldgJetAxis2")    inputs.push_back(TrijetSubldgJetAxis2);
    else if (var == "TrijetLdgJetMult")        inputs.push_back(TrijetLdgJetMult);
    else if (var == "TrijetSubldgJetMult")     inputs.push_back(TrijetSubldgJetMult);
    else std::cout<<"Variable "<<var<<" is not one of the input variables!"<<std::endl;
  }
  if (inputNames.size() > inputs.size()) std::cout<<"Top tagger takes "<<inputNames.size()<<" inputs  ("<<inputs.size()<<" have been given!)"<<std::endl;
  return inputs;
}
