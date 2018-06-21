#include "EventSelection/interface/MVASelection.h"

#include "Framework/interface/ParameterSet.h"
#include "EventSelection/interface/CommonPlots.h"
#include "DataFormat/interface/Event.h"
#include "Framework/interface/HistoWrapper.h"


//TMVA Stuff

MVASelection::Data::Data() 
: bPassedSelection(false) { }

MVASelection::Data::~Data(){ }

MVASelection::MVASelection(const ParameterSet& config, EventCounter& eventCounter, HistoWrapper& histoWrapper, CommonPlots* commonPlots, const std::string& prefix, const std::string& postfix)
: BaseSelection(eventCounter, histoWrapper, commonPlots, prefix+postfix),
  cPassedMVASelection(fEventCounter.addCounter("passed MVA selection "+prefix+" ("+postfix+")")),
  cSubAll(fEventCounter.addSubCounter("MVA selection "+prefix+" ("+postfix+")", "All events"))
{
  initialize(config, postfix);
}


MVASelection::MVASelection(const ParameterSet& config, const std::string& postfix)
: BaseSelection(),
  cPassedMVASelection(fEventCounter.addCounter("passed MVA selection ("+postfix+")")),
  cSubAll(fEventCounter.addSubCounter("MVA selection ("+postfix+")", "All events"))
{
  initialize(config, postfix);
  bookHistograms(new TDirectory());
}

MVASelection::~MVASelection() {
  delete hMVAValueAll;
  delete reader;
}

void MVASelection::initialize(const ParameterSet& config, const std::string& postfix) {
  reader = new TMVA::Reader( "!Color:Silent" );
  reader->AddVariable("Tau_pt",&Tau_pt);
  reader->AddVariable("Bjet_pt",&Bjet_pt);
  reader->AddVariable("MET",&MET);

  reader->AddVariable("DPhi_tau_miss",&DPhi_tau_miss);
  reader->AddVariable("DPhi_bjet_miss",&DPhi_bjet_miss);
  reader->AddVariable("Dist_tau_bjet",&Dist_tau_bjet);

  reader->AddVariable("Upsilon",&Upsilon);
  reader->AddVariable("Transmass",&Transmass);

  reader->BookMVA("BDTG method","../../../data/MyClassification_BDTG.weights.xml");
//  reader->BookMVA("DNN method","MyClassification_DNN.weights.xml");

}

void MVASelection::bookHistograms(TDirectory* dir) {
  TDirectory* subdir = fHistoWrapper.mkdir(HistoLevel::kDebug, dir, "MVASelection_"+sPostfix);
  hMVAValueAll = fHistoWrapper.makeTH<TH1F>(HistoLevel::kDebug, subdir, "MVAValueAll", "MVAValue, all; BDT", 40, -1.0, 1.0);
}

MVASelection::Data MVASelection::silentAnalyze(const Event& event, TMVA::Reader& reader, const BJetSelection::Data& bJetData) {
  ensureSilentAnalyzeAllowed(event.eventID());
  // Disable histogram filling and counter
  disableHistogramsAndCounters();
  Data myData = privateAnalyze(event, bJetData);
  enableHistogramsAndCounters();
  return myData;
}

MVASelection::Data MVASelection::analyze(const Event& event, TMVA::Reader& reader, const BJetSelection::Data& bJetData) {
  ensureAnalyzeAllowed(event.eventID());
  MVASelection::Data data = privateAnalyze(event, bJetData);
  if (fCommonPlotsIsEnabled())
    fCommonPlots->fillControlPlotsAtMVASelection(event, data);
  return data;
}

MVASelection::Data MVASelection::privateAnalyze(const Event& event, const BJetSelection::Data& bJetData) {
  cSubAll.increment();
  bool passedMVA = false;
  Data output;


  Tau_pt=event.taus()[0].pt();
  Tau_eta=event.taus()[0].eta();
  Tau_phi=event.taus()[0].phi();
  Tau_lChTrkPt=event.taus()[0].lChTrkPt();

  auto& bJets = bJetData.getSelectedBJets();
  if (bJets.size()>0){
	  Bjet_pt=bJets[0].pt();
	  Bjet_eta=bJets[0].eta();
	  Bjet_phi=bJets[0].phi();
  } else {
          Bjet_pt=event.jets()[0].pt();
          Bjet_eta=event.jets()[0].eta();
          Bjet_phi=event.jets()[0].phi();
  }

  Float_t met_x,met_y;
  met_x=event.met().x();
  met_y=event.met().y();
  MET=sqrt(pow(met_x,2)+pow(met_y,2));
  double METphi=atan2(met_y,met_x);
  Dist_tau_bjet = sqrt(pow((std::min)(abs(Bjet_phi-Tau_phi),float(2.*TMath::Pi())-abs(Bjet_phi-Tau_phi))+abs(Bjet_eta-Tau_eta),2));

  DPhi_tau_miss = min(abs((Tau_phi)-METphi),2.*TMath::Pi()-abs((Tau_phi)-METphi));
  DPhi_bjet_miss = min(abs((Bjet_phi)-METphi),2.*TMath::Pi()-abs((Bjet_phi)-METphi));

  Upsilon=2.0*Tau_lChTrkPt/Tau_pt-1.0;
  Transmass=sqrt(2.*Tau_pt*MET*(1-cos(METphi+Tau_phi)));


  output.setValue(reader->EvaluateMVA("BDTG method"));
//  std::cout<<output.mvaValue()<<std::endl;
//  output.setValue(reader->EvaluateMVA("DNN method"));
  hMVAValueAll->Fill(output.mvaValue());
//  passedMVA=(output.mvaValue()>-0.3);
  passedMVA=(output.mvaValue()>0.4);
  if(passedMVA){
    output.setTrue();
    cPassedMVASelection.increment();
  }
  return output;
}

