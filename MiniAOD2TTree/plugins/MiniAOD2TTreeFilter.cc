#include "HiggsAnalysis/MiniAOD2TTree/interface/MiniAOD2TTreeFilter.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include <regex>

MiniAOD2TTreeFilter::MiniAOD2TTreeFilter(const edm::ParameterSet& iConfig) :
    //prescaleWeight(iConfig.getParameter<edm::ParameterSet>("PrescaleProvider"), consumesCollector(), this),
    outputFileName(iConfig.getParameter<std::string>("OutputFileName")),
    PUInfoInputFileName(iConfig.getParameter<std::string>("PUInfoInputFileName")),
//    TopPtInputFileName(iConfig.getParameter<std::string>("TopPtInputFileName")),
    codeVersion(iConfig.getParameter<std::string>("CodeVersion")),
    dataVersion(iConfig.getParameter<std::string>("DataVersion")),
    cmEnergy(iConfig.getParameter<int>("CMEnergy")),
    eventInfoCollections(iConfig.getParameter<edm::ParameterSet>("EventInfo"))
{
  if (0) std::cout << "MiniAOD2TTreeFilter::MiniAOD2TTreeFilter() " << std::endl;
  PUInfoPSInputFileName = "";  
  if (iConfig.exists("PUInfoPSInputFileName")) 
    {
      PUInfoPSInputFileName = iConfig.getParameter<std::string>("PUInfoPSInputFileName");
    }
  fOUT = TFile::Open(outputFileName.c_str(),"RECREATE");	
  Events = new TTree("Events","");
  
  eventInfo = new EventInfoDumper(consumesCollector(), eventInfoCollections);
  eventInfo->book(Events);
  
  skimDumper = 0;
  if (iConfig.exists("Skim")) 
    {
      if (0) std::cout << "=== MiniAOD2TTreeFilter (Skim)" << std::endl;
      skim = iConfig.getParameter<edm::ParameterSet>("Skim");
      skimDumper = new SkimDumper(consumesCollector(), skim);
      skimDumper->book();
    } 
  else
    {
      std::cout << "MiniAOD2TTreeFilter::MiniAOD2TTreeFilter() skimDumper ignored, because 'Skim' is missing from config" << std::endl;
    }
  
  trgDumper = 0;
  if (iConfig.exists("Trigger")) 
    {
      if (0) std::cout << "=== MiniAOD2TTreeFilter (Trigger)" << std::endl;
      trigger = iConfig.getParameter<edm::ParameterSet>("Trigger");
      trgDumper = new TriggerDumper(consumesCollector(), trigger);
      trgDumper->book(Events);
      hltProcessName = trigger.getParameter<edm::InputTag>("TriggerResults").process();
    } 
  else
    {
      std::cout << "MiniAOD2TTreeFilter::MiniAOD2TTreeFilter() trgDumper ignored, because 'Trigger' is missing from config" << std::endl;
    }
  
  metNoiseFilterDumper = 0;
  if (iConfig.exists("METNoiseFilter")) 
    {
      if (0) std::cout << "=== MiniAOD2TTreeFilter (METNoiseFilter)" << std::endl;
      metNoiseFilter = iConfig.getParameter<edm::ParameterSet>("METNoiseFilter");
      metNoiseFilterDumper = new METNoiseFilterDumper(consumesCollector(), metNoiseFilter);
      metNoiseFilterDumper->book(Events);
    } else 
    {
      std::cout << "MiniAOD2TTreeFilter::MiniAOD2TTreeFilter() metNoiseFilterDumper ignored, because 'METNoiseFilter' is missing from config" << std::endl;
    }
  
  TopPtInputFileName = "";
  if (iConfig.exists("TopPtInputFileName")) 
    {
      TopPtInputFileName = iConfig.getParameter<std::string>("TopPtInputFileName");
    }
  
  tauDumper = 0;
  if (iConfig.exists("Taus")) 
    {
      if (0) std::cout << "=== MiniAOD2TTreeFilter (Taus)" << std::endl;
      tauCollections = iConfig.getParameter<std::vector<edm::ParameterSet>>("Taus");
      tauDumper = new TauDumper(consumesCollector(), tauCollections);
      tauDumper->book(Events);
    } 
  else
    {
      std::cout << "MiniAOD2TTreeFilter::MiniAOD2TTreeFilter() tauDumper ignored, because 'Taus' is missing from config" << std::endl;
    }
  
  electronDumper = 0;
  if (iConfig.exists("Electrons")) 
    {
      if (0) std::cout << "=== MiniAOD2TTreeFilter (Electrons)" << std::endl;
      electronCollections = iConfig.getParameter<std::vector<edm::ParameterSet>>("Electrons");
      electronDumper = new ElectronDumper(consumesCollector(), electronCollections);
      electronDumper->book(Events);
    } 
  else
    {
      std::cout << "MiniAOD2TTreeFilter::MiniAOD2TTreeFilter() electronDumper ignored, because 'Electrons' is missing from config" << std::endl;
    }
  
  muonDumper = 0;
  if (iConfig.exists("Muons")) 
    {
      if (0) std::cout << "=== MiniAOD2TTreeFilter (Muons)" << std::endl;
      muonCollections = iConfig.getParameter<std::vector<edm::ParameterSet>>("Muons");
      muonDumper = new MuonDumper(consumesCollector(), muonCollections, eventInfoCollections.getParameter<edm::InputTag>("OfflinePrimaryVertexSrc"));
      muonDumper->book(Events);
    }
  else
    {
      std::cout << "MiniAOD2TTreeFilter::MiniAOD2TTreeFilter() muonDumper ignored, because 'Muons' is missing from config" << std::endl;
    }

  
  jetDumper = 0;
  if (iConfig.exists("Jets")) 
    {
      if (0) std::cout << "=== MiniAOD2TTreeFilter (Jets)" << std::endl;
      jetCollections = iConfig.getParameter<std::vector<edm::ParameterSet>>("Jets");
      jetDumper = new JetDumper(consumesCollector(), jetCollections);
      jetDumper->book(Events);
    } 
  else
    {
      std::cout << "MiniAOD2TTreeFilter::MiniAOD2TTreeFilter() jetDumper ignored, because 'Jets' is missing from config" << std::endl;
    }
  
  topDumper = 0;
  if (iConfig.exists("Top")) 
    {
      if (0) std::cout << "=== MiniAOD2TTreeFilter (Top)" << std::endl;
      topCollections = iConfig.getParameter<std::vector<edm::ParameterSet>>("Top");
      topDumper = new TopDumper(consumesCollector(), topCollections);
      topDumper->book(Events);
    } 
  else
    {
      std::cout << "MiniAOD2TTreeFilter::MiniAOD2TTreeFilter() topDumper ignored, because Ttop' is missing from config" << std::endl;
    }
  
  metDumper = 0;
  if (iConfig.exists("METs")) 
    {
      if (0) std::cout << "=== MiniAOD2TTreeFilter (METs)" << std::endl;
      metCollections = iConfig.getParameter<std::vector<edm::ParameterSet>>("METs");
      metDumper = new METDumper(consumesCollector(), metCollections, this->isMC());
      metDumper->book(Events);
    } 
  else
    {
      std::cout << "MiniAOD2TTreeFilter::MiniAOD2TTreeFilter() metDumper ignored, because 'METs' is missing from config" << std::endl;
    }
  
  trackDumper = 0;
  if (iConfig.exists("Tracks")) 
    {
      if (0) std::cout << "=== MiniAOD2TTreeFilter (Tracks)" << std::endl;
      trackCollections = iConfig.getParameter<std::vector<edm::ParameterSet>>("Tracks");
      trackDumper = new TrackDumper(consumesCollector(), trackCollections);
      trackDumper->book(Events);
    } 
  else 
    {
      std::cout << "MiniAOD2TTreeFilter::MiniAOD2TTreeFilter() trackDumper ignored, because 'Tracks' is missing from config" << std::endl;
    }
  
  genMetDumper = 0;
  genWeightDumper = 0;
  genParticleDumper = 0;
  genJetDumper = 0;
  
  if(this->isMC())
    {
      if (iConfig.exists("GenMETs")) 
	{
	  if (0) std::cout << "=== MiniAOD2TTreeFilter (GenMETs)" << std::endl;
	  genMetCollections = iConfig.getParameter<std::vector<edm::ParameterSet>>("GenMETs");
	  genMetDumper = new GenMETDumper(consumesCollector(), genMetCollections);
	  genMetDumper->book(Events);
	} 
      else 
	{
	  std::cout << "MiniAOD2TTreeFilter::MiniAOD2TTreeFilter() genMETDumper ignored, because 'GenMETs' is missing from config" << std::endl;
	}
      
      if (iConfig.exists("GenWeights"))
	{
	  if (0) std::cout << "=== MiniAOD2TTreeFilter (GenWeights)" << std::endl;
	  genWeightCollections = iConfig.getParameter<std::vector<edm::ParameterSet>>("GenWeights");
	  genWeightDumper = new GenWeightDumper(consumesCollector(), genWeightCollections);
	  genWeightDumper->book(Events);
	} 
      else 
	{
	  std::cout << "MiniAOD2TTreeFilter::MiniAOD2TTreeFilter() genWeightDumper ignored, because 'GenWeights' is missing from config" << std::endl;
	}
      
      if (iConfig.exists("GenParticles")) 
	{
	  if (0) std::cout << "=== MiniAOD2TTreeFilter (GenParticles)" << std::endl;
	  genParticleCollections = iConfig.getParameter<std::vector<edm::ParameterSet>>("GenParticles");
	  genParticleDumper = new GenParticleDumper(consumesCollector(), genParticleCollections);
	  genParticleDumper->book(Events);
	} 
      else
	{
	  std::cout << "MiniAOD2TTreeFilter::MiniAOD2TTreeFilter() genParticleDumper ignored, because 'GenParticles' is missing from config" << std::endl;
	}

      if (iConfig.exists("GenJets")) 
	{
	  if (0) std::cout << "=== MiniAOD2TTreeFilter (GenJets)" << std::endl;
	  genJetCollections = iConfig.getParameter<std::vector<edm::ParameterSet>>("GenJets");
	  genJetDumper = new GenJetDumper(consumesCollector(), genJetCollections);
	  genJetDumper->book(Events);
	} 
      else
	{
	  std::cout << "MiniAOD2TTreeFilter::MiniAOD2TTreeFilter() genJetDumper ignored, because 'GenJets' is missing from config" << std::endl;
	}
      }  
}

MiniAOD2TTreeFilter::~MiniAOD2TTreeFilter() {
  system("ls -lt");
}

void MiniAOD2TTreeFilter::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup) {
    bool changed = true;
    if (0) std::cout << "=== MiniAOD2TTreeFilter::beginRun()" << std::endl;
    hltConfig.init(iRun,iSetup,hltProcessName,changed);
    if(trgDumper != 0) trgDumper->book(iRun,hltConfig);

    return;
}

void MiniAOD2TTreeFilter::beginJob(){
  if (0) std::cout << "MiniAOD2TTreeFilter::beginJob() " << std::endl;
  return;
}   

bool MiniAOD2TTreeFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
    if (0) std::cout << "=== MiniAOD2TTreeFilter::filter()" << std::endl;
    reset();

    eventInfo->fill(iEvent,iSetup);

    bool accept = true;

    if (0) std::cout << "=== MiniAOD2TTreeFilter::filter() trgDumper" << std::endl;
    if (trgDumper) accept = accept && trgDumper->fill(iEvent, iSetup);

    if (0) std::cout << "=== MiniAOD2TTreeFilter::filter() metNoiseFilterDumper" << std::endl;
    if (metNoiseFilterDumper) accept = accept && metNoiseFilterDumper->fill(iEvent,iSetup);

    if (tauDumper) 
      {
	if (0) std::cout << "=== MiniAOD2TTreeFilter::filter() tauDumper" << std::endl;
	accept = accept && tauDumper->fill(iEvent,iSetup);
        if (trgDumper) trgDumper->triggerMatch(trigger::TriggerTau,tauDumper->selected());
      }

    if (electronDumper) 
      {
	if (0) std::cout << "=== MiniAOD2TTreeFilter::filter() electronDumper" << std::endl;
	accept = accept && electronDumper->fill(iEvent,iSetup);
	if (trgDumper) trgDumper->triggerMatch(trigger::TriggerCluster,electronDumper->selected());
      }

    if (muonDumper) 
      {
	if (0) std::cout << "=== MiniAOD2TTreeFilter::filter() muonDumper" << std::endl;
	accept = accept && muonDumper->fill(iEvent,iSetup);
	if (trgDumper) trgDumper->triggerMatch(trigger::TriggerMuon,muonDumper->selected());
      }

    if (0) std::cout << "=== MiniAOD2TTreeFilter::filter() jetDumper" << std::endl;
    if (jetDumper) accept = accept && jetDumper->fill(iEvent,iSetup);

    if (0) std::cout << "=== MiniAOD2TTreeFilter::filter() topDumper" << std::endl;
    if (topDumper) accept = accept && topDumper->fill(iEvent,iSetup);

    if (0) std::cout << "=== MiniAOD2TTreeFilter::filter() metDumper" << std::endl;
    if (metDumper) accept = accept && metDumper->fill(iEvent,iSetup);

    if (0) std::cout << "=== MiniAOD2TTreeFilter::filter() genMetDumper" << std::endl;
    if (genMetDumper) accept = accept && genMetDumper->fill(iEvent,iSetup);

    if (0) std::cout << "=== MiniAOD2TTreeFilter::filter() genWeightDumper" << std::endl;
    if (genWeightDumper) accept = accept && genWeightDumper->fill(iEvent,iSetup);

    if (0) std::cout << "=== MiniAOD2TTreeFilter::filter() trackDumper" << std::endl;
    if (trackDumper) accept = accept && trackDumper->fill(iEvent,iSetup);

    if (0) std::cout << "=== MiniAOD2TTreeFilter::filter() genParticleDumper" << std::endl;
    if (genParticleDumper) accept = accept && genParticleDumper->fill(iEvent,iSetup);

    if (0) std::cout << "=== MiniAOD2TTreeFilter::filter() genJetDumper" << std::endl;
    if (genJetDumper) accept = accept && genJetDumper->fill(iEvent,iSetup);

    if(accept) Events->Fill();

    return accept;
}

void MiniAOD2TTreeFilter::reset()
{
  if (0) std::cout << "=== MiniAOD2TTreeFilter::reset()" << std::endl;
  if (skimDumper) skimDumper->reset();
  if (trgDumper) trgDumper->reset();
  if (metNoiseFilterDumper) metNoiseFilterDumper->reset();
  if (tauDumper) tauDumper->reset();
  if (electronDumper) electronDumper->reset();
  if (muonDumper) muonDumper->reset();
  if (jetDumper) jetDumper->reset();
  if (topDumper) topDumper->reset();
  if (metDumper) metDumper->reset();
  if (genMetDumper) genMetDumper->reset();
  if (genWeightDumper) genWeightDumper->reset();
  if (trackDumper) trackDumper->reset();
  if (genParticleDumper) genParticleDumper->reset();
  if (genJetDumper) genJetDumper->reset();
  if (0) std::cout << "=== MiniAOD2TTreeFilter::reset() return" << std::endl;
  return;
}

#include <time.h>
#include "TH1F.h"
void MiniAOD2TTreeFilter::endJob(){
  if (0) std::cout << "=== MiniAOD2TTreeFilter::endJob()" << std::endl;
    fOUT->cd();

// write date
    time_t rawtime;
    time (&rawtime);
    TString dateName = "Generated "+TString(ctime(&rawtime));
    TNamed* dateString = new TNamed("","");
    dateString->Write(dateName);

// write commit string
    TString versionName = "Commit "+codeVersion;
    TNamed* versionString = new TNamed("","");
    versionString->Write(versionName);

// write config info
    TDirectory* infodir = fOUT->mkdir("configInfo");
    infodir->cd();

    TNamed* dv = new TNamed("dataVersion",dataVersion);
    dv->Write();

    int nbins = 2;
    TH1F* cfgInfo = new TH1F("configinfo","",nbins,0,nbins);
    cfgInfo->SetBinContent(1,1.0);
    cfgInfo->GetXaxis()->SetBinLabel(1,"control");
    cfgInfo->SetBinContent(2,cmEnergy);
    cfgInfo->GetXaxis()->SetBinLabel(2,"energy");
    cfgInfo->Write();

    if(skimDumper){
      TH1F* skimCounter = skimDumper->getCounter();
      skimCounter->Write();
    }

    fOUT->cd();

// write TTree
    Events->Write();

    std::cout << std::endl << "List of branches:" << std::endl;
    TObjArray* branches = Events->GetListOfBranches();
    for(int i = 0; i < branches->GetEntries(); ++i){
      int hltCounterAll    = 0;
      int hltCounterPassed = 0;
      if (trgDumper) {
	std::pair<int,int> hltCounters = trgDumper->counters(branches->At(i)->GetName());
	if(hltCounters.first > 0) {
	  hltCounterAll    = hltCounters.first;
	  hltCounterPassed = hltCounters.second;
	}
      }
      if(hltCounterAll > 0){
	std::string name(branches->At(i)->GetName());
	while(name.length() < 70) name += " ";
	
	std::cout << "    " << name << " " << hltCounterAll << " " << hltCounterPassed << std::endl;
      }else{
	std::cout << "    " << branches->At(i)->GetName() << std::endl;
      }
    }
    std::cout << "Number of events saved " << Events->GetEntries() << std::endl << std::endl;

// copy PU histogram from separate file (makes merging of root files so much easier)
    if (PUInfoInputFileName.size()) {
      TFile* fPU = TFile::Open(PUInfoInputFileName.c_str());
      if (fPU) {
        // File open is successful
        TH1F* hPU = dynamic_cast<TH1F*>(fPU->Get("pileup"));
        if (hPU) {
          // Histogram exists
          TH1F* hPUclone = dynamic_cast<TH1F*>(hPU->Clone());
          hPUclone->SetDirectory(fOUT);
	  infodir->cd();
          if(dataVersion.find("data") < dataVersion.length()) hPUclone->SetName("pileupPS"); 
          hPUclone->Write();
        }
      }
      fPU->Close();
    }
// in case there is a PU distribution for a prescaled trigger..
    if (PUInfoPSInputFileName.size()) {
      TFile* fPU_PS = TFile::Open(PUInfoPSInputFileName.c_str());
      if (fPU_PS) {
        // File open is successful
        TH1F* hPU_PS = dynamic_cast<TH1F*>(fPU_PS->Get("pileup"));
        if (hPU_PS) {
          // Histogram exists
          TH1F* hPUclone = dynamic_cast<TH1F*>(hPU_PS->Clone());
          hPUclone->SetDirectory(fOUT);
          infodir->cd();
          hPUclone->SetName("pileupPS");
          hPUclone->Write();
        }
      }
      fPU_PS->Close();
    }

// copy top pt weight histogram from separate file (makes merging of root files so much easier)
    if (TopPtInputFileName.size()) {
      TFile* fTopPt = TFile::Open(TopPtInputFileName.c_str());
      if (fTopPt) {
        // File open is successful
        TH1F* hTopPt = dynamic_cast<TH1F*>(fTopPt->Get("topPtWeightAllEvents"));
        if (hTopPt) {
          // Histogram exists
          TH1F* hTopPtClone = dynamic_cast<TH1F*>(hTopPt->Clone());
          hTopPtClone->SetDirectory(fOUT);
          infodir->cd();
          hTopPtClone->Write();
        }
      }
      fTopPt->Close();
    }    
    
// close output file
    fOUT->Close();
}

void MiniAOD2TTreeFilter::endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup) {
    if(skimDumper) skimDumper->fill(iLumi,iSetup);
}

bool MiniAOD2TTreeFilter::isMC(){
    std::regex data_re("data");
    if(std::regex_search(dataVersion, data_re)) return false;
    return true;
}          
DEFINE_FWK_MODULE(MiniAOD2TTreeFilter);
