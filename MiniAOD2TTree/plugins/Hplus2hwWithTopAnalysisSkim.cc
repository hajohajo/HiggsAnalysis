//================================================================================================   
//
// INSTRUCTIONS:
// EDFilter for skimming MiniAOD files using custom selections on all high-level objects;
// electrons, muons, taus, jets, MET
// This specific file is used to select events for H+ -> H0 W+ semi-leptonic channel
//
// 
// USAGE:
// multicrab.py --create [options[
// 
//
// LAST USED:
// cd CMSSW_8_0_30/src/HiggsAnalysis/MiniAOD2TTree/test
// cmsenv
// source /cvmfs/cms.cern.ch/crab3/crab.csh
// multicrab.py --create -s T3_US_FNALLPC -p miniAOD2TTree_Hplus2hwWithTopAnalysisSkim_cfg.py
// 
//
// LINKS:
// https://gitlab.cern.ch/HPlus/HiggsAnalysis/wikis/home
//
//================================================================================================   

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "HiggsAnalysis/MiniAOD2TTree/interface/MiniIsolation.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/RefToBase.h"

#include "FWCore/Framework/interface/LuminosityBlock.h"

#include <iostream>
#include <regex>
#include <string>
#include <vector>
#include <memory>

class Hplus2hwWithTopAnalysisSkim : public edm::EDFilter {

public:
  explicit Hplus2hwWithTopAnalysisSkim(const edm::ParameterSet&);
  ~Hplus2hwWithTopAnalysisSkim();

  virtual bool filter(edm::Event&, const edm::EventSetup& );

private:
  const bool cfg_verbose;
  edm::EDGetTokenT<edm::TriggerResults> cfg_trgResultsToken;
  std::vector<std::string> cfg_triggerBits;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate> > cfg_pfcandsToken;
  edm::EDGetTokenT<edm::View<reco::Vertex> > cfg_vertexToken;
  //
  edm::EDGetTokenT<pat::ElectronCollection> cfg_electronToken;
  edm::EDGetTokenT<double> cfg_electronRhoToken;
  std::string cfg_electronID;
  edm::EDGetTokenT<edm::ValueMap<float> > cfg_electronMVAToken;
  const double cfg_electronMiniRelIsoEA;
  const bool   cfg_electronUseTightID;
  const double cfg_electronPtCut;
  const double cfg_electronEtaCut;
  const int cfg_electronNCut;
  //
  edm::EDGetTokenT<edm::View<pat::Muon> > cfg_muonToken;
  std::string cfg_muonID;
  const double cfg_muonMiniRelIsoEA;
  const double cfg_muonPtCut;
  const double cfg_muonEtaCut;
  const int cfg_muonNCut;
  //
  const int cfg_leptonNCut;
  //
  edm::EDGetTokenT<edm::View<pat::Tau> > cfg_tauToken;
  std::vector<std::string>  cfg_tauDiscriminators;
  const double cfg_tauPtCut;
  const double cfg_tauEtaCut;
  const int cfg_tauNCut;
  //
  edm::EDGetTokenT<edm::View<pat::Jet>> cfg_jetToken;
  std::vector<std::string> cfg_jetUserFloats;
  const double cfg_jetEtCut;
  const double cfg_jetEtaCut;
  const int cfg_jetsNCut;
  //
  // edm::EDGetTokenT<edm::View<pat::MET>> cfg_metToken;
  // const double cfg_MetCut;

  int nEvents;
  int nSelectedEvents;
};

Hplus2hwWithTopAnalysisSkim::Hplus2hwWithTopAnalysisSkim(const edm::ParameterSet& iConfig)
  : cfg_verbose(iConfig.getParameter<bool>("Verbose")),
    cfg_trgResultsToken(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
    cfg_triggerBits(iConfig.getParameter<std::vector<std::string> >("HLTPaths")),
    cfg_pfcandsToken(consumes<edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("PackedCandidatesCollection"))),        
    cfg_vertexToken(consumes<edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("VertexCollection"))),       
    //
    cfg_electronToken(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("ElectronCollection"))),
    cfg_electronRhoToken(consumes<double>(iConfig.getParameter<edm::InputTag>("ElectronRhoSource"))),
    cfg_electronID(iConfig.getParameter<std::string>("ElectronID")),
    cfg_electronMVAToken(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("ElectronMVA"))),
    cfg_electronMiniRelIsoEA(iConfig.getParameter<double>("ElectronMiniRelIsoEA")),
    cfg_electronUseTightID(iConfig.getParameter<bool>("ElectronUseTightID")),
    cfg_electronPtCut(iConfig.getParameter<double>("ElectronPtCut")),
    cfg_electronEtaCut(iConfig.getParameter<double>("ElectronEtaCut")),
    cfg_electronNCut(iConfig.getParameter<int>("ElectronNCut")),
    //
    cfg_muonToken(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("MuonCollection"))),
    cfg_muonID(iConfig.getParameter<std::string>("MuonID")),
    cfg_muonMiniRelIsoEA(iConfig.getParameter<double>("MuonMiniRelIsoEA")),
    cfg_muonPtCut(iConfig.getParameter<double>("MuonPtCut")),
    cfg_muonEtaCut(iConfig.getParameter<double>("MuonEtaCut")),
    cfg_muonNCut(iConfig.getParameter<int>("MuonNCut")),
    //
    cfg_leptonNCut(iConfig.getParameter<int>("LeptonNCut")), // Muon or Electron
    //
    cfg_tauToken(consumes<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("TauCollection"))),
    cfg_tauDiscriminators(iConfig.getParameter<std::vector<std::string> >("TauDiscriminators")),
    cfg_tauPtCut(iConfig.getParameter<double>("TauPtCut")),
    cfg_tauEtaCut(iConfig.getParameter<double>("TauEtaCut")),
    cfg_tauNCut(iConfig.getParameter<int>("TauNCut")),
    //
    cfg_jetToken(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("JetCollection"))),
    // cfg_jetUserFloats(iConfig.getParameter<std::vector<std::string> >("JetUserFloats")),
    cfg_jetEtCut(iConfig.getParameter<double>("JetEtCut")),
    cfg_jetEtaCut(iConfig.getParameter<double>("JetEtaCut")),
    cfg_jetsNCut(iConfig.getParameter<int>("JetNCut")),
    //
    // cfg_metToken(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("METCollection"))),
    // cfg_MetCut(iConfig.getParameter<double>("MetCut")),

    nEvents(0),
    nSelectedEvents(0)
{
  
}


Hplus2hwWithTopAnalysisSkim::~Hplus2hwWithTopAnalysisSkim(){
    double eff = 0;
    if(nEvents > 0) eff = ((double)nSelectedEvents)/((double) nEvents);
    std::cout <<"Hplus2hwWithTopAnalysisSkim: "<<std::endl; //  	edm::LogVerbatim("Hplus2hwWithTopAnalysisSkim") 
    std::cout <<"  Number of events read                  = "<<nEvents<<std::endl;
    std::cout <<"  Number of events kept                  = "<<nSelectedEvents<<std::endl;
    std::cout <<"  Efficiency                             = "<<eff<<std::endl;
}


bool Hplus2hwWithTopAnalysisSkim::filter(edm::Event& iEvent, const edm::EventSetup& iSetup ){

    nEvents++;

    // Trigger bits
    edm::Handle<edm::TriggerResults> trghandle;
    iEvent.getByToken(cfg_trgResultsToken, trghandle);
    std::vector<std::string> firedTriggers;

    // Sanity check
    if(trghandle.isValid())
      {
        edm::TriggerResults tr = *trghandle;
        bool fromPSetRegistry;
        edm::Service<edm::service::TriggerNamesService> tns;
        std::vector<std::string> hlNames; 
        tns->getTrigPaths(tr, hlNames, fromPSetRegistry);
	bool passed = false;
        bool trgBitFound = false;

	if (cfg_verbose*1) std::cout << "=== Trigger bits:" << std::endl;

	// Check that at least 1 trigger bit is defined in python cfg file
	if(cfg_triggerBits.size() > 0)
	  {

	  // For-loop: All user-defined trigger bits
          for(size_t i = 0; i < cfg_triggerBits.size(); ++i)
	    {
	      std::regex hlt_re(cfg_triggerBits[i]);
	      int n = 0;

	      // For-loop: All available trigger bits
	      for(std::vector<std::string>::const_iterator j = hlNames.begin(); j!= hlNames.end(); ++j)
		{
		 
		  // Search with regular expressions
		  if (std::regex_search(*j, hlt_re)) 
		    {
		      
		      if (cfg_verbose*1) std::cout << "\t" << *j << " = " << trghandle->accept(n) << std::endl;
		      trgBitFound = true;

		      if ( trghandle->accept(n) ) firedTriggers.push_back(*j);
		      
		      // Check the trigger bit. Fire or not?
		      if(trghandle->accept(n)) 
			{
			  passed = true;
			  break;
			}
		    }
		  n++;
		}
	    }

	  // Inform user that none of the triggers was found in the miniAOD files
          if(!trgBitFound) 
	    {
	      std::cout << "Skimming with trigger bit, but none of the triggers was found!" << std::endl;
	      std::cout << "Looked for triggers:" << std::endl;
	      for (auto& p: cfg_triggerBits) {
                std::cout << "    " << p << std::endl;
	      }
	      
	      std::cout << "Available triggers in dataset:" << std::endl;
	      for(std::vector<std::string>::const_iterator j = hlNames.begin(); j!= hlNames.end(); ++j)
		{
		  std::cout << "    " << *j << std::endl;
		}
	      exit(1);
	    }
	  }
	else
	  {
	    // If no trigger bits are defined then all events are accepted
	    passed = true;
	  }
	if(!passed) return false; 
      }
    else
      {
	// trghandle.isValid() == False:
	// If the trigger handle is not valid, The event will always be accepted, regardless of whether any of the skim trigger-bits are true or not.
	// So for example, all non-reHLT samples will NOT be skimmed by any trigger. Instead all events will be regarded as having passed the skim trigger.
      }

    /*
    if (cfg_verbose) 
      {
	std::cout << "=== Passed Trigger:\n" << std::endl;
	for (auto& p: cfg_triggerBits) std::cout << "\t" << p << std::endl;
	std::cout << "\n" << std::endl;
      }	
    */
  
    // Packed Candidates
    edm::Handle<edm::View<pat::PackedCandidate> > pfcandHandle;
    iEvent.getByToken(cfg_pfcandsToken, pfcandHandle);

    // Vertex (for Muon ID)
    edm::Handle<edm::View<reco::Vertex> > vertexHandle;
    iEvent.getByToken(cfg_vertexToken, vertexHandle);

    // Electrons
    edm::Handle<pat::ElectronCollection>  electronHandle;
    iEvent.getByToken(cfg_electronToken, electronHandle);
    int nElectrons = 0;
    edm::Handle<edm::ValueMap<float> > electronMVAHandle;
    iEvent.getByToken(cfg_electronMVAToken, electronMVAHandle);
    edm::Handle<double> rhoHandle;
    iEvent.getByToken(cfg_electronRhoToken, rhoHandle);

    // Sanity check    
    if(electronHandle.isValid())
      {
	
	int iEle = -1;
	// For-loop: All electrons
	for (const pat::Electron &obj: *electronHandle)
	  {
	    
	    iEle++;
	    edm::RefToBase<pat::Electron> ref ( edm::Ref<pat::ElectronCollection >(electronHandle, iEle));
	    
	    // Calculate mini relative isolation (RelIso) for the electron with effective area (EA)
	    double ele_miniRelIsoEA = getMiniIsolation_EffectiveArea(pfcandHandle, dynamic_cast<const reco::Candidate *>(&obj), 0.05, 0.2, 10., false, false, *rhoHandle);
	    float ele_mvaID         = (*electronMVAHandle)[ref];
	    float ele_pt            = obj.p4().pt();
	    float ele_absEta        = fabs(obj.p4().eta());
	    bool ele_mvaLooseID     = false;
	    bool ele_mvaTightID     = false;

	    // Determine Loose ID
	    if (ele_absEta <= 0.8 and ele_mvaID >= -0.041) ele_mvaLooseID = true;
	    if (ele_absEta > 0.8 and ele_absEta < 1.479 and ele_mvaLooseID >= 0.383) ele_mvaLooseID = true;
	    if (ele_absEta >= 1.479 and ele_mvaID >= -0.515) ele_mvaLooseID = true;

	    // Determine Tight ID
	    if (ele_absEta <= 0.8 and ele_mvaID >= -0.674) ele_mvaTightID = true;
	    if (ele_absEta > 0.8 and ele_absEta < 1.479 and ele_mvaTightID >= 0.744) ele_mvaTightID = true;
	    if (ele_absEta >= 1.479 and ele_mvaID >= 0.170) ele_mvaTightID = true;
	
	    // Apply selections
	    if (cfg_electronUseTightID)
	      {
		if (!ele_mvaTightID) continue;
	      }
	    else
	      {
		if (!ele_mvaLooseID) continue;	   
	      }
	    if (ele_miniRelIsoEA  > cfg_electronMiniRelIsoEA) continue;
	    if (ele_pt < cfg_electronPtCut) continue;
	    if (ele_absEta > cfg_electronEtaCut) continue;
	    
	    // Count the number of electrons that pass all selections
	    nElectrons++;
	  }
      } 
    // // Apply electron multiplicity selection
    // if (nElectrons < cfg_electronNCut) return false;
    // if (cfg_verbose) std::cout << "=== nElectrons:\n\t" << nElectrons << " < " << cfg_electronNCut << std::endl;

    
    // Muons    
    edm::Handle<edm::View<pat::Muon> > muonHandle;
    iEvent.getByToken(cfg_muonToken, muonHandle);
    std::vector<pat::Muon> selectedMuons;
    std::vector<float>     selectedMuons_Iso;    
    int nMuons = 0;

    // Sanity check
    if(muonHandle.isValid())
      {

      // For-loop: All muons
      for(size_t i = 0; i < muonHandle->size(); ++i) 
	{
	  const pat::Muon& obj = muonHandle->at(i);
	  
	  // bool isGlobal = obj.isGlobalMuon();
	  float muon_pt      = obj.p4().pt();
	  float muon_absEta  = fabs(obj.p4().eta());
	  bool muon_isLoose  = obj.isLooseMuon();
	  bool muon_isMedium = obj.isMediumMuon();
	  bool muon_isTight  = false;
	  if (vertexHandle->size() == 0)
	    {
	      muon_isTight = false;
	    }
	  else
	    {
	      muon_isTight = obj.isTightMuon(vertexHandle->at(0));
	    }
	  
	  // Apply muon selections
	  double muon_miniRelIsoEA = getMiniIsolation_EffectiveArea(pfcandHandle, dynamic_cast<const reco::Candidate *>(&obj), 0.05, 0.2, 10., false, false, *rhoHandle);
	  
	  if (cfg_muonID == "loose" || cfg_muonID == "Loose")
	    {
	      if (muon_isLoose  == false) continue;
	    }
	  else if (cfg_muonID == "medium" || cfg_muonID == "Medium")
	    {
	      if (muon_isMedium == false) continue;
	    }
	  else if (cfg_muonID == "tight" || cfg_muonID == "Tight")
	    {
	      if (muon_isTight == false) continue;
	    }
	  else {
	    throw cms::Exception("config") << "Invalid muonID option '" << cfg_muonID << "'! Options: 'loose', 'medium', 'tight'";
	  }
	  
	  // Apply selections cuts
	  if (muon_miniRelIsoEA > cfg_muonMiniRelIsoEA) continue;
	  if (muon_pt < cfg_muonPtCut) continue;
	  if (muon_absEta > cfg_muonEtaCut) continue;
	  
	  // Save all selected muons
	  selectedMuons.push_back(obj);
	  selectedMuons_Iso.push_back(muon_miniRelIsoEA);
	  
	  // Count the number of muons that pass all selections	  
	  nMuons++;
	}
      }
    // // Apply muon multiplicity selection
    // if (nMuons < cfg_muonNCut) return false; //fixme
    // if (cfg_verbose) std::cout << "=== Passed Muons:\n\t" << nMuons << " < " << cfg_muonNCut << std::endl;


    // Leptons (Electrons or Muons)
    int nLeptons = 0;
    nLeptons = nElectrons + nMuons;

    // Apply lepton multiplicity selection
    if(nLeptons < cfg_leptonNCut) return false;
    else
      {
	// Apply electron multiplicity selection
	if (nElectrons < cfg_electronNCut) return false;
	if (cfg_verbose*1) std::cout << "=== nElectrons:\n\tSelected = " << nElectrons << ", Requested = " << cfg_electronNCut << std::endl;
	
	// Apply muon multiplicity selection
	if (nMuons < cfg_muonNCut) return false; //fixme
	if (cfg_verbose*1) std::cout << "=== Passed Muons:\n\tSelected = " << nMuons << ", Requested =  " << cfg_muonNCut << std::endl;

      }
    if (cfg_verbose*1) std::cout << "=== Passed Leptons:\n\tSelected = " << nLeptons << ", Requested = " << cfg_leptonNCut << std::endl;
    

    // Taus
    edm::Handle<edm::View<pat::Tau> > tauHandle;
    iEvent.getByToken(cfg_tauToken, tauHandle);      
    int nTaus = 0;

    // Sanity check
    if(tauHandle.isValid())
      {
	
      // For-loop: All taus
      for (const pat::Tau &obj: *tauHandle)
	{
	  float tau_pt     = obj.p4().pt();
	  float tau_absEta = fabs(obj.p4().eta());
	  bool tau_discr   = true;

	  // For-loop: All user-defined discriminators
	  for(size_t j=0; j < cfg_tauDiscriminators.size(); ++j) 
	    {
	      tau_discr = tau_discr && obj.tauID(cfg_tauDiscriminators[j]);
	    }

	  // Apply selections cuts
	  if (tau_pt < cfg_tauPtCut) continue;
	  if (tau_absEta > cfg_tauEtaCut) continue;
	  if (!tau_discr) continue;
	  
	  // Count the number of taus that pass all selections
	  nTaus++;
	}
      }
    // Apply tau multiplicity selection
    if (nTaus < cfg_tauNCut) return false;
    if (cfg_verbose*1) std::cout << "=== Passed Taus:\n\tSelected = " << nTaus << ", Requested = " << cfg_tauNCut << std::endl;
    

    // Jets
    edm::Handle<edm::View<pat::Jet> > jethandle;
    iEvent.getByToken(cfg_jetToken, jethandle);
    int nJets = 0;
    // Sanity check
    if(jethandle.isValid())
      {
	
	// For-loop: All jets
	for(size_t i=0; i<jethandle->size(); ++i) 
	  {
	    const pat::Jet& obj = jethandle->at(i);

	    float jet_pt     = obj.p4().pt();
	    float jet_absEta = fabs(obj.p4().eta());
	    bool jet_floats  = true;
	    
	    // For-loop: All floats
	    for(size_t j = 0; j < cfg_jetUserFloats.size(); ++j)
	      {
		if(obj.userFloat(cfg_jetUserFloats[j]) < 0) 
		  {
		    jet_floats = false;
		    break;
		}
	      }
	    

	  // Apply selection cuts
	  if(jet_pt < cfg_jetEtCut) continue;
	  if(jet_absEta > cfg_jetEtaCut) continue;
	  if(!jet_floats) continue;
	  
	  // Count the number of jets that pass all selections
	  nJets++;
	  }
      }
    // Apply jet multiplicity selection
    if(nJets < cfg_jetsNCut) return false;
    if (cfg_verbose*1) std::cout << "=== Passed Jets:\n\tSelected = " << nJets << ", Requested = " << cfg_jetsNCut << std::endl;
    

    // All selections passed
    if (cfg_verbose)
      {
	std::cout << "=== Passed Selections:\n\tLeptons = " << nLeptons << " (Electrons = " << nElectrons << ", Muons = " << nMuons << "), Taus = " << nTaus << ", Jets  = " << nJets  << std::endl;
	for(std::vector<std::string>::const_iterator j = firedTriggers.begin(); j!= firedTriggers.end(); ++j) std::cout << "\t" << *j << std::endl;
	std::cout << "===" << std::endl;
      }

    nSelectedEvents++;
    return true;
}

DEFINE_FWK_MODULE(Hplus2hwWithTopAnalysisSkim);   

