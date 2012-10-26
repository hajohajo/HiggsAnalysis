#include "HiggsAnalysis/HeavyChHiggsToTauNu/interface/FakeTauIdentifier.h"
#include "HiggsAnalysis/HeavyChHiggsToTauNu/interface/HistoWrapper.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

namespace HPlus {
  
  FakeTauIdentifier::FakeTauIdentifier(const edm::ParameterSet& iConfig, HPlus::HistoWrapper& histoWrapper, std::string label):
    fSFFakeTauBarrelElectron(iConfig.getUntrackedParameter<double>("scalefactorFakeTauBarrelElectron")),
    fSFFakeTauEndcapElectron(iConfig.getUntrackedParameter<double>("scalefactorFakeTauEndcapElectron")),
    fSFFakeTauBarrelMuon(iConfig.getUntrackedParameter<double>("scalefactorFakeTauBarrelMuon")),
    fSFFakeTauEndcapMuon(iConfig.getUntrackedParameter<double>("scalefactorFakeTauEndcapMuon")),
    fSFFakeTauBarrelJet(iConfig.getUntrackedParameter<double>("scalefactorFakeTauBarrelJet")),
    fSFFakeTauEndcapJet(iConfig.getUntrackedParameter<double>("scalefactorFakeTauEndcapJet")),
    fSystematicsFakeTauBarrelElectron(iConfig.getUntrackedParameter<double>("systematicsFakeTauBarrelElectron")),
    fSystematicsFakeTauEndcapElectron(iConfig.getUntrackedParameter<double>("systematicsFakeTauEndcapElectron")),
    fSystematicsFakeTauBarrelMuon(iConfig.getUntrackedParameter<double>("systematicsFakeTauBarrelMuon")),
    fSystematicsFakeTauEndcapMuon(iConfig.getUntrackedParameter<double>("systematicsFakeTauEndcapMuon")),
    fSystematicsFakeTauBarrelJet(iConfig.getUntrackedParameter<double>("systematicsFakeTauBarrelJet")),
    fSystematicsFakeTauEndcapJet(iConfig.getUntrackedParameter<double>("systematicsFakeTauEndcapJet"))
  {
    edm::Service<TFileService> fs;
    // Create histograms
    TFileDirectory myDir = fs->mkdir("FakeTauIdentifier_"+label);
    hTauMatchType = histoWrapper.makeTH<TH1F>(HistoWrapper::kVital, myDir, "TauMatchType", "TauMatchType", kkNumberOfSelectedTauMatchTypes, 0, kkNumberOfSelectedTauMatchTypes);
    if (hTauMatchType->isActive()) {
      hTauMatchType->GetXaxis()->SetBinLabel(1+kkNoMC, "NoMatch");
      hTauMatchType->GetXaxis()->SetBinLabel(1+kkElectronToTau, "e#rightarrow#tau");
      hTauMatchType->GetXaxis()->SetBinLabel(1+kkElectronFromTauDecayToTau, "tau#rightarrowe#rightarrow#tau");
      hTauMatchType->GetXaxis()->SetBinLabel(1+kkMuonToTau, "#mu#rightarrow#tau");
      hTauMatchType->GetXaxis()->SetBinLabel(1+kkMuonFromTauDecayToTau, "tau#rightarrow#mu#rightarrow#tau");
      hTauMatchType->GetXaxis()->SetBinLabel(1+kkTauToTau, "genuine #tau");
      hTauMatchType->GetXaxis()->SetBinLabel(1+kkJetToTau, "jet#rightarrow#tau");
      hTauMatchType->GetXaxis()->SetBinLabel(1+kkElectronToTauAndTauOutsideAcceptance, "e#rightarrow#tau, #tau outside");
      hTauMatchType->GetXaxis()->SetBinLabel(1+kkElectronFromTauDecayToTauAndTauOutsideAcceptance, "#tau#rightarrowe#rightarrow#tau, #tau outside");
      hTauMatchType->GetXaxis()->SetBinLabel(1+kkMuonToTauAndTauOutsideAcceptance, "#mu#rightarrow#tau, #tau outside");
      hTauMatchType->GetXaxis()->SetBinLabel(1+kkMuonFromTauDecayToTauAndTauOutsideAcceptance, "#tau#rightarrow#mu#rightarrow#tau, #tau outside");
      hTauMatchType->GetXaxis()->SetBinLabel(1+kkTauToTauAndTauOutsideAcceptance, "genuine #tau, #tau outside");
      hTauMatchType->GetXaxis()->SetBinLabel(1+kkJetToTauAndTauOutsideAcceptance, "jet#rightarrow#tau, #tau outside");
    }
    hTauOrigin = histoWrapper.makeTH<TH1F>(HistoWrapper::kVital, myDir, "TauOrigin", "TauOrigin", 7, 0, 7);
    if (hTauOrigin->isActive()) {
      hTauOrigin->GetXaxis()->SetBinLabel(1+kkUnknownOrigin, "unknown");
      hTauOrigin->GetXaxis()->SetBinLabel(1+kkFromW, "from W");
      hTauOrigin->GetXaxis()->SetBinLabel(1+kkFromZ, "from Z");
      hTauOrigin->GetXaxis()->SetBinLabel(1+kkFromHplus, "from H+");
      hTauOrigin->GetXaxis()->SetBinLabel(1+kkFromWTau, "from W#rightarrow#tau");
      hTauOrigin->GetXaxis()->SetBinLabel(1+kkFromZTauTau, "from Z#rightarrow#tautau");
      hTauOrigin->GetXaxis()->SetBinLabel(1+kkFromHplusTau, "from H+#rightarrow#tau");
    }
    hMuOrigin = histoWrapper.makeTH<TH1F>(HistoWrapper::kVital, myDir, "MuOrigin", "MuOrigin", 7, 0, 7);
    if (hMuOrigin->isActive()) {
      hMuOrigin->GetXaxis()->SetBinLabel(1+kkUnknownOrigin, "unknown");
      hMuOrigin->GetXaxis()->SetBinLabel(1+kkFromW, "from W");
      hMuOrigin->GetXaxis()->SetBinLabel(1+kkFromZ, "from Z");
      hMuOrigin->GetXaxis()->SetBinLabel(1+kkFromHplus, "from H+");
      hMuOrigin->GetXaxis()->SetBinLabel(1+kkFromWTau, "from W#rightarrow#tau");
      hMuOrigin->GetXaxis()->SetBinLabel(1+kkFromZTauTau, "from Z#rightarrow#tautau");
      hMuOrigin->GetXaxis()->SetBinLabel(1+kkFromHplusTau, "from H+#rightarrow#tau");
    }
    hElectronOrigin = histoWrapper.makeTH<TH1F>(HistoWrapper::kVital, myDir, "ElectronOrigin", "ElectronOrigin", 7, 0, 7);
    if (hElectronOrigin->isActive()) {
      hElectronOrigin->GetXaxis()->SetBinLabel(1+kkUnknownOrigin, "unknown");
      hElectronOrigin->GetXaxis()->SetBinLabel(1+kkFromW, "from W");
      hElectronOrigin->GetXaxis()->SetBinLabel(1+kkFromZ, "from Z");
      hElectronOrigin->GetXaxis()->SetBinLabel(1+kkFromHplus, "from H+");
      hElectronOrigin->GetXaxis()->SetBinLabel(1+kkFromWTau, "from W#rightarrow#tau");
      hElectronOrigin->GetXaxis()->SetBinLabel(1+kkFromZTauTau, "from Z#rightarrow#tautau");
      hElectronOrigin->GetXaxis()->SetBinLabel(1+kkFromHplusTau, "from H+#rightarrow#tau");
    }
  }
  
  FakeTauIdentifier::~FakeTauIdentifier() {
    
  }
  
  FakeTauIdentifier::MCSelectedTauMatchType FakeTauIdentifier::matchTauToMC(const edm::Event& iEvent, const reco::Candidate& tau) {
    FakeTauIdentifier::MCSelectedTauMatchType myMatchType = kkNoMC;
    if (iEvent.isRealData()) return kkNoMC;
    bool foundMCTauOutsideAcceptanceStatus = false;
    bool isMCTau = false;
    bool isMCElectron = false;
    bool isMCMuon = false;
    size_t myTauIndex = 0;
    size_t myElectronIndex = 0;
    size_t myMuIndex = 0;

    edm::Handle <reco::GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
    //std::cout << "matchfinding:" << std::endl;
    for (size_t i=0; i < genParticles->size(); ++i) {
      const reco::Candidate & p = (*genParticles)[i];
      if (std::abs(p.pdgId()) == 11 || std::abs(p.pdgId()) == 13 || std::abs(p.pdgId()) == 15) {
        // Check match with tau
        if (reco::deltaR(p, tau.p4()) < 0.1) {
          if (p.pt() > 10.) {
            //std::cout << "  match found, pid=" << p.pdgId() << " eta=" << std::abs(p.eta()) << " pt=" << p.pt() << std::endl;
            if (std::abs(p.pdgId()) == 15) {
              isMCTau = true;
              myTauIndex = i;
            } else if (std::abs(p.pdgId()) == 11) {
              isMCElectron = true;
              myElectronIndex = i;
            } else if (std::abs(p.pdgId()) == 13) {
              isMCMuon = true;
              myMuIndex = i;
            }
          }
        }
        // Check if there is a tau outside the acceptance in the event
        if (!foundMCTauOutsideAcceptanceStatus && std::abs(p.pdgId()) == 15) {
          if (p.pt() < 40 || abs(p.eta()) > 2.1)
            foundMCTauOutsideAcceptanceStatus = true;
        }
      }
    }
    if (!foundMCTauOutsideAcceptanceStatus) {
      if (isMCElectron) {
        if (isMCTau) {
          myMatchType = kkElectronFromTauDecayToTau;
        } else {
          myMatchType = kkElectronToTau;
        }
      } else if (isMCMuon) {
        if (isMCTau) {
          myMatchType = kkMuonFromTauDecayToTau;
        } else {
          myMatchType = kkMuonToTau;
        }
      } else if (isMCTau) myMatchType = kkTauToTau;
      else myMatchType = kkJetToTau;
    } else {
      if (isMCElectron) {
        if (isMCTau) {
          myMatchType = kkElectronFromTauDecayToTauAndTauOutsideAcceptance;
        } else {
          myMatchType = kkElectronToTauAndTauOutsideAcceptance;
        }
      } else if (isMCMuon) {
        if (isMCTau) {
          myMatchType = kkMuonFromTauDecayToTauAndTauOutsideAcceptance;
        } else {
          myMatchType = kkMuonToTauAndTauOutsideAcceptance;
        }
      }
      else if (isMCTau) myMatchType = kkTauToTauAndTauOutsideAcceptance;
      else myMatchType = kkJetToTauAndTauOutsideAcceptance;
    }
    // Look at ancestor information
    MCSelectedTauOriginType myOriginType = kkUnknownOrigin;
    if (myMatchType != kkJetToTau) {
      size_t myIndex = 0;
      if (isMCElectron)
        myIndex = myElectronIndex;
      else if (isMCMuon)
        myIndex = myMuIndex;
      else if (isMCTau)
        myIndex = myTauIndex;
      bool myWStatus = false;
      bool myZStatus = false;
      bool myHPlusStatus = false;
      bool myTauStatus = false;
      //reco::Candidate* p = const_cast<reco::Candidate*>((*genParticles)[myIndex].mother());
      const reco::Candidate* p = (*genParticles)[myIndex].mother();
      while (p) {
        if (std::abs(p->pdgId()) == 15)
          myTauStatus = true;
        if (std::abs(p->pdgId()) == 24)
          myWStatus = true;
        if (std::abs(p->pdgId()) == 23)
          myZStatus = true;
        if (std::abs(p->pdgId()) == 37)
          myHPlusStatus = true;
        // move to next
        p = p->mother();
      }
      if (!myTauStatus) {
        if (myWStatus) myOriginType = kkFromW;
        else if (myZStatus) myOriginType = kkFromZ;
        else if (myHPlusStatus) myOriginType = kkFromHplus;
      } else {
        if (myWStatus) myOriginType = kkFromWTau;
        else if (myZStatus) myOriginType = kkFromZTauTau;
        else if (myHPlusStatus) myOriginType = kkFromHplusTau; 
      }
    }
    
    // fill histograms
    hTauMatchType->Fill(myMatchType);
    if (isMCElectron)
      hElectronOrigin->Fill(myOriginType);
    else if (isMCMuon)
      hMuOrigin->Fill(myOriginType);
    else if (isMCTau)
      hTauOrigin->Fill(myOriginType);
    
    return myMatchType;
  }

  double FakeTauIdentifier::getFakeTauScaleFactor(FakeTauIdentifier::MCSelectedTauMatchType matchType, double eta) {
    if (matchType == FakeTauIdentifier::kkElectronToTau || matchType == FakeTauIdentifier::kkElectronToTauAndTauOutsideAcceptance ||
        matchType == FakeTauIdentifier::kkElectronFromTauDecayToTau || matchType == FakeTauIdentifier::kkElectronFromTauDecayToTauAndTauOutsideAcceptance) {
      if (std::fabs(eta) < 1.5) {
        return fSFFakeTauBarrelElectron;
      } else {
        return fSFFakeTauEndcapElectron;
      }
    } else if (matchType == FakeTauIdentifier::kkMuonToTau || matchType == FakeTauIdentifier::kkMuonToTauAndTauOutsideAcceptance ||
               matchType == FakeTauIdentifier::kkMuonFromTauDecayToTau || matchType == FakeTauIdentifier::kkMuonFromTauDecayToTauAndTauOutsideAcceptance) {
      if (std::fabs(eta) < 1.5) {
        return fSFFakeTauBarrelMuon;
      } else {
        return fSFFakeTauEndcapMuon;
      }
    } else if (matchType == FakeTauIdentifier::kkJetToTau || matchType == FakeTauIdentifier::kkJetToTauAndTauOutsideAcceptance) {
      if (std::fabs(eta) < 1.5) {
        return fSFFakeTauBarrelJet;
      } else {
        return fSFFakeTauEndcapJet;
      }
    }
    return 1.0;
  }

  double FakeTauIdentifier::getFakeTauSystematics(MCSelectedTauMatchType matchType, double eta) {
    if (matchType == FakeTauIdentifier::kkElectronToTau || matchType == FakeTauIdentifier::kkElectronToTauAndTauOutsideAcceptance ||
        matchType == FakeTauIdentifier::kkElectronFromTauDecayToTau || matchType == FakeTauIdentifier::kkElectronFromTauDecayToTauAndTauOutsideAcceptance) {
      if (std::fabs(eta) < 1.5) {
        return fSystematicsFakeTauBarrelElectron;
      } else {
        return fSystematicsFakeTauEndcapElectron;
      }
    } else if (matchType == FakeTauIdentifier::kkMuonToTau || matchType == FakeTauIdentifier::kkMuonToTauAndTauOutsideAcceptance ||
               matchType == FakeTauIdentifier::kkMuonFromTauDecayToTau || matchType == FakeTauIdentifier::kkMuonFromTauDecayToTauAndTauOutsideAcceptance) {
      if (std::fabs(eta) < 1.5) {
        return fSystematicsFakeTauBarrelMuon;
      } else {
        return fSystematicsFakeTauEndcapMuon;
      }
    } else if (matchType == FakeTauIdentifier::kkJetToTau || matchType == FakeTauIdentifier::kkJetToTauAndTauOutsideAcceptance) {
      if (std::fabs(eta) < 1.5) {
        return fSystematicsFakeTauBarrelJet;
      } else {
        return fSystematicsFakeTauEndcapJet;
      }
    }
    return 0.0;
  }

}
