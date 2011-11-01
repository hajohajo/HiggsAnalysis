#include "HiggsAnalysis/HeavyChHiggsToTauNu/interface/SignalAnalysis.h"
#include "HiggsAnalysis/HeavyChHiggsToTauNu/interface/TransverseMass.h"
#include "HiggsAnalysis/HeavyChHiggsToTauNu/interface/DeltaPhi.h"
#include "HiggsAnalysis/HeavyChHiggsToTauNu/interface/EvtTopology.h"
#include "HiggsAnalysis/HeavyChHiggsToTauNu/interface/MakeTH.h"

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
#include "TNamed.h"

namespace {
  template <typename T>
  void copyPtrToVector(const edm::PtrVector<T>& src, std::vector<T>& dst) {
    dst.reserve(src.size());
    for(typename edm::PtrVector<T>::const_iterator i = src.begin(); i != src.end(); ++i) {
      dst.push_back(**i);
    }
  }
}

namespace HPlus {
  SignalAnalysis::CounterGroup::CounterGroup(EventCounter& eventCounter) :
    fOneTauCounter(eventCounter.addCounter("nonQCDType2:taus == 1")),
    fElectronVetoCounter(eventCounter.addCounter("nonQCDType2:electron veto")),
    fMuonVetoCounter(eventCounter.addCounter("nonQCDType2:muon veto")),
    fMETCounter(eventCounter.addCounter("nonQCDType2:MET")),
    fNJetsCounter(eventCounter.addCounter("nonQCDType2:njets")),
    fBTaggingCounter(eventCounter.addCounter("nonQCDType2:btagging")),
    fDeltaPhiCounter(eventCounter.addCounter("nonQCDType2:deltaphi")),
    fFakeMETVetoCounter(eventCounter.addCounter("nonQCDType2:fake MET veto")),
    fTopSelectionCounter(eventCounter.addCounter("nonQCDType2:Top Selection cut")) { }
  SignalAnalysis::CounterGroup::CounterGroup(EventCounter& eventCounter, std::string prefix) :
    fOneTauCounter(eventCounter.addSubCounter(prefix,":taus == 1")),
    fElectronVetoCounter(eventCounter.addSubCounter(prefix,":electron veto")),
    fMuonVetoCounter(eventCounter.addSubCounter(prefix,":muon veto")),
    fMETCounter(eventCounter.addSubCounter(prefix,":MET")),
    fNJetsCounter(eventCounter.addSubCounter(prefix,":njets")),
    fBTaggingCounter(eventCounter.addSubCounter(prefix,":btagging")),
    fDeltaPhiCounter(eventCounter.addSubCounter(prefix,":deltaphi")),
    fFakeMETVetoCounter(eventCounter.addSubCounter(prefix,":fake MET veto")),
    fTopSelectionCounter(eventCounter.addSubCounter(prefix,":Top Selection cut")) { }
  SignalAnalysis::CounterGroup::~CounterGroup() { }

  SignalAnalysis::SignalAnalysis(const edm::ParameterSet& iConfig, EventCounter& eventCounter, EventWeight& eventWeight):
    fEventWeight(eventWeight),
    //    fmetEmulationCut(iConfig.getUntrackedParameter<double>("metEmulationCut")),
    fAllCounter(eventCounter.addCounter("All events")),
    fTriggerCounter(eventCounter.addCounter("Trigger and HLT_MET cut")),
    fPrimaryVertexCounter(eventCounter.addCounter("primary vertex")),
    fTausExistCounter(eventCounter.addCounter("taus > 0")),
    fOneTauCounter(eventCounter.addCounter("taus == 1")),
    fRtauAfterTauIDCounter(eventCounter.addCounter("RtauAfterTauID")),
    fElectronVetoCounter(eventCounter.addCounter("electron veto")),
    fMuonVetoCounter(eventCounter.addCounter("muon veto")),
    fNJetsCounter(eventCounter.addCounter("njets")),
    fMETCounter(eventCounter.addCounter("MET")),
    fBTaggingCounter(eventCounter.addCounter("btagging")),
    fDeltaPhiTauMETCounter(eventCounter.addCounter("DeltaPhi(Tau,MET) upper limit")),
    fdeltaPhiTauMET10Counter(eventCounter.addCounter("deltaPhiTauMET lower limit")),
    fdeltaPhiTauMET160Counter(eventCounter.addCounter("deltaPhiTauMET upper limit")),
    fFakeMETVetoCounter(eventCounter.addCounter("fake MET veto")),
    fdeltaPhiTauMET160FakeMetCounter(eventCounter.addCounter("deltaPhi160 and fake MET veto")),

    fRtauAfterCutsCounter(eventCounter.addCounter("RtauAfterCuts")),
    fForwardJetVetoCounter(eventCounter.addCounter("forward jet veto")),
    ftransverseMassCut80Counter(eventCounter.addCounter("transverseMass > 60")),
    ftransverseMassCut100Counter(eventCounter.addCounter("transverseMass > 80")),
    ftransverseMassCut80NoRtauCounter(eventCounter.addCounter("transverseMass > 60 no Rtau")),
    ftransverseMassCut100NoRtauCounter(eventCounter.addCounter("transverseMass > 80 no Rtau")),
    fZmassVetoCounter(eventCounter.addCounter("ZmassVetoCounter")),
    fTopSelectionCounter(eventCounter.addCounter("Top Selection cut")),
    ftransverseMassCut100TopCounter(eventCounter.addCounter("transverseMass > 100 top cut")),
    fTriggerSelection(iConfig.getUntrackedParameter<edm::ParameterSet>("trigger"), eventCounter, eventWeight),
    fPrimaryVertexSelection(iConfig.getUntrackedParameter<edm::ParameterSet>("primaryVertexSelection"), eventCounter, eventWeight),
    fGlobalElectronVeto(iConfig.getUntrackedParameter<edm::ParameterSet>("GlobalElectronVeto"), eventCounter, eventWeight),
    fGlobalMuonVeto(iConfig.getUntrackedParameter<edm::ParameterSet>("GlobalMuonVeto"), eventCounter, eventWeight),
    fOneProngTauSelection(iConfig.getUntrackedParameter<edm::ParameterSet>("tauSelection"), eventCounter, eventWeight, 1, "tauID"),
    fJetSelection(iConfig.getUntrackedParameter<edm::ParameterSet>("jetSelection"), eventCounter, eventWeight),
    fMETSelection(iConfig.getUntrackedParameter<edm::ParameterSet>("MET"), eventCounter, eventWeight, "MET"),
    fBTagging(iConfig.getUntrackedParameter<edm::ParameterSet>("bTagging"), eventCounter, eventWeight),
    fFakeMETVeto(iConfig.getUntrackedParameter<edm::ParameterSet>("fakeMETVeto"), eventCounter, eventWeight),
    fJetTauInvMass(iConfig.getUntrackedParameter<edm::ParameterSet>("jetTauInvMass"), eventCounter, eventWeight),
    fTopSelection(iConfig.getUntrackedParameter<edm::ParameterSet>("topSelection"), eventCounter, eventWeight),
    //    ftransverseMassCut(iConfig.getUntrackedParameter<edm::ParameterSet>("transverseMassCut")),
    fGenparticleAnalysis(iConfig.getUntrackedParameter<edm::ParameterSet>("GenParticleAnalysis"), eventCounter, eventWeight),
    fForwardJetVeto(iConfig.getUntrackedParameter<edm::ParameterSet>("forwardJetVeto"), eventCounter, eventWeight),
    fCorrelationAnalysis(eventCounter, eventWeight),
    fEvtTopology(iConfig.getUntrackedParameter<edm::ParameterSet>("EvtTopology"), eventCounter, eventWeight),
    fTriggerEfficiencyScaleFactor(iConfig.getUntrackedParameter<edm::ParameterSet>("triggerEfficiencyScaleFactor"), fEventWeight),
    fVertexWeight(iConfig.getUntrackedParameter<edm::ParameterSet>("vertexWeight")),
    fTree(iConfig.getUntrackedParameter<edm::ParameterSet>("Tree"), fBTagging.getDiscriminator()),
    // Non-QCD Type II related
    fNonQCDTypeIIGroup(eventCounter),
    fAllTausCounterGroup(eventCounter, "All"),
    fElectronToTausCounterGroup(eventCounter, "e->tau"),
    fMuonToTausCounterGroup(eventCounter, "mu->tau"),
    fGenuineToTausCounterGroup(eventCounter, "tau->tau"),
    fJetToTausCounterGroup(eventCounter, "jet->tau"),
    fAllTausAndTauOutsideAcceptanceCounterGroup(eventCounter, "All with tau outside acceptance"),
    fElectronToTausAndTauOutsideAcceptanceCounterGroup(eventCounter, "e->tau with tau outside acceptance"),
    fMuonToTausAndTauOutsideAcceptanceCounterGroup(eventCounter, "mu->tau with tau outside acceptance"),
    fGenuineToTausAndTauOutsideAcceptanceCounterGroup(eventCounter, "tau->tau with tau outside acceptance"),
    fJetToTausAndTauOutsideAcceptanceCounterGroup(eventCounter, "jet->tau with tau outside acceptance"),
    fProduce(iConfig.getUntrackedParameter<bool>("produceCollections", false))
  {
    edm::Service<TFileService> fs;
    // Save the module configuration to the output ROOT file as a TNamed object
    fs->make<TNamed>("parameterSet", iConfig.dump().c_str());

    // Book histograms filled in the analysis body
    
    // Vertex histograms
    TFileDirectory myVertexDir = fs->mkdir("Vertices");
    hVerticesBeforeWeight = makeTH<TH1F>(myVertexDir, "verticesBeforeWeight", "Number of vertices without weighting", 40, 0, 40);
    hVerticesAfterWeight = makeTH<TH1F>(myVertexDir, "verticesAfterWeight", "Number of vertices with weighting", 40, 0, 40);
    hVerticesTriggeredBeforeWeight = makeTH<TH1F>(myVertexDir, "verticesTriggeredBeforeWeight", "Number of vertices without weighting", 40, 0, 40);
    hVerticesTriggeredAfterWeight = makeTH<TH1F>(myVertexDir, "verticesTriggeredAfterWeight", "Number of vertices with weighting", 40, 0, 40);
    //    hmetAfterTrigger = makeTH<TH1F>(*fs, "metAfterTrigger", "metAfterTrigger", 50, 0., 200.);
    
    hTransverseMass = makeTH<TH1F>(*fs, "transverseMass", "transverseMass;m_{T}(tau,MET), GeV/c^{2};N_{events} / 10 GeV/c^{2}", 400, 0., 400.);
    hTransverseMassAfterDeltaPhi = makeTH<TH1F>(*fs, "transverseMass", "transverseMass;m_{T}(tau,MET), GeV/c^{2};N_{events} / 10 GeV/c^{2}", 400, 0., 400.);
    hNonQCDTypeIITransverseMass = makeTH<TH1F>(*fs, "transverseMass", "transverseMass;m_{T}(tau,MET), GeV/c^{2};N_{events} / 10 GeV/c^{2}", 400, 0., 400.);
    hNonQCDTypeIITransverseMassAfterDeltaPhi = makeTH<TH1F>(*fs, "transverseMass", "transverseMass;m_{T}(tau,MET), GeV/c^{2};N_{events} / 10 GeV/c^{2}", 400, 0., 400.);
    
    hDeltaPhi = makeTH<TH1F>(*fs, "deltaPhi", "deltaPhi;#Delta#phi(tau,MET);N_{events} / 10 degrees", 360, 0., 180.);
    hDeltaPhiJetMet = makeTH<TH1F>(*fs, "deltaPhiJetMet", "deltaPhiJetMet", 400, 0., 3.2);  
    hAlphaT = makeTH<TH1F>(*fs, "alphaT", "alphaT", 100, 0.0, 5.0);
    hAlphaTInvMass = makeTH<TH1F>(*fs, "alphaT-InvMass", "alphaT-InvMass", 100, 0.0, 1000.0);    
    hAlphaTVsRtau = makeTH<TH2F>(*fs, "alphaT(y)-Vs-Rtau(x)", "alphaT-Vs-Rtau",  120, 0.0, 1.2, 500, 0.0, 5.0);
    //    hMet_AfterTauSelection = makeTH<TH1F>(*fs, "met_AfterTauSelection", "met_AfterTauSelection", 100, 0.0, 400.0);
    //    hMet_BeforeTauSelection = makeTH<TH1F>(*fs, "met_BeforeTauSelection", "met_BeforeTauSelection", 100, 0.0, 400.0);
    
    TFileDirectory mySelectedTauDir = fs->mkdir("SelectedTau");
    hSelectedTauEt = makeTH<TH1F>(mySelectedTauDir, "SelectedTau_pT_AfterTauID", "SelectedTau_pT_AfterTauID;#tau p_{T}, GeV/c;N_{events} / 10 GeV/c", 400, 0.0, 400.0);
    //    hSelectedTauEtMetCut = makeTH<TH1F>(*fs, "SelectedTau_pT_AfterTauID_MetCut", "SelectedTau_pT_AfterTauID_MetCut;#tau p_{T}, GeV/c;N_{events} / 10 GeV/c", 400, 0.0, 400.0);
    hSelectedTauEta = makeTH<TH1F>(mySelectedTauDir, "SelectedTau_eta_AfterTauID", "SelectedTau_eta_AfterTauID;#tau #eta;N_{events} / 0.1", 300, -3.0, 3.0);
    hSelectedTauEtAfterCuts = makeTH<TH1F>(mySelectedTauDir, "SelectedTau_pT_AfterCuts", "SelectedTau_pT_AfterCuts;#tau p_{T}, GeV/c;N_{events} / 10 GeV/c", 40, 0.0, 400.0);
    hSelectedTauEtaAfterCuts = makeTH<TH1F>(mySelectedTauDir, "SelectedTau_eta_AfterCuts", "SelectedTau_eta_AfterCuts;#tau #eta;N_{events} / 0.1", 30, -3.0, 3.0);
    hSelectedTauPhi = makeTH<TH1F>(mySelectedTauDir, "SelectedTau_phi_AfterTauID", "SelectedTau_eta_AfterTauID;#tau #eta;N_{events} / 0.087", 360, -3.1415926, 3.1415926);
    hSelectedTauRtau = makeTH<TH1F>(mySelectedTauDir, "SelectedTau_Rtau_AfterTauID", "SelectedTau_Rtau_AfterTauID;R_{#tau};N_{events} / 0.1", 360, 0., 1.2);
    hSelectedTauRtauAfterCuts = makeTH<TH1F>(mySelectedTauDir, "SelectedTau_Rtau_AfterCuts", "SelectedTau_Rtau_AfterCuts;R_{#tau};N_{events} / 0.1", 360, 0., 1.2);
    hSelectedTauLeadingTrackPt = makeTH<TH1F>(mySelectedTauDir, "SelectedTau_TauLeadingTrackPt", "SelectedTau_TauLeadingTrackPt;#tau p_{T}, GeV/c;N_{events} / 10 GeV/c", 400, 0.0, 400.0);
    hNonQCDTypeIISelectedTauEtAfterCuts = makeTH<TH1F>(mySelectedTauDir, "NonQCDTypeII_SelectedTau_pT_AfterCuts", "SelectedTau_pT_AfterCuts;#tau p_{T}, GeV/c;N_{events} / 10 GeV/c", 40, 0.0, 400.0);
    hNonQCDTypeIISelectedTauEtaAfterCuts = makeTH<TH1F>(mySelectedTauDir, "NonQCDTypeII_SelectedTau_eta_AfterCuts", "SelectedTau_eta_AfterCuts;#tau #eta;N_{events} / 0.1", 30, -3.0, 3.0);

    hMetAfterCuts = makeTH<TH1F>(*fs, "Met_AfterCuts", "Met_AfterCuts", 500, 0.0, 500.0);

    //    hMetAfterCuts = makeTH<TH1F>(*fs, "Met_AfterCuts", "Met_AfterCuts", 400, 0.0, 400.0);
    
    hSelectionFlow = makeTH<TH1F>(*fs, "SignalSelectionFlow", "SignalSelectionFlow;;N_{events}", 8, 0, 8);
    hSelectionFlow->GetXaxis()->SetBinLabel(1+kSignalOrderTrigger,"Trigger");
    //hSelectionFlow->GetXaxis()->SetBinLabel(1+kSignalOrderVertexSelection,"Vertex");
    hSelectionFlow->GetXaxis()->SetBinLabel(1+kSignalOrderTauID,"#tau ID");
    hSelectionFlow->GetXaxis()->SetBinLabel(1+kSignalOrderElectronVeto,"Isol. e veto");
    hSelectionFlow->GetXaxis()->SetBinLabel(1+kSignalOrderMuonVeto,"Isol. #mu veto");
    hSelectionFlow->GetXaxis()->SetBinLabel(1+kSignalOrderMETSelection,"MET");
    hSelectionFlow->GetXaxis()->SetBinLabel(1+kSignalOrderJetSelection,"#geq 3 jets");
    hSelectionFlow->GetXaxis()->SetBinLabel(1+kSignalOrderBTagSelection,"#geq 1 b jet");
    hSelectionFlow->GetXaxis()->SetBinLabel(1+kSignalOrderDeltaPhiSelection,"#Delta#phi(#tau,MET) cut");
    //hSelectionFlow->GetXaxis()->SetBinLabel(1+kSignalOrderFakeMETVeto,"Further QCD rej.");
    //hSelectionFlow->GetXaxis()->SetBinLabel(1+kSignalOrderTopSelection,"Top mass");

    hEMFractionAll = makeTH<TH1F>(*fs, "NonQCDTypeII_FakeTau_EMFraction_All", "FakeTau_EMFraction_All", 22, 0., 1.1);
    hEMFractionElectrons = makeTH<TH1F>(*fs, "NonQCDTypeII_FakeTau_EMFraction_Electrons", "FakeTau_EMFraction_Electrons", 22, 0., 1.1);   

    fTree.init(*fs);
  }

  SignalAnalysis::~SignalAnalysis() { }

  void SignalAnalysis::produces(edm::EDFilter *producer) const {
    if(fProduce) {
      producer->produces<std::vector<pat::Tau> >("selectedTaus");
      producer->produces<std::vector<pat::Jet> >("selectedJets");
      producer->produces<std::vector<pat::Jet> >("selectedBJets");
      producer->produces<std::vector<pat::Electron> >("selectedVetoElectrons");
      producer->produces<std::vector<pat::Muon> >("selectedVetoMuonsBeforeIsolation");
      producer->produces<std::vector<pat::Muon> >("selectedVetoMuons");
    }
  }


  bool SignalAnalysis::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    fEventWeight.updatePrescale(iEvent); // set prescale
    fTree.setPrescaleWeight(fEventWeight.getWeight());

    // Vertex weight
    std::pair<double, size_t> weightSize = fVertexWeight.getWeightAndSize(iEvent, iSetup);
    if(!iEvent.isRealData()) {
      fEventWeight.multiplyWeight(weightSize.first);
      fTree.setPileupWeight(weightSize.first);
    }
    hVerticesBeforeWeight->Fill(weightSize.second);
    hVerticesAfterWeight->Fill(weightSize.second, fEventWeight.getWeight());
    fTree.setNvertices(weightSize.second);

    increment(fAllCounter);
    
    // Apply trigger and HLT_MET cut or trigger parametrisation
    TriggerSelection::Data triggerData = fTriggerSelection.analyze(iEvent, iSetup);
    if (!triggerData.passedEvent()) return false;
    increment(fTriggerCounter);
    hSelectionFlow->Fill(kSignalOrderTrigger, fEventWeight.getWeight());
    if(triggerData.hasTriggerPath()) // protection if TriggerSelection is disabled
      fTree.setHltTaus(triggerData.getTriggerTaus());

    hVerticesTriggeredBeforeWeight->Fill(weightSize.second);
    hVerticesTriggeredAfterWeight->Fill(weightSize.second, fEventWeight.getWeight());

    // GenParticle analysis (must be done here when we effectively trigger all MC)
    if (!iEvent.isRealData()) {
      GenParticleAnalysis::Data genData = fGenparticleAnalysis.analyze(iEvent, iSetup);
      fTree.setGenMET(genData.getGenMET());
    }

    // Primary vertex
    VertexSelection::Data pvData = fPrimaryVertexSelection.analyze(iEvent, iSetup);
    if(!pvData.passedEvent()) return false;
    increment(fPrimaryVertexCounter);
    //hSelectionFlow->Fill(kSignalOrderVertexSelection, fEventWeight.getWeight());

    // TauID
    TauSelection::Data tauData = fOneProngTauSelection.analyze(iEvent, iSetup);
    if(!tauData.passedEvent()) return false; // Require at least one tau
    // plot leading track without pt cut
    hSelectedTauLeadingTrackPt->Fill(tauData.getSelectedTaus()[0]->leadPFChargedHadrCand()->pt(), fEventWeight.getWeight());
    increment(fTausExistCounter);
    if(tauData.getSelectedTaus().size() != 1) return false; // Require exactly one tau
    // Apply trigger scale factor here, because it depends only on tau
    TriggerEfficiencyScaleFactor::Data triggerWeight = fTriggerEfficiencyScaleFactor.applyEventWeight(*(tauData.getSelectedTaus()[0]));
    fTree.setTriggerWeight(triggerWeight.getEventWeight());
    increment(fOneTauCounter);
    hSelectionFlow->Fill(kSignalOrderTauID, fEventWeight.getWeight());
    if(fProduce) {
      std::auto_ptr<std::vector<pat::Tau> > saveTaus(new std::vector<pat::Tau>());
      copyPtrToVector(tauData.getSelectedTaus(), *saveTaus);
      iEvent.put(saveTaus, "selectedTaus");
    }
    hSelectedTauRtau->Fill(tauData.getRtauOfSelectedTau(), fEventWeight.getWeight());  
    hSelectedTauLeadingTrackPt->Fill(tauData.getSelectedTaus()[0]->leadPFChargedHadrCand()->pt(), fEventWeight.getWeight());
    hSelectedTauEt->Fill(tauData.getSelectedTaus()[0]->pt(), fEventWeight.getWeight());
    hSelectedTauEta->Fill(tauData.getSelectedTaus()[0]->eta(), fEventWeight.getWeight());
    hSelectedTauPhi->Fill(tauData.getSelectedTaus()[0]->phi(), fEventWeight.getWeight());
    // Obtain MC matching - for EWK without genuine taus
    MCSelectedTauMatchType myTauMatch = matchTauToMC(iEvent, tauData.getSelectedTaus()[0]);
    bool myTypeIIStatus = false; // True if the selected tau is a fake
    if (myTauMatch == kkNoMC) {
      if (!(myTauMatch == kkTauToTau || myTauMatch == kkTauToTauAndTauOutsideAcceptance)) {
        myTypeIIStatus = true;
      }
    }
    fAllTausCounterGroup.incrementOneTauCounter();
    fillNonQCDTypeIICounters(myTauMatch, kSignalOrderTauID, tauData);
    if (myTauMatch == kkElectronToTau)
      hEMFractionElectrons->Fill(tauData.getSelectedTaus()[0]->emFraction());
    hEMFractionAll->Fill(tauData.getSelectedTaus()[0]->emFraction());

    
    // Global electron veto
    GlobalElectronVeto::Data electronVetoData = fGlobalElectronVeto.analyze(iEvent, iSetup);
    if (!electronVetoData.passedEvent()) return false;
    increment(fElectronVetoCounter);
    hSelectionFlow->Fill(kSignalOrderElectronVeto, fEventWeight.getWeight());
    fillNonQCDTypeIICounters(myTauMatch, kSignalOrderElectronVeto, tauData);
    if(fProduce) {
      std::auto_ptr<std::vector<pat::Electron> > saveElectrons(new std::vector<pat::Electron>());
      copyPtrToVector(electronVetoData.getSelectedElectrons(), *saveElectrons);
      iEvent.put(saveElectrons, "selectedVetoElectrons");
    }


    // Global muon veto
    GlobalMuonVeto::Data muonVetoData = fGlobalMuonVeto.analyze(iEvent, iSetup, pvData.getSelectedVertex());
    if (!muonVetoData.passedEvent()) return false;
    increment(fMuonVetoCounter);
    hSelectionFlow->Fill(kSignalOrderMuonVeto, fEventWeight.getWeight());
    fillNonQCDTypeIICounters(myTauMatch, kSignalOrderMuonVeto, tauData);
    if(fProduce) {
      std::auto_ptr<std::vector<pat::Muon> > saveMuons(new std::vector<pat::Muon>());
      copyPtrToVector(muonVetoData.getSelectedMuonsBeforeIsolation(), *saveMuons);
      iEvent.put(saveMuons, "selectedVetoMuonsBeforeIsolation");
      saveMuons.reset(new std::vector<pat::Muon>());
      copyPtrToVector(muonVetoData.getSelectedMuons(), *saveMuons);
      iEvent.put(saveMuons, "selectedVetoMuons");
    }


    // Hadronic jet selection
    JetSelection::Data jetData = fJetSelection.analyze(iEvent, iSetup, tauData.getSelectedTaus()[0]); 
    if(!jetData.passedEvent()) return false;
    increment(fNJetsCounter);
    hSelectionFlow->Fill(kSignalOrderJetSelection, fEventWeight.getWeight());
    fillNonQCDTypeIICounters(myTauMatch, kSignalOrderJetSelection, tauData);
    if(fProduce) {
      std::auto_ptr<std::vector<pat::Jet> > saveJets(new std::vector<pat::Jet>());
      copyPtrToVector(jetData.getSelectedJets(), *saveJets);
      iEvent.put(saveJets, "selectedJets");
    }


    // MET, no event cut
    METSelection::Data metData = fMETSelection.analyze(iEvent, iSetup);
    double transverseMass = TransverseMass::reconstruct(*(tauData.getSelectedTaus()[0]), *(metData.getSelectedMET()) );
    // b tagging, no event cut
    BTagging::Data btagData = fBTagging.analyze(iEvent, iSetup, jetData.getSelectedJets());
    // Top reco, no event cut
    TopSelection::Data TopSelectionData = fTopSelection.analyze(iEvent, iSetup, jetData.getSelectedJets(), btagData.getSelectedJets());
    // Calculate alphaT
    EvtTopology::Data evtTopologyData = fEvtTopology.analyze(*(tauData.getSelectedTaus()[0]), jetData.getSelectedJets()); 

    if(metData.getRawMET().isNonnull())
      fTree.setRawMET(metData.getRawMET());
    if(metData.getType1MET().isNonnull())
      fTree.setType1MET(metData.getType1MET());
    if(metData.getType2MET().isNonnull())
      fTree.setType2MET(metData.getType2MET());
    if(metData.getCaloMET().isNonnull())
      fTree.setCaloMET(metData.getCaloMET());
    if(metData.getTcMET().isNonnull())
      fTree.setTcMET(metData.getTcMET());

    // Write the stuff to the tree
    fTree.setFillWeight(fEventWeight.getWeight());
    fTree.setBTagging(btagData.passedEvent(), btagData.getScaleFactor());
    fTree.setTop(TopSelectionData.getTopP4());
    fTree.fill(iEvent, tauData.getSelectedTaus(), jetData.getSelectedJets(),
               evtTopologyData.alphaT().fAlphaT);


    // MET cut
    if(!metData.passedEvent()) return false;
    increment(fMETCounter);
    hSelectionFlow->Fill(kSignalOrderMETSelection, fEventWeight.getWeight());
    fillNonQCDTypeIICounters(myTauMatch, kSignalOrderMETSelection, tauData);


    // b tagging cut
    if(!btagData.passedEvent()) return false;
    // Apply scale factor as weight to event
    btagData.fillScaleFactorHistograms(); // Important!!! Needs to be called before scale factor is applied as weight to the event; Uncertainty is determined from these histograms
    fEventWeight.multiplyWeight(btagData.getScaleFactor());
    increment(fBTaggingCounter);
    hSelectionFlow->Fill(kSignalOrderBTagSelection, fEventWeight.getWeight());
    fillNonQCDTypeIICounters(myTauMatch, kSignalOrderBTagSelection, tauData);
    if(fProduce) {
      std::auto_ptr<std::vector<pat::Jet> > saveBJets(new std::vector<pat::Jet>());
      copyPtrToVector(btagData.getSelectedJets(), *saveBJets);
      iEvent.put(saveBJets, "selectedBJets");
    }

    // Fill transverse mass histograms    
    hTransverseMass->Fill(transverseMass, fEventWeight.getWeight());
    if (myTypeIIStatus) hNonQCDTypeIITransverseMass->Fill(transverseMass, fEventWeight.getWeight());

    
    // Delta phi(tau,MET) cut
    double deltaPhi = DeltaPhi::reconstruct(*(tauData.getSelectedTaus()[0]), *(metData.getSelectedMET())) * 57.3; // converted to degrees
    hDeltaPhi->Fill(deltaPhi, fEventWeight.getWeight());
    if ( deltaPhi > 10) 
      increment(fdeltaPhiTauMET10Counter); 
    if ( deltaPhi < 160)
      increment(fdeltaPhiTauMET160Counter);
    if (deltaPhi > 150) return false;
    increment(fDeltaPhiTauMETCounter);
    fillNonQCDTypeIICounters(myTauMatch, kSignalOrderDeltaPhiSelection, tauData);
 
    // Fill transverse mass histograms after Deltaphi cut
    hTransverseMassAfterDeltaPhi->Fill(transverseMass, fEventWeight.getWeight());
    if (myTypeIIStatus) hNonQCDTypeIITransverseMassAfterDeltaPhi->Fill(transverseMass, fEventWeight.getWeight());

    
    
    // -------------------------------------------
    // Experimental cuts, counters, and histograms
    
    // plot deltaPhi(jet,met)
    for(edm::PtrVector<pat::Jet>::const_iterator iJet = jetData.getSelectedJets().begin(); iJet != jetData.getSelectedJets().end(); ++iJet) {
      deltaPhi = DeltaPhi::reconstruct(**iJet, *(metData.getSelectedMET()));
      hDeltaPhiJetMet->Fill(deltaPhi*57.3, fEventWeight.getWeight());
    }

    hSelectedTauRtauAfterCuts->Fill(tauData.getRtauOfSelectedTau(), fEventWeight.getWeight());
    hSelectedTauEtAfterCuts->Fill(tauData.getSelectedTaus()[0]->pt(), fEventWeight.getWeight());
    hSelectedTauEtaAfterCuts->Fill(tauData.getSelectedTaus()[0]->eta(), fEventWeight.getWeight());
    hMetAfterCuts->Fill(metData.getSelectedMET()->et(), fEventWeight.getWeight());

   // top mass with possible event cuts
    if (TopSelectionData.passedEvent()&& tauData.getRtauOfSelectedTau() > 0.8 ) {
      increment(fTopSelectionCounter);
      //      hSelectionFlow->Fill(kSignalOrderTopSelection, fEventWeight.getWeight());      
      if(transverseMass > 80 ) increment(ftransverseMassCut100TopCounter);   
    } 

    // Fake MET veto a.k.a. further QCD suppression
    FakeMETVeto::Data fakeMETData = fFakeMETVeto.analyze(iEvent, iSetup, tauData.getSelectedTaus(), jetData.getSelectedJets(), metData.getSelectedMET());
    if (fakeMETData.passedEvent() ) {
      increment(fFakeMETVetoCounter);
      
      if ( deltaPhi < 160) {
        increment(fdeltaPhiTauMET160FakeMetCounter);
        hTransverseMassDeltaPhiUpperCutFakeMet->Fill(transverseMass, fEventWeight.getWeight());  
      } 
    }
    //hSelectionFlow->Fill(kSignalOrderFakeMETVeto, fEventWeight.getWeight());
    fillNonQCDTypeIICounters(myTauMatch, kSignalOrderFakeMETVeto, tauData);

    if(transverseMass > 80 )
      increment(ftransverseMassCut80Counter);
    if(transverseMass > 100 )
      increment(ftransverseMassCut100Counter);

    
    // Correlation analysis
    fCorrelationAnalysis.analyze(tauData.getSelectedTaus(), btagData.getSelectedJets());
    // Alpha T
    //if(!evtTopologyData.passedEvent()) return false;
    //    EvtTopology::AlphaStruc sAlphaT = evtTopologyData.alphaT();
    //    hAlphaT->Fill(sAlphaT.fAlphaT, fEventWeight.getWeight()); // FIXME: move this histogramming to evt topology


   
    // Forward jet veto                                                                                                                                                                                                           
    //    ForwardJetVeto::Data forwardJetData = fForwardJetVeto.analyze(iEvent, iSetup);
    //    if (!forwardJetData.passedEvent()) return false;
    //    increment(fForwardJetVetoCounter);

    return true;
  }

  SignalAnalysis::MCSelectedTauMatchType SignalAnalysis::matchTauToMC(const edm::Event& iEvent, const edm::Ptr<pat::Tau> tau) {
    if (iEvent.isRealData()) return kkNoMC;
    bool foundMCTauOutsideAcceptanceStatus = false;
    bool isMCTau = false;
    bool isMCElectron = false;
    bool isMCMuon = false;

    edm::Handle <reco::GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
    //std::cout << "matchfinding:" << std::endl;
    for (size_t i=0; i < genParticles->size(); ++i) {
      const reco::Candidate & p = (*genParticles)[i];
      if (std::abs(p.pdgId()) == 11 || std::abs(p.pdgId()) == 13 || std::abs(p.pdgId()) == 15) {
        // Check match with tau
        if (reco::deltaR(p, tau->p4()) < 0.1) {
          if (p.pt() > 10.) {
            //std::cout << "  match found, pid=" << p.pdgId() << " eta=" << std::abs(p.eta()) << " pt=" << p.pt() << std::endl;
            if (std::abs(p.pdgId()) == 11) isMCElectron = true;
            if (std::abs(p.pdgId()) == 13) isMCMuon = true;
            if (std::abs(p.pdgId()) == 15) isMCTau = true;
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
      if (isMCElectron) return kkElectronToTau;
      if (isMCMuon) return kkMuonToTau;
      if (isMCTau) return kkTauToTau;
      return kkJetToTau;
    }
    if (isMCElectron) return kkElectronToTauAndTauOutsideAcceptance;
    if (isMCMuon) return kkMuonToTauAndTauOutsideAcceptance;
    if (isMCTau) return kkTauToTauAndTauOutsideAcceptance;
    return kkJetToTauAndTauOutsideAcceptance;
  }

  SignalAnalysis::CounterGroup* SignalAnalysis::getCounterGroupByTauMatch(MCSelectedTauMatchType tauMatch) {
    if (tauMatch == kkElectronToTau) return &fElectronToTausCounterGroup;
    else if (tauMatch == kkMuonToTau) return &fMuonToTausCounterGroup;
    else if (tauMatch == kkTauToTau) return &fGenuineToTausCounterGroup;
    else if (tauMatch == kkJetToTau) return &fJetToTausCounterGroup;
    else if (tauMatch == kkElectronToTauAndTauOutsideAcceptance) return &fElectronToTausAndTauOutsideAcceptanceCounterGroup;
    else if (tauMatch == kkMuonToTauAndTauOutsideAcceptance) return &fMuonToTausAndTauOutsideAcceptanceCounterGroup;
    else if (tauMatch == kkTauToTauAndTauOutsideAcceptance) return &fGenuineToTausAndTauOutsideAcceptanceCounterGroup;
    else if (tauMatch == kkJetToTauAndTauOutsideAcceptance) return &fJetToTausAndTauOutsideAcceptanceCounterGroup;
    return 0;
  }
  
  void SignalAnalysis::fillNonQCDTypeIICounters(MCSelectedTauMatchType tauMatch, SignalSelectionOrder selection, const TauSelection::Data& tauData) {
    // Get out if no match has been found
    if (tauMatch == kkNoMC) return;
    // Obtain status for main counter
    bool myTypeIIStatus = true;
    // Define event as type II if no genuine tau was identified as the selected tau
    if (tauMatch == kkTauToTau || tauMatch == kkTauToTauAndTauOutsideAcceptance)
        myTypeIIStatus = false;
    // Fill main and subcounter for the selection
    if (selection == kSignalOrderTauID) {
      if (myTypeIIStatus) fNonQCDTypeIIGroup.incrementOneTauCounter();
      getCounterGroupByTauMatch(tauMatch)->incrementOneTauCounter();
    } else if (selection == kSignalOrderMETSelection) {
      if (myTypeIIStatus) fNonQCDTypeIIGroup.incrementMETCounter();
      getCounterGroupByTauMatch(tauMatch)->incrementMETCounter();
    } else if (selection == kSignalOrderElectronVeto) {
      if (myTypeIIStatus) fNonQCDTypeIIGroup.incrementElectronVetoCounter();
      getCounterGroupByTauMatch(tauMatch)->incrementElectronVetoCounter();
    } else if (selection == kSignalOrderMuonVeto) {
      if (myTypeIIStatus) fNonQCDTypeIIGroup.incrementMuonVetoCounter();
      getCounterGroupByTauMatch(tauMatch)->incrementMuonVetoCounter();
    } else if (selection == kSignalOrderJetSelection) {
      if (myTypeIIStatus) fNonQCDTypeIIGroup.incrementNJetsCounter();
      getCounterGroupByTauMatch(tauMatch)->incrementNJetsCounter();
    } else if (selection == kSignalOrderBTagSelection) {
      if (myTypeIIStatus) {
        fNonQCDTypeIIGroup.incrementBTaggingCounter();
        // Fill histograms
        hNonQCDTypeIISelectedTauEtAfterCuts->Fill(tauData.getSelectedTaus()[0]->pt(), fEventWeight.getWeight());
        hNonQCDTypeIISelectedTauEtaAfterCuts->Fill(tauData.getSelectedTaus()[0]->eta(), fEventWeight.getWeight());
      }
      getCounterGroupByTauMatch(tauMatch)->incrementBTaggingCounter();
    } else if (selection == kSignalOrderDeltaPhiSelection) {
      if (myTypeIIStatus) fNonQCDTypeIIGroup.incrementDeltaPhiCounter();
      getCounterGroupByTauMatch(tauMatch)->incrementDeltaPhiCounter();
    } else if (selection == kSignalOrderFakeMETVeto) {
      if (myTypeIIStatus) fNonQCDTypeIIGroup.incrementFakeMETVetoCounter();
      getCounterGroupByTauMatch(tauMatch)->incrementFakeMETVetoCounter();
    } else if (selection == kSignalOrderTopSelection) {
      if (myTypeIIStatus) fNonQCDTypeIIGroup.incrementTopSelectionCounter();
      getCounterGroupByTauMatch(tauMatch)->incrementTopSelectionCounter();
    }
  }
}
