import FWCore.ParameterSet.Config as cms

skim = cms.EDFilter("Hplus2hwWithTopAnalysisSkim",
                    Verbose        = cms.bool(False),
                    TriggerResults = cms.InputTag("TriggerResults::HLT"),
                    HLTPaths       = cms.vstring(# SingleMuon Primary Dataset (PD)
                                                 "HLT_IsoMu24_v", 
                                                 "HLT_IsoTkMu24_v",
                                                 # SingleElectron Primary Dataset (PD)
                                                 "HLT_Ele27_WPTight_Gsf_v",
                                                 # DoubleMuon Primary Dataset (PD)
                                                 "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
                                                 "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
                                                 "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
                                                 "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
                                                 # DoubleEG Primary Dataset (PD)
                                                 "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
                                                 "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
                                                 "HLT_DoubleEle33_CaloIdL_v",
                                                 # MuonEG Primary Dataset (PD)
                                                 "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
                                                 "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
                                                 "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
                                                 "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
                                                 "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
                                                 "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
                                                 ),
                    ### PF Candidates
                    PackedCandidatesCollection = cms.InputTag("packedPFCandidates"),
                    ### Vertex
                    VertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                    ### Electrons (https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2)
                    ElectronCollection   = cms.InputTag("slimmedElectrons"),
                    ElectronRhoSource    = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
                    ElectronID           = cms.string("cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                    ElectronMVA          = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
                    ElectronMiniRelIsoEA = cms.double(1.0),
                    ElectronUseTightID   = cms.bool(False),     # True = Tight MVA WP, False = Loose MVA WP
                    ElectronPtCut        = cms.double(10.0),
                    ElectronEtaCut       = cms.double(3.0),
                    ElectronNCut         = cms.int32(0),        # This is checked ONLY if the "LeptonNCut" is satisfied
                    ### Muons (https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2)
                    MuonCollection       = cms.InputTag("slimmedMuons"),
                    MuonID               = cms.string("Loose"), # [options: Loose, Medium, Tight]
                    MuonMiniRelIsoEA     = cms.double(1.0),     # [default: 0.4]
                    MuonPtCut            = cms.double(10.0),
                    MuonEtaCut           = cms.double(3.0),
                    MuonNCut             = cms.int32(0),        # This is checked ONLY if the "LeptonNCut" is satisfied
                    ### Leptons (Electrons or Muons)
                    LeptonNCut           = cms.int32(0),        # Concerns the above-defined electrons and muon selections
                    ### Taus
                    TauCollection        = cms.InputTag("slimmedTaus"),
                    TauDiscriminators    = cms.vstring("decayModeFinding",
                                                       "byVLooseIsolationMVArun2v1DBoldDMwLT",),
                    TauPtCut             = cms.double(20),
                    TauEtaCut            = cms.double(3.0),
                    TauNCut              = cms.int32(0),
                    ### Jets (https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data)
                    JetCollection        = cms.InputTag("slimmedJets"),
                    JetUserFloats        = cms.vstring(),       # DISABLED
                    # JetUserFloats        = cms.vstring("pileupJetId:fullDiscriminant",),
                    JetEtCut             = cms.double(30.0),
                    JetEtaCut            = cms.double(2.4),
                    JetNCut              = cms.int32(0),
                    )
