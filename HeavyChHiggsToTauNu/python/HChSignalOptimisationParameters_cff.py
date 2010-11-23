import FWCore.ParameterSet.Config as cms

# WARNING: the trigger path is modified in signalOptimisation_cfg.py depending on
# the data version
trigger = cms.untracked.PSet(
    src = cms.untracked.InputTag("patTriggerEvent"),
    triggers = cms.untracked.vstring("HLT_SingleLooseIsoTau20", #35X/36X MC and Run2010A data (prescaled with values 20,50)
                                     "HLT_SingleIsoTau20_Trk5",
                                     "HLT_SingleIsoTau20_Trk15_MET20",
                                     "HLT_SingleIsoTau20_Trk15_MET25_v3",
                                     "HLT_SingleIsoTau20_Trk15_MET25_v4"
                                     ),
    hltMetCut = cms.untracked.double(30.0),
)
TriggerMETEmulation = cms.untracked.PSet(
    src = cms.untracked.InputTag("patMETs"), # calo MET
    metEmulationCut = cms.untracked.double(0.0)
)

useFactorizedTauID = cms.untracked.bool(False) # only use for QCD. Otherwise set to "False"

import HiggsAnalysis.HeavyChHiggsToTauNu.HChTauIDFactorization_cfi as factorizationParams
tauSelectionBase = cms.untracked.PSet(
    src = cms.untracked.InputTag("selectedPatTausShrinkingConePFTauTauTriggerMatched"),
    selection = cms.untracked.string(""),
    ptCut = cms.untracked.double(30),
    etaCut = cms.untracked.double(2.4),
    leadingTrackPtCut = cms.untracked.double(20),
    rtauCut = cms.untracked.double(0.8), # 0.3 or 0.7. Try 0.8 now
    invMassCut = cms.untracked.double(1.5),
    factorization = factorizationParams.tauIDFactorizationParameters
)

tauSelectionCaloTauCutBased = tauSelectionBase.clone()
tauSelectionCaloTauCutBased.src = cms.untracked.InputTag("selectedPatTausCaloRecoTauTauTriggerMatched")
tauSelectionCaloTauCutBased.selection = cms.untracked.string("CaloTauCutBased")

tauSelectionShrinkingConeCutBased = tauSelectionBase.clone()
tauSelectionShrinkingConeCutBased.src = cms.untracked.InputTag("selectedPatTausShrinkingConePFTauTauTriggerMatched")
tauSelectionShrinkingConeCutBased.selection = cms.untracked.string("ShrinkingConePFTauCutBased")

tauSelectionShrinkingConeTaNCBased = tauSelectionBase.clone()
tauSelectionShrinkingConeTaNCBased.src = cms.untracked.InputTag("selectedPatTausShrinkingConePFTauTauTriggerMatched")
tauSelectionShrinkingConeTaNCBased.selection = cms.untracked.string("ShrinkingConePFTauTaNCBased")

tauSelectionHPSTauBased = tauSelectionBase.clone()
tauSelectionHPSTauBased.src = cms.untracked.InputTag("selectedPatTausHpsPFTauTauTriggerMatched")
tauSelectionHPSTauBased.selection = cms.untracked.string("HPSTauBased")

#tauSelection = tauSelectionShrinkingConeCutBased
#tauSelection = tauSelectionHPSTauBased
#tauSelection = tauSelectionShrinkingConeTaNCBased
tauSelection = tauSelectionCaloTauCutBased


jetSelection = cms.untracked.PSet(
    #src = cms.untracked.InputTag("selectedPatJets"),       # Calo jets
    #src = cms.untracked.InputTag("selectedPatJetsAK5JPT"), # JPT jets 
    src = cms.untracked.InputTag("selectedPatJetsAK5PF"),  # PF jets
    src_met = cms.untracked.InputTag("patMETsPF"), # calo MET 
    cleanTauDR = cms.untracked.double(0.5),
    ptCut = cms.untracked.double(30),
    etaCut = cms.untracked.double(2.4),
    minNumber = cms.untracked.uint32(1),
    METCut = cms.untracked.double(-10.0) # for histogramming purposes of deltaPhi
)

MET = cms.untracked.PSet(
    # src = cms.untracked.InputTag("patMETs"), # calo MET
    src = cms.untracked.InputTag("patMETsPF"), # PF MET
    #src = cms.untracked.InputTag("patMETsTC"), # tc MET
    METCut = cms.untracked.double(-10.0)
)

bTagging = cms.untracked.PSet(
    discriminator = cms.untracked.string("trackCountingHighEffBJetTags"),
    discriminatorCut = cms.untracked.double(2.0),
    ptCut = cms.untracked.double(30),
    etaCut = cms.untracked.double(2.4),
    minNumber = cms.untracked.uint32(0)
)

transverseMassCut = cms.untracked.double(-10)

EvtTopology = cms.untracked.PSet(
    #discriminator = cms.untracked.string("test"),
    #discriminatorCut = cms.untracked.double(0.0),
    #alphaT = cms.untracked.double(-5.00)
    alphaT = cms.untracked.double(-10.0)
)

GlobalElectronVeto = cms.untracked.PSet(
    ElectronCollectionName = cms.untracked.InputTag("selectedPatElectrons"),
    ElectronSelection = cms.untracked.string("simpleEleId90relIso"),
    ElectronPtCut = cms.untracked.double(20.0),
    ElectronEtaCut = cms.untracked.double(2.5)
)

GlobalMuonVeto = cms.untracked.PSet(
    MuonCollectionName = cms.untracked.InputTag("selectedPatMuons"),
    MuonSelection = cms.untracked.string("GlobalMuonPromptTight"),
    MuonPtCut = cms.untracked.double(20.0),
    MuonEtaCut = cms.untracked.double(2.5)
)

fakeMETVeto = cms.untracked.PSet(
  src = MET.src,
  maxDeltaPhi = cms.untracked.double(999.)
)
