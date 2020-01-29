'''
For miniAOD instructions see: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015 
'''
#================================================================================================  
# Import Modules
#================================================================================================  
import FWCore.ParameterSet.Config as cms
import HiggsAnalysis.MiniAOD2TTree.tools.git as git
from HiggsAnalysis.MiniAOD2TTree.tools.HChOptions import getOptionsDataVersion


#================================================================================================  
# Options
#================================================================================================  
maxEvents    = 1000
maxWarnings  = 100
reportEvery  = 1
testWithData = True

# Define ROOT files for local testing
testData = {}
testData["SingleElectron"] = [
    '/store/data/Run2016B/SingleElectron/MINIAOD/03Feb2017_ver2-v2/110001/FE93EA39-A2EB-E611-A1E5-0CC47AD99176.root',                                                                                   
    '/store/data/Run2016B/SingleElectron/MINIAOD/03Feb2017_ver2-v2/110001/FE23EF13-84EB-E611-B99A-0025905A6070.root',                                                                                   
    '/store/data/Run2016B/SingleElectron/MINIAOD/03Feb2017_ver2-v2/110001/FAEF3272-79EB-E611-84C3-0CC47AACFCDE.root',   
    # '/store/data/Run2016G/SingleElectron/MINIAOD/03Feb2017-v1/50001/FEB9169C-39EB-E611-9E94-0CC47AD99116.root',
    # '/store/data/Run2016G/SingleElectron/MINIAOD/03Feb2017-v1/50001/FE9CA978-3CEB-E611-8617-90B11C27E14D.root',
    # '/store/data/Run2016G/SingleElectron/MINIAOD/03Feb2017-v1/50001/FCB0EDF4-28EB-E611-ABAF-A0000420FE80.root',
    # '/store/data/Run2016H/SingleElectron/MINIAOD/PromptReco-v2/000/284/035/00000/F063D9E7-449F-E611-9C38-02163E014291.root',
    # '/store/data/Run2016H/SingleElectron/MINIAOD/PromptReco-v2/000/284/035/00000/EA3918D4-449F-E611-8B15-FA163E8769C6.root',
    # '/store/data/Run2016H/SingleElectron/MINIAOD/PromptReco-v2/000/284/035/00000/DEBB0DF5-449F-E611-909F-02163E01425C.root',
    ]
testData["SingleMuon"] = [
    '/store/data/Run2016G/SingleMuon/MINIAOD/03Feb2017-v1/820000/DA8E43DC-61F0-E611-8411-70106F4A94F0.root',
    '/store/data/Run2016G/SingleMuon/MINIAOD/03Feb2017-v1/820000/C037802E-62F0-E611-94AA-70106F48BBEE.root',
    '/store/data/Run2016G/SingleMuon/MINIAOD/03Feb2017-v1/820000/AE7F0E4C-62F0-E611-A831-0CC47A7FC7B8.root',
    #'/store/data/Run2016H/SingleMuon/MINIAOD/PromptReco-v2/000/284/035/00000/F2E41147-4D9F-E611-A17C-FA163E8E25D7.root',
    #'/store/data/Run2016H/SingleMuon/MINIAOD/PromptReco-v2/000/284/035/00000/F29156D1-4E9F-E611-BD9C-FA163EDC8AC8.root',
    #'/store/data/Run2016H/SingleMuon/MINIAOD/PromptReco-v2/000/284/035/00000/E6006CCE-4C9F-E611-9197-02163E01391D.root',
    ]
testData["DoubleEG"] = [
    '/store/data/Run2016G/DoubleEG/MINIAOD/03Feb2017-v1/80000/FEA99FAC-F6EA-E611-9BE7-003048F5B69C.root',
    '/store/data/Run2016G/DoubleEG/MINIAOD/03Feb2017-v1/80000/FE5CB43C-EAEA-E611-BF19-008CFA582BF4.root',
    '/store/data/Run2016G/DoubleEG/MINIAOD/03Feb2017-v1/80000/FE5134DE-BAEA-E611-A9F6-0CC47A4D7678.root',
    #'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/281/707/00000/2C258D6C-C386-E611-866E-02163E014732.root',
    #'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/281/707/00000/1E4E52DD-CC86-E611-920F-02163E0142E7.root',
    #'/store/data/Run2016H/DoubleEG/MINIAOD/PromptReco-v2/000/281/707/00000/0E695BBC-C586-E611-A229-02163E01249F.root',
    ]
testData["DoubleMuon"] = [
    '/store/data/Run2016G/DoubleMuon/MINIAOD/03Feb2017-v1/100000/FE8491B3-F5EA-E611-9708-001E67A3AEB8.root',
    '/store/data/Run2016G/DoubleMuon/MINIAOD/03Feb2017-v1/100000/FE0B09B4-78EB-E611-956C-ECB1D7B67E10.root',
    '/store/data/Run2016G/DoubleMuon/MINIAOD/03Feb2017-v1/100000/FCD1A6F4-58EB-E611-8A1D-002481CFE5C4.root',
    '/store/data/Run2016G/DoubleMuon/MINIAOD/03Feb2017-v1/100000/FAB6B333-EEEA-E611-B18E-90B11C066D31.root',
    #'/store/data/Run2016H/DoubleMuon/MINIAOD/PromptReco-v2/000/284/035/00000/EEA9D66F-449F-E611-A9B5-02163E0144C7.root',
    #'/store/data/Run2016H/DoubleMuon/MINIAOD/PromptReco-v2/000/284/035/00000/B4ABD959-449F-E611-8766-FA163EE30F6B.root',
    #'/store/data/Run2016H/DoubleMuon/MINIAOD/PromptReco-v2/000/284/035/00000/94DAB68C-439F-E611-BDF9-FA163E7CCCBE.root',
    ]
testData["MuonEG"] = [    
    '/store/data/Run2016G/MuonEG/MINIAOD/03Feb2017-v1/80000/FA796197-7DEB-E611-957E-0025904C66EC.root',
    '/store/data/Run2016G/MuonEG/MINIAOD/03Feb2017-v1/80000/F65CA7FE-6EEB-E611-8BB8-0025905C3D6C.root',
    '/store/data/Run2016G/MuonEG/MINIAOD/03Feb2017-v1/80000/E023AA96-7DEB-E611-9003-0025905C3D98.root',
    #'/store/data/Run2016H/MuonEG/MINIAOD/PromptReco-v2/000/281/727/00000/FA84435B-FA86-E611-A578-02163E0139CB.root',
    #'/store/data/Run2016H/MuonEG/MINIAOD/PromptReco-v2/000/281/727/00000/D4CCAFBD-E686-E611-B6CA-FA163E5B52C5.root',
    #'/store/data/Run2016H/MuonEG/MINIAOD/PromptReco-v2/000/281/727/00000/86A5CBD5-FF86-E611-AE99-FA163EF1D56D.root',
    ]
testData["TT"] = [
        '/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root',
        '/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0806AB92-99BE-E611-9ECD-0025905A6138.root',
        '/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/165F54A0-A3BE-E611-B3F7-0025905A606A.root',
        '/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/18E31463-B3BE-E611-B6A3-0CC47A4D7678.root',
        ]
testData["QCD"] = [
     '/store/mc/RunIISummer16MiniAODv2/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/110000/0055B499-54B6-E611-9F86-FA163E1F94C5.root',
     '/store/mc/RunIISummer16MiniAODv2/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/110000/02B77462-7CB5-E611-A061-0025905B8568.root',
     '/store/mc/RunIISummer16MiniAODv2/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/110000/02CDB360-51B5-E611-A568-002590747E0E.root',
     '/store/mc/RunIISummer16MiniAODv2/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/110000/0486046D-0BB6-E611-B533-002590D9D8C0.root',
     ]
testData["ZJetsToQQ"] = [
    '/store/mc/RunIISummer16MiniAODv2/ZJetsToQQ_HT600toInf_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/16DC0526-F4FA-E611-938E-6CC2173BBA40.root',
    '/store/mc/RunIISummer16MiniAODv2/ZJetsToQQ_HT600toInf_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/18DFBF7B-B5FC-E611-80D2-002481DE49B6.root',
    '/store/mc/RunIISummer16MiniAODv2/ZJetsToQQ_HT600toInf_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/20C84FF5-11FB-E611-AA3C-C4346BC70B58.root'
    ]

if testWithData:
    dataVersion  = "80Xdata"
    datasetFiles = testData["SingleElectron"]
else:
    dataVersion  = "80Xmc" 
    datasetFiles = testData["TT"] # "TT", "ZJetsToQQ"
    
# For debugging purposes
debug       = False
RunNum_1    = 1
LumiBlock_1 = 2
EvtNum_1    = 240
RunNum_2    = 1
LumiBlock_2 = 2
EvtNum_2    = 260


#================================================================================================  
# Setup the Process
#================================================================================================  
process = cms.Process("TTreeDump")
process.options = cms.untracked.PSet(
    #SkipEvent         = cms.untracked.vstring('ProductNotFound'),
    wantSummary       = cms.untracked.bool(debug),
    printDependencies = cms.untracked.bool(debug),
)
process.maxEvents = cms.untracked.PSet(
    input         = cms.untracked.int32(maxEvents)
    )


#================================================================================================  
# Tracer service is for debugging purposes (Tells user what cmsRun is accessing)
#================================================================================================  
if debug:
    process.Tracer = cms.Service("Tracer")


#================================================================================================  
# Setup the Process
#================================================================================================  
process.load("FWCore/MessageService/MessageLogger_cfi")
process.MessageLogger.categories.append("TriggerBitCounter")
process.MessageLogger.cerr.FwkReport.reportEvery = reportEvery # print the event number for every 100th event
process.MessageLogger.cerr.TriggerBitCounter = cms.untracked.PSet(limit = cms.untracked.int32(maxWarnings)) # print max 100 warnings

#================================================================================================  
# Set the process options -- Display summary at the end, enable unscheduled execution
#================================================================================================  
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
    )

#================================================================================================  
# Define the input files 
#================================================================================================  
process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring('/store/data/Run2016B/JetHT/MINIAOD/PromptReco-v2/000/273/150/00000/66051AAF-D819-E611-BD3D-02163E011D55.root',)
                            fileNames = cms.untracked.vstring(datasetFiles)
                            #eventsToProcess = cms.untracked.VEventRange('%s:%s:%s-%s:%s:%s' % (RunNum_1, LumiBlock_1, EvtNum_1, RunNum_2, LumiBlock_2, EvtNum_2) ), 
                            )


#================================================================================================  
# Get Dataset version and Global tag
#================================================================================================  
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
options, dataVersion = getOptionsDataVersion(dataVersion)
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, str(dataVersion.getGlobalTag()), '')


#================================================================================================  
# Print Settings
#================================================================================================  
msgAlign = "{:<10} {:<55} {:<25} {:<25}"
title    =  msgAlign.format("Data", "Global Tag", "Trigger Source", "Trigger Tag")
print "="*len(title)
print title
print "="*len(title)
print msgAlign.format(dataVersion.version, dataVersion.getGlobalTag(), dataVersion.getMETFilteringProcess(), dataVersion.getTriggerProcess())
print 


#================================================================================================
# Load processes
#================================================================================================
process.load("HiggsAnalysis/MiniAOD2TTree/PUInfo_cfi")
process.load("HiggsAnalysis/MiniAOD2TTree/TopPt_cfi")
process.load("HiggsAnalysis/MiniAOD2TTree/Tau_cfi")
process.load("HiggsAnalysis/MiniAOD2TTree/Electron_cfi")
process.load("HiggsAnalysis/MiniAOD2TTree/Muon_cfi")
process.load("HiggsAnalysis/MiniAOD2TTree/Jet_cfi")
process.load("HiggsAnalysis/MiniAOD2TTree/FatJet_cfi")
process.load("HiggsAnalysis/MiniAOD2TTree/Top_cfi")
process.load("HiggsAnalysis/MiniAOD2TTree/MET_cfi")
process.load("HiggsAnalysis/MiniAOD2TTree/METNoiseFilter_cfi")

process.METNoiseFilter.triggerResults = cms.InputTag("TriggerResults::"+str(dataVersion.getMETFilteringProcess())) 

process.dump = cms.EDFilter('MiniAOD2TTreeFilter',
    OutputFileName      = cms.string("miniaod2tree.root"),
    PUInfoInputFileName = process.PUInfo.OutputFileName,
    TopPtInputFileName  = process.TopPtProducer.OutputFileName,
    CodeVersion         = cms.string(git.getCommitId()),
    DataVersion         = cms.string(str(dataVersion.version)),
    CMEnergy            = cms.int32(13),
    Skim = cms.PSet(
	Counters = cms.VInputTag(
	    "skimCounterAll",
            "skimCounterPassed"
        ),
    ),
    EventInfo = cms.PSet(
        PileupSummaryInfoSrc        = process.PUInfo.PileupSummaryInfoSrc, 
	LHESrc                      = cms.untracked.InputTag("externalLHEProducer"),
	OfflinePrimaryVertexSrc     = cms.InputTag("offlineSlimmedPrimaryVertices"),
	TopPtProducer               = cms.InputTag("TopPtProducer"),
    ),
    Trigger = cms.PSet(
	TriggerResults = cms.InputTag("TriggerResults::"+str(dataVersion.getTriggerProcess())),
	TriggerBits    = cms.vstring(
            # SingleMuon Primary Dataset (PD)
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
	L1Extra        = cms.InputTag("l1extraParticles:MET"),
        L1EtSumObjects = cms.InputTag("caloStage2Digis:EtSum"),
	TriggerObjects = cms.InputTag("selectedPatTrigger"),
        TriggerMatch   = cms.untracked.vstring(
            # SingleMuon Primary Dataset (PD)
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
        TriggerPrescales = cms.untracked.PSet(
            src   = cms.InputTag("patTrigger",""),
            paths = cms.vstring(
                )
            ),
	filter = cms.untracked.bool(False)
    ),
    METNoiseFilter = process.METNoiseFilter,
    #Taus      = process.Taus_TauPOGRecommendation,
    Taus      = process.Taus,
    Electrons = process.Electrons,
    Muons     = process.Muons,
    Jets      = process.Jets,
    FatJets   = process.FatJets,
    Top       = process.Top,
    METs      = process.METs,
    GenWeights = cms.VPSet(
        cms.PSet(
            branchname = cms.untracked.string("GenWeights"),
            src        = cms.InputTag("generator"),
            filter     = cms.untracked.bool(False)
        )
    ),
    GenMETs = cms.VPSet(
        cms.PSet(
            branchname = cms.untracked.string("GenMET"),
            src        = cms.InputTag("genMetTrue"),
            filter     = cms.untracked.bool(False)
        )
    ),
    GenJets = cms.VPSet(      
        cms.PSet(
            branchname = cms.untracked.string("GenJets"),
            src        = cms.InputTag("slimmedGenJets"), # ak4
        )
    ),
    GenParticles = cms.VPSet(      
        cms.PSet(
            branchname          = cms.untracked.string("genParticles"),
            src                 = cms.InputTag("prunedGenParticles"),
            saveAllGenParticles = cms.untracked.bool(True),
            saveGenBooleans     = cms.untracked.bool(True),
            saveGenStatusFlags  = cms.untracked.bool(True),
            saveGenElectrons    = cms.untracked.bool(False),
            saveGenMuons        = cms.untracked.bool(False),
            saveGenTaus         = cms.untracked.bool(False),
            saveGenNeutrinos    = cms.untracked.bool(False),
            saveTopInfo         = cms.untracked.bool(False),
            saveWInfo           = cms.untracked.bool(False),
            saveHplusInfo       = cms.untracked.bool(False),
        )
    ),

)

#================================================================================================ 
# Setup skim counters
#================================================================================================ 
process.load("HiggsAnalysis.MiniAOD2TTree.Hplus2hwWithTopAnalysisSkim_cfi")
process.skimCounterAll        = cms.EDProducer("HplusEventCountProducer")
process.skimCounterPassed     = cms.EDProducer("HplusEventCountProducer")
process.skim.TriggerResults = cms.InputTag("TriggerResults::"+str(dataVersion.getTriggerProcess()))

#================================================================================================ 
# Setup customizations
#================================================================================================ 
from HiggsAnalysis.MiniAOD2TTree.CommonFragments import produceCustomisations
produceCustomisations(process, dataVersion.isData()) # This produces process.CustomisationsSequence which needs to be included to path

from HiggsAnalysis.MiniAOD2TTree.CommonFragments import produceAK8Customisations
produceAK8Customisations(process, dataVersion.isData())   # This produces process.AK8CustomisationsSequence which needs to be included to path


# Set up electron MVA ID for Skimming
# https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

print "=== Adding Electron MVA: ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

for idmod in ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff']:
    setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)

#================================================================================================ 
# Module execution
#================================================================================================ 
process.runEDFilter = cms.Path(process.PUInfo*
                               process.TopPtProducer*
                               process.egmGsfElectronIDSequence* # needed?
                               process.skimCounterAll*
                               process.skim*
                               process.skimCounterPassed*
                               process.CustomisationsSequence*
                               process.AK8CustomisationsSequence*
                               process.dump)


#process.output = cms.OutputModule("PoolOutputModule",
#   outputCommands = cms.untracked.vstring(
#       "keep *_*AK8*_*_*",
#       "keep *_*AK4*_*_*",
#       "keep *_selected*_*_*",
#       "keep *_updated*_*_*",
##      "keep *",
#   ),
#   fileName = cms.untracked.string("CMSSW.root")
#)
#process.out_step = cms.EndPath(process.output)
