#!/usr/bin/env python

import os
import sys
import re
import datetime
import subprocess

from optparse import OptionParser
from collections import OrderedDict

from CRABClient.UserUtilities import setConsoleLogLevel
from CRABClient.UserUtilities import getUsernameFromSiteDB
from CRABClient.ClientUtilities import LOGLEVEL_MUTE
from CRABClient.UserUtilities import getConsoleLogLevel

from CRABAPI.RawCommand import crabCommand


class Dataset:
    def __init__(self, url, dbs="global", dataVersion="94Xmc", lumiMask="", name=""):
        self.URL = url
        self.DBS = dbs
        self.dataVersion = dataVersion
        self.lumiMask = lumiMask
        self.name = name

    def isData(self):
        if "data" in self.dataVersion:
            return True
        return False

    def getName(self):
        return self.name

    def getYear(self):
        year_re = re.compile("/Run(?P<year>201\d)\S")
        match = year_re.search(self.URL)
        if match:
            return match.group("year")
        else:
            return None

######################################################
# DATASETS
######################################################
datasets = []
#lumimask="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt"
#lumimask="/afs/cern.ch/work/s/slehti/hplusAnalysis/94x/CMSSW_9_4_9/src/HiggsAnalysis/MiniAOD2TTree/data/Cert_297050-299329_13TeV_EOY2017ReReco_Collisions17_JSON_v1_Run2017B.txt"



#2018
datasets2018 = []
datasets2018.append(Dataset('/Tau/Run2018A-Nano1June2019-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_315257-316995_13TeV_PromptReco_Collisions18_JSON_Run2018A.txt"))
datasets2018.append(Dataset('/DoubleMuon/Run2018B-Nano1June2019-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_317080-319310_13TeV_PromptReco_Collisions18_JSON_Run2018B.txt"))
datasets2018.append(Dataset('/DoubleMuon/Run2018C-Nano1June2019-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_319337-320065_13TeV_PromptReco_Collisions18_JSON_Run2018C.txt"))
datasets2018.append(Dataset('/DoubleMuon/Run2018D-Nano1June2019_ver2-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_320413-325172_13TeV_PromptReco_Collisions18_JSON_Run2018D.txt"))

datasets2018.append(Dataset('/EGamma/Run2018A-Nano1June2019-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_315257-316995_13TeV_PromptReco_Collisions18_JSON_Ru\
n2018A.txt"))
datasets2018.append(Dataset('/EGamma/Run2018B-Nano1June2019-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_317080-319310_13TeV_PromptReco_Collisions18_JSON_Ru\
n2018B.txt"))
datasets2018.append(Dataset('/EGamma/Run2018C-Nano1June2019-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_319337-320065_13TeV_PromptReco_Collisions18_JSON_Ru\
n2018C.txt"))
datasets2018.append(Dataset('/EGamma/Run2018D-Nano1June2019-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_320413-325172_13TeV_PromptReco_Collisions18_JS\
ON_Run2018D.txt"))

datasets2018.append(Dataset('/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM', dataVersion="mc"))
datasets2018.append(Dataset('/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM', dataVersion="mc"))
datasets2018.append(Dataset('/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM', dataVersion="mc"))
datasets2018.append(Dataset('/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM', dataVersion="mc"))

datasets2018.append(Dataset('/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM', dataVersion="mc"))

#2017
datasets2017 = []
datasets2017.append(Dataset('/DoubleMuon/Run2017B-Nano1June2019-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_297050-299329_13TeV_EOY2017ReReco_Collisions17_JSON_v1_Run2017B.txt"))
datasets2017.append(Dataset('/DoubleMuon/Run2017C-Nano1June2019-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_299368-302029_13TeV_EOY2017ReReco_Collisions17_JSON_v1_Run2017C.txt"))
datasets2017.append(Dataset('/DoubleMuon/Run2017D-Nano1June2019-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_302031-302663_13TeV_EOY2017ReReco_Collisions17_JSON_v1_Run2017D.txt"))
datasets2017.append(Dataset('/DoubleMuon/Run2017E-Nano1June2019-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_303824-304797_13TeV_EOY2017ReReco_Collisions17_JSON_v1_Run2017E.txt"))
datasets2017.append(Dataset('/DoubleMuon/Run2017F-Nano1June2019-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_305040-306460_13TeV_EOY2017ReReco_Collisions17_JSON_v1_Run2017F.txt"))

datasets2017.append(Dataset('/DoubleEG/Run2017B-Nano1June2019-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_297050-299329_13TeV_EOY2017ReReco_Collisio\
ns17_JSON_v1_Run2017B.txt"))
datasets2017.append(Dataset('/DoubleEG/Run2017C-Nano1June2019-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_299368-302029_13TeV_EOY2017ReReco_Collisio\
ns17_JSON_v1_Run2017C.txt"))
datasets2017.append(Dataset('/DoubleEG/Run2017D-Nano1June2019-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_302031-302663_13TeV_EOY2017ReReco_Collisio\
ns17_JSON_v1_Run2017D.txt"))
datasets2017.append(Dataset('/DoubleEG/Run2017E-Nano1June2019-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_303824-304797_13TeV_EOY2017ReReco_Collisio\
ns17_JSON_v1_Run2017E.txt"))
datasets2017.append(Dataset('/DoubleEG/Run2017F-Nano1June2019-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_305040-306460_13TeV_EOY2017ReReco_Collisio\
ns17_JSON_v1_Run2017F.txt"))

datasets2017.append(Dataset('/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv5-PU2017_12Apr2018_Nano1June2019_v3_102X_mc2017_realistic_v7_ext1-v1/NANOAODSIM', dataVersion="mc"))
datasets2017.append(Dataset('/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv5-PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/NANOAODSIM', dataVersion="mc"))
datasets2017.append(Dataset('/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv5-PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7_ext1-v1/NANOAODSIM', dataVersion="mc"))
datasets2017.append(Dataset('/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv5-PU2017_12Apr2018_Nano1June2019_v2_102X_mc2017_realistic_v7-v1/NANOAODSIM', dataVersion="mc"))

datasets2017.append(Dataset('/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv5-PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/NANOAODSIM', dataVersion="mc"))

# 2016
datasets2016 = []
datasets2016.append(Dataset('/Tau/Run2016B_ver2-Nano1June2019_ver2-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_273301-275376_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016B.txt"))

#datasets2016.append(Dataset('/DY1JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv5-PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM', dataVersion="mc"))
datasets2016.append(Dataset('/DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv5-PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM', dataVersion="mc"))
#datasets2016.append(Dataset('/DY2JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv5-PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM', dataVersion="mc"))
#datasets2016.append(Dataset('/DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv5-PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM', dataVersion="mc"))
#datasets2016.append(Dataset('/DY3JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv5-PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM', dataVersion="mc"))
#datasets2016.append(Dataset('/DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv5-PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM', dataVersion="mc"))
#datasets2016.append(Dataset('/DY4JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv5-PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM', dataVersion="mc"))
#datasets2016.append(Dataset('/DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv5-PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM', dataVersion="mc"))

#datasets2016.append(Dataset('/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv5-PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7-v1/NANOAODSIM', dataVersion="mc"))

#datasets2016.append(Dataset('/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv5-PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7_ext1-v1/NANOAODSIM', dataVersion="mc"))
#datasets2016.append(Dataset('/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv5-PUMoriond17_Nano1June2019_102X_mcRun2_asymptotic_v7_ext2-v1/NANOAODSIM', dataVersion="mc"))

datasets = datasets2016
#datasets = datasets2017
#datasets = datasets2018

alldatasets = []
alldatasets.extend(datasets2016)
alldatasets.extend(datasets2017)
alldatasets.extend(datasets2018)


"""
datasets.append(Dataset('/DoubleMuon/Run2016G-Nano14Dec2018-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_278820-280385_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016G.txt"))
datasets.append(Dataset('/DoubleMuon/Run2016H-Nano14Dec2018-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_281613-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016H.txt"))

datasets.append(Dataset('/DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/NANOAODSIM', dataVersion="mc"))
datasets.append(Dataset('/DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/NANOAODSIM', dataVersion="mc"))
datasets.append(Dataset('/DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/NANOAODSIM', dataVersion="mc"))
datasets.append(Dataset('/DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/NANOAODSIM', dataVersion="mc"))

datasets.append(Dataset('/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-herwigpp_30M/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/NANOAODSIM', dataVersion="mc", name="DYJetsToLL_M_50_herwig"))

datasets.append(Dataset('/DYBJetsToLL_M-50_Zpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/NANOAODSIM', dataVersion="mc", name="DYBJetsToLL_M_50_ZPT_100to200"))
datasets.append(Dataset('/DYBJetsToLL_M-50_Zpt-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/NANOAODSIM', dataVersion="mc", name="DYBJetsToLL_M_50_ZPT_100to200_ext"))
datasets.append(Dataset('/DYBJetsToLL_M-50_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/NANOAODSIM', dataVersion="mc", name="DYBJetsToLL_M_50_ZPT_200toInf"))
datasets.append(Dataset('/DYBJetsToLL_M-50_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext1-v1/NANOAODSIM', dataVersion="mc", name="DYBJetsToLL_M_50_ZPT_200toInf_ext"))



datasets.append(Dataset('/DoubleMuon/Run2017B-Nano14Dec2018-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_297050-299329_13TeV_EOY2017ReReco_Collisions17_JSON_v1_Run2017B.txt"))
datasets.append(Dataset('/DoubleMuon/Run2017C-Nano14Dec2018-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_299368-302029_13TeV_EOY2017ReReco_Collisions17_JSON_v1_Run2017C.txt"))
datasets.append(Dataset('/DoubleMuon/Run2017D-Nano14Dec2018-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_302031-302663_13TeV_EOY2017ReReco_Collisions17_JSON_v1_Run2017D.txt"))
datasets.append(Dataset('/DoubleMuon/Run2017E-Nano14Dec2018-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_303824-304797_13TeV_EOY2017ReReco_Collisions17_JSON_v1_Run2017E.txt"))
datasets.append(Dataset('/DoubleMuon/Run2017F-Nano14Dec2018-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_305040-306460_13TeV_EOY2017ReReco_Collisions17_JSON_v1_Run2017F.txt"))

#datasets.append(Dataset('/DoubleEG/Run2017B-Nano14Dec2018-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_297050-299329_13TeV_EOY2017ReReco_Collisions17_JSON_v1_Run2017B.txt"))
#datasets.append(Dataset('/DoubleEG/Run2017C-Nano14Dec2018-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_299368-302029_13TeV_EOY2017ReReco_Collisions17_JSON_v1_Run2017C.txt"))
#datasets.append(Dataset('/DoubleEG/Run2017D-Nano14Dec2018-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_302031-302663_13TeV_EOY2017ReReco_Collisions17_JSON_v1_Run2017D.txt"))
#datasets.append(Dataset('/DoubleEG/Run2017E-Nano14Dec2018-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_303824-304797_13TeV_EOY2017ReReco_Collisions17_JSON_v1_Run2017E.txt"))
#datasets.append(Dataset('/DoubleEG/Run2017F-Nano14Dec2018-v1/NANOAOD', dataVersion="data",lumiMask="/afs/cern.ch/work/s/slehti/JME/skim/CMSSW_10_2_14/src/PhysicsTools/NanoAODTools/crab/Cert_305040-306460_13TeV_EOY2017ReReco_Collisions17_JSON_v1_Run2017F.txt"))

datasets.append(Dataset('/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_v3_102X_mc2017_realistic_v6_ext1-v1/NANOAODSIM', dataVersion="mc"))
datasets.append(Dataset('/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM', dataVersion="mc"))
datasets.append(Dataset('/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6_ext1-v1/NANOAODSIM', dataVersion="mc"))
datasets.append(Dataset('/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM', dataVersion="mc"))
datasets.append(Dataset('/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6_ext1-v1/NANOAODSIM', dataVersion="mc"))
datasets.append(Dataset('/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_v2_102X_mc2017_realistic_v6-v1/NANOAODSIM', dataVersion="mc"))


#datasets.append(Dataset('/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017RECOSIMstep_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv4-PU2017RECOSIMstep_12Apr2018_Nano14Dec2018_102X_mc2017_realistic_v6_ext1-v1/NANOAODSIM', dataVersion="mc"))

####datasets.append(Dataset('/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-herwigpp_30M/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/NANOAODSIM', dataVersion="mc", name="DYJetsToLL_M_50_herwig"))
#/DYJetsToLL_M-50_TuneCUETHS1_13TeV-madgraphMLM-herwigpp/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/NANOAODSIM

datasets.append(Dataset('/DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/NANOAODSIM', dataVersion="mc"))
datasets.append(Dataset('/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/NANOAODSIM', dataVersion="mc"))
datasets.append(Dataset('/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/NANOAODSIM', dataVersion="mc"))
datasets.append(Dataset('/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/NANOAODSIM', dataVersion="mc"))
datasets.append(Dataset('/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/NANOAODSIM', dataVersion="mc"))
datasets.append(Dataset('/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/NANOAODSIM', dataVersion="mc"))
datasets.append(Dataset('/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/NANOAODSIM', dataVersion="mc"))
datasets.append(Dataset('/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAODv4-PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6-v1/NANOAODSIM', dataVersion="mc"))


#datasets.append(Dataset('/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/NANOAODSIM', dataVersion="mc"))

#datasets.append(Dataset('/DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))


#datasets.append(Dataset('/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#
#datasets.append(Dataset('/ST_s-channel_4f_InclusiveDecays_13TeV-amcatnlo-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM', dataVersion="mc"))

#datasets.append(Dataset('/WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v2/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v2/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/NANOAODSIM', dataVersion="mc"))
#datasets.append(Dataset('/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16NanoAOD-PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2_ext1-v1/NANOAODSIM', dataVersion="mc"))
"""

######################################################

#================================================================================================
# Class Definition
#================================================================================================
class colors:
    '''
    \033[  Escape code, this is always the same
    1 = Style, 1 for normal.
    32 = Text colour, 32 for bright green.
    40m = Background colour, 40 is for black.
        
    WARNING:
    Python doesn't distinguish between 'normal' characters and ANSI colour codes, which are also characters that the terminal interprets.
    In other words, printing '\x1b[92m' to a terminal may change the terminal text colour, Python doesn't see that as anything but a set of 5 characters.
    If you use print repr(line) instead, python will print the string literal form instead, including using escape codes for non-ASCII printable characters
    (so the ESC ASCII code, 27, is displayed as \x1b) to see how many have been added.
        
    You'll need to adjust your column alignments manually to allow for those extra characters.
    Without your actual code, that's hard for us to help you with though.
    
    Useful Links:
    http://ozzmaker.com/add-colour-to-text-in-python/
    http://stackoverflow.com/questions/15580303/python-output-complex-line-with-floats-colored-by-value
    '''
    colordict = { 
                'RED'     :'\033[91m',
                'GREEN'   :'\033[92m',
                'BLUE'    :'\033[34m',
                'GRAY'    :'\033[90m',
                'WHITE'   :'\033[00m',
                'ORANGE'  :'\033[33m',
                'CYAN'    :'\033[36m',
                'PURPLE'  :'\033[35m',
                'LIGHTRED':'\033[91m',
                'PINK'    :'\033[95m',
                'YELLOW'  :'\033[93m',
                }
    if sys.stdout.isatty(): 
        RED      = colordict['RED']
        GREEN    = colordict['GREEN']
        BLUE     = colordict['BLUE']
        GRAY     = colordict['GRAY']
        WHITE    = colordict['WHITE']
        ORANGE   = colordict['ORANGE']
        CYAN     = colordict['CYAN']
        PURPLE   = colordict['PURPLE']
        LIGHTRED = colordict['LIGHTRED']
        PINK     = colordict['PINK']
        YELLOW   = colordict['YELLOW']
    else:
        RED, GREEN, BLUE, GRAY, WHITE, ORANGE, CYAN, PURPLE, LIGHTRED, PINK, YELLOW = '', '', '', '', '', '', '', '', '', '', ''

#================================================================================================
# Class Definition
#================================================================================================
class Report:
    def __init__(self, name, allJobs, idle, retrieved, running, finished, failed, transferring, retrievedLog, retrievedOut, eosLog, eosOut, status, dashboardURL):
        '''
        Constructor
        '''
        self.name            = name   
        self.allJobs         = str(allJobs)
        self.retrieved       = str(retrieved)
        self.running         = str(running)
        self.dataset         = self.name.split("/")[-1]
        self.dashboardURL    = dashboardURL
        self.status          = self.GetTaskStatusStyle(status)
        self.finished        = str(len(finished))
        self.failed          = failed
        self.idle            = idle 
        self.transferring    = transferring
        self.retrievedLog    = str(len(retrievedLog))
        self.retrievedOut    = str(len(retrievedOut))
        self.eosLog          = eosLog 
        self.eosOut          = eosOut   
        return

    def Print(self, printHeader=True):
        '''
        Simple function to print report.
        '''
        name = os.path.basename(self.name)
        while len(name) < 30:
            name += " "
    
        fName = GetSelfName()
        cName = self.__class__.__name__
        name  = fName + ": " + cName
        if printHeader:
            print "=== ", name
        msg  = '{:<20} {:<40}'.format("\t %sDataset"           % (colors.WHITE) , ": " + self.dataset)
        msg += '\n {:<20} {:<40}'.format("\t %sRetrieved Jobs" % (colors.WHITE) , ": " + self.retrieved + " / " + self.allJobs)
        msg += '\n {:<20} {:<40}'.format("\t %sStatus"         % (colors.WHITE) , ": " + self.status)
        msg += '\n {:<20} {:<40}'.format("\t %sDashboard"      % (colors.WHITE) , ": " + self.dashboardURL)
        print msg
        return
        
        
    def GetURL():
        return self.dashboardURL

    def GetTaskStatusStyle(self, status):
        '''   
        NEW, RESUBMIT, KILL: Temporary statuses to indicate the action ('submit', 'resubmit' or 'kill') that has to be applied to the task.
        QUEUED: An action ('submit', 'resubmit' or 'kill') affecting the task is queued in the CRAB3 system.
        SUBMITTED: The task was submitted to HTCondor as a DAG task. The DAG task is currently running.
        SUBMITFAILED: The 'submit' action has failed (CRAB3 was unable to create a DAG task).
        FAILED: The DAG task completed all nodes and at least one is a permanent failure.
        COMPLETED: All nodes have been completed
        KILLED: The user killed the task. 
        KILLFAILED: The 'kill' action has failed.
        RESUBMITFAILED: The 'resubmit' action has failed.
        '''
        # Remove all whitespace characters (space, tab, newline, etc.)
        status = ''.join(status.split())
        if status == "NEW":
            status = "%s%s%s" % (colors.BLUE, status, colors.WHITE)
        elif status == "RESUBMIT":
            status = "%s%s%s" % (colors.BLUE, status, colors.WHITE)
        elif status == "QUEUED":
            status = "%s%s%s" % (colors.GRAY, status, colors.WHITE)
        elif status == "SUBMITTED":
            status = "%s%s%s" % (colors.BLUE, status, colors.WHITE)
        elif status == "SUBMITFAILED":
            status = "%s%s%s" % (colors.RED, status, colors.WHITE)
        elif status == "FAILED":
            status = "%s%s%s" % (colors.RED, status, colors.WHITE)
        elif status == "COMPLETED":
            status = "%s%s%s" % (colors.GREEN, status, colors.WHITE)
        elif status == "KILLED":
            status = "%s%s%s" % (colors.ORANGE, status, colors.WHITE)
        elif status == "KILLFAILED":
            status = "%s%s%s" % (colors.ORANGE, status, colors.WHITE)
        elif status == "RESUBMITFAILED":
            status = "%s%s%s" % (colors.ORANGE, status, colors.WHITE)
        elif status == "?":
            status = "%s%s%s" % (colors.PINK, status, colors.WHITE)
        elif status == "UNDETERMINED":
            status = "%s%s%s" % (colors.CYAN, status, colors.WHITE)
        elif status == "UNKNOWN":
            status = "%s%s%s" % (colors.LIGHTRED, status, colors.WHITE)
        else:
            raise Exception("Unexpected task status %s." % (status) )
        return status

def EnsurePathDoesNotExist(taskDirName, requestName):
    '''                                                                                                                                     
    Ensures that file does not already exist                                                                                                
    '''
    filePath = os.path.join(taskDirName, requestName)

    if not os.path.exists(filePath):
        return
    else:
        msg = "File '%s' already exists!" % (filePath)
        Print(msg + "\n\tProceeding to overwrite file.")
    return

def GetCMSSW():
    '''                                                                                                                                     
    Get a command-line-friendly format of the CMSSW version currently use.                                                                  
    https://docs.python.org/2/howto/regex.html                                                                                              
    '''

    # Get the current working directory                                                                                                     
    pwd = os.getcwd()

    # Create a compiled regular expression object                                                                                           
    cmssw_re = re.compile("/CMSSW_(?P<version>\S+?)/")

    # Scan through the string 'pwd' & look for any location where the compiled RE 'cmssw_re' matches                                        
    match = cmssw_re.search(pwd)

    # Return the string matched by the RE. Convert to desirable format                                                                      
    version = ""
    if match:
        version = match.group("version")
        version = version.replace("_","")
        version = version.replace("pre","p")
        version = version.replace("patch","p")
    return version

def CreateTaskDir(pset,dirName):
    '''                                                                                                                                     
    Create the CRAB task directory and copy inside it the PSET to be used for the CRAB job.                                                 
    '''
    # Copy file to be used (and others to be tracked) to the task directory                                                                 
    cmd = "cp %s %s" %(pset, dirName)

    if not os.path.exists(dirName):
        os.mkdir(dirName)

        os.system(cmd)
    else:
        pass

####    # Write the commit id, "git status", "git diff" command output the directory created for the multicrab task                             
####    gitFileList = git.writeCodeGitInfo(dirName, False)
    return

def GetRequestName(dataset):
    '''                                                                                                                                     
    Return the file name and path to an (empty) crabConfig_*.py file where "*"                                                              
    contains the dataset name and other information such as tune, COM, Run number etc..                                                     
    of the Data or MC sample used                                                                                                           
    '''
    #print "check GetRequestName",dataset.getName()
    if len(dataset.getName()) > 0:
        return dataset.getName()

    # Create compiled regular expression objects                                                                                            
    datadataset_re = re.compile("^/(?P<name>\S+?)/(?P<run>Run\S+?)/")
    mcdataset_re   = re.compile("^/(?P<name>\S+?)/")
    tune_re        = re.compile("(?P<name>\S+)_Tune")
    tev_re         = re.compile("(?P<name>\S+)_13TeV")
    ext_re         = re.compile("(?P<name>_ext\d+)-")
    runRange_re    = re.compile("Cert_(?P<RunRange>\d+-\d+)_")
    # runRange_re    = re.compile("Cert_(?P<RunRange>\d+-\d+)_13TeV_PromptReco_Collisions15(?P<BunchSpacing>\S*)_JSON(?P<Silver>(_\S+|))\."\)                                                                                                                                           
    # runRange_re    = re.compile("Cert_(?P<RunRange>\d+-\d+)_13TeV_PromptReco_Collisions15(?P<BunchSpacing>\S*)_JSON")                     
    # runRange_re    = re.compile("Cert_(?P<RunRange>\d+-\d+)_13TeV_PromptReco_Collisions15_(?P<BunchSpacing>\d+ns)_JSON_v")                

    # Scan through the string 'dataset.URL' & look for any location where the compiled RE 'mcdataset_re' matches                            
    match = mcdataset_re.search(dataset.URL)
    if dataset.isData():
        match = datadataset_re.search(dataset.URL)

    # Append the dataset name                                                                                                               
    if match:
        requestName = match.group("name")
    firstName = requestName

    # Append the Run number (for Data samples only)                                                                                         
    if dataset.isData():
        requestName+= "_"
        requestName+= match.group("run")

    # Append the MC-tune (for MC samples only)                                                                                              
    tune_match = tune_re.search(requestName)
    if tune_match:
        requestName = tune_match.group("name")

    # Append the COM Energy (for MC samples only)                                                                                           
    tev_match = tev_re.search(requestName)
    if tev_match:
        requestName = tev_match.group("name")

    # Simple hack to prevent overwrite of special TT samples                                                                                
####    requestName = GetTTbarSystematicsName(firstName, requestName)

    # Append the Ext                                                                                                                        
    ext_match = ext_re.search(dataset.URL)
    if ext_match:
        requestName+=ext_match.group("name")

    # Append the Run Range (for Data samples only)                                                                                          
    if dataset.isData():
        runRangeMatch = runRange_re.search(dataset.lumiMask)
        if runRangeMatch:
            runRange= runRangeMatch.group("RunRange")
            runRange = runRange.replace("-","_")
            #bunchSpace = runRangeMatch.group("BunchSpacing")                                                                               
            requestName += "_" + runRange #+ bunchSpace                                                                                     
            #Ag = runRangeMatch.group("Silver")                                                                                             
            #if Ag == "_Silver": # Use  chemical element of silver (Ag)                                                                     
            #    requestName += Ag

    # Finally, replace dashes with underscores                                                                                              
    requestName = requestName.replace("-","_")
    return requestName

def CreateCfgFile(dataset, taskDirName, requestName, infilePath, opts):
    '''                                                                                                                                     
    Creates a CRAB-specific configuration file which will be used in the submission                                                         
    of a job. The function uses as input a generic cfg file which is then customised                                                        
    based on the dataset type used.                                                                                                         
                                                                                                                                            
    infilePath = "crabConfig.py"                                                                                                            
    '''

    outfilePath = os.path.join(taskDirName, "crabConfig_" + requestName + ".py")

    # Check that file does not already exist                                                                                                
    EnsurePathDoesNotExist(taskDirName, outfilePath)

    # Open input file (read mode) and output file (write mode)                                                                              
    fIN  = open(infilePath , "r")
    fOUT = open(outfilePath, "w")

    # Create compiled regular expression objects                                                                                            
    crab_requestName_re     = re.compile("config.General.requestName")
    crab_workArea_re        = re.compile("config.General.workArea")
    crab_transferOutputs_re = re.compile("config.General.transferOutputs")
    crab_transferLogs_re    = re.compile("config.General.transferLogs")
    crab_pset_re            = re.compile("config.JobType.psetName")
    crab_psetParams_re      = re.compile("config.JobType.pyCfgParams")
    crab_dataset_re         = re.compile("config.Data.inputDataset")
    crab_split_re           = re.compile("config.Data.splitting")
    crab_splitunits_re      = re.compile("config.Data.unitsPerJob")
    crab_dbs_re             = re.compile("config.Data.inputDBS")
    crab_storageSite_re     = re.compile("config.Site.storageSite")
    crab_outLFNDirBase_re   = re.compile("config.Data.outLFNDirBase")

    # For-loop: All line of input fine                                                                                                      
    for line in fIN:

        # Skip lines which are commented out                                                                                                
        if line[0] == "#":
            continue

        # Set the "inputDataset" field which specifies the name of the dataset. Can be official CMS dataset or a dataset produced by a user.                                                                                                                                           
        match = crab_dataset_re.search(line)
        if match:
            line = "config.Data.inputDataset = '" + dataset.URL + "'\n"

        # Set the "requestName" field which specifies the request/task name. Used by CRAB to create a project directory (named crab_<requestName>)                                                                                                                                     
        match = crab_requestName_re.search(line)
        if match:
            line = "config.General.requestName = '" + requestName + "'\n"

        # Set the "workArea" field which specifies the (full or relative path) where to create the CRAB project directory.                  
        match = crab_workArea_re.search(line)
        if match:
            line = "config.General.workArea = '" + taskDirName + "'\n"

        # Set the "psetName" field which specifies the name of the CMSSW pset_cfg.py file that will be run via cmsRun.                      
        match = crab_pset_re.search(line)
        if match:
            line = "config.JobType.psetName = '" + opts.pset +"'\n"

def CreateCfgFile(dataset, taskDirName, requestName, infilePath, opts):
    '''                                                                                                                                     
    Creates a CRAB-specific configuration file which will be used in the submission                                                         
    of a job. The function uses as input a generic cfg file which is then customised                                                        
    based on the dataset type used.                                                                                                         
                                                                                                                                            
    infilePath = "crabConfig.py"                                                                                                            
    '''

    outfilePath = os.path.join(taskDirName, "crabConfig_" + requestName + ".py")

    # Check that file does not already exist                                                                                                
    EnsurePathDoesNotExist(taskDirName, outfilePath)

    # Open input file (read mode) and output file (write mode)                                                                              
    fIN  = open(infilePath , "r")
    fOUT = open(outfilePath, "w")

    # Create compiled regular expression objects                                                                                            
    crab_requestName_re     = re.compile("config.General.requestName")
    crab_workArea_re        = re.compile("config.General.workArea")
    crab_transferOutputs_re = re.compile("config.General.transferOutputs")
    crab_transferLogs_re    = re.compile("config.General.transferLogs")
    crab_pset_re            = re.compile("config.JobType.psetName")
    crab_psetParams_re      = re.compile("config.JobType.pyCfgParams")
    crab_dataset_re         = re.compile("config.Data.inputDataset")
    crab_split_re           = re.compile("config.Data.splitting")
    crab_splitunits_re      = re.compile("config.Data.unitsPerJob")
    crab_dbs_re             = re.compile("config.Data.inputDBS")
    crab_storageSite_re     = re.compile("config.Site.storageSite")
    crab_outLFNDirBase_re   = re.compile("config.Data.outLFNDirBase")

    # For-loop: All line of input fine                                                                                                      
    for line in fIN:

        # Skip lines which are commented out                                                                                                
        if line[0] == "#":
            continue

        # Set the "inputDataset" field which specifies the name of the dataset. Can be official CMS dataset or a dataset produced by a user.                                                                                                                                           
        match = crab_dataset_re.search(line)
        if match:
            line = "config.Data.inputDataset = '" + dataset.URL + "'\n"

        # Set the "requestName" field which specifies the request/task name. Used by CRAB to create a project directory (named crab_<requestName>)                                                                                                                                     
        match = crab_requestName_re.search(line)
        if match:
            line = "config.General.requestName = '" + requestName + "'\n"

        # Set the "workArea" field which specifies the (full or relative path) where to create the CRAB project directory.                  
        match = crab_workArea_re.search(line)
        if match:
            line = "config.General.workArea = '" + taskDirName + "'\n"

        # Set the "psetName" field which specifies the name of the CMSSW pset_cfg.py file that will be run via cmsRun.                      
        match = crab_pset_re.search(line)
        if match:
            line = "config.JobType.psetName = '" + opts.pset +"'\n"

        # Set the "pyCfgParams" field which contains list of parameters to pass to the pset_cfg.py file.                                    
        match = crab_psetParams_re.search(line)
        if match:
            line = "config.JobType.pyCfgParams = ['dataVersion=" + dataset.dataVersion +"']\n"

        # Set the "inputDBS" field which specifies the URL of the DBS reader instance where the input dataset is published                  
        match = crab_dbs_re.search(line)
        if match:
            line = "config.Data.inputDBS = '" + dataset.DBS + "'\n"

        # Set the "storageSite" field which specifies the destination site for submission [User MUST have write access to destination site!]                                                                                                                                           
        match = crab_storageSite_re.search(line)
        if match:
            line = "config.Site.storageSite = '" + opts.storageSite + "'\n"

        match = crab_outLFNDirBase_re.search(line)
        if match:
            if opts.dirName.endswith("/"):
                mcrabDir = os.path.basename(opts.dirName[:-1])  # exclude last "/", either-wise fails                                       
            else:
                mcrabDir = os.path.basename(opts.dirName)
            fullDir  = "/store/user/%s/CRAB3_TransferData/%s" % (getUsernameFromSiteDB(), mcrabDir) # NOT getpass.getuser()                 
            line     = "config.Data.outLFNDirBase = '" + fullDir + "'\n"

        # Only if dataset is real data                                                                                                      
        if dataset.isData():

            # Set the "splitting" field which specifies the mode to use to split the task in jobs ('FileBased', 'LumiBased', or 'EventAwareLumiBased')                                                                                                                                 
            match = crab_split_re.search(line)
####            if match:
####                line = "config.Data.splitting = 'LumiBased'\n"
####                line+= "config.Data.lumiMask = '"+ dataset.lumiMask + "'\n"

            # Set the "unitsPerJob" field which suggests (but not impose) how many files, lumi sections or events to include in each job.   
            match = crab_splitunits_re.search(line)
            if match:
####                line = "config.Data.unitsPerJob = 25\n"
                line = "config.Data.unitsPerJob = 1\n"
        else:
            pass

        # Write line to the output file                                                                                                     
        fOUT.write(line)

    # Close input and output files                                                                                                          
    fOUT.close()
    fIN.close()

    return

def SubmitTaskDir(taskDirName, requestName):
    '''                                                                                                                                     
    Submit a given CRAB task using the specific cfg file.                                                                                   
    '''

    outfilePath = os.path.join(taskDirName, "crabConfig_" + requestName + ".py")

    # Submit the CRAB task                                                                                                                  
    cmd_submit = "crab submit " + outfilePath
    os.system(cmd_submit)

    # Rename the CRAB task directory (remove "crab_" from its name)                                                                         
    cmd_mv = "mv " + os.path.join(taskDirName, "crab_" + requestName) + " " + os.path.join(taskDirName, requestName)
    #print "check cmd_mv",cmd_mv
    os.system(cmd_mv)
    return

def GetTaskDirName(analysis, version, datasets, opts):
    '''                                                                                                                                     
    Get the name of the CRAB task directory to be created. For the user's benefit this                                                      
    will include the CMSSW version and possibly important information from                                                                  
    the dataset used, such as the bunch-crossing time.                                                                                      
    '''

    # Constuct basic task directory name                                                                                                    
    dirName = "multicrab"
    dirName+= "_"  + analysis
    dirName+= "_v" + version

    # Add dataset-specific info, like bunch-crossing info                                                                                   
    bx_re = re.compile("\S+(?P<bx>\d\dns)_\S+")
    match = bx_re.search(datasets[0].URL)
    if match:
        dirName+= "_"+match.group("bx")

    run_re = re.compile("(?P<run>Run201\d)(?P<letter>\S)")
    runs = ""
    runletters = []
    for d in datasets:
        if not d.isData():
            continue
        match = run_re.search(d.URL)
        if match and len(runs) == 0:
            runs = match.group("run")
        if match and match.group("letter") not in runletters:
            runletters.append(match.group("letter"))
    for l in runletters:
        runs+=l
    dirName+= "_"+runs

    # Append the creation time to the task directory name                                                                                   
    # time = datetime.datetime.now().strftime("%d%b%Y_%Hh%Mm%Ss")                                                                           
    time = datetime.datetime.now().strftime("%Y%m%dT%H%M")
    dirName+= "_" + time

    # If directory already exists (resubmission)                                                                                            
    if os.path.exists(opts.dirName) and os.path.isdir(opts.dirName):
        dirName = opts.dirName

    return dirName

def CreateJob(opts, args):
    '''                                                                                                                                     
    Create & submit a CRAB task, using the user-defined PSET and list of datasets.                                                          
    '''

    # Get general info                                                                                                                      
    version      = GetCMSSW()
    analysis     = opts.name.replace('.py','')
    #datasets     = datasets #DatasetGroup(analysis).GetDatasetList()
    taskDirName  = GetTaskDirName(analysis, version, datasets, opts)
    opts.dirName = taskDirName

    # Create CRAB task diractory                                                                                                            
    CreateTaskDir(opts.pset,taskDirName)

    # For-loop: All datasets                                                                                                                
    for dataset in datasets:
        requestName = GetRequestName(dataset)
        fullDir     = taskDirName + "/" + requestName

        if os.path.exists(fullDir) and os.path.isdir(fullDir):
            continue

        CreateCfgFile(dataset, taskDirName, requestName, "crab_cfg.py", opts)
        SubmitTaskDir(taskDirName, requestName)

    print taskDirName
    return 0

def CheckJob(opts, args):
    '''
    Check status, retrieve, resubmit, kill CRAB tasks.
    '''
    
    # Force crabCommand to stay quite
    setConsoleLogLevel(LOGLEVEL_MUTE)
    
    # Retrieve the current crabCommand console log level:
    crabConsoleLogLevel = getConsoleLogLevel()
                      
    # Get the paths for the datasets (absolute paths)
    datasets = GetDatasetsPaths(opts)
    if len(datasets) < 1:
        Print("Found %s CRAB tasks under %s! Exit .." % (opts.dirName) )
        return
        
    # Create a dictionary to map TaskName <-> CRAB Report
    reportDict = GetCrabReportDictionary(datasets,opts)
    
    # Print a summary table with information on each CRAB Task
    PrintTaskSummary(reportDict)
    return

def GetDatasetsPaths(opts):
    '''
    Return the absolute path for each task/dataset located inside
    the working multi-CRAB directory. If the --inlcudeTask or the
    --excludeTask options are used, they are taken into consideration
    accordingly.
    '''

    # Get the multi-CRAB working dir
    multicrabDirPath = [opts.dirName]
    
    # Get the absolute path for each task(=dataset)
    datasetsDirPaths = GetDatasetAbsolutePaths(multicrabDirPath)

    # Check include/exclude options to get final datasets list
    datasets = GetIncludeExcludeDatasets(datasetsDirPaths, opts)
    
    return datasets

def GetDatasetAbsolutePaths(datasetdirs):

    datasets = []
    # For-loop: All CRAB dirs (absolute paths)
    for d in datasetdirs:

        if os.path.exists( os.path.join(d, "results") ):
            datasets.append(d)
    
        # Get the contents of this directory
        cands = Execute("ls -tr %s"%d)

        # For-loop: All directory contents
        for c in cands:
            path = os.path.join(d, c)
            # Get all dataset directories
            if os.path.exists( os.path.join(path, "results") ):
                datasets.append(path)
    return datasets

def Execute(cmd):
    '''
    Executes a given command and return the output.
    '''
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    
    stdin  = p.stdout
    stdout = p.stdout
    ret    = []
    for line in stdout:
        ret.append(line.replace("\n", ""))
    
    stdout.close()
    return ret

def GetIncludeExcludeDatasets(datasets, opts):
    '''
    Does nothing by default, unless the user specifies a dataset to include (--includeTasks <datasetNames>) or
    to exclude (--excludeTasks <datasetNames>) when executing the script. This function filters for the inlcude/exclude
    datasets and returns the lists of datasets and baseNames to be used further in the program.
    '''

    # Initialise lists
    newDatasets  = []
    
    # Exclude datasets
    if opts.excludeTasks != "None":
        tmp = []
        exclude = GetRegularExpression(opts.excludeTasks)
    
        for d in datasets:
            task  = GetBasename(d)
            found = False 
    
            for e_re in exclude:
                if e_re.search(task):
                    found = True
                    break
            if found:
                continue
            newDatasets.append(d)
        return newDatasets

    # Include datasets
    if opts.includeTasks != "None":
        tmp = []
        include = GetRegularExpression(opts.includeTasks)
    
        for d in datasets:
            task  = GetBasename(d)
            found = False

            for i_re in include:
                if i_re.search(task):
                    found = True
                    break
            if found:
                newDatasets.append(d)
        return newDatasets

    return datasets

def GetCrabReportDictionary(datasets,opts):
    '''
    Loops over all datasets paths.
    Retrieves the report object for the given task
    and saves it into a dictionary, thus mapping the
    task name (basename of dataset path) to the CRAB
    report for that task.
    '''

    reportDict = {}
    # For-loop: All (absolute) paths of the datasets
    for index, d in enumerate(datasets):

        # Check if task is in "DONE" state
        if GetTaskStatusBool(d):
            continue  
    
        # Get the CRAB task report & add to dictionary (retrieves job output!)
        report = GetTaskReports(d, opts)
        reportDict[d.split("/")[-1]] = report
    
    return reportDict

def PrintTaskSummary(reportDict):
    '''                                                                                                                                                    
    Print a summary table of all submitted tasks with their information.                                                                                   
    The purpose it to easily determine which jobs are done, running and failed.                                                                            
    '''
    Verbose("PrintTaskSummary()")

    reports  = []
    msgAlign = "{:<3} {:<50} {:^16} {:^16} {:^16} {:^16} {:^16} {:^16} {:^16} {:^16} {:^16} {:^16}"
    header   = msgAlign.format("#", "Task",
                               "%s%s" % (colors.GRAY  , "Idle"    ),
                               "%s%s" % (colors.RED   , "Failed"  ),
                               "%s%s" % (colors.ORANGE, "Running" ),
                               "%s%s" % (colors.ORANGE, "Transfer"),
                               "%s%s" % (colors.WHITE , "Done"    ),
                               "%s%s" % (colors.PURPLE, "Logs"    ),
                               "%s%s" % (colors.BLUE  , "Out"     ),
                               "%s%s" % (colors.CYAN  , "Logs"    ),
                               "%s%s" % (colors.CYAN  , "Out"     ),
                               "%s%s" % (colors.WHITE , "Status"  ),
                               )
    hLine = colors.WHITE + "="*175
    reports.append(hLine)
    reports.append(header)
    reports.append(hLine)

    # Alphabetical sorting of tasks                                                                                                                        
    ReportDict = OrderedDict(sorted(reportDict.items(), key=lambda t: t[0]))
    # For-loop: All datasets (key) and corresponding status (value)                                                                                        
    for i, dataset in enumerate(ReportDict):
        report     = reportDict[dataset]
        index      = i+1
        task       = dataset
        status     = report.status
        idle       = '{0: >3}'.format(report.idle)
        allJobs    = '{0: <3}'.format(report.allJobs)
        running    = '{0: >3}'.format(report.running)
        finished   = '{0: >3}'.format(report.finished)
        transfer   = '{0: >3}'.format(report.transferring)
        failed     = '{0: >3}'.format(len(report.failed))
        rLogs      = '{0: >3}'.format(report.retrievedLog)
        rOutput    = '{0: >3}'.format(report.retrievedOut)
        rLogsEOS   = '{0: >3}'.format(report.eosLog)
        rOutputEOS = '{0: >3}'.format(report.eosOut)
        line = msgAlign.format(index, task,
                               "%s%s/%s" % (colors.GRAY  , idle      , allJobs),
                               "%s%s/%s" % (colors.RED   , failed    , allJobs),
                               "%s%s/%s" % (colors.ORANGE, running   , allJobs),
                               "%s%s/%s" % (colors.ORANGE, transfer  , allJobs),
                               "%s%s/%s" % (colors.WHITE , finished  , allJobs),
                               "%s%s/%s" % (colors.PURPLE, rLogs     , allJobs),
                               "%s%s/%s" % (colors.BLUE  , rOutput   , allJobs),
                               "%s%s/%s" % (colors.CYAN  , rLogsEOS  , allJobs),
                               "%s%s/%s" % (colors.CYAN  , rOutputEOS, allJobs),
                               "%s"   % (status), #already with colour                                                                                     
                               )
        reports.append(line)
    reports.append(hLine)

    # For-loop: All lines in report table                                                                                                                  
    print
    for r in reports:
        print r
    print
    return

def GetTaskStatusBool(datasetPath):
    '''
    Check the crab.log for the given task to determine the status.
    If the the string "Done" is found inside skip it.
    '''
    crabLog      = os.path.join(datasetPath,"crab.log")   
    stringToGrep = "Done"
    cmd          = "grep '%s' %s" % (stringToGrep, crabLog)
    
    if os.system(cmd) == 0:
        return True
    return False

def GetTaskDashboardURL(datasetPath):                                                                                                                  
    '''                                                                                                                                                
    Call the "grep" command to look for the dashboard URL from the crab.log file                                                                       
    of a given dataset. It uses as input parameter the absolute path of the task dir (datasetPath)                                                     
    '''                                                                                                                                                
                                                                                                                                                       
    # Variable Declaration                                                                                                                             
    crabLog      = os.path.join(datasetPath, "crab.log")                                                                                               
    grepFile     = os.path.join(datasetPath, "grep.tmp")                                                                                               
    stringToGrep = "Dashboard monitoring URL"                                                                                                          
    cmd          = "grep '%s' %s > %s" % (stringToGrep, crabLog, grepFile )                                                                            
    dashboardURL = "UNKNOWN"                                                                                                                           
                                                                                                                                                       
    # Execute the command                                                                                                                              
    if os.system(cmd) == 0:                                                                                                                            
                                                                                                                                                       
        if os.path.exists( grepFile ):                                                                                                                 
            results      = [i for i in open(grepFile, 'r').readlines()]                                                                                
            dashboardURL = FindBetween( results[0], "URL:\t", "\n" )                                                                                   
            os.system("rm -f %s " % (grepFile) )                                                                                                       
        else:                                                                                                                                          
            raise Exception("File %s not found!" % (grepFile) )                                                                                        
    else:                                                                                                                                              
        raise Exception("Could not execute command %s" % (cmd) )                                                                                       
    return dashboardURL

def GetTaskStatus(datasetPath):
    '''
    Call the "grep" command to look for the "Task status" from the crab.log file
    of a given dataset. It uses as input parameter the absolute path of the task dir (datasetPath)
    '''
    
    # Variable Declaration
    crabLog      = os.path.join(datasetPath, "crab.log")
    grepFile     = os.path.join(datasetPath, "grep.tmp")
    #stringToGrep = "Task status:"
    #stringToGrep = "Status on the CRAB server:"
    #stringToGrep = "Jobs status:"
    stringToGrep  = "Status on the scheduler:"
    cmd          = "grep '%s' %s > %s" % (stringToGrep, crabLog, grepFile )
    status       = "UNKNOWN"
            
    if not os.path.exists( crabLog ):
        raise Exception("File %s not found!" % (crabLog) )
            
    # Execute the command
    if os.system(cmd) == 0:
    
        if os.path.exists( grepFile ):
            results = [i for i in open(grepFile, 'r').readlines()]
            status  = FindBetween( results[-1], stringToGrep, "\n" )
            os.system("rm -f %s " % (grepFile) )
        else:
            raise Exception("File %s not found!" % (grepFile) )
    else:
        raise Exception("Could not execute command %s" % (cmd) )
    return status

def GetTaskReports(datasetPath, opts):
    '''
    Execute "crab status", get task logs and output.
    Resubmit or kill task according to user options. 
    '''
    report = None
    
    # Get all files under <dataset_dir>/results/
    files = Execute("ls %s" % os.path.join( datasetPath, "results") )
    print "check datasetPath",datasetPath    
    try:
        d = GetBasename(datasetPath)

        # Execute "crab status --dir=datasetPath"
        result = crabCommand('status', dir=datasetPath)
    
        # Get CRAB task status
        status = GetTaskStatus(d).replace("\t", "")

        # Get CRAB task dashboard URL
        dashboardURL = GetTaskDashboardURL(d)

        # Assess JOB success/failure for task
        idle, running, finished, transferring, failed, retrievedLog, retrievedOut, eosLog, eosOut = RetrievedFiles(datasetPath, result, dashboardURL, True, opts)

        # Get the task logs & output ?
        GetTaskLogs(datasetPath, retrievedLog, finished)
       
        # Get the task output
        GetTaskOutput(datasetPath, retrievedOut, finished)
    
        # Resubmit task if failed jobs found
        ResubmitTask(datasetPath, failed)
    
        # Kill task which are active
        KillTask(datasetPath)
        # Count retrieved/all jobs
        retrieved = min(finished, retrievedLog, retrievedOut)
        alljobs   = GetTotalJobsFromStatus(result)
        
        # Append the report
        report = Report(datasetPath, alljobs, idle, retrieved, running, finished, failed, transferring, retrievedLog, retrievedOut, eosLog, eosOut, status, dashboardURL)

        # Determine if task is DONE or not
        if retrieved == alljobs and retrieved > 0:
            absolutePath = os.path.join(datasetPath, "crab.log")
            os.system("sed -i -e '$a\DONE! (Written by multicrabCheck.py)' %s" % absolutePath )

    # Catch exceptions (Errors detected during execution which may not be "fatal")
    except:
        msg = sys.exc_info()[1]
        report = Report(datasetPath, "?", "?", "?", "?", "?", "?", "?", "?", "?", "?", "?", "?", "?")
        Print("crab status failed with message %s. Skipping ..." % ( msg ), True)
    return report

def Verbose(msg, printHeader=False):
    '''                                                                                                                                                    
    Calls Print() only if verbose options is set to true.                                                                                                  
    '''
    if not opts.verbose:
        return
    Print(msg, printHeader)
    return

def Print(msg, printHeader=True):
    '''
    Simple print function. If verbose option is enabled prints, otherwise does nothing.
    '''
    fName = __file__.split("/")[-1]
    if printHeader:
        print "=== ", fName
    print "\t", msg
    return

def GetBasename(fullPath):
    return os.path.basename(fullPath)

def FindBetween(myString, first, last ):
    try:
        start = myString.index( first ) + len( first )
        end   = myString.index( last, start )
        return myString[start:end]
    except ValueError:
        return ""

if __name__ == "__main__":

    # Default Values                                                                                                                        
    VERBOSE = False
    PSET    = "PSet.py"
    SITE    = "T2_FI_HIP"
    DIRNAME = ""

    parser = OptionParser(usage="Usage: %prog [options]")
    parser.add_option("--create", dest="create", default=False, action="store_true",
                      help="Flag to create a CRAB job [default: False")

    parser.add_option("--status", dest="status", default=False, action="store_true",
                      help="Flag to check the status of all CRAB jobs [default: False")

    parser.add_option("--get", dest="get", default=False, action="store_true",
                      help="Get output and log files of finished jobs [defaut: False]")

    parser.add_option("--log", dest="log", default=False, action="store_true",
                      help="Get log files of finished jobs [defaut: False]")

    parser.add_option("--out", dest="out", default=False, action="store_true",
                      help="Get output files of finished jobs [defaut: False]")

    parser.add_option("--resubmit", dest="resubmit", default=False, action="store_true",
                      help="Resubmit all failed jobs [defaut: False]")

    parser.add_option("--kill", dest="kill", default=False, action="store_true",
                      help="Kill all submitted jobs [defaut: False]")

    parser.add_option("-v", "--verbose", dest="verbose", default=VERBOSE, action="store_true",
                      help="Verbose mode for debugging purposes [default: %s]" % (VERBOSE))

    parser.add_option("-p", "--pset", dest="pset", default=PSET, type="string",
                      help="The python cfg file to be used by cmsRun [default: %s]" % (PSET))

    parser.add_option("-n", "--name", dest="name", default=None, type="string",
                      help="The python cfg file to be used by cmsRun [default: \"\"]")

    parser.add_option("-d", "--dir", dest="dirName", default=DIRNAME, type="string",
                      help="Custom name for CRAB directory name [default: %s]" % (DIRNAME))

    parser.add_option("-s", "--site", dest="storageSite", default=SITE, type="string",
                      help="Site where the output will be copied to [default: %s]" % (SITE))

    parser.add_option("-i", "--includeTasks", dest="includeTasks", default="None", type="string",
                      help="Only perform action for this dataset(s) [default: \"\"]")

    parser.add_option("-e", "--excludeTasks", dest="excludeTasks", default="None", type="string",
                      help="Exclude this dataset(s) from action [default: \"\"]")

    (opts, args) = parser.parse_args()

    if opts.create == False and opts.dirName == "":
        opts.dirName = os.getcwd()

    if opts.create == True and opts.status == True:
        raise Exception("Cannot both create and check a CRAB job!")

    if opts.create == True:
        sys.exit( CreateJob(opts, args) )
    elif opts.status == True or opts.get == True or opts.out == True or opts.log == True or opts.resubmit == True or opts.kill == True:
        if opts.dirName == "":
            raise Exception("Must provide a multiCRAB dir with the -d option!")
        else:
            sys.exit( CheckJob(opts, args) )
    else:
        raise Exception("Must either create or check a CRAB job!")


