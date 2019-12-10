#!/usr/bin/env python
import os
import sys
import math
import PSet
import re

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 

#this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import ROOT

from array import array


class Counter():
    def __init__(self,name):
        self.name  = name
        self.count = 0
        
    def increment(self):
        self.count += 1
        
    def Print(self):
        self.name += " "
        while len(self.name) < 39:
            self.name += "."
        print self.name,self.count

class Skim(Module):
    def __init__(self):
        self.cControl       = Counter("Skim: control")
	self.cControl.increment()
        self.cAllEvents     = Counter("Skim: All events")
        self.cTrigger       = Counter("Skim: Trigger selection")
	self.cMETCleaning   = Counter("Skim: METCleaning")
        self.cPassedEvents  = Counter("Skim: Passed events")

	self.objs = []


    def __del__(self):
        self.cAllEvents.Print()
        self.cTrigger.Print()
	self.cMETCleaning.Print()
        self.cPassedEvents.Print()

    def beginJob(self):

        self.h_pileup = ROOT.TH1F('pileup','',100,0,100)
        self.addObject(self.h_pileup)

	self.h_skimcounter = ROOT.TH1F("skimCounter","",5,0,5)
        self.addObject(self.h_skimcounter)


    def endJob(self):
	pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.dir = outputFile.mkdir("configInfo")


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):

	outputFile.cd()

	self.dir.cd()
        self.h_pileup.Write()

        
        self.h_skimcounter.SetBinContent(1,self.cControl.count)
        self.h_skimcounter.GetXaxis().SetBinLabel(1,self.cControl.name)
        self.h_skimcounter.SetBinContent(2,self.cAllEvents.count)
        self.h_skimcounter.GetXaxis().SetBinLabel(2,self.cAllEvents.name)
	self.h_skimcounter.SetBinContent(3,self.cTrigger.count)
        self.h_skimcounter.GetXaxis().SetBinLabel(3,self.cTrigger.name)
        self.h_skimcounter.SetBinContent(4,self.cMETCleaning.count)
        self.h_skimcounter.GetXaxis().SetBinLabel(4,self.cMETCleaning.name)
        self.h_skimcounter.SetBinContent(5,self.cPassedEvents.count)
        self.h_skimcounter.GetXaxis().SetBinLabel(5,self.cPassedEvents.name)
	self.h_skimcounter.Write()


    def analyze(self, event):

        if event._tree.GetListOfBranches().FindObject("Pileup_nTrueInt"):
            self.h_pileup.Fill(event.Pileup_nTrueInt)
        
        self.cAllEvents.increment()
        
        # selection
        # 2016 trigger
        triggerDecision = False
	if hasattr(event._tree, 'HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90'):
            triggerDecision = triggerDecision or event.HLT_LooseIsoPFTau50_Trk30_eta2p1

        # 2017 trigger
        if hasattr(event._tree, 'MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100'):
            triggerDecision = triggerDecision or event.MediumChargedIsoPFTau50_Trk30_eta2p1_1pr

        # 2018 trigger
        if hasattr(event._tree, 'MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100'):
            triggerDecision = triggerDecision or event.MediumChargedIsoPFTau50_Trk30_eta2p1_1pr

        if not triggerDecision:
            return False
	self.cTrigger.increment()

        # MET cleaning
	cleaningDecision = event.Flag_HBHENoiseFilter and event.Flag_HBHENoiseIsoFilter and event.Flag_EcalDeadCellTriggerPrimitiveFilter and event.Flag_eeBadScFilter and event.Flag_goodVertices and event.Flag_globalTightHalo2016Filter and event.Flag_BadPFMuonFilter and event.Flag_BadChargedCandidateFilter
	if not cleaningDecision:
            return False
	self.cMETCleaning.increment()

        self.cPassedEvents.increment()

	return True

SkimModule = lambda : Skim()



files=["root://xrootd-cms.infn.it//store/data/Run2016H/Tau/NANOAOD/Nano1June2019-v1/240000/FB1AF208-9FEC-0445-AB21-772D27660951.root"]

if not "CMSSITE" in os.environ:
    p=PostProcessor(".",files,"",modules=[SkimModule()],provenance=True,fwkJobReport=True,haddFileName="events.root",jsonInput=runsAndLumis(),maxEntries=20000)
else:
    p=PostProcessor(".",inputFiles(),"",modules=[SkimModule()],provenance=True,fwkJobReport=True,haddFileName="events.root",jsonInput=runsAndLumis())

p.run()

os.system("mv *_Skim.root events.root")

print "DONE"

