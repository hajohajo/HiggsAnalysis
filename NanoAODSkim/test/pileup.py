#!/usr/bin/env python

import os
import sys
import re

from multicrab import Dataset,alldatasets,GetRequestName

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData
calcMode       = "true"
maxPileupBin   = "100" 
numPileupBins  = "100"
pileupHistName = "pileup"
PileUpJSON_2016 = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt"
PileUpJSON_2017 = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt"
PileUpJSON_2018 = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt"
#PileUpJSON = PileUpJSON_2016
#PileUpJSON = PileUpJSON_2017
#PileUpJSON = PileUpJSON_2018
#
# Recommended minimum bias xsection                                                                                                                                         
minBiasXsecNominal = 69200 #from https://twiki.cern.ch/twiki/bin/viewauth/CMS/POGRecipesICHEP2016
minBiasXsec = minBiasXsecNominal
puUncert    = 0.05 

def usage():
    print
    print "### Usage:  ",sys.argv[0]," <multicrabdir>"
    print

def CallPileupCalc(inputFile,inputLumiJSON,minBiasXsec,fOUT,pileupHistName):
    cmd = ["./pileupCalc.py", "-i", inputFile, "--inputLumiJSON", inputLumiJSON, "--calcMode", calcMode,
           "--minBiasXsec", minBiasXsec, "--maxPileupBin", maxPileupBin, "--numPileupBins", numPileupBins,
           fOUT, "--pileupHistName", pileupHistName]
    sys_cmd = " ".join([str(c) for c in cmd])
    print sys_cmd
    os.system(sys_cmd)

def main():

    if len(sys.argv) == 1:
        usage()
        sys.exit()

    multicrabdir = sys.argv[1]
    for dataset in alldatasets:
        if not dataset.isData():
            continue

        dsetname = GetRequestName(dataset)
        path = os.path.join(multicrabdir,dsetname)
        #print path
        if os.path.exists(path):
            print dsetname
            fOUT = os.path.join(multicrabdir,dsetname,"results","PileUp.root")
            if dataset.getYear() == "2016":
                PileUpJSON = PileUpJSON_2016
            if dataset.getYear() == "2017":
                PileUpJSON = PileUpJSON_2017
            if dataset.getYear() == "2018":
                PileUpJSON = PileUpJSON_2018

            inputLumiJSON = PileUpJSON
            inputFile = dataset.lumiMask # crab report not working for nanoaod yet, assuming 100% jobs successfull
            hName = pileupHistName
            minBiasXsec = minBiasXsecNominal
            CallPileupCalc(inputFile,inputLumiJSON,minBiasXsec,fOUT,hName)

            minBiasXsec_up = minBiasXsec*(1+puUncert)
            fOUT_up        = fOUT.replace(".root","_up.root")
            hName_up       = pileupHistName+"_up"
            CallPileupCalc(inputFile,inputLumiJSON,minBiasXsec_up,fOUT_up,hName_up)

            minBiasXsec_down = minBiasXsec*(1-puUncert)
            fOUT_down        = fOUT.replace(".root","_down.root")
            hName_down       = pileupHistName+"_down"
            CallPileupCalc(inputFile,inputLumiJSON,minBiasXsec_down,fOUT_down,hName_down)

if __name__ == "__main__":
    main()
