#!/usr/bin/env python
'''
DESCRIPTIONM:
Script that plots Data/MC for all histograms under a given folder (passsed as option to the script)
Good for sanity checks for key points in the cut-flow


USAGE:
./plotDataMC.py -m <pseudo_mcrab_directory> [opts]


EXAMPLES:
./plotDataMC.py --folder ForDataDrivenCtrlPlots --ratio -s png -m
./plotDataMC.py --folder PUDependency --ratio -s png -m 
./plotDataMC.py --folder counters/weighted --ratio -s png -m 
./plotDataMC.py --folder muSelection_ --ratio -s png -m
./plotDataMC.py --folder tauSelection_ --ratio -s png -m
./plotDataMC.py --folder jetSelection_ --ratio -s png -m


LAST USED:
./plotDataMC.py --folder ForDataDrivenCtrlPlots --ratio -s png -m
./plotDataMC.py --folder PUDependency --ratio -s png -m 
./plotDataMC.py --folder counters/weighted --ratio -s png -m 
./plotDataMC.py --folder muSelection_ --ratio -s png -m
./plotDataMC.py --folder tauSelection_ --ratio -s png -m
./plotDataMC.py --folder jetSelection_ --ratio -s png -m


'''
#================================================================================================ 
# Imports
#================================================================================================ 
import sys
import math
import copy
import os
import array
import re
from optparse import OptionParser

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *

import HiggsAnalysis.NtupleAnalysis.tools.dataset as dataset
import HiggsAnalysis.NtupleAnalysis.tools.histograms as histograms
import HiggsAnalysis.NtupleAnalysis.tools.counter as counter
import HiggsAnalysis.NtupleAnalysis.tools.tdrstyle as tdrstyle
import HiggsAnalysis.NtupleAnalysis.tools.styles as styles
import HiggsAnalysis.NtupleAnalysis.tools.plots as plots
import HiggsAnalysis.NtupleAnalysis.tools.aux as aux
import HiggsAnalysis.NtupleAnalysis.tools.systematics as systematics
import HiggsAnalysis.NtupleAnalysis.tools.crosssection as xsect
import HiggsAnalysis.NtupleAnalysis.tools.multicrabConsistencyCheck as consistencyCheck
import HiggsAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles


#================================================================================================ 
# Variable definition                                                                                                                                                                                                             
#================================================================================================ 
ss = ShellStyles.SuccessStyle()
ns = ShellStyles.NormalStyle()
ts = ShellStyles.NoteStyle()
hs = ShellStyles.HighlightAltStyle()
ls = ShellStyles.HighlightStyle()
es = ShellStyles.ErrorStyle()
cs = ShellStyles.CaptionStyle()

#================================================================================================ 
# Function Definition
#================================================================================================ 
def Print(msg, printHeader=False):
    fName = __file__.split("/")[-1]
    if printHeader==True:
        print "=== ", fName
        print "\t", msg
    else:
        print "\t", msg
    return

def Verbose(msg, printHeader=True, verbose=False):
    if not opts.verbose:
        return
    Print(msg, printHeader)
    return

def GetLumi(datasetsMgr):
    lumi = 0.0
    for d in datasetsMgr.getAllDatasets():
        if d.isMC():
            continue
        else:
            lumi += d.getLuminosity()
    Verbose("Luminosity = %s (pb)" % (lumi), True)
    return lumi

def GetDatasetsFromDir(opts):
    Verbose("Getting datasets")
    
    if (not opts.includeOnlyTasks and not opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode, 
                                                        analysisName=opts.analysisName,
                                                        optimizationMode=opts.optMode)
    elif (opts.includeOnlyTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode,
                                                        analysisName=opts.analysisName,
                                                        includeOnlyTasks=opts.includeOnlyTasks,
                                                        optimizationMode=opts.optMode)
    elif (opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode,
                                                        analysisName=opts.analysisName,
                                                        excludeTasks=opts.excludeTasks,
                                                        optimizationMode=opts.optMode)
    else:
        raise Exception("This should never be reached")
    return datasets
    
def main(opts):

    # Apply TDR style
    style = tdrstyle.TDRStyle()
    style.setOptStat(False)
    style.setGridX(opts.gridX)
    style.setGridY(opts.gridY)

    optModes = [""]
    if opts.optMode != None:
        optModes = [opts.optMode]
        
    # For-loop: All opt Mode
    for opt in optModes:
        opts.optMode = opt

        # Setup & configure the dataset manager 
        datasetsMgr = GetDatasetsFromDir(opts)
        datasetsMgr.updateNAllEventsToPUWeighted()
        datasetsMgr.loadLuminosities() # from lumi.json


        # Merge histograms (see NtupleAnalysis/python/tools/plots.py) 
        plots.mergeRenameReorderForDataMC(datasetsMgr) 

        # Set signal cross-section
        newOrder = datasetsMgr.getAllDatasetNames()
        for d in datasetsMgr.getAllDatasets():
            if "ChargedHiggs" in d.getName():
                datasetsMgr.getDataset(d.getName()).setCrossSection(1.0) # ATLAS 13 TeV H->tb exclusion limits
                newOrder.remove(d.getName())
                newOrder.insert(1, d.getName())

        if opts.verbose:
            datasetsMgr.PrintCrossSections()
            datasetsMgr.PrintLuminosities()
            datasetsMgr.PrintInfo()

        # Ensure that signal dataset is plotted last
        if 1:
            datasetsMgr.selectAndReorder(newOrder)
          
        # Merge EWK samples
        if opts.mergeEWK:
            datasetsMgr.merge("EWK", aux.GetListOfEwkDatasets())
            plots._plotStyles["EWK"] = styles.getAltEWKStyle()
            
        # Print dataset information
        datasetsMgr.PrintInfo()

        # Do Data-MC histograms
        folder     = opts.folder
        histoList  = datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDirectoryContent(folder)        
        histoPaths = [os.path.join(folder, h) for h in histoList]
        myHistos   = []
        skipList   = [
            "counters",
            "Weighting",
            "SplittedBinInfo",
            "ForDataDrivenCtrlPlots",
            "ForDataDrivenCtrlPlotsEWKFakeB",
            "ForDataDrivenCtrlPlotsEWKGenuineB",
            "ForDataDrivenCtrlPlotsEWKFakeTaus",
            "ForDataDrivenCtrlPlotsEWKGenuineTaus",
            "AngularCuts_Collinear",
            "AngularCuts_BackToBack",
            "PUDependency",
            "jetSelection_",
            "bjetSelection_",
            "eSelection_Veto",
            "muSelection_",
            "tauSelection_",
            "metSelection_",
            "topSelectionBDT_",
            "config",
            "NSelectedVsRunNumber",
            ]

        # For-loop: All histograms
        for h in histoPaths:
            myHistos.append(h)

        # For-loop: All histos
        for i, h in enumerate(myHistos, 1):
            msg   = "{:<9} {:>3} {:<1} {:<3} {:<50}".format("Histogram", "%i" % i, "/", "%s:" % (len(myHistos)), h)
            Print(ss + msg + ns, i==1)
            
            if "_Vs_" in h or "etaphi" in h:
                Print 
                continue

            if opts.folder == "":
                if h in skipList:
                    continue
            DataMCHistograms(datasetsMgr, h)
        
    Print("All plots saved under directory %s" % (ts + aux.convertToURL(opts.saveDir, opts.url) + ns), True)    
    return

def GetHistoKwargs(h, opts):

    # Common bin settings
    _legRM  = {"dx": -10000.23, "dy": -10000.01, "dh": -0.1}
    _legNE  = {"dx": -0.23, "dy": -0.01, "dh": 0.1}
    _legNW  = {"dx": -0.52, "dy": -0.01, "dh": 0.1}
    _legSE  = {"dx": -0.23, "dy": -0.50, "dh": 0.1}
    _legSW  = {"dx": -0.52, "dy": -0.50, "dh": 0.1}
    _legCE  = {"dx": -0.08, "dy": -0.44, "dh": 0.1}
    _logY   = True
    _yLabel = "Events / %.0f "
    if "_AfterAllSelections" in h or "_AfterTopSelection" in h or "_AfterMetSelection" in h:
        _yMin = 1e-2
    else:
        _yMin = 1e-1

    if _logY:
        _yMaxF = 20
    else:
        _yMaxF = 1.2

    _kwargs = {
        "ylabel"           : _yLabel,
        "rebinX"           : 1,
        "rebinY"           : None,
        "ratioType"        : "errorScale",
        "ratioErrorOptions": {"numeratorStatSyst": False},
        "ratioYlabel"      : "Data/MC ",
        "ratio"            : opts.ratio, 
        "stackMCHistograms": True,
        "ratioInvert"      : False, 
        "addMCUncertainty" : True, 
        "addLuminosityText": True,
        "addCmText"        : True,
        "cmsExtraText"     : "Very Preliminary",
        "opts"             : {"ymin": _yMin, "ymaxfactor": _yMaxF},
        "opts2"            : {"ymin": 0.0, "ymax": 3.0},
        "divideByBinWidth" : False,
        "log"              : _logY,
        "moveLegend"       : _legNE,
        "moveBlindedText"  : {"dx": -0.23, "dy": +0.08, "dh": 0.0},
        }

    kwargs = copy.deepcopy(_kwargs)
    
    if "pt" in h.lower():
        units                = "GeV"
        kwargs["rebinX"]     = 5
        kwargs["ylabel"]     = _yLabel + units
        kwargs["moveLegend"] = _legNE
        kwargs["opts"]   = {"xmax": +800.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "DeltaEta" in h:
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["opts"]   = {"xmax": +6.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "CollinearAngularCuts" in h:
        units = "#circ" #^{#circ}"
        if "Jet" in h:
            kwargs["xlabel"] = "R_{coll}^{jet} (%s)" % units
            kwargs["moveLegend"] = _legNW
        if "Minimum" in h:
            kwargs["xlabel"] = "R_{coll}^{min} (%s)" % units                
            kwargs["moveLegend"] = _legNE
        kwargs["ylabel"] = _yLabel + units
        #kwargs["cutBox"] = {"cutValue": 400.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["rebinX"] = 4
        kwargs["opts"]   = {"ymin": _yMin, "ymaxfactor": _yMaxF*2.0}

    if "BackToBackAngularCuts" in h:
        units = "#circ" #^{#circ}"
        if "Jet" in h:
            kwargs["xlabel"] = "R_{bb}^{jet} (%s)" % units
            kwargs["moveLegend"] = _legNW
        if "Minimum" in h:
            kwargs["xlabel"] = "R_{bb}^{min} (%s)" % units                
            kwargs["moveLegend"] = _legNE
        kwargs["ylabel"] = _yLabel + units
        #kwargs["cutBox"] = {"cutValue": 400.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["rebinX"] = 4
        kwargs["opts"]   = {"ymin": _yMin, "ymaxfactor": _yMaxF*2.0}

    if "DeltaR" in h or "DeltaY" in h or "DR" in h:
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["opts"]   = {"xmax": +8.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if h == "TopCandMass": #after all cuts
        units            = "GeV" #"GeV/c^{2}"
        kwargs["rebinX"] = 2
        kwargs["xlabel"] = "m_{jjb}^{BDT} (%s)" % units
        kwargs["ylabel"] = _yLabel + units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +500.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 173.21, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if h == "BDTmultiplicity":
        kwargs["rebinX"] = 1
        kwargs["xlabel"] = "top multiplicity"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +40.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 1.5, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "BDT" in h or "topbdt" in h.lower():
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["cutBox"] = {"cutValue": -0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -1.0, "xmax": +1.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["rebinX"] = 2        
        kwargs["xlabel"] = "BDTG scrore"
        if "BDT_Selected" in h:
            kwargs["moveLegend"] = _legNW

    if h == "BDTGresponse":
        kwargs["rebinX"] = 1
        kwargs["xlabel"] = "BDT discriminant"
        kwargs["ylabel"] = "Events / %.1f "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +1.0, "ymin": _yMin, "ymax": 3e5} #_yMaxF}
        kwargs["cutBox"] = {"cutValue": 0.9, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["moveLegend"] = {"dx": -0.49, "dy": -0.5, "dh": 0.0}

    if h == "LdgBjetPt":
        ROOT.gStyle.SetNdivisions(8, "X")
        units            = "GeV" #"GeV/c"
        #kwargs["rebinX"] = 2
        kwargs["xlabel"] = "p_{T} (%s)" % units
        kwargs["ylabel"] = _yLabel + units
        kwargs["cutBox"] = {"cutValue": 40.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
   
    if "bdisc" in h.lower():
        kwargs["rebinX"] = 2
        kwargs["xlabel"] = "b-jet discriminant"
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["cutBox"] = {"cutValue": 0.8484, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.4, "xmax": +1.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["moveLegend"] = {"dx": -0.49, "dy": -0.0, "dh": 0.1}
 
    if "dgTrijetBJetPt" in h:
        units            = "GeV" #"GeV/c"
        #kwargs["rebinX"] = 2
        kwargs["xlabel"] = "p_{T} (%s)" % units
        kwargs["ylabel"] = _yLabel + units
        kwargs["cutBox"] = {"cutValue": 40.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +600.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "DiJetEta" in h:
        kwargs["rebinX"] = 1
        kwargs["xlabel"] = "#eta"
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -2.5, "xmax": +2.5, "ymin": _yMin, "ymaxfactor": _yMaxF}
        
    if "DiJetPt" in h:
        units            = "GeV" #"GeV/c"
        #kwargs["rebinX"] = 2
        kwargs["xlabel"] = "p_{T,jj}^{BDT} (%s)" % units
        kwargs["ylabel"] = _yLabel + units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +800.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if  "TopJet1Eta" in h or "Jet2Eta" in h:
        kwargs["rebinX"] = 1
        kwargs["xlabel"] = "#eta"
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": -2.5, "xmax": +2.5, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if  "TopJet1Pt" in h or "ldgTrijetJet2Pt" in h:
        units            = "GeV" #"GeV" #"GeV/c"
        #kwargs["rebinX"] = 2
        kwargs["xlabel"] = "p_{T} (%s)" % units
        kwargs["ylabel"] = _yLabel + units
        kwargs["cutBox"] = {"cutValue": 40.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +600.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "TopPt" in h:
        ROOT.gStyle.SetNdivisions(8, "X")
        units            = "GeV" #"GeV" #"GeV/c"
        kwargs["rebinX"] = systematics._dataDrivenCtrlPlotBinning["LdgTrijetPt_AfterAllSelections"]
        kwargs["xlabel"] = "p_{T} (%s)" % (units)
        kwargs["ylabel"] = _yLabel + units

    if "nselectedcleaned" in h.lower():
        kwargs["rebinX"] = 1
        kwargs["xlabel"] = "top multiplicity"
        kwargs["ylabel"] = _yLabel
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +8.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 1.5, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "topmultiplicity" in h.lower():
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +20.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 1.5, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    #if "topbdt_allcandidates" in h.lower():
    #    kwargs["rebinX"] = 10
      
    if "nalltops" in h.lower():
        kwargs["rebinX"] = 1
        kwargs["xlabel"] = "top multiplicity"
        #kwargs["ylabel"] = _yLabel
        kwargs["ylabel"] = "Events / %.1f "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +15.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 1.5, "fillColor": 16, "box": False, "line": False, "greaterThan": True}

    if "nallcleanedtops" in h.lower():
        kwargs["rebinX"] = 1
        kwargs["xlabel"] = "top multiplicity"
        kwargs["ylabel"] = _yLabel
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +8.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 1.5, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
                
    if "counters" in opts.folder:
        ROOT.gStyle.SetLabelSize(16.0, "X")
        kwargs["moveLegend"] = {"dx": -0.08, "dy": 0.0, "dh": 0.1}
        kwargs["moveLegend"] = _legNE

    if "fatjetNPassed" in h:
        kwargs["opts"]   = {"xmax": 4.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        #kwargs["cutBox"] = {"cutValue": opts.signalMass, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "fatjetPt" in h:
        units            = "GeV" #"GeV/c"
        kwargs["rebinX"] = 5
        kwargs["xlabel"] = "p_{T} (%s)" % units
        #kwargs["cutBox"] = {"cutValue": 0.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 450.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "IsolPt" in h:
        ROOT.gStyle.SetNdivisions(8, "X")
        units            = "GeV" #"GeV/c"
        kwargs["rebinX"] = 2
        kwargs["xlabel"] = "p_{T} (%s)" % units
        kwargs["ylabel"] = _yLabel + units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +400.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        
    if "RelIsoBefore" in h:
        kwargs["ylabel"] = "Events / %.0f "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +60.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["rebinX"] = 5

    if "RelIsoAfter" in h or "RelIsoPassed" in h:
        kwargs["ylabel"] = "Events / %.0f "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +20.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["rebinX"] = 5

    if "RelIsoAll" in h:
        kwargs["ylabel"] = "Events / %.0f "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +200.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["rebinX"] = 10

    if "IsolMiniIso" in h or "MiniIso" in h:
        kwargs["rebinX"]     = 10
        kwargs["ylabel"]     = "Events / %.2f "
        kwargs["xlabel"]     = "relative mini-isolation"
        kwargs["opts"]       = {"xmin": 0.0, "xmax": 60.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["moveLegend"] = _legNE

        if "after" in h.lower() or "passed" in h.lower():
            kwargs["cutBox"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
            kwargs["rebinX"] = 1
            kwargs["ylabel"] = "Events / %.2f "
            kwargs["opts"]   = {"xmin": 0.0, "xmax": +5.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "electronPt" in h:
        units            = "GeV" #"GeV/c"
        kwargs["xlabel"] = "p_{T} (%s)" % units
        kwargs["ylabel"] = _yLabel + units

    #if "jetPtAll" in h or "jetPtPassed" in h:
    if "jetpt" in h.lower():
        units            = "GeV" #"GeV/c"
        kwargs["xlabel"] = "p_{T} (%s)" % units
        kwargs["ylabel"] = _yLabel + units
        kwargs["opts"]    = {"xmin": 0.0, "xmax": 700.0, "ymin": 1e-2, "ymaxfactor": _yMaxF}

    if h == "tauNpassed":
        units            = "GeV" #"GeV/c"
        kwargs["xlabel"] = "#tau-jet multiplicity"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 7.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 2.5, "fillColor": 16, "box": False, "line": True, "greaterThan": False}

    if "HT" in h or "JT" in h:
        ROOT.gStyle.SetNdivisions(8, "X")
        kwargs["rebinX"] = 5 #systematics._dataDrivenCtrlPlotBinning["HT_AfterAllSelections"]
        units            = "GeV"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +1000.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        #kwargs["ylabel"] = "<%s>" % (_yLabel + units)
        kwargs["ylabel"] = "%s" % (_yLabel + units)
        kwargs["xlabel"] = "H_{T} (%s)" % (units)
        if "JT" in h:            
            kwargs["xlabel"] = "J_{T} (%s)" % (units)

    regex = re.compile('selectedJets.*JetPt')
    if(regex.search(h)):
        ROOT.gStyle.SetNdivisions(8, "X")
        units            = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)" % units
        kwargs["ylabel"] = _yLabel + units
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["rebinX"] = 2
        if "second" in h.lower():
            kwargs["opts"] = {"xmin": 0.0, "xmax": +500.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        if "third" in h.lower() or "fourth" in h.lower():
            kwargs["opts"] = {"xmin": 0.0, "xmax": +300.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        if "fifth" in h.lower():
            kwargs["opts"] = {"xmin": 0.0, "xmax": +300.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        if "sixth" in h.lower():
            kwargs["opts"] = {"xmin": 0.0, "xmax": +100.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        
    regex = re.compile('selectedBJets.*JetPt')
    if(regex.search(h)):
        #ROOT.gStyle.SetNdivisions(8, "X")
        kwargs["rebinX"] = 2
        kwargs["cutBox"] = {"cutValue": 40.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    regex = re.compile('selectedBJets.*BDisc')
    if(regex.search(h)):
        kwargs["rebinX"] = 2
        ROOT.gStyle.SetNdivisions(8, "X")
        kwargs["xlabel"] = "b-tag discriminator"
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["cutBox"] = {"cutValue": 0.8484, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.4, "xmax": +1.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["moveLegend"] = {"dx": -0.49, "dy": -0.00, "dh": 0.1}

    if "btagdiscriminator" in h.lower():
        kwargs["rebinX"] = 2
        kwargs["xlabel"] = "b-jet discriminant"
        kwargs["ylabel"] = "Events / %.2f "
        ROOT.gStyle.SetNdivisions(8, "X")
        kwargs["cutBox"] = {"cutValue": 0.8484, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +1.0, "ymin": 1e0, "ymaxfactor": _yMaxF}
        kwargs["moveLegend"] = {"dx": +1000.0, "dy": +1000.00, "dh": -10000.1}

    if "toptagsf" in h.lower():
        kwargs["rebinX"] = 5
        kwargs["xlabel"] = "top-tag SF"
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 5.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if h == "btagSF":
        kwargs["xlabel"] = "b-jet SF"
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["opts"]   = {"xmin": 0.5, "xmax": 2.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["rebinX"] = 5

    if "met_" in h.lower():
        units            = "GeV"
        kwargs["rebinX"] =  systematics._dataDrivenCtrlPlotBinning["MET_AfterAllSelections"] #4
        kwargs["xlabel"] = "E_{T}^{miss} (%s)" % (units)
        kwargs["ylabel"] = _yLabel + units
        kwargs["opts"]   = {"ymin": _yMin, "xmax": 500.0, "ymaxfactor": _yMaxF}
        if "METFilter" in h:
            kwargs["rebinX"] = 1
            kwargs["xlabel"] = ""
            kwargs["ylabel"] = "Events / %.0f "
            kwargs["opts"]   = {"xmin": 0.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
            
        if "metphi" in h.lower():
            units            = "rads"
            kwargs["rebinX"] =  2
            kwargs["xlabel"] = "E_{T}^{miss} #phi (%s)" % (units)
            kwargs["ylabel"] = _yLabel + units
            kwargs["opts"]   = {"xmin": -3.2, "xmax": 3.2, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "LdgTrijetBjetPt" in h:
        units            = "GeV" #"GeV/c"
        #kwargs["rebinX"] = 2
        kwargs["xlabel"] = "p_{T} (%s)" % units
        kwargs["ylabel"] = _yLabel + units
        # kwargs["cutBox"] = {"cutValue": 40.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        if "AllSelections" in h:
            kwargs["opts"]   = {"xmin": 0.0, "xmax": +600.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        else:
            kwargs["opts"]   = {"xmin": 0.0, "xmax": +1000.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "MHT" in h:
        kwargs["rebinX"] = 5 #systematics._dataDrivenCtrlPlotBinning["MET_AfterAllSelections"] #2
        units            = "GeV"
        kwargs["xlabel"] = "MHT (%s)" % units
        kwargs["ylabel"] = _yLabel + units
        kwargs["opts"]   = {"xmin": 0.0, "xmax": +500.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "NBjets" in h:
        kwargs["xlabel"] = "b-jet multiplicity"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 8.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 1.5, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "Njets" in h:
        kwargs["xlabel"] = "jet multiplicity"
        kwargs["opts"]   = {"xmin": 3.0, "xmax": 12.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 3.5, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "vtx" in h.lower():
        units            = ""
        kwargs["xlabel"] = "vertex multiplicity"
        kwargs["rebinX"] = 5
        kwargs["divideByBinWidth"] = False
        kwargs["ylabel"] = "%s" % (_yLabel + units)

    if "NVertices" in h:
        units = ""
        kwargs["rebinX"] = 5
        kwargs["ylabel"] = "%s" % (_yLabel + units)
        kwargs["xlabel"] = "vertex multiplicity"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 65.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "counters" in opts.folder.lower():
        kwargs["moveLegend"] = _legSW
        kwargs["opts"]       = {"ymin": 1e0, "ymax": 1e9}
        if h == "counter":
            xMin =  8.0
            xMax = 25.0 # 27.0 is actual last. using 1 before because blinding fails there
            #xMin =  8.0
            #xMax = 25.0 # 26.0 is actual last. using 1 before because blinding fails there
            kwargs["opts"]   = {"xmin": xMin, "xmax": xMax, "ymin": 1e0, "ymax": 5e9}
            kwargs["moveLegend"] = _legNE
            kwargs["blindingRangeString"] = "0-100"
            #kwargs["moveBlindedText"]     = {"dx": 0.0, "dy": 0.04, "dh": 0.0}
        else:
            kwargs["opts"]   = {"ymin": 1e0, "ymaxfactor": 10}

    if "DR" in h:
        kwargs["rebinX"] = 4

    if "Mass" in h:
        units = "GeV" #"GeV/c^{2}"
        kwargs["ylabel"] = _yLabel + units

        if "TransverseMass_" in h:
            kwargs["rebinX"] = 5
            kwargs["opts"]   = {"xmax": 1000.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
            kwargs["xlabel"] = "m_{T} (%s)" % units
        elif "TauTau" in h:
            kwargs["rebinX"] = 5
            kwargs["opts"]   = {"xmax": 600.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
            kwargs["xlabel"] = "m_{vis} (%s)" % units
        elif "TopMassWMassRatio" in h:
            units            = ""
            kwargs["rebinX"] = 2
            kwargs["ylabel"] = "Events / %.2f"
            kwargs["opts"]   = {"xmin": 0.0, "xmax": +5.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
            kwargs["cutBox"] = {"cutValue": (173.21/80.385), "fillColor": 16, "box": False, "line": True, "greaterThan": True} 
            #if "AfterStandardSelections" in h:
            #    kwargs["opts"]   = {"xmin": 0.0, "xmax": +10.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        elif "TopMass" in h:
            kwargs["rebinX"] = 10 #systematics._dataDrivenCtrlPlotBinning["LdgTrijetMass_AfterAllSelections"] #2
            kwargs["xlabel"] = "m_{jjb} (%s)" % units
            kwargs["cutBox"] = {"cutValue": 173.21, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
            kwargs["opts"]   = {"xmin": 0.0, "xmax": 1000.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
            # kwargs["rebinX"] = 4
        elif "dijetmass" in h.lower():
            kwargs["rebinX"] = 4
            #kwargs["rebinX"] = systematics._dataDrivenCtrlPlotBinning["LdgTrijetDijetMass_AfterAllSelections"] #2
            kwargs["xlabel"] = "m_{jj} (%s)" % units
            kwargs["opts"]   = {"xmax": +300.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
            kwargs["cutBox"] = {"cutValue": 80.385, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        else:
            kwargs["rebinX"] = 4
            kwargs["opts"]   = {"xmax": 400.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
            units            = "GeV"
            kwargs["xlabel"] = "m_{jjb} (%s)" % units
            kwargs["ylabel"] = _yLabel + units
            #Print("Unexpected histogram \"%s\"." % h, True)

    if "resolution" in h.lower():
        kwargs["moveLegend"] = _legCE
        if "eta" in h.lower():
            kwargs["opts"] = {"xmin": -1.2, "xmax": 1.2, "ymin": _yMin, "ymaxfactor": _yMaxF}
        if "pt" in h.lower() and "tau" in opts.folder.lower():
            kwargs["opts"] = {"xmin": -1.0, "xmax":1.0, "ymin": _yMin, "ymaxfactor": _yMaxF}


    if "eta" in h.lower():
        kwargs["rebinX"] = 2
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["opts"]   = {"xmin": -2.5, "xmax": +2.5, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["moveLegend"] = _legRM #_legNE

    if "_AfterAllSelections" in h or "_AfterTopSelection" in h or "_AfterMetSelection" in h:
        kwargs["blindingRangeString"] = "0-%s" % (5000)
        kwargs["ratio"] = False

    if kwargs["divideByBinWidth"]:
        kwargs["ylabel"] = "<" + ["ylabel"] + ">"
    
    if  opts.folder == "ForDataDrivenCtrlPlotsEWKFakeB" or opts.folder == "ForDataDrivenCtrlPlotsEWKGenuineB":
        kwargs["ratio"] = False

    if "resolution" in h.lower() or h == "btagSF" or "counters" in opts.folder:
        kwargs["ratio"]      = False

    return kwargs
    

def GetBinWidthMinMax(binList):
    if not isinstance(binList, list):
        raise Exception("Argument is not a list instance!")

    minWidth = +1e6
    maxWidth = -1e6
    # For-loop: All bin values (centre)
    for i in range(0, len(binList)-1):
        j = i + 1
        iBin = binList[i]
        jBin = binList[j]
        wBin = jBin-iBin
        if wBin < minWidth:
            minWidth = wBin

        if wBin > maxWidth:
            maxWidth = wBin
    return minWidth, maxWidth

def getHistos(datasetsMgr, histoName):

    h1 = datasetsMgr.getDataset("Data").getDatasetRootHisto(histoName)
    h1.setName("Data")

    h2 = datasetsMgr.getDataset("EWK").getDatasetRootHisto(histoName)
    h2.setName("EWK")
    return [h1, h2]

def DataMCHistograms(datasetsMgr, histoName):
    Verbose("Plotting Data-MC Histograms")

    # Skip 2-D plots
    skipStrings = []
    if opts.folder == "topbdtSelection_":
        skipStrings = ["_Vs_", "Vs", "Matched", "MCtruth", "TopQuark", 
                       "RealSelected", "DeltaMVAgt1", "SelectedTop", 
                       "LdgTrijetFake", "LdgTrijetFakeJJB", "TrijetFake",
                       "FakeInTopDir", "LdgTrijetFakeJJB_BDT", "LdgTrijetFake_BDT"]

    if opts.folder == "counters":
        skipStrings = ["weighted"]
    if opts.folder == "eSelection_Veto":
        skipStrings = ["Resolution"]
    if opts.folder == "muSelection":
        skipStrings = ["Resolution"]
    if opts.folder == "tauSelection_":
        skipStrings = ["riggerMatch", "NprongsMatrix", "Resolution"]
    if opts.folder == "PUDependency":
        skipStrings = ["WithProbabilisticBtag", "AngularCuts", "AntiIsolatedTau", "NvtxTau", "METSelection", "NvtxAllSelections"] #latter is empty!
    if opts.folder == "jetSelection_":
        skipStrings = ["JetMatching", "SeventhJet"]
    if opts.folder == "bjetSelection_":
        skipStrings = ["MatchDeltaR", "btagSFRelUncert", "_dEta", "_dPhi", "_dPt", "_dR"]
    if opts.folder == "metSelection_":
        skipStrings = [""]
    if opts.folder == "topologySelection_":
        skipStrings = ["_Vs_"]
    if opts.folder == "topSelectionBDT_":
        skipStrings = ["RelUncert"]

    if "ForDataDrivenCtrlPlots" in opts.folder:
        skipStrings = ["_Vs_", "JetEtaPhi", "MinDeltaPhiJet", "MaxDeltaPhiJet", "MinDeltaRJet"]

    # Skip histograms if they contain a given string
    for keyword in skipStrings:
        if keyword in histoName:
            return

    # Get Histogram name and its kwargs
    saveName = histoName.rsplit("/")[-1]
    kwargs_  = GetHistoKwargs(saveName, opts)

    # Create the plotting object
    p = plots.DataMCPlot(datasetsMgr, histoName, saveFormats=[])

    # Overwite signal style?
    if 0:
        if opts.signalMass != 0:
            p.histoMgr.forHisto(opts.signal, styles.getSignalStyleHToTB_M(opts.signalMass))

    # Apply blinding of signal region
    if "blindingRangeString" in kwargs_:
        startBlind = float(kwargs_["blindingRangeString"].split("-")[1])
        endBlind   = float(kwargs_["blindingRangeString"].split("-")[0])
        plots.partiallyBlind(p, maxShownValue=startBlind, minShownValue=endBlind, invert=True, moveBlindedText=kwargs_["moveBlindedText"])

    # Draw and save the plot
    plots.drawPlot(p, saveName, **kwargs_) #the "**" unpacks the kwargs_ dictionary

    # Replace bin labels
    if "counter" in opts.folder:
        p.getFrame().GetXaxis().LabelsOption("v") #vertical orientation of bin labels
        if histoName.split("/")[-1] == "counter":
            p.getFrame().GetXaxis().LabelsOption("v") #vertical orientation of bin labels
            replaceBinLabels(p, saveName)
        pass

    # Save the plots in custom list of saveFormats
    SavePlot(p, saveName, os.path.join(opts.saveDir, opts.optMode, opts.folder), opts.saveFormats)
    return

def replaceBinLabels(p, histoName):
    '''
    https://root.cern.ch/doc/master/classTAttText.html#T5
    '''
    myBinList = []
    if histoName == "counter" or histoName == "weighted/counter":
        #myBinList = ["trigger", "filter", "pv", "e veto", "1 #mu", "2 #tau jets", "genuine #tau",
        #             "#tau N", "#tau OS", "#tau SF", "fake #tau", "#geq 3 jets", "#geq 1 b jets",
        #             "b jets SF", "1 top" ]
        #myBinList = ["trigger", "filter", "pv", "e veto", "1 #mu", "2 #tau jets", "genuine #tau",
        #             "#tau N", "#tau OS", "#tau SF", "fake #tau", "#geq 3 jets", "#geq 1 b jets",
        #             "b jets SF", "E_{T}^{miss}","1 top", "top SF", "selected"]
        #myBinList = ["trigger", "filter", "pv", "e veto", "1 #mu", "#geq 1 #tau jets", "genuine #tau",
        #             "2 #tau jets", "#tau OS", "#tau SF", "fake #tau SF", "#geq 3 jets", "R_{coll}^{min}", "#geq 1 b jets",
        #             "b jets SF", "E_{T}^{miss}", "R_{bb}^{min}", "1 top", "selected"] # "top SF" = "selected"
        myBinList = ["trigger", "filter", "PV", "e veto", "1 #mu", "#geq 1 #tau jets", "genuine #tau",
                     "2 #tau jets", "#tau OS", "#tau SF", "fake #tau SF", "#geq 3 jets", "#geq 1 b jets",
                     "b jets SF", "E_{T}^{miss}", "1 top", "top SF", "selected"] # "top SF" = "selected"

    elif "bjet" in histoName:
        myBinList = ["All", "#eta", "p_{T}", "CSVv2 (M)", "Trg Match", "#geq 3"]
    elif "jet" in histoName:
        myBinList = ["All", "jet ID", "PU ID", "#tau match", "#eta", "p_{T}", "#geq 7", "H_{T}", "J_{T}", "MHT"]
    else:
        pass
    for i in range(0, len(myBinList)):
        p.getFrame().GetXaxis().SetBinLabel(i+1, myBinList[i])
    return

def SavePlot(plot, plotName, saveDir, saveFormats = [".C", ".png", ".pdf"]):
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    # Create the name under which plot will be saved
    saveName = os.path.join(saveDir, plotName.replace("/", "_"))
    saveName = saveName.replace(" ", "_")
    saveName = saveName.replace(")", "")
    saveName = saveName.replace("(", "")

    # For-loop: All save formats
    for i, ext in enumerate(saveFormats):
        saveNameURL = saveName + ext
        saveNameURL = aux.convertToURL(saveNameURL, opts.url)
        Verbose(saveNameURL, i==0)
        plot.saveAs(saveName, formats=saveFormats)
    return

#================================================================================================ 
# Main
#================================================================================================ 
if __name__ == "__main__":
    '''
    https://docs.python.org/3/library/argparse.html
 
    name or flags...: Either a name or a list of option strings, e.g. foo or -f, --foo.
    action..........: The basic type of action to be taken when this argument is encountered at the command line.
    nargs...........: The number of command-line arguments that should be consumed.
    const...........: A constant value required by some action and nargs selections.
    default.........: The value produced if the argument is absent from the command line.
    type............: The type to which the command-line argument should be converted.
    choices.........: A container of the allowable values for the argument.
    required........: Whether or not the command-line option may be omitted (optionals only).
    help............: A brief description of what the argument does.
    metavar.........: A name for the argument in usage messages.
    dest............: The name of the attribute to be added to the object returned by parse_args().
    '''
    
    # Default Settings
    ANALYSISNAME = "Hplus2hwAnalysisWithTop"
    SEARCHMODE   = "80to1000"
    DATAERA      = "Run2016"
    GRIDX        = False
    GRIDY        = False
    OPTMODE      = None
    BATCHMODE    = True
    INTLUMI      = -1.0
    SIGNALMASS   = 500
    MERGEEWK     = False
    URL          = False
    SAVEDIR      = None
    VERBOSE      = False
    RATIO        = False
    HISTOLEVEL   = "Vital" # 'Vital' , 'Informative' , 'Debug' 
    FOLDER       = "ForDataDrivenCtrlPlots" # "topSelectionBDT_" #"ForDataDrivenCtrlPlots" #jetSelection_
    SAVEFORMATS  = "pdf,png,C"
    
    # Define the available script options
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("-m", "--mcrab", dest="mcrab", action="store", 
                      help="Path to the multicrab directory for input")

    parser.add_option("-o", "--optMode", dest="optMode", type="string", default=OPTMODE, 
                      help="The optimization mode when analysis variation is enabled  [default: %s]" % OPTMODE)

    parser.add_option("-b", "--batchMode", dest="batchMode", action="store_false", default=BATCHMODE, 
                      help="Enables batch mode (canvas creation does NOT generate a window) [default: %s]" % BATCHMODE)

    parser.add_option("--analysisName", dest="analysisName", type="string", default=ANALYSISNAME,
                      help="Override default analysisName [default: %s]" % ANALYSISNAME)

    parser.add_option("--ratio", dest="ratio", action="store_true", default=RATIO,
                      help="Enable ratio pad for Data/Bkg comparison [default: %s]" % RATIO)

    parser.add_option("--intLumi", dest="intLumi", type=float, default=INTLUMI,
                      help="Override the integrated lumi [default: %s]" % INTLUMI)

    parser.add_option("--searchMode", dest="searchMode", type="string", default=SEARCHMODE,
                      help="Override default searchMode [default: %s]" % SEARCHMODE)

    parser.add_option("--dataEra", dest="dataEra", type="string", default=DATAERA, 
                      help="Override default dataEra [default: %s]" % DATAERA)

    parser.add_option("--mergeEWK", dest="mergeEWK", action="store_true", default=MERGEEWK, 
                      help="Merge all EWK samples into a single sample called \"EWK\" [default: %s]" % MERGEEWK)

    parser.add_option("--gridX", dest="gridX", action="store_true", default=GRIDX, 
                      help="Enable the x-axis grid lines [default: %s]" % GRIDX)

    parser.add_option("--gridY", dest="gridY", action="store_true", default=GRIDY, 
                      help="Enable the y-axis grid lines [default: %s]" % GRIDY)

    parser.add_option("--signalMass", dest="signalMass", type=int, default=SIGNALMASS, 
                     help="Mass value of signal to use [default: %s]" % SIGNALMASS)

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR, 
                      help="Directory where all pltos will be saved [default: %s]" % SAVEDIR)

    parser.add_option("--url", dest="url", action="store_true", default=URL, 
                      help="Don't print the actual save path the histogram is saved, but print the URL instead [default: %s]" % URL)
    
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=VERBOSE, 
                      help="Enables verbose mode (for debugging purposes) [default: %s]" % VERBOSE)

    parser.add_option("-i", "--includeOnlyTasks", dest="includeOnlyTasks", action="store", 
                      help="List of datasets in mcrab to include")

    parser.add_option("-e", "--excludeTasks", dest="excludeTasks", action="store", 
                      help="List of datasets in mcrab to exclude")

    parser.add_option("--folder", dest="folder", type="string", default = FOLDER,
                      help="ROOT file folder under which all histograms to be plotted are located [default: %s]" % (FOLDER) )

    parser.add_option("-s", "--saveFormats", dest="saveFormats", default = SAVEFORMATS,
                      help="Save formats for all plots [default: %s]" % SAVEFORMATS)

    (opts, parseArgs) = parser.parse_args()

    # Require at least two arguments (script-name, path to multicrab)
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    if opts.mcrab == None:
        Print("Not enough arguments passed to script execution. Printing docstring & EXIT.")
        parser.print_help()
        #print __doc__
        sys.exit(1)

    if opts.saveDir == None:
        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="", postfix="DataMC")

    # Sanity check
    allowedMass = [180, 200, 220, 250, 300, 350, 400, 500, 650, 800, 1000, 1500, 2000, 3000]
    if opts.signalMass!=0 and opts.signalMass not in allowedMass:
        Print("Invalid signal mass point (=%.0f) selected! Please select one of the following:" % (opts.signalMass), True)
        for m in allowedMass:
            Print(m, False)
        sys.exit()
    else:
        opts.signal = "ChargedHiggs_HplusTB_HplusToTB_M_%.0f" % opts.signalMass

    # Create save formats
    if "," in opts.saveFormats:
        opts.saveFormats = opts.saveFormats.split(",")
    else:
        opts.saveFormats = [opts.saveFormats]
    opts.saveFormats = ["." + s for s in opts.saveFormats]          

    # Call the main function
    main(opts)

    if not opts.batchMode:
        raw_input("=== plotDataMC.py: Press any key to quit ROOT ...")
