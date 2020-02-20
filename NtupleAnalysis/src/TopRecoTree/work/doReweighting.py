#!/usr/bin/env python
'''
DESCRIPTION:
Script for plotting TH2 histograms only.


USAGE:
./doReweighting.py -m <pseudo-mcrab> [--options]


EXAMPLES:
./doReweighting.py -m TopRecoTree_multFloat_AnalysisSelections --normalizeToOne --url

LAST USED:
./doReweighting.py -m TopRecoTree_multFloat_AnalysisSelections --normalizeToOne --ref TrijetPtDR --url
./doReweighting.py -m TopRecoTree_200217_114252_FirstNewBranches_hadronic --normalizeToOne --ref trijetPt --url --treeName treeS --logZ
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
import math
import getpass
import socket
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
import HiggsAnalysis.NtupleAnalysis.tools.crosssection as xsect
import HiggsAnalysis.NtupleAnalysis.tools.aux as aux
import HiggsAnalysis.NtupleAnalysis.tools.multicrabConsistencyCheck as consistencyCheck
import HiggsAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles


#import warnings
#warnings.filterwarnings("ignore")

#ROOT.gErrorIgnoreLevel = ROOT.kError

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

    optModes = [""]

    if opts.optMode != None:
        optModes = [opts.optMode]
        
    # Apply TDR style
    style = tdrstyle.TDRStyle()
    style.setOptStat(False)
    style.setGridX(opts.gridX)
    style.setGridY(opts.gridY)
    style.setLogX(opts.logX)
    style.setLogY(opts.logY)
    style.setLogZ(opts.logZ)


    # For-loop: All opt Mode
    for opt in optModes:
        opts.optMode = opt

        # Setup & configure the dataset manager 
        datasetsMgr = GetDatasetsFromDir(opts)
        datasetsMgr.updateNAllEventsToPUWeighted()
        datasetsMgr.loadLuminosities() # from lumi.json
        if opts.verbose:
            datasetsMgr.PrintCrossSections()
            datasetsMgr.PrintLuminosities()
            datasetsMgr.PrintInfo()

        # Merge histograms (see NtupleAnalysis/python/tools/plots.py) 
        plots.mergeRenameReorderForDataMC(datasetsMgr) 
        
        # Print merged datasets and MC samples
        if 0:
            datasetsMgr.PrintInfo()

        # Get Luminosity
        if opts.intLumi < 0:
            if "Data" in datasetsMgr.getAllDatasetNames():
                opts.intLumi = datasetsMgr.getDataset("Data").getLuminosity()
            else:
                opts.intLumi = 1.0 # Why?
         
        # Set/Overwrite cross-sections. Remove all but 1 signal mass
        removeList = []
        for d in datasetsMgr.getAllDatasets():
            if "ChargedHiggs" in d.getName():
                datasetsMgr.getDataset(d.getName()).setCrossSection(1.0)
                if d.getName() != opts.signal:
                    removeList.append(d.getName())

        # Custom filtering of datasets
        for i, d in enumerate(removeList, 1):
            msg = "Removing datasets %s from dataset manager" % (ShellStyles.NoteStyle() + d + ShellStyles.NormalStyle())
            Verbose(msg, i==1)
            datasetsMgr.remove(filter(lambda name: d == name, datasetsMgr.getAllDatasetNames()))


        # Print dataset information
        if 1:
            datasetsMgr.PrintInfo()
                    
        # Get reference branch (all the variables will be plotted as a function of the reference) 
        ref = opts.refHisto
        doReweighting(datasetsMgr, ref, opts.treeName)
    Print("All plots saved under directory %s" % (ShellStyles.NoteStyle() + aux.convertToURL(opts.saveDir, opts.url) + ShellStyles.NormalStyle()), True)    
    return


def GetHistoKwargs(h, opts):

    # Defaults
    xMin   = 0
    xMax   = 800
    xBins  = 10

    zMin   =   0
    zMax   = None
    yLabel = "y-axis"

    if opts.logX:
        xMin =  1
    if opts.logY:
        yMin =  1

    if opts.normalizeToLumi:
        yLabel  = "Events"
        yMin    = 1e0
        yMaxF   = 1
    elif opts.normalizeByCrossSection:
        yLabel  = "#sigma (pb)"
        #zMin    = 0 #1e-3
        yMaxF   = 1
    elif opts.normalizeToOne:
        yLabel  = "Arbitrary Units"
        yMin    = None #5e-4
        yMax    = None #1e-1
        yMaxF   = 1.1
    else:
        yLabel = "Unknown"

    xLabel = "xlabel"

    cutBox      = {"cutValue": 400.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True} #box = True works
    cutBoxY     = {"cutValue": 200.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
                       
    _format = "%0.0f "

    kwargs = {
        "stackMCHistograms": False,
        "addMCUncertainty" : False,
        "addLuminosityText": opts.normalizeToLumi,
        "addCmsText"       : True,
        "cmsExtraText"     : " Preliminary",
        "xmin"             : xMin,
        "xmax"             : xMax,
        "xbins"            : xBins,
        "ymin"             : yMin,        
        "ymax"             : yMax,
        "cutBox"           : cutBox,
        "cutBoxY"          : cutBoxY,
        #"moveLegend"       : {"dx": -2.0, "dy": 0.0, "dh": 0}, #hack to remove legend (tmp)
        "createLegend"     : {"x1": 0.58, "y1": 0.65, "x2": 0.92, "y2": 0.82},
        "xlabel"           : xLabel,
        "ylabel"           : yLabel,
        "opts"             : {"ymin": yMin, "ymaxfactor": yMaxF},
        "format"           : _format
        }
    
    if "mass" in h.lower():
        _units  = "GeV/c^{2}"
        kwargs['format'] = "%0.0f " + _units
        kwargs["xlabel"] = "M (%s)" % _units

    if "pt" in h.lower():
        _units  = "GeV/c"
        kwargs['format'] = "%0.0f " + _units
        kwargs["xlabel"] = "p_{T} (%s)" % _units
        if "trijet" in h.lower():
            kwargs["xlabel"] = "p_{T,t} (%s)" % _units
        if "dijet" in h.lower():
            kwargs["xlabel"] = "p_{T,w} (%s)" % _units
        if "ldgjetpt" in h.lower():
            kwargs["xmax"] = 200
            kwargs["xlabel"] = "leading jet p_{T} (%s)" % _units

    if "trijetmass" in h.lower():
        _units  = "GeV/c^{2}"
        kwargs['format'] = "%0.0f " + _units
        kwargs["xlabel"] = "m_{top} (%s)" % _units
        kwargs["xmax"] = 805 #1005

    if "jet_mass" in h.lower():
        kwargs["xmax"] = 750

    if "ldgjetmass" in h.lower():
        kwargs["xlabel"] = "m_{leading jet} (%s)" % _units
        kwargs["xmax"] = 50
        if "subldg" in h.lower():
            kwargs["xlabel"] = "m_{subleading jet} (%s)" % _units

    if "bjetmass" in h.lower():
        kwargs["xlabel"] = "m_{b-tagged jet} (%s)" % _units
        kwargs["xmax"] = 50

    if "bjetldgjet_mass" in h.lower():
        kwargs["xlabel"]= "m_{b, ldg jet} (%s)" % _units
        kwargs["xmax"] = 705
    if "bjetsubldgjet_mass" in h.lower():
        kwargs["xlabel"] = "m_{b-tagged, subldg jet} (%s)" % _units
        kwargs["xmax"] = 705

    if "bjetldgjetmass" in h.lower():
        kwargs["xlabel"]= "m_{b, ldg jet} (%s)" % _units
        kwargs["xmax"] = 705
    if "bjetsubldgjetmass" in h.lower():
        kwargs["xlabel"] = "m_{b-tagged, subldg jet} (%s)" % _units
        kwargs["xmax"] = 705

    if "mult" in h.lower():
        kwargs['format'] = "%0.0f"        
        kwargs["xmax"] = 50
        kwargs["xlabel"] = "mult"
        if "ldg" in h.lower():
            kwargs["xlabel"] = "Leading jet mult"
        if "subldg" in h.lower():
            kwargs["xlabel"] = "Subleading jet mult"
        if "avg" in h.lower():
            kwargs["xlabel"] = "avg CvsL"

    if "cvsl" in h.lower():
        kwargs['format'] = "%0.2f"
        kwargs["xmax"] = 1
        kwargs["xlabel"] = "CvsL"
        if "ldg" in h.lower():
            kwargs["xlabel"] = "Leading jet CvsL"
        if "subldg" in h.lower():
             kwargs["xlabel"] = "Subleading jet CvsL"
        if "avg" in h.lower():
            kwargs["xlabel"] = "avg CvsL"

    if "dijetmass" in h.lower():
        _units  = "GeV/c^{2}"
        kwargs['format'] = "%0.0f " + _units
        kwargs["xlabel"] = "m_{W} (%s)" % _units
        kwargs["xmax"] = 600

    if "bdisc" in h.lower():
        kwargs['format'] = "%0.2f"
        kwargs["xlabel"] = "CSV"
        kwargs["xmax"] = 1
        if "subldg" in h.lower():
            kwargs["xlabel"] = "Subleading jet CSV"
        elif "ldg" in h.lower():
            _xlabel = "Leading jet CSV"
        else:
            _xlabel = "b-tagged jet CSV"

    if "over" in h.lower():
         kwargs['format'] = "%0.2f "
         kwargs["xlabel"] = "m_{W}/m_{t}"
         kwargs["xmax"] = 1
         kwargs["xmin"] = 0

    if "eta" in h.lower():
         kwargs['format'] = "%0.2f "
         kwargs["xlabel"] = "#eta"
         kwargs["xmax"] = 5
         kwargs["xmin"] = -5         
         if "DEta" in h:
             kwargs["xlabel"] = "#Delta#eta"
             kwargs["xmin"] = 0

    if "phi" in h.lower():
         kwargs['format'] = "%0.2f "
         kwargs["xlabel"] = "#phi"
         kwargs["xmax"] = 3.5
         kwargs["xmin"] = -3.5         
         if "DPhi" in h:
             kwargs["xlabel"] = "#Delta#phi"
             #kwargs["xmin"] = 0

    if "DR" in h:
         kwargs['format'] = "%0.2f "
         kwargs["xlabel"] = "#Delta R"
         kwargs["xmax"] = 6
         kwargs["xmin"] = 0

    if "ptd" in h.lower():
        kwargs['format'] = "%0.2f"
        kwargs["xlabel"] = "p_{T}D"
        kwargs["xmax"] = 1
        if "ldg" in h.lower():
            kwargs["xlabel"] = "Leading jet p_{T}D"
        elif "subldg" in h.lower():
            kwargs["xlabel"]= "Subleading jet p_{T}D"

    if "ptdr" in h.lower():
        kwargs["xmax"] =800
        kwargs['format'] = "%0.0f"
        kwargs["xlabel"] = "p_{T}#Delta R_{t}"        
        if "dijetptdr" in h.lower():
            kwargs["xlabel"] = "p_{T}#Delta R_{W}"

    if "axis2" in h.lower():
        kwargs['format'] = "%0.3f"
        kwargs["xmax"] = 0.3
        kwargs["xmin"] = 0
        kwargs["xlabel"] = "axis2"
        if "ldg" in h.lower():
            kwargs["xlabel"] = "Leading jet axis2"
        if "subldg" in h.lower():
            kwargs["xlabel"] = "Subleading jet axis2"
        if "avg" in h.lower():
            kwargs["xlabel"] = "avg axis2"

    if "likelihood" in h.lower():
        kwargs['format'] = "%0.2f"
        kwargs["xmax"] = 1
        kwargs["xlabel"] = "QGL"
        if "ldg" in h.lower():
            kwargs["xlabel"] = "Leading jet QGL"
        if "subldg" in h.lower():
            kwargs["xlabel"] = "Subleading jet QGL"
        if "avg" in h.lower():
            kwargs["xlabel"] = "avg QGL"

    if "softdrop" in h.lower():
        kwargs["xlabel"] = "Soft Drop n_{2}"
        kwargs["xmax"] = 2

    if "chisquared" in h.lower():
        kwargs["xlabel"] = "#chi^{2}"
        kwargs["xmax"] = 10


    # Get Bins
    binWidth = 10
    Xmax = [0.001, 0.01, 0.1, 1., 10., 100., 1000.]

    for x in Xmax:
        if kwargs["xmax"] <= x:
            binWidth = x/100.
            break

    if (int(kwargs["xmax"] - kwargs["xmin"])/binWidth < 10):
            binWidth = binWidth/10.

    kwargs["xbins"] = int((kwargs["xmax"] - kwargs["xmin"])/binWidth)
    #

    if 0:
        ROOT.gStyle.SetNdivisions(8, "X")
        ROOT.gStyle.SetNdivisions(8, "Y")
        
    return kwargs
    
def CreateHistograms(inputList, kwargs, txt):
    histoMap = {}
    Print("List of branches:", True)

    for i, brName in enumerate(inputList):

        kwargs   = GetHistoKwargs(brName, opts)

        Print("%s. %s" % (i, brName), False)
        
        xbins = kwargs["xbins"]
        xmin  = kwargs["xmin"]
        xmax  = kwargs["xmax"]

        hName = "%s_%s" % (brName, txt)
        # Create histogram
        histo = ROOT.TH1F(hName, hName, xbins, xmin, xmax)
        histo.GetXaxis().SetTitle(kwargs["xlabel"])
        histo.GetYaxis().SetTitle(kwargs["ylabel"])

        histoMap.update({brName : histo})
        
    return histoMap


def GetListOfBranches(tree):
    branchList = []
    brsList = tree.GetListOfBranches()
    nBranches = len(brsList)
    for i in range(nBranches):
        if "Weight" in brsList[i].GetName():
            continue
        # Branches with old names
        if "Trijet" not in brsList[i].GetName():
            continue
        if "over" not in brsList[i].GetName().lower():
            continue
        branchList.append(brsList[i].GetName())        
    return branchList
def GetWeights(sigTree, bkgTree, refBranch):
    kwargs   = GetHistoKwargs(refBranch, opts)    
            
    xmin  = kwargs["xmin"]
    xmax  = kwargs["xmax"]
    xbins = 2*kwargs["xbins"] #(xmax - xmin)/10 #kwargs["xbins"]

    # Create histogram
    sigName = "%s_S" % refBranch
    bkgName = "%s_B" % refBranch
    Print("Creating reference histograms %s %s" % (sigName, bkgName), True)

    hSig = ROOT.TH1F(sigName, sigName, xbins, xmin, xmax)
    hBkg = ROOT.TH1F(bkgName, bkgName, xbins, xmin, xmax)

    # Number of entries in trees
    nEntries_s = sigTree.GetEntries();
    nEntries_b = nEntries_s #bkgTree.GetEntries();
       
    if (nEntries_b == nEntries_s):
        Print("Warning! using nEntries_b = nEntries_s", True)

    # loop over signal entries and fill hSig
    for ientry in range(nEntries_s):    
        sigTree.GetEntry(ientry)
        xvalue = getattr(sigTree, refBranch)
        hSig.Fill(xvalue)

    # loop over background entries and fill hBkg
    for ientry in range(nEntries_b):    
        bkgTree.GetEntry(ientry)
        xvalue = getattr(bkgTree, refBranch)
        hBkg.Fill(xvalue)

    # Normalize histograms to unity
    hSig.Scale(1./hSig.Integral())
    hBkg.Scale(1./hBkg.Integral())

    # Get weights
    hWeights = hSig.Clone("weights")
    hWeights.Divide(hBkg)
        
    # 
    hSig.Scale(hWeights.Integral())
    hBkg.Scale(hWeights.Integral())
    
    hSig.SetFillColorAlpha(ROOT.kBlue, 0.2)
    hSig.SetLineColorAlpha(ROOT.kBlue, 0.2)

    hBkg.SetFillColorAlpha(ROOT.kRed, 0.2)
    hBkg.SetLineColorAlpha(ROOT.kRed, 0.2)

    hWeights.SetLineWidth(3)
    kwargs["ymax"] = hWeights.GetMaximum() + 1
    
    binwidth = (xmax - xmin)/xbins
    kwargs["ylabel"] = "Weights (a.u.) / %s GeV" % str(binwidth)

    kwargs["opts"]   = {"ymin": 0.0, "ymaxfactor": 1.1}

    p = plots.ComparisonManyPlot(histograms.Histo(hWeights,"weights", "l", "L"), 
                                 [histograms.Histo(hBkg,"background", "l", "L"),
                                  histograms.Histo(hSig,"signal", "l", "L")], saveFormats=[])

    p.histoMgr.setHistoLegendLabelMany({"weights" : "Signal/Bkg", "signal" : "Signal", "background" : "Bkg"})
    p.histoMgr.setHistoLegendStyle("weights", "L")
    p.histoMgr.setHistoLegendStyle("signal", "F")
    p.histoMgr.setHistoLegendStyle("background", "F")
    kwargs["createLegend"] = {"x1": 0.58, "y1": 0.65, "x2": 0.92, "y2": 0.82}

    # Save plot in all formats    
    saveName = "weights_%s" % (refBranch)
    plots.drawPlot(p, saveName, **kwargs)
    SavePlot(p, saveName, os.path.join(opts.saveDir), opts.saveFormats )

    return hWeights

def PlotOverlay(s,b, kw, saveName):
    
    #Define Signal style
    signalStyle     = styles.StyleCompound([styles.StyleMarker(markerSize=1.2, markerColor=ROOT.kAzure-3, markerSizes=None, markerStyle=ROOT.kFullDiamond),
                                            styles.StyleLine(lineColor=ROOT.kAzure-3, lineStyle=ROOT.kSolid, lineWidth=3),
                                            styles.StyleFill(fillColor=ROOT.kAzure-3)])
    #Define Background style
    backgroundStyle = styles.StyleCompound([styles.StyleMarker(markerSize=1.2, markerColor=ROOT.kRed-4, markerSizes=None, markerStyle=ROOT.kFullDiamond),
                                            styles.StyleLine(lineColor=ROOT.kRed-4, lineStyle=ROOT.kSolid, lineWidth=3)])
                                            #styles.StyleFill(fillColor=ROOT.kRed+2, fillStyle=3001)])
    
    p = plots.ComparisonPlot(histograms.Histo(b,"background", "pl", "PL"),histograms.Histo(s,"signal", "pl", "PL"),) 
    dName = opts.dataset
    dName = dName.replace("TT", "t#bar{t}")
    if "ChargedHiggs" in dName:
        dName = dName.replace("ChargedHiggs_HplusTB_HplusToTB_M_", "m_{H^{^#pm}} = ")
        dName = dName+"GeV"
        
    p.histoMgr.setHistoLegendLabelMany({"signal": "Truth-matched", "background": "Unmatched"}) 
    p.histoMgr.forHisto("background", backgroundStyle ) 
    p.histoMgr.setHistoDrawStyle("background", "HIST")
    p.histoMgr.setHistoLegendStyle("background", "L") #F
    
    p.histoMgr.forHisto("signal", signalStyle) 
    p.histoMgr.setHistoDrawStyle("signal" , "HIST")
    p.histoMgr.setHistoLegendStyle("signal", "L") #L
    
    
    _kwargs = {
        "xlabel"           : kw['xlabel'],
        "ylabel"           : "Arbitrary Units / %s" % (kw['format']),
        "ratioYlabel"      : "Ratio ",
        "ratio"            : False,
        "ratioInvert"      : False,
        "stackMCHistograms": False,
        "addMCUncertainty" : False,
        "addLuminosityText": False,
        "addCmsText"       : True,
        "cmsExtraText"     : "Preliminary",
        "xmin"             : kw['xmin'],
        "xmax"             : kw['xmax'],
        #"opts"             : {"ymin": 0.0, "ymaxfactor": 1.1},
        "opts"             : {"ymin": 0.0, "ymaxfactor": 1.1},
        "opts2"            : {"ymin": 0.6, "ymax": 1.5},
        "log"              : False,
        "createLegend"     : {"x1": 0.58, "y1": 0.65, "x2": 0.92, "y2": 0.82}, 
        }

    # Save plot in all formats    
    #saveName = "%s_%s" % (brName, txt)
    plots.drawPlot(p, saveName, **_kwargs)
    SavePlot(p, saveName, os.path.join(opts.saveDir), opts.saveFormats )
    
    return
    #return p, _kwargs

def doReweighting(datasetsMgr, refBranch, treeName):
    # https://gitlab.cern.ch/HPlus/HiggsAnalysis/blob/master/NtupleAnalysis/src/Keras_ANN/work/func.py#L49
    # Get weights
    #hWeights = GetWeights(ref_sig, ref_bkg) # Do it with branches !
    
    # Get ROOT file
    f = ROOT.TFile(opts.rootFileName)

    # Get signal and background trees
    sigTree = f.Get("treeS")
    bkgTree = f.Get("treeB")

    # Get weights in "refBranch" bins
    hWeights = GetWeights(sigTree, bkgTree, refBranch)
    
    # Number of entries in trees
    nEntries_s = sigTree.GetEntries();
    nEntries_b = nEntries_s #bkgTree.GetEntries();

    if (nEntries_b == nEntries_s):
        Print("Warning! using nEntries_b = nEntries_s", True)

    # GetListOfBranches
    inputList = GetListOfBranches(sigTree)

    kwargs = {}
    # Create map with branch names and corresponding TH1 histograms
    histoMap_s = CreateHistograms(inputList, kwargs, "sig")
    histoMap_b = CreateHistograms(inputList, kwargs, "bkg")
    
    # loop over signal entries
    for ientry in range(nEntries_s):
        sigTree.GetEntry(ientry)
        # loop over all branches
        for brName in inputList:
            # Get branch value
            # https://www.programiz.com/python-programming/methods/built-in/getattr
            xvalue = getattr(sigTree, brName)
            histoMap_s[brName].Fill(xvalue)

    # loop over bkg entries
    for ientry in range(nEntries_b):
        bkgTree.GetEntry(ientry)
        refValue = getattr(bkgTree, refBranch)
        bin      = hWeights.FindBin(refValue)
        weight   = hWeights.GetBinContent(bin)
        # loop over all branches
        for brName in inputList:
            # Get branch value
            # https://www.programiz.com/python-programming/methods/built-in/getattr
            xvalue   = getattr(bkgTree, brName)
            # Fill weighted bkg histogram
            histoMap_b[brName].Fill(xvalue, weight)

    txt = "Rew"
    if (refBranch == "trijetMass"):
        txt = "topMassRew"
    if (refBranch == "trijetPt"):
        txt = "topPtRew"

    #Define Signal style
    signalStyle     = styles.StyleCompound([styles.StyleMarker(markerSize=1.2, markerColor=ROOT.kAzure-3, markerSizes=None, markerStyle=ROOT.kFullDiamond),
                                            styles.StyleLine(lineColor=ROOT.kAzure-3, lineStyle=ROOT.kSolid, lineWidth=3),
                                            styles.StyleFill(fillColor=ROOT.kAzure-3)])
    #Define Background style
    backgroundStyle = styles.StyleCompound([styles.StyleMarker(markerSize=1.2, markerColor=ROOT.kRed-4, markerSizes=None, markerStyle=ROOT.kFullDiamond),
                                            styles.StyleLine(lineColor=ROOT.kRed-4, lineStyle=ROOT.kSolid, lineWidth=3)])
                                            #styles.StyleFill(fillColor=ROOT.kRed+2, fillStyle=3001)])
    for brName in inputList:
        s = histoMap_s[brName]
        b = histoMap_b[brName]
        # Normalize histograms to unity
        s.Scale(1./s.Integral())
        b.Scale(1./b.Integral())

        # plot signal and background on same canvas        
        #p, _kwargs = PlotOverlay(s,b, GetHistoKwargs(brName, opts))
        PlotOverlay(s,b, GetHistoKwargs(brName, opts), "%s_%s" % (brName, txt))
            
    return

def SavePlot(plot, plotName, saveDir, saveFormats = [".C", ".png", ".pdf"]):
    Verbose("Saving the plot in %s formats: %s" % (len(saveFormats), ", ".join(saveFormats) ) )

    # Check that path exists
    #if not os.path.exists(saveDir):
    #    os.makedirs(saveDir)
        
    # Create the name under which plot will be saved
    saveName = os.path.join(saveDir, plotName.replace("/", "_"))
    #print "soti", saveName
    # For-loop: All save formats
    for i, ext in enumerate(saveFormats):
        saveNameURL = saveName + ext
        saveNameURL = aux.convertToURL(saveNameURL, opts.url)
        Print(saveNameURL, i==0)
        plot.saveAs(saveName, formats=saveFormats)
    return

def SaveInRootFile(inputList, histoMap_s, histoMap_b):

#Save reweighted histograms under input rootFile

    # Open input root file
    outFileName = "%s.root" % opts.mcrab
    outputFile = ROOT.TFile.Open(outFileName,"RECREATE")
    outputFile.cd()
    
    # Check if path exists and create
    outDirName_sig = "reweighted_%s_Genuine" % opts.refHisto
    outDirName_bkg = "reweighted_%s_Fake" % opts.refHisto
    
    nFolders = len(outputFile.GetListOfKeys())
    makeDir = True
    for out in [outDirName_sig, outDirName_bkg]:
        # For-loop: All folders  in ROOT file
        for k, key in enumerate(outputFile.GetListOfKeys(), 1):
            kName = key.GetName()
            if (kName == out):
                makeDir = False
        if (makeDir):
            outDir = outputFile.mkdir(out)
            Print("Directory %s has been created." % out)
            outDir.cd()
        else:
            outDir = outputFile.Get(out)
            outDir.cd()

        for brName in inputList:
            if "Genuine" in out:
                h = histoMap_s[brName]
                h.Write()
            elif "Fake" in out:
                h = histoMap_b[brName]                
                h.Write()
    outputFile.Write()
    outputFile.Close()
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
    ANALYSISNAME = "TopRecoTree"
    SEARCHMODE   = "80to1000"
    DATAERA      = "Run2016"
    OPTMODE      = None
    BATCHMODE    = True
    GRIDX        = False
    GRIDY        = False
    LOGX         = False
    LOGY         = False
    LOGZ         = False
    MERGEEWK     = False
    URL          = False
    SAVEDIR      = None
    VERBOSE      = False
    FOLDER       = "ForDataDrivenCtrlPlots" #"topbdtSelection_" #jetSelection_
    SIGNALMASS   = 800
    SIGNAL       = None
    INTLUMI      = -1.0
    NORM2ONE     = False
    NORM2XSEC    = False
    NORM2LUMI    = False
    REFHISTO     = "trijetMass"
    TREENAME     = "treeS"
    DATASET      = "TT"
    ROOTFILENAME = None
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

    parser.add_option("--searchMode", dest="searchMode", type="string", default=SEARCHMODE,
                      help="Override default searchMode [default: %s]" % SEARCHMODE)

    parser.add_option("--dataEra", dest="dataEra", type="string", default=DATAERA, 
                      help="Override default dataEra [default: %s]" % DATAERA)

    parser.add_option("--signalMass", dest="signalMass", type=int, default=SIGNALMASS,
                     help="Mass value of signal to use [default: %s]" % SIGNALMASS)

    parser.add_option("--gridX", dest="gridX", action="store_true", default=GRIDX, 
                      help="Enable the x-axis grid lines [default: %s]" % GRIDX)

    parser.add_option("--gridY", dest="gridY", action="store_true", default=GRIDY, 
                      help="Enable the y-axis grid lines [default: %s]" % GRIDY)

    parser.add_option("--logX", dest="logX", action="store_true", default=LOGX, 
                      help="Set x-axis to logarithm scale [default: %s]" % LOGX)

    parser.add_option("--logY", dest="logY", action="store_true", default=LOGY,
                      help="Set y-axis to logarithm scale [default: %s]" % LOGY)

    parser.add_option("--logZ", dest="logZ", action="store_true", default=LOGZ,
                      help="Set z-axis to logarithm scale [default: %s]" % LOGZ)

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

    parser.add_option("--mergeEWK", dest="mergeEWK", action="store_true", default = MERGEEWK,
                      help="Merge EWK datasets into a single dataset? [default: %s]" % (MERGEEWK) )

    parser.add_option("--intLumi", dest="intLumi", type=float, default=INTLUMI,
                      help="Override the integrated lumi [default: %s]" % INTLUMI)

    parser.add_option("--normalizeToOne", dest="normalizeToOne", action="store_true", default=NORM2ONE,
                      help="Normalise plot to one [default: %s]" % NORM2ONE)

    parser.add_option("--normalizeByCrossSection", dest="normalizeByCrossSection", action="store_true", default=NORM2XSEC,
                      help="Normalise plot by cross-section [default: %s]" % NORM2XSEC)

    parser.add_option("--normalizeToLumi", dest="normalizeToLumi", action="store_true", default=NORM2LUMI,
                      help="Normalise plot to luminosity [default: %s]" % NORM2LUMI)

    parser.add_option("--refHisto", dest="refHisto", type="string", default=REFHISTO, 
                      help="Plot all the histograms as a function of the reference histogram [default: %s]" % REFHISTO)

    parser.add_option("--treeName", dest="treeName", type="string", default=TREENAME, 
                      help="Plot all the histograms as a function of the reference histogram [default: %s]" % TREENAME)

    parser.add_option("--dataset", dest="dataset", type="string", default=DATASET,
                      help="Dataset to draw (only 1 allowed in 2D) [default: %s]" % (DATASET) )

    parser.add_option("--rootFileName", dest="rootFileName", type="string", default=ROOTFILENAME,
                      help="ROOT file folder under which all trees and branches to be plotted are located [default: %s]" % (ROOTFILENAME) )

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

    # Create save formats
    if "," in opts.saveFormats:
        opts.saveFormats = opts.saveFormats.split(",")
    else:
        opts.saveFormats = [opts.saveFormats]
    opts.saveFormats = ["." + s for s in opts.saveFormats]          
    

    if opts.saveDir == None:
        #opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="", postfix="reweighted_%s" % opts.refHisto)
        user    = getpass.getuser()
        initial = getpass.getuser()[0]
        opts.saveDir = "/publicweb/%s/%s/%s/reweighted_%s/" % (initial, user, opts.mcrab, opts.refHisto)
    if not os.path.exists(opts.saveDir):
        os.makedirs(opts.saveDir)

    # Get root file name
    if (opts.rootFileName == None):
            opts.rootFileName = "%s/%s/res/histograms-%s.root" % (opts.mcrab, opts.dataset, opts.dataset)

    if opts.normalizeToOne == False and opts.normalizeByCrossSection == False and opts.normalizeToLumi == False:
        raise Exception("One of the options --normalizeToOne, --normalizeByCrossSection, --normalizeToLumi must be enabled (set to \"True\").")


    # Sanity check
    allowedFolders = ["counters", "counters/weighted", "PUDependency", "Weighting", 
                      "eSelection_Veto", "muSelection_Veto", "tauSelection_Veto",
                      "ForDataDrivenCtrlPlotsEWKFakeB", "ForDataDrivenCtrlPlotsEWKGenuineB",
                      "jetSelection_", "bjetSelection_", "metSelection_", 
                      "topologySelection_", "topbdtSelection_", "ForDataDrivenCtrlPlots", "TrijetCandidate"]

    if opts.folder not in allowedFolders:
        Print("Invalid folder \"%s\"! Please select one of the following:" % (opts.folder), True)
        for m in allowedFolders:
            Print(m, False)
        sys.exit()

     # Sanity check
    allowedMass = [180, 200, 220, 250, 300, 350, 400, 500, 800, 1000, 1500, 2000, 2500, 3000, 5000, 7000]
    if opts.signalMass!=0 and opts.signalMass not in allowedMass:
        Print("Invalid signal mass point (=%.0f) selected! Please select one of the following:" % (opts.signalMass), True)
        for m in allowedMass:
            Print(m, False)
        sys.exit()
    else:
        opts.signal = "ChargedHiggs_HplusTB_HplusToTB_M_%i" % opts.signalMass

    # Call the main function
    main(opts)

    if not opts.batchMode:
        raw_input("=== plotTH1.py: Press any key to quit ROOT ...")
