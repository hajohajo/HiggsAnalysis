#!/usr/bin/env python
'''
DESCRIPTION:
Script for comparing exlusion limits from same channel but different cuts or measurements, 
or even for comparing exclusion limits of different channels.


USAGE:
./plotOutputs.py  [opts]


EXAMPLES:
./plotOutputs.py -s png --plotType significance --yMin 0.0 --yMaxFactor 1.1 --dirs new5,new6,new7,new8 --saveDir /publicweb/a/aattikis/new5_6_7_8
./plotOutputs.py -s png --plotType efficiency --yMin 0.0 --yMax 1.0 --dirs new1,new2,new3 --refIndex 3 --saveDir /publicweb/a/aattikis/Test
./plotOutputs.py -s png --plotType significance --yMin 0.0 --yMaxFactor 1.1 --refIndex 3 --dirs new1,new2,new3 --saveDir /publicweb/a/aattikis/Test 
./plotOutputs.py -s png --logY --plotType output --dirs new1,new2,new3,new4 --refIndex 8 --saveDir /publicweb/a/aattikis/Test && ./plotOutputs.py -s png --plotType efficiency --yMin 0.0 --yMax 1.0 --dirs new1,new2,new3,new4 --refIndex 8 --saveDir /publicweb/a/aattikis/Test && ./plotOutputs.py -s png --plotType significance --yMin 0.0 --yMaxFactor 1.1 --refIndex 8 --dirs new1,new2,new3,new4 --saveDir /publicweb/a/aattikis/Test
./plotOutputs.py -s png --logY --plotType output --dirs new1,new8,new9 --refIndex 0 --saveDir /publicweb/a/aattikis/Test && ./plotOutputs.py -s png --plotType efficiency --dirs new1,new8,new9 --refIndex 8 --saveDir /publicweb/a/aattikis/Test && ./plotOutputs.py -s png --plotType significance --yMin 0.0 --yMaxFactor 1.1 --refIndex 0 --dirs new1,new8,new9 --saveDir /publicweb/a/aattikis/Test && ./plotOutputs.py -s png --plotType roc --logY --dirs new1,new8,new9 --saveDir /publicweb/a/aattikis/Test
./plotOutputs.py -s png --plotType roc --logY -saveDir /publicweb/a/aattikis/Test --logX --xMin 0.6 --yMin 6e-3 --yMax 1e-1 --dirs new2,new8,new9,new10,new11
./plotOutputs.py -s png --plotType roc --xMin 0.01 --logY --dirs new9,new10,new11 --saveDir /publicweb/a/aattikis/BDT --refName BDT --cutLineX 0.93
./plotOutputs.py -s png --plotType roc --logY --saveDir /publicweb/a/aattikis/Test  --dirs 500k_sample,Keras_BDTG2018_Aug2018 --refName BDT --cutLineX 0.93
./plotOutputs.py -s png --logY --plotType output --dirs new10,new11 --refIndex 0 --saveDir /publicweb/a/aattikis/Test && ./plotOutputs.py -s png --plotType efficiency --dirs new10,new11 --refIndex 0 --saveDir /publicweb/a/aattikis/Test && ./plotOutputs.py -s png --plotType significance --yMin 0.0 --yMaxFactor 1.1 --refIndex 0 --dirs new10,new11 --saveDir /publicweb/a/aattikis/Test && ./plotOutputs.py -s png --plotType roc --logY --dirs new10,new11 --saveDir /publicweb/a/aattikis/Test


LAST USED:
./plotOutputs.py -s png --plotType roc --logY --saveDir /publicweb/a/aattikis/Test --dirs 500k_sample,Keras_BDTG2018_Aug2018 --cutLineX 0.93 --xMin 0.0
./plotOutputs.py -s png --plotType var --logY --saveDir /publicweb/a/aattikis/SotiSoti --dirs Keras_Test --yMin 0.5e-3 --yMax 1.0
./plotOutputs.py -s pdf,C,png --plotType var --saveDir /publicweb/a/aattikis/DNN --dirs Keras_3Layers_19relu_190relu_1sigmoid_500Epochs_5000BatchSize_21-Nov-2019_07h21m34s --yMin 0.5e-3 --yMaxFactor 2.0 --boldText
./plotOutputs.py --plotType roc --logY --saveDir /publicweb/a/aattikis/Keras --cutLineX 0.93 --xMin 0.0 --dirs 500k_sample,Keras_BDTG2018_Aug2018 --refName BDTG

./plotOutputs.py --plotType significance --logY --saveDir /publicweb/a/aattikis/Keras --cutLineX 0.93 --xMin 0.0 --dirs 500k_sample,Keras_BDTG2018_Aug2018 --refName BDT && ./plotOutputs.py --plotType roc --logY --saveDir /publicweb/a/aattikis/Keras --cutLineX 0.93 --xMin 0.0 --dirs 500k_sample,Keras_BDTG2018_Aug2018 --refName BDT && ./plotOutputs.py --plotType efficiency --logY --saveDir /publicweb/a/aattikis/Keras --cutLineX 0.93 --xMin 0.0 --dirs 500k_sample,Keras_BDTG2018_Aug2018 --refName BDT

'''
#================================================================================================
# Import modules
#================================================================================================
import os
import getpass
import sys
import glob
import json
import array
import copy
from optparse import OptionParser

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

import HiggsAnalysis.NtupleAnalysis.tools.histograms as histograms
import HiggsAnalysis.NtupleAnalysis.tools.tdrstyle as tdrstyle
import HiggsAnalysis.NtupleAnalysis.tools.plots as plots
import HiggsAnalysis.NtupleAnalysis.tools.styles as styles
import HiggsAnalysis.Keras_ANN.results as _results
import HiggsAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles
import HiggsAnalysis.NtupleAnalysis.tools.aux as aux

#================================================================================================
# Shell Types
#================================================================================================
sh_e = ShellStyles.ErrorStyle()
sh_s = ShellStyles.SuccessStyle()
sh_h = ShellStyles.HighlightStyle()
sh_a = ShellStyles.HighlightAltStyle()
sh_t = ShellStyles.NoteStyle()
sh_n = ShellStyles.NormalStyle()
sh_w = ShellStyles.WarningStyle()

#================================================================================================
# Function definition
#================================================================================================
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

def main():

    # Apply TDR style
    style = tdrstyle.TDRStyle()

    # Enable/Disable grids for axes
    style.setGridX(opts.gridX)
    style.setGridY(opts.gridY)

    # Enable/Disable logscale for axes 
    style.setLogX(opts.logX)
    style.setLogY(opts.logY)

    # Definitions
    savePath = opts.saveDir
    if opts.url:
        savePath = opts.saveDir.replace("/afs/cern.ch/user/a/attikis/public/html", "https://cmsdoc.cern.ch/~%s" % getpass.getuser())

    Verbose("Load module for reading the BR limits from JSON file produced by Combine", True)

    # Definitions
    resultsList = []
    for i,d in enumerate(opts.dirList, 1):
        
        dirs = glob.glob(d)
        dirs.sort()
        if len(dirs) == 0:
            raise Exception("No directories for pattern '%s'" % path)
        directory = dirs[-1]

        Verbose("Picked %s" % directory, True)
        resultsList.append(_results.Output(directory, excludePoints=[]))
        if opts.refName != None:
            if opts.refName in d:
                opts.refIndex = i-1 #iro

    # Do comparison plot
    msg  = "Creating comparison plots (%d) using the following results directories:%s\n\t%s" % (len(opts.dirList), sh_t, "\n\t".join([os.path.basename(d) for d in opts.dirList]) )
    Print(msg + sh_n, True)
    if "output" in opts.plotType.lower():
        doOutput(opts.saveName, resultsList) 
    elif "efficiency" in opts.plotType.lower():
        doEfficiency(opts.saveName, resultsList)
    elif "significance" in opts.plotType.lower():
        doSignificance(opts.saveName, resultsList) 
    elif "roc" in opts.plotType.lower():
        doROC(opts.saveName, resultsList)
    elif "var" in opts.plotType.lower():
        for var in ["TrijetPtDR", "TrijetDijetPtDR", "TrijetBjetMass", "TrijetLdgJetBDisc",
                    "TrijetSubldgJetBDisc", "TrijetBJetLdgJetMass", "TrijetBJetSubldgJetMass",
                    "TrijetMass", "TrijetDijetMass", "TrijetBJetBDisc","TrijetSoftDrop_n2",
                    "TrijetLdgJetCvsL", "TrijetSubldgJetCvsL", "TrijetLdgJetPtD", "TrijetSubldgJetPtD",
                    "TrijetLdgJetAxis2", "TrijetSubldgJetAxis2", "TrijetLdgJetMult","TrijetSubldgJetMult"]:
            doVariables(opts.saveName, var, resultsList)
    elif "metric" in opts.plotType.lower():
        for metric in ["TrainAccuracy", "ValAccuracy", "TrainLoss", "ValLoss"]:
            doMetrics(opts.saveName, metric, resultsList)
    else:
        pass
    
    # inform user of output location
    Print("Plots saved under directory %s"% (sh_s + aux.convertToURL(opts.saveDir, opts.url) + sh_n), True)
    return

def doOutput(name, resultsList):

    # Do the comparison plot
    Verbose("Creating the expected plots", True)
    gSigList  = []
    gBkgList  = []
    legList   = []

    # For-loop: All Output-class objects
    for r in resultsList:
        print r.getDirectory()
        gSig, lSig = r.getGraphs(opts.plotName + "_signal")
        gBkg, lBkg = r.getGraphs(opts.plotName + "_background")
        if opts.yMin == None:
            opts.yMin = r.getYMin()
        if opts.yMax == None and opts.yMax == None:
            opts.yMax = r.getYMax()*1.10
        gSigList.extend(gSig)
        gBkgList.extend(gBkg)
        legList.extend(lSig)
        #legList.extend(lBkg) # same as "lSig"

    # Re-arrange legend
    legList.insert(0, legList.pop(opts.refIndex))
    gSigList.insert(0, gSigList.pop(opts.refIndex))
    gBkgList.insert(0, gBkgList.pop(opts.refIndex))
    # Plot the graphs
    kwargs = GetKwargs(opts)
    doPlot(legList, gSigList, opts.saveName + "_Signal"    , **kwargs)
    doPlot(legList, gBkgList, opts.saveName + "_Background", **kwargs)

    # Do the relative plot
    kwargs = GetKwargsRatio(kwargs)
    gSigList = GetRelativeGraphs(gSigList, opts.refIndex)
    gBkgList = GetRelativeGraphs(gBkgList, opts.refIndex)
    doPlot(legList, gSigList, opts.saveName + "_SignalRel"    , **kwargs)
    doPlot(legList, gBkgList, opts.saveName + "_BackgroundRel", **kwargs)
    return

def doEfficiency(name, resultsList):

    # Do the comparison plot
    Verbose("Creating the expected plots", True)
    gSigList = []
    gBkgList = []
    legList  = []

    # For-loop: All Output-class objects
    for r in resultsList:
        gSig, lSig = r.getGraphs(opts.plotName + "Sig")
        gBkg, lBkg = r.getGraphs(opts.plotName + "Bkg")
        if opts.yMin == None:
            opts.yMin = r.getYMin()
        if opts.yMax == None and opts.yMax == None:
            opts.yMax = r.getYMax()*1.10
        gSigList.extend(gSig)
        gBkgList.extend(gBkg)
        legList.extend(lSig)
        
    # Re-arrange legend
    legList.insert(0, legList.pop(opts.refIndex))
    gSigList.insert(0, gSigList.pop(opts.refIndex))
    gBkgList.insert(0, gBkgList.pop(opts.refIndex))
    # Plot the graph
    kwargs = GetKwargs(opts)
    doPlot(legList, gSigList, opts.saveName + "_Sig", **kwargs)
    doPlot(legList, gBkgList, opts.saveName + "_Bkg", **kwargs)
        
    # Do the relative plot
    kwargs = GetKwargsRatio(kwargs)
    gSigList = GetRelativeGraphs(gSigList, opts.refIndex)
    gBkgList = GetRelativeGraphs(gBkgList, opts.refIndex)
    doPlot(legList, gSigList, opts.saveName + "_SigRel", **kwargs)
    doPlot(legList, gBkgList, opts.saveName + "_BkgRel", **kwargs)
    return

def doSignificance(name, resultsList):

    # Do the comparison plot
    Verbose("Creating the expected plots", True)
    gSigList = []
    gBkgList = []
    legList  = []

    # For-loop: All Output-class objects
    for r in resultsList:
        gSig, lSig = r.getGraphs(opts.plotName + "A")
        gBkg, lBkg = r.getGraphs(opts.plotName + "B")
        if opts.yMin == None:
            opts.yMin = r.getYMin()
        if opts.yMax == None and opts.yMax == None:
            opts.yMax = r.getYMax()*1.10
        gSigList.extend(gSig)
        gBkgList.extend(gBkg)
        legList.extend(lSig)
        
    # Re-arrange lists
    if opts.refIndex >= len(legList):
        opts.refIndex = opts.refIndex-1
    legList.insert(0, legList.pop(opts.refIndex))
    gSigList.insert(0, gSigList.pop(opts.refIndex))
    gBkgList.insert(0, gBkgList.pop(opts.refIndex))

    # Plot the graph
    kwargs = GetKwargs(opts)

    # Default definition
    kwargs["ylabel"] = "S/#sqrt{S+B}"
    doPlot(legList, gSigList, opts.saveName + "_A", **kwargs)

    # Alternative definition
    kwargs["ylabel"] = "2(#sqrt{S+B} - #sqrt{B})"
    doPlot(legList, gBkgList, opts.saveName + "_B", **kwargs)
        
    # Do the relative plot
    kwargs = GetKwargsRatio(kwargs)
    gSigList = GetRelativeGraphs(gSigList, opts.refIndex)
    gBkgList = GetRelativeGraphs(gBkgList, opts.refIndex)
    doPlot(legList, gSigList, opts.saveName + "_ARel", **kwargs)
    doPlot(legList, gBkgList, opts.saveName + "_BRel", **kwargs)
    return

def doROC(name, resultsList):
    
    # Do the comparison plot
    Verbose("Creating the expected plots", True)
    gSigList = []
    gBkgList = []
    legList  = []
    gSvBList = []

    # For-loop: All Output-class objects
    for r in resultsList:
        gSig, lSig = r.getGraphs("EfficiencySig")
        gBkg, lBkg = r.getGraphs("EfficiencyBkg")
        if opts.yMin == None:
            opts.yMin = r.getYMin()
        if opts.yMax == None and opts.yMax == None:
            opts.yMax = r.getYMax()*1.10
        gSigList.extend(gSig)
        gBkgList.extend(gBkg)
        if "BDT" in r.getDirectory():
            lSig = ["BDTG"]
        legList.extend(lSig)
        
    # Re-arrange legend
    legList.insert(0, legList.pop(opts.refIndex))
    gSigList.insert(0, gSigList.pop(opts.refIndex))
    gBkgList.insert(0, gBkgList.pop(opts.refIndex))

    # Plot the graph
    kwargs = GetKwargs(opts)
    # For-loop: All TGraphs
    for i, g in enumerate(gSigList, 0):
        x  = []
        y  = []
        N  = gSigList[i].GetN()
        nS = gSigList[i].GetName()
        nB = gBkgList[i].GetName()
        Verbose("nS = %s, nB = %s" % (nS, nB), True)

        # For-loop: All points in graph
        for j in range(0, N):
            xVal = gSigList[i].GetY()[j]
            yVal = gBkgList[i].GetY()[j]

            if xVal <= opts.xMin:
                continue

            x.append(xVal)
            y.append(yVal)
            Verbose("%s) x = %s, y = %s" % (g.GetName(), xVal, yVal), False)

        # Create the TGraph
        g = ROOT.TGraph(N, array.array('d', x), array.array('d',y))
        gSvBList.append(g)

    # Make the plot
    doPlot(legList, gSvBList, opts.saveName, **kwargs)
    return

def doVariables(name, variable, resultsList):
    
    # Do the comparison plot
    Verbose("Creating the expected plots", True)
    gSigList = []
    gBkgList = []
    legList  = []
    gSvBList = []
    opts.saveName = variable

    # For-loop: All Output-class objects
    for r in resultsList:
        gSig, lSig = r.getGraphs(variable + "_signal")
        gBkg, lBkg = r.getGraphs(variable + "_background")
        if opts.yMin == None:
            opts.yMin = r.getYMin()
        if opts.yMax == None and opts.yMax == None:
            opts.yMax = r.getYMax()*1.10
        gSigList.extend(gSig)
        gBkgList.extend(gBkg)
        legList.extend(["signal"]) 
        legList.extend(["background"])
        
    # Plot the graph
    kwargs = GetVarKwargs(variable, opts)
    gList = [gSigList[0], gBkgList[0]]

    # For-loop: All TGraphs
    #for i, g in enumerate(gSigList, 0):
    for i, g in enumerate(gList, 0):
        x  = []
        y  = []
        N  = g.GetN()
        nS = g.GetName()
        nB = g.GetName()
        Verbose("nS = %s, nB = %s" % (nS, nB), True)

        # For-loop: All points in graph
        for j in range(0, N):
            xVal = g.GetX()[j]
            yVal = g.GetY()[j]

            if xVal <= opts.xMin:
                continue

            # Protection: don't include (x,y) = (0.0, 0.0)
            if xVal !=0.0 and yVal !=0.0:
                x.append(xVal)
                y.append(yVal)
            Verbose("%d) x = %0.3f, y = %0.3f" % (j, xVal, yVal), j==0)

        # Create the TGraph
        Verbose("Creating the TGraph with len(x) = %d, len(x) = %d" % (len(x), len(y)), True)
        g = ROOT.TGraph(len(x), array.array('d', x), array.array('d',y))
        gSvBList.append(g)

    # Make the plot
    Verbose("Making the plot", True)
    doPlot(legList, gSvBList, opts.saveName, **kwargs)
    return

def doMetrics(name, metric, resultsList):
    
    # Do the comparison plot
    Verbose("Creating the expected plots", True)
    grList  = []
    legList = []
    opts.saveName = metric

    # For-loop: All Output-class objects
    for r in resultsList:
        gr, leg = r.getGraphs(metric)
        if opts.yMin == None:
            opts.yMin = r.getYMin()
        if opts.yMax == None and opts.yMax == None:
            opts.yMax = r.getYMax()*1.10
        grList.extend(gr)
        legList.extend(leg)
        
    # Plot the graph
    kwargs = GetMetricKwargs(metric, opts)
    gList  = []

    # For-loop: All TGraphs
    for i, g in enumerate(grList, 0):
        x  = []
        y  = []
        N  = g.GetN()

        # For-loop: All points in graph
        for j in range(0, N):
            xVal = g.GetX()[j]
            yVal = g.GetY()[j]

            if xVal <= opts.xMin:
                continue

            # Protection: don't include (x,y) = (0.0, 0.0)
            if xVal !=0.0 and yVal !=0.0:
                x.append(xVal)
                y.append(yVal)
            Verbose("%d) x = %0.3f, y = %0.3f" % (j, xVal, yVal), j==0)

        # Create the TGraph
        Verbose("Creating the TGraph with len(x) = %d, len(x) = %d" % (len(x), len(y)), True)
        g = ROOT.TGraph(len(x), array.array('d', x), array.array('d',y))
        gList.append(g)

    # Make the plot
    Verbose("Making the plot", True)
    doPlot(legList, gList, opts.saveName, **kwargs)
    return

def doPlot(legList, graphList, saveName, **kwargs):
    
    # Definitions
    hgList = []
    lList  = {}

    # For-loop: All TGraphs
    for i, g in enumerate(graphList, 0):
        if opts.boldText:
            gName = legList[i]
        else:
            gName = "#font[42]{%s}" % legList[i]
        hg = histograms.HistoGraph(graphList[i], gName, drawStyle="L", legendStyle="l")
        hgList.append(hg)

    # Create a plot-base object
    Verbose("Creating the plot-base object", True)
    # plot = plots.PlotBase(hgList, saveFormats=[])
    #hgList.insert(0, hgList.pop(opts.refIndex))
    plot = plots.ComparisonManyPlot(hgList[0], hgList[1:], saveFormats=[])
    #plot = plots.ComparisonManyPlot(hgList[-1], hgList[:-1], saveFormats=[])

    # Apply histo style
    Verbose("Applying the histogram styles (generator)", True)
    plot.histoMgr.forEachHisto(styles.generator())
    # plot.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetMarkerSize(1.2))
    # plot.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetLineWidth(3))

    def sty(h):
        r = h.getRootHisto()
        r.SetLineWidth(3)
        r.SetMarkerSize(1.2)
        return

    # Apply style and set label
    Verbose("Applying the histogram styles (forEachHisto)", True)
    plot.histoMgr.forEachHisto(sty)
    if opts.plotType != "var":
        plot.setLegendHeader("Sequential Model (Keras)")

    # Draw the plot
    Verbose("Drawing the plot", True)
    plots.drawPlot(plot, saveName, **kwargs)
    if opts.removeLegend:
        plot.removeLegend()

    # Save plots and return
    Verbose("Saving the plot as %s" % (saveName), True)
    SavePlot(plot, opts.saveDir, saveName, opts.saveFormats)
    
    Verbose("Plots saved under directory %s"% (sh_s + aux.convertToURL(opts.saveDir, opts.url) + sh_n), True)
    return


def GetKwargsRatio(kwargs):
    kwargs["ylabel"] = "Ratio"
    kwargs["opts"]["ymin"] = 0.0
    kwargs["opts"]["ymax"] = 3.0
    kwargs["log"] = False
    tdrstyle.TDRStyle().setLogY(False)
    return kwargs

def GetMetricKwargs(metric, opts):
    dh = -0.12 + (len(opts.dirList)-3)*0.04
    legNW  = {"dx": -0.53, "dy": -0.05, "dh": dh}
    legSW  = {"dx": -0.53, "dy": -0.30, "dh": dh}
    legNE  = {"dx": -0.53, "dy": -0.05, "dh": dh}
    legSE  = {"dx": -0.53, "dy": -0.45, "dh": dh}
    if opts.paper:
        #histograms.cmsTextMode = histograms.CMSMode.PAPER
        cmsExtraText = ""
    else:
        #histograms.cmsTextMode = histograms.CMSMode.PRELIMINARY
        #histograms.cmsTextMode = histograms.CMSMode.UNPUBLISHED
        cmsExtraText = "Simulation" #"Preliminary"

    kwargs  = {
        "xlabel"           : "DNN Output",
        "ylabel"           : "Entries",
        "addMCUncertainty" : False,
        "addLuminosityText": False,
        "addCmsText"       : True,
        "cmsExtraText"     : cmsExtraText,
        "cmsTextPosition"  : "outframe",
        "opts"             : {"xmin": opts.xMin, "xmax": opts.xMax},#, "ymin": opts.yMin, "ymax": opts.yMax},
        "log"              : opts.logY,
        "moveLegend"       : legNE,
        "xtitlesize"       : 0.1,
        "ytitlesize"       : 0.1,
        "cutBox"           : {"cutValue":  0.0, "fillColor": 16, "box": False, "line": False , "cutGreaterThan": False},
        #"cutBoxY"          : {"cutValue": 1000.0, "fillColor": 16, "box": False, "line": True, "cutGreaterThan": False}  # does not work!
        }

    # General settings
    if opts.yMaxFactor:
        kwargs["opts"]["ymaxfactor"] = opts.yMaxFactor
        #del kwargs["opts"]["ymax"]
    if opts.logY and kwargs["opts"]["ymin"] <= 0.0:
        msg = "Cannot have log-y enabled and y-min set to 0.0. Setting to y-min to 1.0."
        Print(sh_h + msg + sh_n, True)
        kwargs["opts"]["ymin"] = 1e0
    if opts.logX and kwargs["opts"]["xmin"] <= 0.0:
        msg = "Cannot have log-x enabled and y-min set to <=0.0. Setting to x-min to 1e-2"
        Print(sh_h + msg + sh_n, True)
        kwargs["opts"]["xmin"] = 1e-2
    if opts.cutLineX != None:
        kwargs["cutBox"]  = {"cutValue": opts.cutLineX, "fillColor": 16, "box": False, "line": True , "cutGreaterThan": False}
    if opts.cutLineY != None:
        Print("This does not work! Bug-fixing required!", True)
        kwargs["cutBoxY"] = {"cutValue": opts.cutLineY, "fillColor": 16, "box": True, "line": True, "cutGreaterThan": False} 

    if "loss" in metric.lower():
        #kwargs["opts"]["xmin"]   =  0.0
        #kwargs["opts"]["xmax"]   = 800.0
        kwargs["xlabel"] = "epoch"
        kwargs["yMin"]   = 0.0 #1e-3
        kwargs["ylabel"] = "loss"

    if "accuracy" in metric.lower():
        #kwargs["opts"]["xmin"]   =   0.0
        #kwargs["opts"]["xmax"]   = 800.0
        kwargs["xlabel"] = "epoch"
        kwargs["yMin"]   = 0.0 # 1e-3
        kwargs["ylabel"] = "accuracy"
        kwargs["moveLegend"] = legSE

    # No need for this since we plot TGraphs
    if opts.xMin == None:
        del kwargs["opts"]["xmin"]
    else:
        kwargs["opts"]["xmin"] = opts.xMin
    if opts.xMax == None:
        del kwargs["opts"]["xmax"]
    else:
        kwargs["opts"]["xmax"] = opts.xMax

    return kwargs

def GetVarKwargs(var, opts):
    dh = -0.12 + (len(opts.dirList)-3)*0.04
    legNW  = {"dx": -0.53, "dy": -0.05, "dh": dh}
    legSW  = {"dx": -0.53, "dy": -0.30, "dh": dh}
    legNE  = {"dx": -0.05, "dy": -0.05, "dh": dh}
    legSE  = {"dx": -0.18, "dy": -0.45, "dh": dh}
    if opts.paper:
        #histograms.cmsTextMode = histograms.CMSMode.PAPER
        cmsExtraText = ""
    else:
        #histograms.cmsTextMode = histograms.CMSMode.PRELIMINARY
        #histograms.cmsTextMode = histograms.CMSMode.UNPUBLISHED
        cmsExtraText = "Simulation" #"Preliminary"

    kwargs  = {
        "xlabel"           : "DNN Output",
        "ylabel"           : "Entries",
        "addMCUncertainty" : False,
        "addLuminosityText": False,
        "addCmsText"       : True,
        "cmsExtraText"     : cmsExtraText,
        "cmsTextPosition"  : "outframe",
        "opts"             : {"xmin": opts.xMin, "xmax": opts.xMax, "ymin": opts.yMin, "ymax": opts.yMax},
        "opts2"            : {"ymin": 0.59, "ymax": 1.41},
        "log"              : opts.logY,
        "moveLegend"       : legNE,
        "xtitlesize"       : 0.1,
        "ytitlesize"       : 0.1,
        "cutBox"           : {"cutValue":  0.0, "fillColor": 16, "box": False, "line": False , "cutGreaterThan": False},
        #"cutBoxY"          : {"cutValue": 1000.0, "fillColor": 16, "box": False, "line": True, "cutGreaterThan": False}  # does not work!
        }

    # General settings
    if opts.yMaxFactor:
        kwargs["opts"]["ymaxfactor"] = opts.yMaxFactor
        del kwargs["opts"]["ymax"]
    if opts.logY and kwargs["opts"]["ymin"] <= 0.0:
        msg = "Cannot have log-y enabled and y-min set to 0.0. Setting to y-min to 1.0."
        Print(sh_h + msg + sh_n, True)
        kwargs["opts"]["ymin"] = 1e0
    if opts.logX and kwargs["opts"]["xmin"] <= 0.0:
        msg = "Cannot have log-x enabled and y-min set to <=0.0. Setting to x-min to 1e-2"
        Print(sh_h + msg + sh_n, True)
        kwargs["opts"]["xmin"] = 1e-2
    if opts.cutLineX != None:
        kwargs["cutBox"]  = {"cutValue": opts.cutLineX, "fillColor": 16, "box": False, "line": True , "cutGreaterThan": False}
    if opts.cutLineY != None:
        Print("This does not work! Bug-fixing required!", True)
        kwargs["cutBoxY"] = {"cutValue": opts.cutLineY, "fillColor": 16, "box": True, "line": True, "cutGreaterThan": False} 

    if var == "TrijetPtDR":
        kwargs["opts"]["xmin"]   =  0.0
        kwargs["opts"]["xmax"]   = 800.0
        kwargs["xlabel"] = "p_{T}#DeltaR_{t}"
        kwargs["yMin"]   = 1e-3
        #kwargs["ylabel"] = "a.u. / %0.0f"

    if var == "TrijetDijetPtDR":
        kwargs["opts"]["xmin"]   =  0.0
        kwargs["opts"]["xmax"]   = 800.0
        kwargs["xlabel"] = "p_{T}#DeltaR_{W}"
        kwargs["yMin"]   = 1e-3
        #kwargs["ylabel"] = "a.u. / %0.f"

    if var == "TrijetBjetMass":
        kwargs["opts"]["xmin"]   =   0.0
        kwargs["opts"]["xmax"]   = 100.0
        kwargs["xlabel"] = "m_{b} [GeV]"
        #kwargs["ylabel"] = "a.u. / %0.f GeV"

    if var == "TrijetLdgJetBDisc":
        kwargs["opts"]["xmin"]   =   0.0
        kwargs["opts"]["xmax"]   =   1.0
        kwargs["xlabel"] = "b-tag discr."
        #kwargs["ylabel"] = "a.u. / %0.2f"

    if var == "TrijetSubldgJetBDisc":
        kwargs["opts"]["xmin"]   =   0.0
        kwargs["opts"]["xmax"]   =   1.0
        kwargs["xlabel"] = "b-tag discr."
        #kwargs["ylabel"] = "a.u. / %0.2f"

    if var == "TrijetBJetLdgJetMass":
        kwargs["opts"]["xmin"]   =  0.0
        kwargs["opts"]["xmax"]   = 600.0
        kwargs["xlabel"] = "m_{b,j_{1}} [GeV]"
        kwargs["yMin"]   = 1e-3

    if var == "TrijetBJetSubldgJetMass":
        kwargs["opts"]["xmin"]   =   0.0
        kwargs["opts"]["xmax"]   = 600.0
        kwargs["xlabel"] = "m_{b,j_{2}} [GeV]"
        kwargs["yMin"]   = 1e-3
        #kwargs["ylabel"] = "a.u. / %0.0f GeV"

    if var == "TrijetMass":
        kwargs["opts"]["xmin"]   =   0.0
        kwargs["opts"]["xmax"]   = 800.0
        kwargs["xlabel"] = "m_{W} [GeV]"
        kwargs["yMin"]   = 1e-3
        #kwargs["ylabel"] = "a.u. / %0.0f GeV"

    if var == "TrijetDijetMass":
        kwargs["opts"]["xmin"]   =   0.0
        kwargs["opts"]["xmax"]   = 500.0
        kwargs["xlabel"] = "m_{W} [GeV]"
        kwargs["yMin"]   = 1e-3
        #kwargs["ylabel"] = "a.u. / %0.0f GeV"

    if var == "TrijetBJetBDisc":
        kwargs["opts"]["xmin"]   =  0.0
        kwargs["opts"]["xmax"]   =  1.0
        kwargs["xlabel"] = "b-tag discr."
        #kwargs["ylabel"] = "a.u. / %0.2f"

    if var == "TrijetSoftDrop_n2":
        kwargs["opts"]["xmin"]   =   0.0
        kwargs["opts"]["xmax"]   =   2.0
        kwargs["xlabel"] = "soft-drop n_{2}"
        kwargs["yMin"]   = 1e-3
        #kwargs["ylabel"] = "a.u. / %0.2f"

    if var == "TrijetLdgJetCvsL":
        kwargs["opts"]["xmin"]   =  -1.0
        kwargs["opts"]["xmax"]   =  +1.0
        kwargs["xlabel"] = "CvsL discr."
        kwargs["yMin"]   = 1e-3
        kwargs["yMax"]   = 1e-1
        #kwargs["ylabel"] = "a.u. / %0.2f"

    if var == "TrijetSubldgJetCvsL":
        kwargs["opts"]["xmin"]   =  -1.0
        kwargs["opts"]["xmax"]   =  +1.0
        kwargs["xlabel"] = "CvsL discr."
        kwargs["yMin"]   = 1e-3
        kwargs["yMax"]   = 1e-1
        #kwargs["ylabel"] = "a.u. / %0.2f"

    if var == "TrijetLdgJetPtD":
        kwargs["opts"]["xmin"]   =   0.0
        kwargs["opts"]["xmax"]   =   1.0
        kwargs["xlabel"] = "p_{T}D"
        #kwargs["ylabel"] = "a.u. / %0.2f"

    if var == "TrijetSubldgJetPtD":
        kwargs["opts"]["xmin"]   =   0.0
        kwargs["opts"]["xmax"]   =   1.0
        kwargs["xlabel"] = "p_{T}D"
        #kwargs["ylabel"] = "a.u. / %0.2f"

    if var == "TrijetLdgJetAxis2":
        kwargs["opts"]["xmin"]   =   0.0
        kwargs["opts"]["xmax"]   =   0.16
        kwargs["xlabel"] = "axis2"
        #kwargs["ylabel"] = "a.u. / %0.3f"

    if var == "TrijetSubldgJetAxis2":
        kwargs["opts"]["xmin"]   =   0.0
        kwargs["opts"]["xmax"]   =   0.16
        kwargs["xlabel"] = "axis2"
        #kwargs["ylabel"] = "a.u. / %0.3f"

    if var == "TrijetLdgJetMult":
        kwargs["opts"]["xmin"]  =  0.0
        #kwargs["opts"]["xmax"]   = 40.0
        kwargs["xlabel"] = "constituent multiplicity"
        #kwargs["ylabel"] = "a.u. / %0.0f"

    if var == "TrijetSubldgJetMult":
        kwargs["opts"]["xmin"]   =  0.0
        kwargs["opts"]["xmax"]   = 40.0
        kwargs["xlabel"] = "constituent multiplicity"
        #kwargs["ylabel"] = "a.u. / %0.0f"

    # No need for this since we plot TGraphs
    if opts.xMin == None:
        del kwargs["opts"]["xmin"]
    else:
        kwargs["opts"]["xmin"] = opts.xMin
    if opts.xMax == None:
        del kwargs["opts"]["xmax"]
    else:
        kwargs["opts"]["xmax"] = opts.xMax

    return kwargs

def GetKwargs(opts):
    dh = -0.12 + (len(opts.dirList)-3)*0.04
    legNW  = {"dx": -0.53, "dy": -0.05, "dh": dh}
    legSW  = {"dx": -0.53, "dy": -0.30, "dh": dh}
    legNE  = {"dx": -0.05, "dy": -0.05, "dh": dh}
    legSE  = {"dx": -0.18, "dy": -0.45, "dh": dh}
    if opts.paper:
        #histograms.cmsTextMode = histograms.CMSMode.PAPER
        cmsExtraText = ""
    else:
        #histograms.cmsTextMode = histograms.CMSMode.PRELIMINARY
        #histograms.cmsTextMode = histograms.CMSMode.UNPUBLISHED
        cmsExtraText = "Preliminary"

    kwargs  = {
        "xlabel"           : "DNN Output",
        "ylabel"           : "Entries",
        "addMCUncertainty" : False,
        "addLuminosityText": False,
        "addCmsText"       : True,
        "cmsExtraText"     : cmsExtraText,
        "cmsTextPosition"  : "outframe",
        "opts"             : {"xmin": opts.xMin, "xmax": opts.xMax, "ymin": opts.yMin, "ymax": opts.yMax},
        "opts2"            : {"ymin": 0.59, "ymax": 1.41},
        "log"              : opts.logY,
        "moveLegend"       : legNW,
        "xtitlesize"       : 0.1,
        "ytitlesize"       : 0.1,
        "cutBox"           : {"cutValue":  0.0, "fillColor": 16, "box": False, "line": False , "cutGreaterThan": False},
        #"cutBoxY"          : {"cutValue": 1000.0, "fillColor": 16, "box": False, "line": True, "cutGreaterThan": False}  # does not work!
        }

    # General settings
    if opts.yMaxFactor:
        kwargs["opts"]["ymaxfactor"] = opts.yMaxFactor
        del kwargs["opts"]["ymax"]
    if opts.logY and kwargs["opts"]["ymin"] <= 0.0:
        msg = "Cannot have log-y enabled and y-min set to 0.0. Setting to y-min to 1.0."
        Print(sh_h + msg + sh_n, True)
        kwargs["opts"]["ymin"] = 1e0
    if opts.logX and kwargs["opts"]["xmin"] <= 0.0:
        msg = "Cannot have log-x enabled and y-min set to <=0.0. Setting to x-min to 1e-2"
        Print(sh_h + msg + sh_n, True)
        kwargs["opts"]["xmin"] = 1e-2
    if opts.cutLineX != None:
        kwargs["cutBox"]  = {"cutValue": opts.cutLineX, "fillColor": 16, "box": False, "line": True , "cutGreaterThan": False}
    if opts.cutLineY != None:
        Print("This does not work! Bug-fixing required!", True)
        kwargs["cutBoxY"] = {"cutValue": opts.cutLineY, "fillColor": 16, "box": True, "line": True, "cutGreaterThan": False} 

    # Plot-type specific
    if opts.plotType == "efficiency":
        kwargs["ylabel"] = "Efficiency"
    if opts.plotType == "significance":
        kwargs["ylabel"] = "Significance"
    if opts.plotType == "roc":
        kwargs["ylabel"] = "background efficiency"
        kwargs["xlabel"] = "signal efficiency"

    return kwargs

def GetRelativeGraphs(grList, refIndex=0):
    
    # Sanity
    if refIndex > len(grList)-1:
        raise Exception("Index for reference histogram (=%d) is out of range (max is %d)" % (refIndex, len(grList)-1))

    newList = []
    # For-loop: all graphs in list
    for i, g in enumerate(grList, 0):
        newList.append(_results.divideGraph(grList[i], grList[refIndex]))
        
    # Get the number of points
    grRef = grList[refIndex]
    nPoints = grRef.GetN()

    # For-loop: All Results
    for i in range(0, nPoints):
        grRef.SetPoint(i, grRef.GetX()[i], 1.0)
    newList[refIndex] = grRef
    return newList

def SavePlot(plot, saveDir, plotName, saveFormats = [".png", ".pdf"]):
    # Check that path exists
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    # Create the name under which plot will be saved
    saveName = os.path.join(saveDir, plotName.replace("/", "_"))

    # For-loop: All save formats
    for i, ext in enumerate(saveFormats, 1):
        saveNameURL = saveName + ext
        saveNameURL = saveNameURL.replace("/afs/cern.ch/user/a/attikis/public/html", "https://cmsdoc.cern.ch/~%s" % getpass.getuser())
        if opts.url:
            Verbose(saveNameURL, i==0)
        else:
            Verbose(saveName + ext, i==0)
        plot.saveAs(saveName, formats=saveFormats)
    return

if __name__ == "__main__":

    # Default options
    ANALYSISTYPE = "HToHW"
    BATCHMODE    = True    
    CFGJSON      = "config.json"
    CUTLINEX     = None
    CUTLINEY     = None
    REFINDEX     = 0
    REFNAME      = None
    EXCLUDE      = ""
    REMOVELEGEND = False
    GRIDX        = False
    GRIDY        = False
    RESULTSJSON  = "results.json"
    LOGX         = False
    LOGY         = False
    XMIN         = None
    XMAX         = None
    YMIN         = None
    YMAX         = None
    YMAXFACTOR   = None
    PAPER        = False
    SAVENAME     = None
    SAVEDIR      = None
    SAVEFORMATS  = "pdf,png" # pdf,C,png" # is .C format problematic? (hangs)
    URL          = False
    VERBOSE      = False
    BOLDTEXT     = False
    DIRS         = None
    PLOTTYPE     = "output"

    parser = OptionParser(usage="Usage: %prog [options]", add_help_option=True, conflict_handler="resolve")
    parser.add_option("-v", "--verbose", dest="verbose", default=VERBOSE, action="store_true",
                      help="Verbose mode for debugging purposes [default: %s]" % (VERBOSE) )
    
    parser.add_option("--boldText", dest="boldText", default=BOLDTEXT, action="store_true",
                      help="Use bold text printed on canvas? [default: %s]" % (BOLDTEXT))

    parser.add_option("--dirs", dest="dirs", default=DIRS,
                      help="List for Keras directories output (comma separated WITHOUT space) [default: %s]" % (DIRS))

    parser.add_option("--resultsJSON", dest="resultsJSON", default=RESULTSJSON,
                      help="JSON file containing the results [default: %s]" % (RESULTSJSON))
    
    parser.add_option("--cfgJSON", dest="cfgJSON", default=CFGJSON,
                      help="JSON file containing the configurations [default: %s]" % (CFGJSON))

    parser.add_option("--exclude", dest="exclude", default=EXCLUDE,
                      help="List for points to exclude (comma separated WITHOUT space) [default: %s]" % (EXCLUDE))
  
    parser.add_option("--analysisType", dest="analysisType", default=ANALYSISTYPE,
                      help="Flag to indicate the analysis type (e.g. \"HToTauNu\", \"HToTB\", \"HToHW\") [default: %s]" % (ANALYSISTYPE) )
 
    parser.add_option("--paper", dest="paper", default=PAPER, action="store_true",
                      help="Paper mode [default: %s]" % (PAPER) )
    
    parser.add_option("--saveName", dest="saveName", type="string", default=SAVENAME,
                      help="Name of the output plot [default = %s]" % (SAVENAME))

    parser.add_option("--url", dest="url", action="store_true", default=URL,
                      help="Don't print the actual save path the plots are saved, but print the URL instead [default: %s]" % URL)

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR,
                      help="Directory where all plots will be saved [default: %s]" % SAVEDIR)

    parser.add_option("--plotType", dest="plotType", type="string", default=PLOTTYPE,
                      help="Type of plot to produce (poutput, pred, test, train)  saved [default: %s]" % PLOTTYPE)

    parser.add_option("--logX", dest="logX", action="store_true", default=LOGX,
                      help="Plot x-axis (mass) as logarithmic [default: %s]" % (LOGX) )
    
    parser.add_option("--logY", dest="logY", action="store_true", default=LOGY,
                      help="Plot y-axis (exlusion limit) as logarithmic [default: %s]" % (LOGY) )
    
    parser.add_option("--removeLegend", dest="removeLegend", default=REMOVELEGEND, action="store_true",
                      help="Remove the legend from the canvas? (cleaner canvas) [default: %s]" % (REMOVELEGEND) )

    parser.add_option("--gridX", dest="gridX", default=GRIDX, action="store_true",
                      help="Enable the grid for the x-axis [default: %s]" % (GRIDX) )

    parser.add_option("--gridY", dest="gridY", default=GRIDY, action="store_true",
                      help="Enable the grid for the y-axis [default: %s]" % (GRIDY) )

    parser.add_option("--yMin", dest="yMin", default=YMIN, type="float",
                      help="Overwrite automaticly calculated minimum value of y-axis [default: %s]" % (YMIN) )
    
    parser.add_option("--yMax", dest="yMax", default=YMAX, type="float",
                      help="Overwrite automaticly calculated maximum value of y-axis [default: %s]" % (YMAX) )

    parser.add_option("--yMaxFactor", dest="yMaxFactor", default=YMAXFACTOR, type="float",
                      help="Overwrite automaticly calculated maximum value of y-axis [default: %s]" % (YMAXFACTOR) )

    parser.add_option("--xMin", dest="xMin", default=XMIN, type="float",
                      help="Overwrite minimum value of x-axis [default: %s]" % (XMIN) )
    
    parser.add_option("--xMax", dest="xMax", default=XMAX, type="float",
                      help="Overwrite maximum value of x-axis [default: %s]" % (XMAX) )

    parser.add_option("-s", "--saveFormats", dest="saveFormats", default = SAVEFORMATS,
                      help="Save formats for all plots [default: %s]" % SAVEFORMATS)

    parser.add_option("-b", "--batchMode", dest="batchMode", action="store_false", default=BATCHMODE,
                      help="Enables batch mode (canvas creation does NOT generate a window) [default: %s]" % BATCHMODE)

    parser.add_option("--cutLineX", dest="cutLineX", default=CUTLINEX, type=float,
                      help="Value for x-axis line to be drawn on canvas [default: %s]" % (CUTLINEX))

    parser.add_option("--cutLineY", dest="cutLineY", default=CUTLINEY, type=float,
                      help="Value for y-axis line to be drawn on canvas [default: %s]" % (CUTLINEY))

    parser.add_option("--refIndex", dest="refIndex", type="int", default=REFINDEX,
                      help="Index to use fo the reference TGraph when doing relative plots [default: %s]" % (REFINDEX) )

    parser.add_option("--refName", dest="refName", type="string", default=REFNAME,
                      help="String in TGraph directory to use as reference when doing relative plots [default: %s]" % (REFNAME) )

    (opts, args) = parser.parse_args()


    # Sanity checks
    opts.dirList = []
    if opts.dirs != None:
        opts.dirList = opts.dirs.split(",")
        if "Keras_" not in opts.dirList[0]:
            opts.dirList = []
            for d in opts.dirs.split(","):
                opts.dirList.extend(sorted([w[0] for w in os.walk(d) if "Keras_" in w[0]]))
            #print opts.dirList
    else:
        msg = "No directories provided!"
        raise Exception(sh_e + msg + sh_n)

    if len(opts.dirList) > 1:
        msg = "Directories considered:%s\n\t%s" % (sh_t, "\n\t".join(opts.dirList))
        Verbose(msg + sh_n, True)
    elif len(opts.dirList) ==  1:
        if opts.plotType != "var":
            msg = "At least 2 directories required. Only %d passed with --dirs argument!" % len(opts.dirList)
            raise Exception(sh_e + msg + sh_n)
    else:
        msg = "This should never be reached"
        raise Exception(sh_e + msg + sh_n)

    # Sanity check (analysis type)
    myAnalyses = ["HToTauNu", "HToTB", "HToHW"]
    if opts.analysisType not in myAnalyses:
        msg = "Invalid analysis type \"%s\". Please selected one of the following: \"%s" % (opts.analysisType, "\", \"".join(myAnalyses) + "\"")
        raise Exception(sh_e + msg + sh_n)
    else:
        msg = "Analysis type is %s" % (sh_t + opts.analysisType + sh_n)
        Print(msg, True)        

    # Sanity checks
    opts.excludePoints = []
    if len(opts.exclude) > 0:
        opts.excludePoints = opts.exclude.split(",")
    if len(opts.excludePoints) > 0:
        msg = "Excluding points \"%s\"" % (", ".join(opts.excludePoints))
        Print(sh_t + msg + sh_n, True)
    else:
        pass    

    # Sanity check
    for i, d in enumerate(opts.dirList, 0):
        if not os.path.isdir(d):
            msg = "Directory \"%s\" does not exist" % (d)
            raise Exception(sh_e + msg + sh_n)
        else:
            d2 = os.path.join(os.path.join(os.getcwd(), d))
            if os.path.isdir(d2):
                opts.dirList[i] = d2
            else:
                msg = "Directory \"%s\" does not exist" % (os.path.join(d))
                raise Exception(sh_e + msg + sh_n)

    for i, d in enumerate(opts.dirList, 0):
        Verbose("index = %d, dir = %s" % (i, d), i==0)
        
    # Check plot type
    plotTypes = ["output", "prediction", "training", "testing", "efficiency", "significance", "roc", "var", "metric"]
    plotNames = {}
    plotNames["output"]       = "Output"
    plotNames["prediction"]   = "OutputPred"
    plotNames["training"]     = "OutputTrain"
    plotNames["testing"]      = "OutputTest"
    plotNames["efficiency"]   = "Efficiency"
    plotNames["significance"] = "Significance"
    plotNames["roc"]          = "ROC"
    plotNames["var"]          = "TrijetLdgJetBDisc"
    plotNames["metric"]       = "metric"

    if opts.plotType.lower() not in plotTypes:
        msg = "Unsupported plot type  \"%s\". Please select from the following:" % (", ".join(plotTypes))
        raise Exception(sh_e + msg + sh_n)
    else:
        opts.plotName = plotNames[opts.plotType]
        if opts.plotType in ["roc", "efficiency", "significance"]:
            if opts.xMin == None:
                opts.xMin = 0.0
            if opts.xMax == None:
                opts.xMax = 1.0

    # Define directory name for saving output
    if opts.saveDir == None:
        #opts.saveDir = opts.dirList[0]
        opts.saveDir = aux.getSaveDirPath(opts.dirList[0], prefix="", postfix=opts.plotType)
    if opts.saveName == None:
        opts.saveName = opts.plotType

    # Create save formats
    if "," in opts.saveFormats:
        opts.saveFormats = opts.saveFormats.split(",")
    else:
        opts.saveFormats = [opts.saveFormats]
    opts.saveFormats = ["." + s for s in opts.saveFormats]

    # Call the main function
    main()
    
    if not opts.batchMode:
        raw_input("=== plotOutputs.py: Press any key to quit ROOT ...")
