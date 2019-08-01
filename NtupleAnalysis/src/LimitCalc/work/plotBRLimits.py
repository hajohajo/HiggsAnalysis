#!/usr/bin/env python
'''
DESCRIPTION:
Script for comparing exlusion limits from same channel but different cuts or measurements, 
or even for comparing exclusion limits of different channels.


USAGE:
./plotBRLimits.py  [opts]


EXAMPLES:
./plotBRLimits.py --logY --relative --yMax 10 --yMin 1e-2 --url --bandValue 5 --saveDir /publicweb/a/aattikis/Combine/ --analysisType "HToHW" --dirs datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_MediumTauID_NoBDTGm1p0_tauLdgTrkPt20_24July2019_autoMCStats,datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_MediumTauID_NoBDTGm1p0_23July2019_autoMCStats,datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_LooseTauID_NoBDTGm1p0_tauLdgTrkPt20_22July2019_autoMCStats,datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_LooseTauID_NoBDTGm1p0_19July2019_autoMCStats

./plotBRLimits.py -s png,pdf --cutLineX 400 --cutLineY 1 --logY --relative --yMax 10 --yMin 1e-2 --url --bandValue 5 --saveDir /publicweb/a/aattikis/Combine/ --analysisType "HToHW" --dirs datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_MediumTauID_NoBDTGm1p0_tauLdgTrkPt20_24July2019_autoMCStats,datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_MediumTauID_NoBDTGm1p0_23July2019_autoMCStats,datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_LooseTauID_NoBDTGm1p0_tauLdgTrkPt20_22July2019_autoMCStats,datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_LooseTauID_NoBDTGm1p0_19July2019_autoMCStats


LAST USED:
./plotBRLimits.py -s png --logY --relative --yMax 10 --yMin 1e-2 --url --saveDir /publicweb/a/aattikis/Combine/ --analysisType "HToHW" --dirs datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_LooseTauID_NoBDTGm1p0_19July2019_autoMCStats,datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_MediumTauID_NoBDTGm1p0_tauLdgTrkPt20_24July2019_autoMCStats,datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_LooseTauID_BDTG0p40_tauLdgTrkPt20_30July2019_autoMCStats,datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_MediumTauID_NoBDTGm1p0_23July2019_autoMCStats,datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_MediumTauID_NoBDTGm1p0_tauLdgTrkPt20_24July2019_autoMCStats,datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_MediumTauID_BDTG0p40_tauLdgTrkPt20_30July2019_autoMCStats,datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_TightTauID_NoBDTG0p40_tauLdgTrkPt20_31July2019_autoMCStats,datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_TightTauID_BDTG0p40_tauLdgTrkPt20_31July2019_autoMCStats 

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
import HiggsAnalysis.LimitCalc.limit as limit
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
# Global definitions
#================================================================================================
styleList = [styles.Style(24, ROOT.kBlack)] + styles.getStyles()


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
    style.setGridX(opts.gridX)
    style.setGridY(opts.gridY)
    # style.setGridX(opts.logX)
    # style.setGridY(opts.logY)

    # Enable OpenGL
    if opts.opaqueLegend:
        ROOT.gEnv.SetValue("OpenGL.CanvasPreferGL", 1)

    # Set text mode
    if opts.paper:
        histograms.cmsTextMode = histograms.CMSMode.PAPER
    else:
        histograms.cmsTextMode = histograms.CMSMode.UNPUBLISHED
        #histograms.cmsTextMode = histograms.CMSMode.PRELIMINARY

    # Definitions
    savePath = opts.saveDir
    if opts.url:
        savePath = opts.saveDir.replace("/afs/cern.ch/user/a/attikis/public/html", "https://cmsdoc.cern.ch/~%s" % getpass.getuser())


    Verbose("Load module for reading the BR limits from JSON file produced by Combine", True)
    opts.directory= os.getcwd() +"/datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_MediumTauID_NoBDTGm1p0_tauLdgTrkPt20_24July2019_autoMCStats/CombineResults_taujets_190729_134022"
    limits = limit.BRLimits(opts.directory, opts.analysisType, opts.excludeMassPoints, opts.limitsFile, opts.cfgFile)
    opts.mHplus       = limit.mHplus()
    opts.sigmaBRlimit = limits.getSigmaBRlimitText()
    opts.hardProcess  = limits.getHardProcessText()
    opts.finalState   = limits.getFinalStateText()
    opts.brAssumption = limits.getBRassumptionText()
    opts.intLumi      = limits.getLuminosity()

    # Create mapping of datacard path and label in a tuple
    compareList = []
    # For-loop: All datacard directories
    for i,d in enumerate(opts.dirList, 1):
        label = GetLabel(d)
        path  = os.path.join(d, "CombineResults*")
        compareList.append( (label, path) )

    # Do comparison plot
    msg  = "Creating 95%% upper limit comparison plots (%d) using the following datacard directories:%s\n\t%s" % (len(opts.dirList), sh_a, "\n\t".join([os.path.basename(d) for d in opts.dirList]) )
    Print(msg + sh_n, True)
    doCompare(opts.name, compareList) 

    # Make overlap plot (two limits on top of each other)
    oList = [compareList[0], compareList[1]]
    msg   = "Creating overlap plot using the first 2 datacard directories:%s\n\t%s" % (sh_t, "\n\t".join([os.path.basename(d) for d in opts.dirList[:2]]) )
    Print(msg + sh_n, True)
    doOverlap(opts.name, oList)

    # inform user of output location
    Print("Plots saved under directory %s"% (sh_s + aux.convertToURL(opts.saveDir, opts.url) + sh_n), True)
    return

def GetLabel(dirName):
    '''
    Predefined labels for datacard directories
    '''
    dirToLabelDict = {}
    dirToLabelDict["datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_LooseTauID_NoBDTGm1p0_19July2019_autoMCStats"] = "HPS-L"
    dirToLabelDict["datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_LooseTauID_NoBDTGm1p0_tauLdgTrkPt20_22July2019_autoMCStats"] = "HPS-L, p_{T}^{tk}>20"
    dirToLabelDict["datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_LooseTauID_BDTG0p40_tauLdgTrkPt20_30July2019_autoMCStats"] =  "HPS-L, p_{T}^{tk}>20, 1 top"
    dirToLabelDict["datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_MediumTauID_NoBDTGm1p0_23July2019_autoMCStats"] = "HPS-M"
    dirToLabelDict["datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_MediumTauID_NoBDTGm1p0_tauLdgTrkPt20_24July2019_autoMCStats"] = "HPS-M, p_{T}^{tk}>20"
    dirToLabelDict["datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_MediumTauID_BDTG0p40_tauLdgTrkPt20_30July2019_autoMCStats"] = "HPS-M, p_{T}^{tk}>20, 1 top"
    dirToLabelDict["datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_TightTauID_BDTG0p40_tauLdgTrkPt20_31July2019_autoMCStats"] = "HPS-T, p_{T}^{tk}>20, 1 top"
    dirToLabelDict["datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_TightTauID_NoBDTG0p40_tauLdgTrkPt20_31July2019_autoMCStats"] = "HPS-T, p_{T}^{tk}>20"
    dirToLabelDict["datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_VLooseID_BDT0p40_tauLdgTrkPt20_31July2019_autoMCStats"] = "HPS-VL, p_{T}^{tk}>20, 1 top"

    dName = os.path.basename(dirName)
    if dName in dirToLabelDict.keys():
        # print "dirToLabelDict[%s] = %s" % (dName, dirToLabelDict[dName])
        return dirToLabelDict[dName]
    else:
        return os.path.basename(dName)

        
def _ifNotNone(value, default):
    if value is None:
        return default
    return value

def doCompare(name, compareList, **kwargs):

    # Definitions
    legendLabels = []
    limits       = []

    # For-loop: All label-path pairs
    for label, path in compareList:
        legendLabels.append(label)
        dirs = glob.glob(path)
        dirs.sort()

        if len(dirs) == 0:
            raise Exception("No directories for pattern '%s'" % path)
        directory = dirs[-1]
        
        Verbose("Picked %s" % directory, True)
        limits.append(limit.BRLimits(directory, excludeMassPoints=["155"]))

    # ================================================================================================================
    # Do the sigma bands
    # ================================================================================================================
    Verbose("Creating the sigma-bands plots", True)
    _opts = copy.deepcopy(opts)
    doPlotSigmaBands(limits, legendLabels, opts.name + "_sigmaBands", _opts)

    # ================================================================================================================
    # Do expected
    # ================================================================================================================
    Verbose("Creating the expected plots", True)
    observedList = [l.observedGraph() for l in limits]
    expectedList = [l.expectedGraph() for l in limits]
    
    # 1) Expected: Median +/- 1,2 sigma
    doPlot(limits, legendLabels, expectedList, opts.name + "_median", opts.sigmaBRlimit, _opts, yTitle="Expected Sigma Bands")

    # 2) Expected: +/- 1 sigma
    list1 = [l.expectedGraph(sigma=+1) for l in limits]
    list2 = [l.expectedGraph(sigma=-1) for l in limits]
    exp1sigmaList = list1 + list2
    legendLabels2 = legendLabels + [None]*len(legendLabels)
    doPlot(limits, legendLabels2, exp1sigmaList, opts.name + "_sigma1", "Expected #pm1#sigma", _opts, yTitle="Expected #pm1sigma")

    # 3) Expected: +/- 2 sigma
    list1 = [l.expectedGraph(sigma=+2) for l in limits]
    list2 = [l.expectedGraph(sigma=-2) for l in limits]
    exp2sigmaList = list1 + list2
    doPlot(limits, legendLabels2, exp2sigmaList, opts.name + "_sigma2", "Expected #pm2#sigma", _opts, yTitle="Expected #pm2sigma")

    # ================================================================================================================
    # Do the observed plots
    # ================================================================================================================
    Verbose("Creating the observed plots", True)
    if opts.unblinded:
        doPlot(limits, legendLabels, observedList, opts.name, opts.sigmaBRlimit, _opts, yTitle="Observed")

    # ================================================================================================================
    # Do the relative plots
    # ================================================================================================================
    Verbose("Creating the relative plots", True)
    if not opts.relative:
        return
    # Overwrite some settings
    _opts.logY = False
    if _opts.logY:
        _opts.yMin    = 1e-1
        _opts.yMax    = 1e+1
    else:
        _opts.yMin    = +0.65
        _opts.yMax    = +2.05
        #_opts.yMin    = +0.85
        #_opts.yMax    = +1.25

    # 1) Relative: median
    relLimits    = GetRelativeLimits(limits)
    expectedList = [l.expectedGraph() for l in relLimits]
    doPlot(relLimits, legendLabels, expectedList, opts.name + "_medianRelative", opts.relativeYlabel, _opts, yTitle="Expected median")

    # 2) Relative: (expected 1 sigma) / (median)
    list1 = [limit.divideGraph(l.expectedGraph(sigma=+1), l.expectedGraph()) for l in limits]
    list2 = [limit.divideGraph(l.expectedGraph(sigma=-1), l.expectedGraph()) for l in limits]
    sigma1List = list1 + list2
    doPlot(limits, legendLabels2, sigma1List, opts.name + "_sigma1Relative", "Expected #pm1#sigma / median", _opts, yTitle="Expected #pm1#sigma / median")

    # 3) Relative: (expected 2 sigma) / (median)
    list1 = [limit.divideGraph(l.expectedGraph(sigma=+2), l.expectedGraph()) for l in limits]
    list2 = [limit.divideGraph(l.expectedGraph(sigma=-2), l.expectedGraph()) for l in limits]
    sigma2List = list1 + list2
    doPlot(limits, legendLabels2, sigma2List, opts.name + "_sigma2Relative", "Expected #pm2#sigma / median", _opts, yTitle="Expected #pm2#sigma / median")
    return


def doOverlap(name, compareList, **kwargs):

    # Define lists
    labels = []
    paths  = []
    limits = []

    # For-loop: All label-path pairs
    for label, path in compareList:
        labels.append(label)
        dirs = glob.glob(path)
        dirs.sort()

        if len(dirs) == 0:
            raise Exception("No directories for pattern '%s'" % path)
        directory = dirs[-1]
        
        Verbose("Picked %s" % directory, True)
        limits.append( limit.BRLimits(directory, opts.analysisType, opts.excludeMassPoints, opts.limitsFile, opts.cfgFile) )
        paths.append(path)

    Print("Creating the sigma-bands plots", True)
    _opts = copy.deepcopy(opts)
    doPlotOverlap(limits, paths, labels, opts.name + "_overlap", _opts)
    return

def doPlot(limits, legendLabels, graphs, name, ylabel, _opts={}, yTitle=None):
    
    # Definitions
    hg = []
    ll = {}
    nGraphs = len(graphs)

    # For-loop: All HistoGraphs
    for i in xrange(nGraphs):
        hg.append(histograms.HistoGraph(graphs[i], "Graph%d"%i, drawStyle="PL", legendStyle="lp"))
        if opts.boldText:
            ll["Graph%d" % (i) ] = legendLabels[i]
        else:
            ll["Graph%d" % (i) ] = "#font[42]{%s}" % legendLabels[i]

    # Create a plot-base object
    plot = plots.PlotBase(hg)
    plot.histoMgr.forEachHisto(styles.Generator(styleList[0:len(limits)]))
    def sty(h):
        r = h.getRootHisto()
        r.SetLineWidth(3)
        r.SetLineStyle(ROOT.kSolid)
        return

    # Apply style and set label
    plot.histoMgr.forEachHisto(sty)
    plot.histoMgr.setHistoLegendLabelMany(ll)

    # Create & set legend
    nGraphs = len(graphs)
    # If sigma bands are drawn each legend entry is plotted twice. Correct this in the count
    if "Sigma1" in name or "Sigma2" in name:
        nGraphs = nGraphs/2.0
    legend = getLegend(nGraphs, limit)
    plot.setLegend(legend)

    # Determine save name, minimum and maximum of y-axis
    ymin, ymax, saveName = getYMinMaxAndName(limits, name)
    if _opts.yMin == -1:
        _opts.yMin = ymin
    if _opts.yMax == -1:
        _opts.yMax = ymax

    # Create the frame and set axes titles
    if _opts.xMax != -1:
        plot.createFrame(saveName, opts={"xmax": _opts.xMax, "ymin": _opts.yMin, "ymax": _opts.yMax})
        if _opts.xMin != -1:
            plot.createFrame(saveName, opts={"xmin": _opts.xMin, "xmax": _opts.xMax, "ymin": _opts.yMin, "ymax": _opts.yMax})
    else:
        plot.createFrame(saveName, opts={"ymin": _opts.yMin, "ymax": _opts.yMax})
    
    # Add cut lines?
    for x in opts.cutLinesX:
        kwargs = {"greaterThan": True}
        plot.addCutBoxAndLine(cutValue=int(x), fillColor=ROOT.kRed, box=False, line=True, **kwargs)

    for y in opts.cutLinesY:
        kwargs = {"greaterThan": True, "mainCanvas": True, "ratioCanvas": False}
        plot.addCutBoxAndLineY(cutValue=int(y), fillColor=ROOT.kRed, box=False, line=True, **kwargs)

    if opts.bandValue != 0:
        if "relative" in name.lower():
            # https://root.cern.ch/doc/master/classTAttFill.html
            kwargs = {"cutValue": 1.0 + float(opts.bandValue)/100.0, "fillColor": ROOT.kGray, "fillStyle": 3001, "box": False, "line": True, "greaterThan": True, "mainCanvas": True, "ratioCanvas": False, "mirror": True}
            plot.addCutBoxAndLineY(**kwargs)    

    # Set axes titles
    plot.frame.GetXaxis().SetTitle(limit.mHplus())
    plot.frame.GetYaxis().SetTitle(ylabel)

    # Enable/Disable logscale for axes 
    ROOT.gPad.SetLogy(_opts.logY)
    ROOT.gPad.SetLogx(_opts.logX)

    # Draw and add text
    plot.draw()
    plot.setLuminosity(limits[0].getLuminosity())
    plot.addStandardTexts(cmsTextPosition="outframe")
    addPhysicsText(histograms, limit, x=0.53)

    # Save plots and return
    SavePlot(plot, _opts.saveDir, saveName, opts.saveFormats)
    return

def GetRelativeLimits(limits):
    
    # Create relative limits by copying limits
    relLimits  = copy.deepcopy(limits)
    massPoints = len(relLimits[0].expectedMedian)

    # For-loop: All limits
    for i in range(1, len(relLimits) ):
        relLimits[i].divideByLimit(relLimits[0])
        
    # For-loop: All mass points 
    for j in range(0, massPoints):
        relLimits[0].expectedMedian[j] = 1.0
        relLimits[0].expectedMinus2[j] = 1.0
        relLimits[0].expectedMinus1[j] = 1.0
        relLimits[0].expectedPlus2[j]  = 1.0
        relLimits[0].expectedPlus1[j]  = 1.0
        relLimits[0].observed[j]       = 1.0
    return relLimits

        
def addPhysicsText(histograms, limit, x=0.45, y=0.84, size=20):
    '''
    Add physics-process text on canvas
    '''
    # Add process text on canvas
    histograms.addText(x, y+0.04, opts.hardProcess, size=size, bold=opts.boldText)

    # Add final-state text
    histograms.addText(x, y, opts.finalState, size=size, bold=opts.boldText)

    # Add BR assumption (if applicaple)
    if len(opts.brAssumption)>0:
        histograms.addText(x, y-0.05, opts.brAssumption, size=size, bold=opts.boldText)
    return

def getLegend(nPlots, limit, xLeg1=0.53):
    dy = - (nPlots-3)*0.03
    
    # Create customised legend
    xLeg2 = 0.93
    yLeg1 = 0.60 + dy
    yLeg2 = 0.80 #+ dy
    if opts.unblinded:
        yLeg2 = 0.92 + dy

    # Adjust legend slightly to visually align with text
    legend = histograms.createLegend(xLeg1*.98, yLeg1, xLeg2, yLeg2)
    legend.SetMargin(0.17)

    # Make room for the final state text
    if opts.opaqueLegend:
        legend.SetFillStyle(1001)
    return legend


def doPlotSigmaBands(limits, legendLabels, saveName, _opts={}):

    # Define graphs to be used
    graphs = [
        histograms.HistoGraph(limits[0].expectedGraph(), "Expected", drawStyle="L"),
        histograms.HistoGraph(limits[0].expectedBandGraph(sigma=1), "Expected1", drawStyle="F", legendStyle="fl"),
        histograms.HistoGraph(limits[0].expectedBandGraph(sigma=2), "Expected2", drawStyle="F", legendStyle="fl"),
        ]

    # Set line style
    graphs[0].getRootHisto().SetLineStyle(ROOT.kSolid)

    # Create plot base object
    plot = plots.PlotBase(graphs)
    plot.setLuminosity(opts.intLumi)
    lExp1 = "%s #pm 1#sigma" % legendLabels[0],
    lExp2 = "%s #pm 2#sigma" % legendLabels[0],
    if not opts.boldText:
        lExp1 = "#font[42]{%s}" % lExp1
        lExp2 = "#font[42]{%s}" % lExp2
        
    ll = {
        "Expected" : None,
        "Expected1": lExp1,
        "Expected2": lExp2,
        }

    stGen   = styles.generator()
    nLimits = len(limits)

    # For-loop: All limits
    for i in xrange(1, nLimits):
        name = "Exp%d" % i
        gr   = histograms.HistoGraph(limits[i].expectedGraph(), name, drawStyle="L")
        stGen(gr)
        gr.getRootHisto().SetLineWidth(3)
        gr.getRootHisto().SetLineStyle(1)
        plot.histoMgr.insertHisto(len(plot.histoMgr)-2, gr, legendIndex=len(plot.histoMgr))
        label = legendLabels[i]
        if not opts.boldText:
            label = "#font[42]{%s}" % label
        ll[name] = label

    # Set histo labels
    plot.histoMgr.setHistoLegendLabelMany(ll)

    # Create & set legend
    nGraphs = len(graphs)

    # If sigma bands are drawn each legend entry is plotted twice. Correct this in the count
    if "Sigma1" in name or "Sigma2" in name:
        nGraphs = nGraphs/2.0
    legend = getLegend(nGraphs+4, limit)
    plot.setLegend(legend)

    # Determine save name, minimum and maximum of y-axis
    ymin, ymax, saveName = getYMinMaxAndName(limits, saveName)
    if _opts.yMin == -1:
        _opts.yMin = ymin
    if _opts.yMax == -1:
        _opts.yMax = ymax

    # Create the frame and set axes titles
    plot.createFrame(saveName, opts={"ymin": _opts.yMin, "ymax": _opts.yMax})

    # Add cut lines?
    for x in opts.cutLinesX:
        kwargs = {"greaterThan": True}
        plot.addCutBoxAndLine(cutValue=int(x), fillColor=ROOT.kRed, box=False, line=True, **kwargs)

    for y in opts.cutLinesY:
        kwargs = {"greaterThan": True, "mainCanvas": True, "ratioCanvas": False}
        plot.addCutBoxAndLineY(cutValue=int(y), fillColor=ROOT.kRed, box=False, line=True, **kwargs)

    # Set axes titles
    plot.frame.GetXaxis().SetTitle(opts.mHplus)
    plot.frame.GetYaxis().SetTitle(opts.sigmaBRlimit)

    # Enable/Disable logscale for axes 
    ROOT.gPad.SetLogy(_opts.logY)
    ROOT.gPad.SetLogx(_opts.logX)

    # Draw the plot with standard texts
    plot.draw()
    plot.addStandardTexts(addLuminosityText=True) #cmsTextPosition="outframe")
    plot.setLuminosity(limits[0].getLuminosity())
    addPhysicsText(histograms, limit, x=0.53)

    # Save the plots & return
    SavePlot(plot, _opts.saveDir, saveName, opts.saveFormats)
    return


def doPlotOverlap(limits, paths, labels, saveName, _opts={}):

    analyses = paths
    graphs   = [
        histograms.HistoGraph(limits[1].expectedGraph(), "%s-Median" % (analyses[1]), drawStyle="L"),
        histograms.HistoGraph(limits[0].expectedGraph(), "%s-Median" % (analyses[0]), drawStyle="L"),
        histograms.HistoGraph(limits[1].expectedBandGraph(sigma=1), "%s-1Sigma" % (analyses[1]), drawStyle="F", legendStyle="F"),
        histograms.HistoGraph(limits[0].expectedBandGraph(sigma=1), "%s-1Sigma" % (analyses[0]), drawStyle="F", legendStyle="F"),
        histograms.HistoGraph(limits[1].expectedBandGraph(sigma=2), "%s-2Sigma" % (analyses[1]), drawStyle="F", legendStyle="F"),
        histograms.HistoGraph(limits[0].expectedBandGraph(sigma=2), "%s-2Sigma" % (analyses[0]), drawStyle="F", legendStyle="F"),
        ]
    
    # Create plot base object
    if 0:
        plot = plots.PlotBase(graphs)
    else:
        plot = plots.ComparisonManyPlot(graphs[0], graphs[1:])

    # Customise legend
    ll = {
        "%s-Median" % analyses[0]: labels[0] + " Median",
        "%s-Median" % analyses[1]: labels[1] + " Median",
        "%s-1Sigma" % analyses[0]: labels[0] + "#pm 1#sigma",
        "%s-2Sigma" % analyses[0]: labels[0] + "#pm 2#sigma",
        "%s-1Sigma" % analyses[1]: None, #labels[1] + "#pm 1#sigma", #None,
        "%s-2Sigma" % analyses[1]: None, #labels[1] + "#pm 2#sigma", #None,
         }

    if not opts.boldText:
        for k in ll:
            if ll[k] == None:
                continue
            else:
                ll[k] = "#font[42]{%s}" % ll[k]

    # https://root.cern.ch/doc/master/classTStyle.html#a94c9ae34f921c1efaa09348c2d16d60a
    sList = [3, 5, 6, 8] # [ROOT.kDashed, ROOT.kDotted]
    cList = [ROOT.kRed, ROOT.kBlue]
    fList = [1001, 1001, 1001, 1001] #[1001, 3005, 1001, 3005]
    count1 = 0
    count2 = 0
    # For-loop: All limits
    for i, gr in enumerate(graphs, 1):
        name = gr.getName().lower()
        gr.getRootHisto().SetLineWidth(3)

        if "median" in name:
            gr.getRootHisto().SetLineStyle(sList[count1])
            gr.getRootHisto().SetLineColor(cList[count1])
            count1+=1
        elif "sigma" in name:
            gr.getRootHisto().SetLineStyle(ROOT.kSolid)
            gr.getRootHisto().SetLineWidth(4)
            gr.getRootHisto().SetFillStyle(fList[count2])
            count2+=1
        else:
            pass

    # Set histo labels
    plot.histoMgr.setHistoLegendLabelMany(ll)

    # If sigma bands are drawn each legend entry is plotted twice. Correct this in the count
    legend = getLegend(0, limit)
    plot.setLegend(legend)
    histograms.moveLegend(legend, dx=+0.0, dy=0.0, dw=0.0, dh=+0.1) #dh=+0.2)
    
    # Determine save name, minimum and maximum of y-axis
    ymin, ymax, saveName = getYMinMaxAndName(limits, saveName)
    if _opts.yMin == -1:
        _opts.yMin = ymin
    if _opts.yMax == -1:
        _opts.yMax = ymax

    # Create the frame and set axes titles
    plot.createFrame(saveName, opts={"ymin": _opts.yMin, "ymax": _opts.yMax})
    
    # Add cut line?
    for x in opts.cutLinesX:
        kwargs = {"greaterThan": True}
        plot.addCutBoxAndLine(cutValue=int(x), fillColor=ROOT.kRed, box=False, line=True, **kwargs)

    for y in opts.cutLinesY:
        kwargs = {"greaterThan": True, "mainCanvas": True, "ratioCanvas": False}
        plot.addCutBoxAndLineY(cutValue=int(y), fillColor=ROOT.kRed, box=False, line=True, **kwargs)

    # Set axes titles
    plot.frame.GetXaxis().SetTitle(opts.mHplus)
    plot.frame.GetYaxis().SetTitle(opts.sigmaBRlimit)

    # Enable/Disable logscale for axes 
    ROOT.gPad.SetLogy(_opts.logY)
    ROOT.gPad.SetLogx(_opts.logX)

    # Draw the plot with standard texts
    plot.draw()
 
    plot.setLuminosity(limits[0].getLuminosity())
    plot.addStandardTexts(addLuminosityText=True, cmsTextPosition="outframe")
    addPhysicsText(histograms, limit, x=0.53)

    # Save the plots & return
    SavePlot(plot, _opts.saveDir, saveName, opts.saveFormats)
    return

def getYMinMaxAndName(limits, name, minIsMedian=False):
    ymin = 1e6
    ymax = -1e6

    # For-loop: all limits
    for l in limits:
        if minIsMedian:
            _ymin = l.getYMinMedian()
        else:
            _ymin = l.getYMin()
        _ymax = l.getYMax()
        if _ymin < ymin:
            ymin = _ymin
        if _ymax > ymax:
            ymax = _ymax
        
    if opts.logY: #fixme
        name += "_logY"
        ymax *= 2
    else:
        ymin =  0.0
        ymax *= 1.2

    if opts.logX: #fixme
        name += "_logX"
    return ymin, ymax, name

def SavePlot(plot, saveDir, plotName, saveFormats = [".png", ".pdf"]):
    # Check that path exists
    if not os.path.exists(saveDir):#fixme
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
    ANALYSISTYPE = "HToTauNu"
    BANDVALUE    = 0
    BATCHMODE    = True    
    CFGFILE      = "configuration.json"
    CUTLINE      = 999.9
    CUTLINEX     = ""
    CUTLINEY     = ""
    EXCLUDE      = ""
    GRIDX        = False
    GRIDY        = False
    LIMITSFILE   = "limits.json"
    LOGX         = False
    LOGY         = False
    MAXX         = -1
    MAXY         = -1
    MINX         = -1
    MINY         = -1
    NAME         = "limitsBr"
    OPAQUELEGEND = False
    PAPER        = False
    RELATIVE     = False
    RELPAIRS     = False
    SAVEDIR      = None
    SAVEFORMATS  = "png" #pdf,png,C"
    UNBLINDED    = False
    URL          = False
    VERBOSE      = False
    BOLDTEXT     = False
    DIRS         = None

    parser = OptionParser(usage="Usage: %prog [options]", add_help_option=True, conflict_handler="resolve")

    parser.add_option("-v", "--verbose", dest="verbose", default=VERBOSE, action="store_true",
                      help="Verbose mode for debugging purposes [default: %s]" % (VERBOSE) )
    
    parser.add_option("--boldText", dest="boldText", default=BOLDTEXT, action="store_true",
                      help="Use bold text printed on canvas? [default: %s]" % (BOLDTEXT))

    parser.add_option("--dirs", dest="dirs", default=DIRS,
                      help="List for datacard directories draw in comparison (comma separated WITHOUT space) [default: %s]" % (DIRS))
    
    parser.add_option("--limitsFile", dest="limitsFile", default=LIMITSFILE,
                      help="JSON file containing the limits [default: %s]" % (LIMITSFILE))
    
    parser.add_option("--cfgFile", dest="cfgFile", default=CFGFILE,
                      help="JSON file containing the configurations [default: %s]" % (CFGFILE))

    parser.add_option("--exclude", dest="exclude", default=EXCLUDE,
                      help="List for mass points to exclude (comma separated WITHOUT space) [default: %s]" % (EXCLUDE))
  
    parser.add_option("--unblinded", dest="unblinded", default=UNBLINDED, action="store_true",
                    help="Enable unblined mode and thus also produced observed limits [default: %s]" % (UNBLINDED) )

    parser.add_option("--analysisType", dest="analysisType", default=ANALYSISTYPE,
                      help="Flag to indicate the analysis type (e.g. \"HToTauNu\", \"HToTB\", \"HToHW\") [default: %s]" % (ANALYSISTYPE) )
 
    parser.add_option("--paper", dest="paper", default=PAPER, action="store_true",
                      help="Paper mode [default: %s]" % (PAPER) )
    
    parser.add_option("--relative", dest="relative", action="store_true", default=RELATIVE, 
                      help="Do comparison relative to the first item [default: %s]" % (RELATIVE) )

    parser.add_option("--relativePairs", dest="relativePairs", action="store_true", default=RELPAIRS, 
                      help="Do multiple relative comparisons. The list of input directories is halved, the first half is the denominator and the second half is the numerator [default: %s]" % (RELPAIRS) )

    parser.add_option("--name", dest="name", type="string", default=NAME,
                      help="Name of the output plot [default = %s]" % (NAME))

    parser.add_option("--relativeYmax", dest="relativeYmax", type="float", default=None, 
                      help="Maximum y-value for relative plots")
    
    parser.add_option("--relativeYlabel", dest="relativeYlabel", default="95% CL upper limit ratio", 
                      help="Y-axis title for relative plots")

    parser.add_option("--url", dest="url", action="store_true", default=URL,
                      help="Don't print the actual save path the plots are saved, but print the URL instead [default: %s]" % URL)

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR,
                      help="Directory where all plots will be saved [default: %s]" % SAVEDIR)

    parser.add_option("--opaqueLegend", dest="opaqueLegend", default=OPAQUELEGEND, action="store_true",
                      help="Add excluded area as in MSSM exclusion plots [default: %s]" % (OPAQUELEGEND) )

    parser.add_option("--logX", dest="logX", action="store_true", default=LOGX,
                      help="Plot x-axis (mass) as logarithmic [default: %s]" % (LOGX) )
    
    parser.add_option("--logY", dest="logY", action="store_true", default=LOGY,
                      help="Plot y-axis (exlusion limit) as logarithmic [default: %s]" % (LOGY) )
    
    parser.add_option("--gridX", dest="gridX", default=GRIDX, action="store_true",
                      help="Enable the grid for the x-axis [default: %s]" % (GRIDX) )

    parser.add_option("--gridY", dest="gridY", default=GRIDY, action="store_true",
                      help="Enable the grid for the y-axis [default: %s]" % (GRIDY) )

    parser.add_option("--yMin", dest="yMin", default=MINY, type="float",
                      help="Overwrite automaticly calculated minimum value of y-axis [default: %s]" % (MINY) )
    
    parser.add_option("--yMax", dest="yMax", default=MAXY, type="float",
                      help="Overwrite automaticly calculated maximum value of y-axis [default: %s]" % (MAXY) )

    parser.add_option("--xMin", dest="xMin", default=MINX, type="float",
                      help="Overwrite minimum value of x-axis [default: %s]" % (MINX) )
    
    parser.add_option("--xMax", dest="xMax", default=MAXX, type="float",
                      help="Overwrite maximum value of x-axis [default: %s]" % (MAXX) )

    parser.add_option("-s", "--saveFormats", dest="saveFormats", default = SAVEFORMATS,
                      help="Save formats for all plots [default: %s]" % SAVEFORMATS)

    parser.add_option("-b", "--batchMode", dest="batchMode", action="store_false", default=BATCHMODE,
                      help="Enables batch mode (canvas creation does NOT generate a window) [default: %s]" % BATCHMODE)

    parser.add_option("--cutLineX", dest="cutLineX", default=CUTLINEX,
                      help="List for values for x-axis lines to be drawn on canvas (comma separated WITHOUT space) [default: %s]" % (CUTLINEX))

    parser.add_option("--cutLineY", dest="cutLineY", default=CUTLINEY,
                      help="List for values for y-axis lines to be drawn on canvas (comma separated WITHOUT space) [default: %s]" % (CUTLINEY))

    parser.add_option("--bandValue", dest="bandValue", type="int", default=BANDVALUE,
                      help="Add a symmetric band around 1.0. Value passed should be the percentage (e.g 10 or 5)  [default: %s]" % (BANDVALUE) )

    (opts, args) = parser.parse_args()


    # Sanity checks
    opts.dirList = []
    if opts.dirs == None:
        msg = "No datacard directories provided!"
        raise Exception(sh_e + msg + sh_n)
    else:
        opts.dirList = opts.dirs.split(",")

    if len(opts.dirList) > 1:
        msg = "Datacard directories considered:%s\n\t%s" % (sh_t, "\n\t".join(opts.dirList))
        Verbose(msg + sh_n, True)
    else:
        msg = "At least 2 datacard directories required. Only %d passed with --dirs argument!" % len(opts.dirsList)
        raise Exception(sh_e + msg + sh_n)


    # Save in current working directory?
    if opts.saveDir =="":
        opts.saveDir = os.getcwd()

    # Sanity check (analysis type)
    myAnalyses = ["HToTauNu", "HToTB", "HToHW"]
    if opts.analysisType not in myAnalyses:
        msg = "Invalid analysis type \"%s\". Please selected one of the following: \"%s" % (opts.analysisType, "\", \"".join(myAnalyses) + "\"")
        raise Exception(sh_e + msg + sh_n)
    else:
        msg = "Analysis type is %s" % (sh_t + opts.analysisType + sh_n)
        Print(msg, True)        

    # Sanity checks
    opts.excludeMassPoints = []
    if len(opts.exclude) > 0:
        opts.excludeMassPoints = opts.exclude.split(",")
    if len(opts.excludeMassPoints) > 0:
        msg = "Excluding mass points \"%s\"" % (", ".join(opts.excludeMassPoints))
        Print(sh_t + msg + sh_n, True)
    else:
        pass

    # Cut lines
    opts.cutLinesX = []
    if len(opts.cutLineX) > 0:
        opts.cutLinesX = opts.cutLineX.split(",")
    if len(opts.cutLinesX) > 0:
        msg = "Adding x-axis lines \"%s\"" % (", ".join(opts.cutLinesX))
        Print(sh_t +  msg + sh_n, True)
    else:
        pass
    opts.cutLinesY = []
    if len(opts.cutLineY) > 0:
        opts.cutLinesY = opts.cutLineY.split(",")
    if len(opts.cutLinesY) > 0:
        msg = "Addingyx-axis lines \"%s\"" % (", ".join(opts.cutLinesY))
        Print(sh_t +  msg + sh_n, True)
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
                               
    if opts.saveDir == None:
        opts.saveDir = opts.dirList[0]

    # Create save formats
    if "," in opts.saveFormats:
        opts.saveFormats = opts.saveFormats.split(",")
    else:
        opts.saveFormats = [opts.saveFormats]
    opts.saveFormats = ["." + s for s in opts.saveFormats]

    # Call the main function
    main()
    
    if not opts.batchMode:
        raw_input("=== plotBRLimits.py: Press any key to quit ROOT ...")
