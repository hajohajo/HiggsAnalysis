#!/usr/bin/env python
'''
DESCRIPTION:


USAGE:
./plotTauFakeRates.py -m <pseudo_mcrab_directory> [opts]


EXAMPLES:
./plotTauFakeRates.py -m TauFakeRate_Attempt4_MuonPt40_AtLeast2Jets_08Feb2020 --numHisto "tauPt_num_1pr" --denHisto "tauPt_den_1pr" -i "SingleMu" --logY --yMin 1e-3 --yMax 1
./plotTauFakeRates.py -m Hplus2hwAnalysis_fake_191016_180629 --numHisto "tauPt_num_1pr" --denHisto "tauPt_den_1pr" --analysisName Hplus2hwAnalysis_fake --dataEra "Run2016" --searchMode 350to3000
./plotTauFakeRates.py -m TauFakeRate_Attempt4_MuonPt40_AtLeast2Jets_08Feb2020 --numHisto "tauPt_num_1pr" --denHisto "tauPt_den_1pr"
./plotTauFakeRates.py -m TauFakeRate_Attempt4_MuonPt40_AtLeast2Jets_08Feb2020 --numHisto "tauPt_num_1pr" --denHisto "tauPt_den_1pr" -s png --cutLineX 45.0 --cutLineY 0.4
./plotTauFakeRates.py -m TauFakeRate_Attempt4_MuonPt40_AtLeast2Jets_08Feb2020 --numHisto "tauPt_num_1pr" --denHisto "tauPt_den_1pr" -s png --gridX --gridY --yMin 0.0 --yMax 0.6
./plotTauFakeRates.py -m TauFakeRate_Attempt4_MuonPt40_AtLeast2Jets_08Feb2020 --numHisto "tauPt_num_1pr" --denHisto "tauPt_den_1pr" -s png --gridX --gridY --yMin 0.0 --yMax 0.6 -e "WJets" --individualMC


LAST USED (HToHW_withTop):
./plotTauFakeRates.py -m TauFakeRate_Attempt4_MuonPt40_AtLeast2Jets_08Feb2020 --numHisto "tauPt_num_1pr" --denHisto "tauPt_den_1pr" --gridX --gridY --yMin 0.0 --yMax 1.0 -e "WJets" --individualMC


LAST USED (HToHW):
./plotTauFakeRates.py -m Hplus2hwAnalysis_fake_191016_180629 --numHisto "tauPt_num_1pr" --denHisto "tauPt_den_1pr" -s png --gridX --gridY --yMin 0.0 --yMax 0.6 --analysisName Hplus2hwAnalysis_fake --dataEra "Run2016" --searchMode 350to3000

'''

#================================================================================================ 
# Imports
#================================================================================================ 
import sys
import math
import copy
import os
import array
from optparse import OptionParser

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import *

import HiggsAnalysis.NtupleAnalysis.tools.dataset as dataset
import HiggsAnalysis.NtupleAnalysis.tools.systematics as systematics
import HiggsAnalysis.NtupleAnalysis.tools.histograms as histograms
import HiggsAnalysis.NtupleAnalysis.tools.counter as counter
import HiggsAnalysis.NtupleAnalysis.tools.tdrstyle as tdrstyle
import HiggsAnalysis.NtupleAnalysis.tools.styles as styles
import HiggsAnalysis.NtupleAnalysis.tools.plots as plots
import HiggsAnalysis.NtupleAnalysis.tools.crosssection as xsect
import HiggsAnalysis.NtupleAnalysis.tools.aux as aux
import HiggsAnalysis.NtupleAnalysis.tools.multicrabConsistencyCheck as consistencyCheck
import HiggsAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles
import HiggsAnalysis.NtupleAnalysis.tools.analysisModuleSelector as analysisModuleSelector
import HiggsAnalysis.NtupleAnalysis.tools.errorPropagation as errorPropagation


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

    # Apply TDR style
    style = tdrstyle.TDRStyle()
    style.setGridX(opts.gridX)
    style.setGridY(opts.gridY)
    #style.setLogY(opts.logY) # problem if kwargs["ratio"] = True
    style.setLogX(opts.logX)

    # Obtain dsetMgrCreator and register it to module selector
    dsetMgrCreator = dataset.readFromMulticrabCfg(directory=opts.mcrab)

    # Get list of eras, modes, and optimisation modes
    erasList      = dsetMgrCreator.getDataEras()
    modesList     = dsetMgrCreator.getSearchModes()
    optList       = dsetMgrCreator.getOptimizationModes()
    sysVarList    = dsetMgrCreator.getSystematicVariations()
    sysVarSrcList = dsetMgrCreator.getSystematicVariationSources()

    # If user does not define optimisation mode do all of them
    if opts.optMode == None:
        if len(optList) < 1:
            optList.append("")
        else:
            pass
        optModes = optList
    else:
        optModes = [opts.optMode]


    # For-loop: All optimisation modes
    for opt in optModes:
        opts.optMode = opt
        mcList = []

        # Setup & configure the dataset manager 
        datasetsMgr = GetDatasetsFromDir(opts)
        datasetsMgr.updateNAllEventsToPUWeighted()
        datasetsMgr.loadLuminosities() # from lumi.json
        
        # Print PSets used for FakeBMeasurement
        if opts.verbose:
            datasetsMgr.printSelections()

        if opts.verbose:
            datasetsMgr.PrintCrossSections()
            datasetsMgr.PrintLuminosities()
            datasetsMgr.PrintInfo()
               
        # Merge histograms (see NtupleAnalysis/python/tools/plots.py) 
        plots.mergeRenameReorderForDataMC(datasetsMgr, keepSourcesMC=False, analysisType="HToHW_withTop") 

        if opts.intLumi < 0:
            if "Data" in datasetsMgr.getAllDatasetNames():
                opts.intLumi = datasetsMgr.getDataset("Data").getLuminosity()
                Verbose("Automatically determined the integrated luminosity to be %.2f 1/fb" % (opts.intLumi), True)            
            else:
                opts.intLumi = 1.0
                Verbose("Could not automatically determine the integrated luminosity. Set value to %.2f 1/fb" % (opts.intLumi), True)

        global dataName
        global mcName

        foundData = False
        foundMC   = False
        # For-loop: All dataset names to set/overwrite cross-sections
        for d in datasetsMgr.getAllDatasetNames():
            if "ChargedHiggs" in d:
                xs = 1.0
                Verbose("Setting cross-section for dataset to %.2f" % xs, True)
                datasetsMgr.getDataset(d).setCrossSection(xs) # HIG-18-014: mH=200-3000 GeV sigmaBR = 0.65->0.005 (pb)
                datasetsMgr.remove(d)
            elif "Data" in d:
                foundData = True
                dataName  = d
            else:
                foundMC = True
                mcList.append(d)
                mcName = "Simulation"

        # Merge all MC datasets into a a MC-soup dataset called "Simulation"  ?
        if 0:
            datasetsMgr.merge(mcName, mcList) #iro

        if not foundData:
            raise Exception("Data samples are required for this script! (foundData = %s)" % (foundData) )
        if not foundMC:
            raise Exception("MC samples are required for this script! (foundMC = %s)" % (foundMC) )

        # Print final dataset summary
        datasetsMgr.PrintInfo()

        # Get the efficiency graphs (Efficiency Vs Cut-flow)
        hGraphList, kwargs = GetHistoGraphs(datasetsMgr)

        # Plot the histo graphs
        PlotHistoGraphs(hGraphList, kwargs)

    Print("All plots saved under directory %s" % (ts + aux.convertToURL(opts.saveDir, opts.url) + ns), True)
    return


def GetHistoGraphs(datasetsMgr):

    # Get histogram customisations
    _kwargs  = GetHistoKwargs(opts)

    # Other definitions
    opts.numHistoG = opts.numHisto.replace("num_","num_g_")
    opts.denHistoG = opts.denHisto.replace("den_","den_g_")

    # Get histos
    p_num   = plots.DataMCPlot(datasetsMgr, opts.numHisto, saveFormats=[])
    p_den   = plots.DataMCPlot(datasetsMgr, opts.denHisto, saveFormats=[])
    pMC_num = plots.MCPlot(datasetsMgr, opts.numHistoG, normalizeToLumi=opts.intLumi, saveFormats=[])
    pMC_den = plots.MCPlot(datasetsMgr, opts.denHistoG, normalizeToLumi=opts.intLumi, saveFormats=[])

    # Make plots-object list
    pList = []
    pList.append(p_num)
    pList.append(pMC_num)
    pList.append(p_den)
    pList.append(pMC_den)

    # Customise histos (unstacked histos need to be beautified)
    for p in pList:
        p.histoMgr.setHistoLegendStyleAll("L")
        p.histoMgr.setHistoDrawStyleAll("HIST")

    # Draw the histograms
    plots.drawPlot(p_num  , opts.numHisto , **_kwargs)
    plots.drawPlot(p_den  , opts.numHistoG, **_kwargs)
    plots.drawPlot(pMC_num, opts.denHisto , **_kwargs)
    plots.drawPlot(pMC_den, opts.denHistoG, **_kwargs)

    # Save the histos
    saveDir = os.path.join(opts.saveDir, "Src")
    SavePlot(p_num  , opts.numHisto , saveDir, saveFormats = opts.saveFormats)
    SavePlot(p_den  , opts.denHisto , saveDir, saveFormats = opts.saveFormats)
    SavePlot(pMC_num, opts.numHistoG, saveDir, saveFormats = opts.saveFormats)
    SavePlot(pMC_den, opts.denHistoG, saveDir, saveFormats = opts.saveFormats)


    # Clone histograms 
    opts.datasets = datasetsMgr.getAllDatasetNames()

    # Convert histos to TGraph
    hgList = []
    gList  = convertHisto2TGraph(datasetsMgr, p_num, pMC_num, p_den, pMC_den, printValues=False)
    dIndex = None
    sIndex = None

    # For-loop: All graphs in the list
    for i, g in enumerate(gList, 0):
        name = g.GetName().replace("Eff_", "")
        Print("%d/%d: Processing graph \"%s\"" % (i+1, len(gList), name), i==0)    
        
        # Apply unique style
        plots._plotStyles[name].apply(g)
        
        # Store indices for "Data" and "Simulation" histograms
        if name == "Data":
            dIndex = i
        if name == "Simulation":
            sIndex = i
            g.SetFillColor(ROOT.kPink)
            g.SetMarkerStyle(ROOT.kFullSquare)
            g.SetFillStyle(3002)

        # Create histoGraph object
        if i == 0:
            hg = histograms.HistoGraph(g, g.GetTitle(), legendStyle="LP", drawStyle="P")
        else:
            hg = histograms.HistoGraph(g, g.GetTitle(), legendStyle="LP", drawStyle="P")
                
        # Append to list
        hgList.append(hg)

    # Move "Data - Genuine #tau_{h}" histo to first in list
    hgList.insert(0, hgList.pop(dIndex))
    hgList.insert(1, hgList.pop(sIndex))
    
    if not opts.individualMC:
        hgList = hgList[:2]
    opts.nDatasets = len(hgList)
    return hgList, _kwargs


def PlotHistoGraphs(hGraphList, kwargs):

    _kwargs = GetGraphKwargs(kwargs)

    # For-loop: All efficiencie TGraphs    
    for g in hGraphList:
        # g.getRootHisto().GetXaxis().SetTitle(_kwargs["xlabel"]) #iro
        g.getRootHisto().GetYaxis().SetTitle(_kwargs["ylabel"])
        if opts.yMin != None:
            g.getRootHisto().SetMinimum(_kwargs["opts"]["ymin"])
        if opts.yMax != None:
            g.getRootHisto().SetMaximum(_kwargs["opts"]["ymax"])

    # Create & draw the plot
    p = plots.ComparisonManyPlot(hGraphList[0], hGraphList[1:], saveFormats=[])
    p.setLuminosity(opts.intLumi)

    if 0:
        p.histoMgr.setHistoLegendStyleAll("LP")
        p.histoMgr.setHistoDrawStyleAll("LP")
        p.histoMgr.setHistoDrawStyle(hGraphList[0].getRootHisto().GetName(), "HIST")

    # Draw the plot
    if 0:
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetLineColor(ROOT.kRed))
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetLineStyle(ROOT.kDotted))
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetLineWidth(3))
        #
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetMarkerColor(ROOT.kBlue))
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetMarkerSize(1.3))
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetMarkerStyle(ROOT.kFullCircle))
     
    # https://root.cern.ch/doc/master/classTGraphPainter.html
    p.histoMgr.setHistoDrawStyle("Simulation", "P2")
    p.histoMgr.setHistoLegendStyle("Simulation", "LP")

    # Draw the comparison plot
    plots.drawPlot(p, opts.saveName, **_kwargs)

    # Save the plot
    SavePlot(p, opts.saveName, os.path.join(opts.saveDir, opts.optMode), saveFormats = opts.saveFormats )
    return


def GetGraphKwargs(kwargs):

    _kwargs = copy.deepcopy(kwargs)

    if opts.yMin != None:
        _kwargs["opts"]["ymin"] = opts.yMin
    if opts.yMax != None:
        _kwargs["opts"]["ymax"] = opts.yMax

    _kwargs["ratio"] = True
    _kwargs["ratioYlabel"] = "Ratio "
    _kwargs["stackMCHistograms"] = False
    _kwargs["addMCUncertainty" ] = False
    _kwargs["moveLegend"] =  {"dx": -0.15, "dy": -0.01, "dh": -0.18 + (opts.nDatasets*0.025)}
    #_kwargs["moveLegend"] =  {"dx": -0.52, "dy": -0.63, "dh": -0.18 + (opts.nDatasets*0.02)}
    _kwargs["log"] = opts.logY
    del _kwargs["rebinX"]
    del _kwargs["rebinY"]
    return _kwargs 


def GetHistoKwargs(opts):

    # Common bin settings
    _legRM  = {"dx": -10000.23, "dy": -10000.01, "dh": -0.1}
    _legNE  = {"dx": -0.08, "dy": -0.01, "dh": 0.1}
    _legNW  = {"dx": -0.52, "dy": -0.01, "dh": 0.1}
    _legSE  = {"dx": -0.23, "dy": -0.50, "dh": 0.1}
    _legSW  = {"dx": -0.52, "dy": -0.48, "dh": 0.1}
    _legCE  = {"dx": -0.08, "dy": -0.44, "dh": 0.1}

    kwargs = {
        "xlabel"           : "#tau_{h} p_{T} (GeV)",
        "ylabel"           : "Fake Rate (jet #rightarrow #tau_{h})",
        "rebinX"           : 1,
        #"rebinX"           : [0, 30, 60, 90, 120], #use with caution!
        "rebinY"           : None,
        "ratioType"        : "errorScale",
        "ratioErrorOptions": {"numeratorStatSyst": False},
        "ratioYlabel"      : "Data/MC ",
        "ratio"            : False, #True,
        "ratioInvert"      : False,
        "stackMCHistograms": False, #True,
        "addMCUncertainty" : True,
        "addLuminosityText": True,
        "addCmsText"       : True,
        "cmsExtraText"     : "Very Preliminary",
        "opts"             : {"ymaxfactor": 10.0},
        "opts2"            : {"ymin": 0.0, "ymax": 2.0},
        "divideByBinWidth" : False,
        "log"              : True,
        "moveLegend"       : _legNE, 
        "moveBlindedText"  : {"dx": -0.23, "dy": +0.08, "dh": 0.0},
        "errorBarsX"       : True,
        #"cutBox"           : {},
        #"cutBoxY"          : {},
        }

        
    if opts.cutLineX != None:
        Verbose("WARNING!? For TGraph better to also use the --yMin and --yMax options so that the line has definite range in y!", True)
        kwargs["cutBox"] = {"cutValue": opts.cutLineX, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
    if opts.cutLineY != None:
        Print("WARNING!? The option --cutLineY seems to be not working. Debugging needed!", True)
        kwargs["cutBoxY"] = {"cutValue": opts.cutLineY, "fillColor": 16, "box": False, "line": True, "greaterThan": True, "mainCanvas": True, "ratioCanvas": True}

    # Set x-axis divisions
    n1 = 10 # primary divisions (8)
    n2 = 0  # second order divisions (5)
    n3 = 0  # third order divisions (@)
    nDivs = n1 + 100*n2 + 10000*n3
    ROOT.gStyle.SetNdivisions(nDivs, "X")
    return kwargs

def SavePlot(plot, plotName, saveDir, saveFormats = [".C", ".png", ".pdf"]):
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    # Create the name under which plot will be saved
    saveName = os.path.join(saveDir, plotName.replace("/", "_"))

    # For-loop: All save formats
    for i, ext in enumerate(saveFormats):
        saveNameURL = saveName + ext
        saveNameURL = aux.convertToURL(saveNameURL, opts.url)
        Verbose(saveNameURL, i==0)
        plot.saveAs(saveName, formats=saveFormats)
    return


def convertHisto2TGraph(datasetsMgr, p_num, pMC_num, p_den, pMC_den, printValues=False):

    # Lists for values
    gList  = []
    hList  = []
    nBinsX = 0

    # Get the numerator histograms (tight taus)
    numHisto_Data    = p_num.histoMgr.getHisto(dataName).getRootHisto().Clone("%s_num" % dataName) 
    numHisto_Genuine = numHisto_Data.Clone("Genuine_num")
    numHisto_Genuine.Reset()
    # For-loop: All dataset names
    for d in datasetsMgr.getAllDatasetNames():
        numHisto_Genuine.Add(pMC_num.histoMgr.getHisto(d).getRootHisto(), +1)
    numHisto = numHisto_Data.Clone("Num")
    numHisto.Add(numHisto_Genuine, -1)

    # Get the denominator histograms (loose taus)
    denHisto_Data    = p_den.histoMgr.getHisto(dataName).getRootHisto().Clone("%s_den" % dataName)
    denHisto_Genuine = denHisto_Data.Clone("Genuine_num")
    denHisto_Genuine.Reset()
    # For-loop: All dataset names
    for d in datasetsMgr.getAllDatasetNames():
        denHisto_Genuine.Add(pMC_den.histoMgr.getHisto(d).getRootHisto(), +1)
    denHisto = denHisto_Data.Clone("Den")
    denHisto.Add(denHisto_Genuine, -1)
    
    # Get the efficiency histograms
    dataLabel = "%s-Genuine #tau_{h}" % (dataName)
    hEff_Data = GetEfficiencyHisto(numHisto, denHisto, "Eff_%s" % dataName, dataLabel, printValues=True) 

    # Create the "Simulation" histogram (MC-soup)
    mcName      = "Simulation"
    mcLabel     = "Simulation"
    numHisto_MC = numHisto.Clone("%s_num" % mcName)
    denHisto_MC = denHisto.Clone("%s_den" % mcName)
    numHisto_MC.Reset()
    denHisto_MC.Reset()
    hEff_MC = GetEfficiencyHisto(numHisto_MC, denHisto_MC, "Eff_%s" % mcName, mcLabel, printValues=True) 

    # For-loop: All MC datasets
    for d in datasetsMgr.getAllDatasets():
        if d.isData():
            continue
        if "Charged" in d.getName():
            continue

        # Get information
        dsetName      = d.getName()
        dsetLabel     = plots._legendLabels[dsetName]
        numHisto_tmp  = p_num.histoMgr.getHisto(dsetName).getRootHisto().Clone("%s_num" % dsetName)
        denHisto_tmp  = p_den.histoMgr.getHisto(dsetName).getRootHisto().Clone("%s_num" % dsetName)
        numHisto_tmpG = pMC_num.histoMgr.getHisto(dsetName).getRootHisto().Clone("%s_numG" % dsetName)
        denHisto_tmpG = pMC_den.histoMgr.getHisto(dsetName).getRootHisto().Clone("%s_numG" % dsetName)

        # Subtract the Genuine tau histograms
        numHisto_tmp.Add(numHisto_tmpG, -1)
        denHisto_tmp.Add(denHisto_tmpG, -1)
        hEff_tmp = GetEfficiencyHisto(numHisto_tmp, denHisto_tmp, "Eff_%s" % dsetName, dsetLabel, printValues=True) 
        hList.append(hEff_tmp)

        # Simulation histogram (Addition of all MC samples)
        numHisto_MC.Add(numHisto_tmp, +1)
        denHisto_MC.Add(denHisto_tmp, +1)

    # Subtract the Genuine tau histograms from the "Simulation" histograms
    #numHisto_MC.Add(numHisto_Genuine, -1)
    #denHisto_MC.Add(denHisto_Genuine, -1)

    # Get the "Simulation" efficinency histo
    hEff_MC = GetEfficiencyHisto(numHisto_MC, denHisto_MC, "Eff_%s" % mcName, mcLabel, printValues=True) 

    # Append efficiency histograms to list
    hList.append(hEff_Data)
    hList.append(hEff_MC)

    # For-loop: All efficinecy histos
    for hEff in hList:
        
        xVal   = []
        yVal   = []
        xErrL  = []
        xErrH  = []
        yErrL  = []
        yErrH  = []
        
        # For-loop: All histo bins
        for j in range (1, hEff_Data.GetNbinsX()+2):
            
            # Calculate the values
            x     = hEff.GetBinCenter(j)
            y     = hEff.GetBinContent(j)
            xLow  = x-hEff.GetXaxis().GetBinLowEdge(j)
            yLow  = hEff.GetBinError(j)
            xHigh = hEff.GetXaxis().GetBinUpEdge(j)-x
            yHigh = hEff.GetBinError(j)

            # Store the values
            xVal.append(x)
            yVal.append(y)
            xErrL.append(xLow)
            yErrL.append(yLow)            
            xErrH.append(xHigh)
            yErrH.append(yHigh)
            #Verbose("xLow = %.1f, x = %.1f, xHigh = %.1f | yLow = %.1f, y = %.1f, yHigh = %.1f" % (xLow, x, xHigh, yLow, y, yHigh), j==0)
        hName  = hEff.GetName()
        hTitle = hEff.GetTitle()
        del hEff

        # Construct info table (debugging)
        table  = []
        align  = "{:>6} {:^20} {:>15} {:>15} {:>10} {:^3} {:<10}"
        header = align.format("Bin", "Range", "Numerator", "Denominator", "Eff Value", "+/-", "Eff Error")
        hLine  = "="*100
        nBinsX = denHisto.GetNbinsX()
        table.append("{:^100}".format(hName))
        table.append(hLine)
        table.append(header)
        table.append(hLine)    
        
        # Create the TGraph with asymmetric errors
        tgraph = ROOT.TGraphAsymmErrors(nBinsX,
                                        array.array("d", xVal),
                                        array.array("d", yVal),
                                        array.array("d", xErrL),
                                        array.array("d", xErrH),
                                        array.array("d", yErrL),
                                        array.array("d", yErrH)
                                        )
        # Set the correct name
        tgraph.SetName(hName)
        tgraph.SetTitle(hTitle)

        # Append to list
        gList.append(tgraph)

        # Construct info table (debugging)
        table  = []
        align  = "{:>3} {:>10} {:^10} {:<10} {:>10} {:^10} {:<10}"
        header = align.format("#", "pT-low", "pT", "pT-high",  "eff-low", "eff", "eff-high")
        hLine  = "="*70
        table.append("")
        table.append(hLine)
        table.append("{:^70}".format(hName))
        table.append(header)
        table.append(hLine)
        
        # For-loop: All values x-y and their errors
        for i, x in enumerate(xVal, 0):
            xAvg = xVal[i]
            xMin = xAvg - xErrL[i]
            xMax = xAvg + xErrH[i]
            yAvg = yVal[i]
            yMin = yAvg - yErrL[i]
            yMax = yAvg + yErrH[i]
            row = align.format(i+1, "%.1f" % xMin,  "%.1f" %  xAvg, "%.1f" %  xMax, "%.3f" %  yMin, "%.3f" %  yAvg, "%.3f" %  yMax)
            table.append(row)
        table.append(hLine)

        if printValues:
            for i, line in enumerate(table, 1):
                Print(line, False)    

    return gList


def GetEfficiencyHisto(numHisto, denHisto, hName, hTitle, printValues=False):
    
    # Define histos here
    hEff = numHisto.Clone(hName)
    hEff.Reset()
    hEff.SetTitle(hTitle)
    hNum = ROOT.TH1D('Numerator'  ,'Numerator'  , 1, 0, 1)
    hDen = ROOT.TH1D('Denominator','Denominator', 1, 0, 1)

    # Construct info table (debugging)
    table  = []
    align  = "{:>6} {:>10} {:>10} {:>10} {:>12} {:>12} {:>10} {:^3} {:<10}"
    header = align.format("#", "pT-low", "pT", "pT-high", "Numerator", "Denominator", "Eff", "+/-", "Error")
    hLine  = "="*100
    nBinsX = denHisto.GetNbinsX()
    table.append("{:^100}".format(hName))
    table.append(hLine)
    table.append(header)
    table.append(hLine)

    # First determine the bin number
    for j in range(0, nBinsX+2):
        
        # Get the denominator
        denValue   = ROOT.Double(0.0)
        denError   = ROOT.Double(0.0)
        denValue   = denHisto.IntegralAndError(j, j, denError)
        ROOT.gErrorIgnoreLevel = ROOT.kFatal #kPrint = 0,  kInfo = 1000, kWarning = 2000, kError = 3000, kBreak = 4000, kSysError = 5000, kFatal = 6000   
        
        # Get the denominator
        numValue    = ROOT.Double(0.0)
        numError    = ROOT.Double(0.0)
        numValue    = numHisto.IntegralAndError(j, j, numError)
        effValue    = 0
        effError    = 0

        # Sanity
        if numValue < 0.0:
            msg = "Integral is less than zero (numValue = %.3f, denValue = %.3f). Setting nunValue to zero!" % (numValue, denValue)
            numValue = 0.0
            Print(es + msg + ns, True)
        if numError < 0.0:
            raise Exception("Error is less than zero (numError = %.3f, denError = %.3f)" % (numError, denError) )
        
        # Numerator and Denominator histos
        Verbose("Evaluating efficiency for \"%s\"" % ("test"), False)
        hNum.SetBinError(1, numError)
        hDen.SetBinError(1, denError)
        hNum.SetBinContent(1, numValue)
        hDen.SetBinContent(1, denValue)

        # Calculate Efficiency
        teff = ROOT.TEfficiency(hNum, hDen)
        teff.SetStatisticOption(ROOT.TEfficiency.kFNormal) #statistic option is 'normal approximation'
        effValue = teff.GetEfficiency(1)
        effError = teff.GetEfficiencyErrorLow(1)

        # Bin-range or overflow bin?
        ptLow  = "%.1f" % (denHisto.GetXaxis().GetBinLowEdge(j))
        ptVal  = "%.1f" % (denHisto.GetXaxis().GetBinCenter(j))
        ptHigh = "%.1f" % (denHisto.GetXaxis().GetBinUpEdge(j) )

        if j >= nBinsX+1:
            ptHigh = "infty"

        # Fill histogram
        hEff.SetBinContent(j, effValue)
        if math.isnan(effError):
            effError = 0.0
        hEff.SetBinError(j, effError)

        # Save information in table
        row = align.format(j, ptLow, ptVal, ptHigh, "%.1f" % numValue, "%.1f" % denValue, "%.2f" % effValue, "+/-", "%.2f" % effError)
        table.append(row)

        # Reset histos 
        hNum.Reset()
        hDen.Reset()               

    # Finalise table
    table.append(hLine)

    # Print purity as function of final shape bins
    if printValues:
        for i, line in enumerate(table):
            Print(line, i==0)

    ROOT.gErrorIgnoreLevel = ROOT.kWarning #kPrint = 0,  kInfo = 1000, kWarning = 2000, kError = 3000, kBreak = 4000, kSysError = 5000, kFatal = 6000   
    return hEff


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
    ANALYSISNAME = "TauFakeRate" #Hplus2hwAnalysisWithTop"
    SEARCHMODE   = "350to3000"
    DATAERA      = "Run2016"
    OPTMODE      = ""
    BATCHMODE    = True
    INTLUMI      = -1.0
    URL          = False
    SAVEDIR      = None
    SAVENAME     = None
    VERBOSE      = False
    FOLDER       = "/"
    NUMHISTO     = None
    DENHISTO     = None
    GRIDX        = False
    GRIDY        = False
    LOGX         = False
    LOGY         = False
    YMIN         = None
    YMAX         = None
    CUTLINEX     = None
    CUTLINEY     = None
    SAVEFORMATS  = "pdf,png,C"
    INDIVIDUALMC = False
    
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

    parser.add_option("--intLumi", dest="intLumi", type=float, default=INTLUMI,
                      help="Override the integrated lumi [default: %s]" % INTLUMI)

    parser.add_option("--denHisto", dest="denHisto", type="string", default=DENHISTO,
                      help="Counter name to use as reference when calculating efficiencies (i.e. start point in the selections) [default: %s]" % DENHISTO)

    parser.add_option("--numHisto", dest="numHisto", type="string", default=NUMHISTO,
                      help="Counter names to use when calculating efficiencies (Use comma-separated values) [default: %s]" % NUMHISTO)

    parser.add_option("--searchMode", dest="searchMode", type="string", default=SEARCHMODE,
                      help="Override default searchMode [default: %s]" % SEARCHMODE)

    parser.add_option("--yMin", dest="yMin", type=float, default=YMIN,
                      help="Minimum value of y-axis [default: %s]" % YMIN)

    parser.add_option("--yMax", dest="yMax", type=float, default=YMAX,
                      help="Maxmimum value of y-axis [default: %s]" % YMAX)

    parser.add_option("--gridX", dest="gridX", action="store_true", default=GRIDX,
                      help="Enable the grid in x-axis [default: %s]" % GRIDX)

    parser.add_option("--gridY", dest="gridY", action="store_true", default=GRIDY,
                      help="Enable the grid in y-axis [default: %s]" % GRIDY)

    parser.add_option("--logX", dest="logX", action="store_true", default=LOGX,
                      help="Enable the log scale in x-axis [default: %s]" % LOGX)

    parser.add_option("--logY", dest="logY", action="store_true", default=LOGY,
                      help="Enable the log scale in y-axis [default: %s]" % LOGY)

    parser.add_option("-s", "--saveFormats", dest="saveFormats", default = SAVEFORMATS,
                      help="Save formats for all plots [default: %s]" % SAVEFORMATS)

    parser.add_option("--dataEra", dest="dataEra", type="string", default=DATAERA, 
                      help="Override default dataEra [default: %s]" % DATAERA)

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR, 
                      help="Directory where all pltos will be saved [default: %s]" % SAVEDIR)

    parser.add_option("--saveName", dest="saveName", type="string", default=SAVENAME, 
                      help="The Name of the histogram as it will be saved [default: %s]" % SAVENAME)

    parser.add_option("--url", dest="url", action="store_true", default=URL, 
                      help="Don't print the actual save path the histogram is saved, but print the URL instead [default: %s]" % URL)
    
    parser.add_option("--cutLineX", dest="cutLineX", type="float", default=CUTLINEX,
                      help="X-axis value for vertical cut line [default: %s]" % (CUTLINEX) )

    parser.add_option("--cutLineY", dest="cutLineY", type="float", default=CUTLINEY,
                      help="Y-axis value for horizontal cut line [default: %s]" % (CUTLINEY) )

    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=VERBOSE, 
                      help="Enables verbose mode (for debugging purposes) [default: %s]" % VERBOSE)

    parser.add_option("-i", "--includeOnlyTasks", dest="includeOnlyTasks", action="store", 
                      help="List of datasets in mcrab to include")

    parser.add_option("-e", "--excludeTasks", dest="excludeTasks", action="store", 
                      help="List of datasets in mcrab to exclude")

    parser.add_option("--folder", dest="folder", type="string", default = FOLDER,
                      help="ROOT file folder under which all histograms to be plotted are located [default: %s]" % (FOLDER) )

    parser.add_option("--individualMC", dest="individualMC", action="store_true", default=INDIVIDUALMC, 
                      help="Include in the Fake Rate comparison plot the individual MC samples (in addition to \"Simulation\"\)? [default: %s]" % INDIVIDUALMC)


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
        
    # Define save directory
    if opts.saveDir == None:
        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="", postfix="")
        
    if opts.numHisto == None:
        msg = "Please provide counters with the --counters option."
        raise Exception(es + msg + ns)

    if opts.denHisto == None:
        msg = "The reference counter \"%s\" is not supported. If you are certain it exists in the ROOT files then add it to the list and re-rerun." % opts.denHisto
        raise Exception(es + msg + ns)

    # Set name for saving
    opts.saveName = "TauFR_%s" % (opts.denHisto.split("_")[-1])

    # Create save formats
    if "," in opts.saveFormats:
        opts.saveFormats = opts.saveFormats.split(",")
    else:
        opts.saveFormats = [opts.saveFormats]
    opts.saveFormats = ["." + s for s in opts.saveFormats]          

    # Call the main function
    main(opts)

    if not opts.batchMode:
        raw_input("=== plotTauFakeRates.py: Press any key to quit ROOT ...")
