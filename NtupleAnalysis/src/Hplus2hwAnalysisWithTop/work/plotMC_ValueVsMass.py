#!/usr/bin/env python
'''
DESCRIPTION:


USAGE:
./plotMC_ValueVsMass.py -m <pseudo_mcrab_directory> [opts]


EXAMPLES:
./plotMC_ValueVsMass.py --logY --yMin 1e-3 --yMax 1.2 --refCounter "ttree: skimCounterAll" -m Hplus2hwAnalysisWithTop_MediumTauID_BDTG0p40_tauLdgTrkPt20_30July2019 --url
./plotMC_ValueVsMass.py --logY --yMin 0.3e-2 --yMax 1.2 --refCounter "passed trigger" -m Hplus2hwAnalysisWithTop_MediumTauID_BDTG0p40_tauLdgTrkPt20_30July2019 --url
./plotMC_ValueVsMass.py --logY --yMin 1e-4 --yMax 1.2 --refCounter "ttree: skimCounterAll" -m Hplus2hwAnalysisWithTop_MediumTauID_BDTG0p40_tauLdgTrkPt20_30July2019 


LAST USED:
./plotMC_ValueVsMass.py --logY --yMin 1e-3 --yMax 1.2 --refCounter "passed trigger" -m Hplus2hwAnalysisWithTop_MediumTauID_BDTG0p40_tauLdgTrkPt20_30July2019
./plotMC_ValueVsMass.py --logY --yMin 1e-1 --yMax 1.2 --refCounter "ttree: skimCounterAll" -m Hplus2hwAnalysisWithTop_MediumTauID_BDTG0p40_tauLdgTrkPt20_30July2019 --counters "ttree: skimCounterAll,ttree: skimCounterPassed,Base::AllEvents,Base::PUReweighting,Base::Prescale,Base::Weighted events with top pT,Base::Weighted events for exclusive samples,all events,passed trigger" --refCounter "ttree: skimCounterAll"

./plotMC_ValueVsMass.py --yMin 0.0 --yMax 1.05 --refCounter "All events" --counters "Passed trigger matching,Passed e discr,Passed mu discr,Passed pt cut,Passed eta cut,Passed ldg.trk pt cut,Passed nprongs,Passed isolation" --counterName "tau selection ()" -m Hplus2hwAnalysisWithTop_MediumTauID_BDTG0p40_tauLdgTrkPt20_30July2019

#./plotMC_ValueVsMass.py --yMin 0.0 --yMax 1.05 --refCounter "All events" --counters "Passed trigger matching,Passed decay mode,Passed generic discriminators,Passed e discr,Passed mu discr,Passed pt cut,Passed eta cut,Passed ldg.trk pt cut,Passed nprongs,Passed isolation" --counterName "tau selection ()" -m Hplus2hwAnalysisWithTop_MediumTauID_BDTG0p40_tauLdgTrkPt20_30July2019


USE WITH:
hplusPrintCounters.py --mainCounterOnly --dataEra "Run2016" --mergeData  --includeTasks "2016|M_500" Hplus2tbAnalysis_TopMassLE400_BDT0p40_Binning4Eta5Pt_Syst_NoTopPtReweightCorrXML_10Jan2019

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
    style.setLogY(opts.logY)
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
        plots.mergeRenameReorderForDataMC(datasetsMgr) 

        if opts.intLumi < 0:
            if "Data" in datasetsMgr.getAllDatasetNames():
                opts.intLumi = datasetsMgr.getDataset("Data").getLuminosity()
                Verbose("Automatically determined the integrated luminosity to be %.2f 1/fb" % (opts.intLumi), True)
            else:
                opts.intLumi = 1.0
                Verbose("Could not automatically determine the integrated luminosity. Set value to %.2f 1/fb" % (opts.intLumi), True)

        # Set/Overwrite cross-sections
        for d in datasetsMgr.getAllDatasetNames():
            if "ChargedHiggs" in d:
                xs = 1.0
                Verbose("Setting cross-section for dataset to %.2f" % xs, True)
                datasetsMgr.getDataset(d).setCrossSection(xs) # HIG-18-014: mH=200-3000 GeV sigmaBR = 0.65->0.005 (pb)
            else:
                Verbose("Removing dataset %s" % d, True)
                datasetsMgr.remove(d)
            
        # Print final dataset summary
        if 0:
            datasetsMgr.PrintInfo()

        # Get the efficiency graphs (Efficiency Vs Cut-flow)
        hGraphList = []
        histoName  = os.path.join(opts.folder, opts.counterName)
        hGraphList, _kwargs = GetHistoGraphs(datasetsMgr, opts.folder, histoName)

        # Plot the histo graphs
        PlotHistoGraphs(hGraphList, _kwargs)

    Print("All plots saved under directory %s" % (ts + aux.convertToURL(opts.saveDir, opts.url) + ns), True)
    return


def GetHistoGraphs(datasetsMgr, folder, hName):

    # Get histogram customisations
    _kwargs  = GetHistoKwargs(opts)

    # Get histos (Data, EWK) for Inclusive
    p1 = plots.MCPlot(datasetsMgr, hName, normalizeToLumi=True, saveFormats=[])
    kwargs = copy.deepcopy(_kwargs)
    kwargs["opts"]       = {"ymin": 1.0, "ymaxfactor": 10.0}
    kwargs["xlabelsize"] = 10
    kwargs["ratio"]      = False
    plots.drawPlot(p1, hName, **kwargs) #countes/weighted
    if 0:
        SavePlot(p1, hName, os.path.join(opts.saveDir, opts.optMode), saveFormats = [".C", ".png", ".pdf"] )

    # Clone histograms 
    opts.datasets = opts.signal
    graphList = []
    hgList    = []

    # For-loop: All bins (i.e. counters) to consider
    for index, binLabel in enumerate(opts.myCounters, 0):

        # Convert histos to TGraph
        tgraph = convertHisto2TGraph(p1, binLabel, printValues=True)
        
        # Apply random histogram styles
        s = GetStyle(index)
        #s = styles.styles[index] 
        s.apply(tgraph)
            
        # Append to list
        graphList.append(tgraph)
        
        # Create histoGraph object
        htgraph = histograms.HistoGraph( tgraph, GetLegendNameForBinLabel(binLabel), "LP", "LP")
                
        # Append to list
        hgList.append(htgraph)

    return hgList, _kwargs


def GetStyle(index):

    myStyles = []
    myStyles.append(styles.wwStyle)
    myStyles.append(styles.dyStyle)
    myStyles.append(styles.ewkStyle)
    myStyles.append(styles.qcdBEnrichedStyle)
    myStyles.append(styles.dibStyle)
    myStyles.append(styles.singleTopStyle)
    #myStyles.append(styles.wjetsStyle)
    myStyles.append(styles.Style(ROOT.kFullCircle, ROOT.kAzure+8))
    myStyles.append(styles.ttStyle)
    myStyles.append(styles.qcdStyle)
    myStyles.append(styles.stsStyle)
    myStyles.append(styles.ttwStyle)
    myStyles.append(styles.sttwStyle)
    myStyles.append(styles.wzStyle)
    myStyles.append(styles.ttjetsStyle)
    myStyles.append(styles.zjetsStyle)
    myStyles.append(styles.stStyle)
    myStyles.append(styles.ttHTobbStyle)
    myStyles.append(styles.ttHStyle)
    myStyles.append(styles.ttXStyle)
    myStyles.append(styles.ttbbStyle)
    myStyles.append(styles.ttttStyle)
    myStyles.append(styles.sttStyle)
    myStyles.append(styles.ttHToGGStyle)
    myStyles.append(styles.ttHToNonbbStyle)
    myStyles.append(styles.ttHToTTStyle)
    myStyles.append(styles.ttzStyle)
    myStyles.append(styles.wStyle)
    myStyles.append(styles.zzStyle)
    #myStyles.reverse()
    return myStyles[index]

def GetLegendNameForBinLabel(binLabel):
    binDict   = {
        "passed trigger"                     : "trigger", 
        "passed PV"                          : "PV", 
        "passed e selection (Veto)"          : "e veto", 
        "passed mu selection (Veto)"         : "#mu veto",
        "Passed tau selection (Veto)"        : "#tau veto", 
        "passed jet selection ()"            : "#geq 7 jets", 
        "passed b-jet selection ()"          : "#geq 3 b-jets",
        "b-tag SF"                           : "b-tag SF",
        "top-tag SF"                         : "top SF", 
        "passed top selection ()"            : "2 tops",
        "ttree: skimCounterAll"              : "skim (all)",
        "ttree: skimCounterPassed"           : "skim (pass)",
        "Base::AllEvents"                    : "base (all)",
        "Base::PUReweighting"                : "base (PU rew.)",
        "Base::Prescale"                     : "base (prescale)",
        "Base::Weighted events with top pT"  : "base (top p_{T})",
        "Base::Weighted events for exclusive samples": "base (excl.)",
        "all events"                         : "all",
        "passed METFilter selection ()"      : "filters",
        "passed mu selection ()"             : "1 #mu",
        "Passed tau selection ()"            : "#geq 1 #tau jets", 
        "Passed tau selection and genuine ()": "genuine #tau",
        "#tau N"                             : "2 #tau jets",
        "passed jet selection ()"            : "#geq 3 jets",
        "passed angular cuts Collinear ()"   : "R_{coll}",
        "passed b-jet selection ()"          : "#geq 1 b jets",
        "b-tag SF"                           : "b jets SF",
        "passed angular cuts BackToBack ()"  : "R_{bb}",
        "passed top selection ()"            : "1 top",
        "top-tag SF"                         : "top SF",
        "Selected Events"                    : "selected",
        #"#tau OS"
        #"#tau SF":
        #"Fake #tau SF":
        #"#tau DM": 
        "All events"                         : "all", 
        "Passed trigger matching"            : "trg-m", 
        "Passed decay mode"                  : "DM", 
        "Passed generic discriminators"      : "discr.", 
        "Passed e discr"                     : "e discr.", 
        "Passed mu discr"                    : "#mu discr.", 
        "Passed pt cut"                      : "p_{T}", 
        "Passed eta cut"                     : "#eta", 
        "Passed ldg.trk pt cut"              : "p^{ldg tk}_{T}", 
        "Passed nprongs"                     : "n-prongs", 
        "Passed isolation"                   : "isolation", 
        "Passed Rtau"                        : "R_{#tau}", 
        "Passed anti-isolation"              : "anti-isolation", 
        #"Passed anti-isolated Rtau"          : "", 
        #"Passed tau selection and genuine"   : "", 
        #"multiple selected taus"             : "", 
        #"Passed anti-isolated tau selection" : "", 
        #"multiple anti-isolated taus"        : "", 
        }
    if binLabel in binDict.keys():
        return binDict[binLabel]
    else:
        return binLabel


def PlotHistoGraphs(hGraphList, _kwargs):

    # Create & draw the plot    
    p = plots.ComparisonManyPlot(hGraphList[0], hGraphList[1:], saveFormats=[])
    p.setLuminosity(opts.intLumi)

    if 0:
        p.histoMgr.setHistoLegendStyleAll("P")
        p.histoMgr.setHistoDrawStyleAll("LP")
        p.histoMgr.setHistoDrawStyle(hGraphList[0].getRootHisto().GetName(), "HIST")

    if opts.saveName == None:
        opts.saveName = hGraphList[0].getRootHisto().GetName().split("_")[0] + "VsMass"

    # Draw the plot
    if 0:
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetLineColor(ROOT.kRed))
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetLineStyle(ROOT.kDotted))
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetLineWidth(3))
        #
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetMarkerColor(ROOT.kBlue))
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetMarkerSize(1.3))
        p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetMarkerStyle(ROOT.kFullCircle))
        
    for i, g in enumerate(hGraphList, 1):
        if 0: #i == 0:
            g.getRootHisto().SetMarkerStyle(0)
        g.getRootHisto().SetLineStyle(i)

    #p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetLineStyle(ROOT.kDotted))
    p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetLineWidth(3))
    p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetMarkerSize(1.3))
    plots.drawPlot(p, opts.saveName, **_kwargs)

    # Save the plot
    SavePlot(p, opts.saveName, os.path.join(opts.saveDir, opts.optMode), saveFormats = [".C", ".png", ".pdf"] )
    return


def GetHistoKwargs(opts):
    nCounters = len(opts.myCounters)
    if nCounters < 8:
        dh = 0
    else:
        dh = 0.02*(nCounters-7)

    _kwargs     = {
        "xlabel"           : "m_{H^{#pm}} (GeV)",
        "ylabel"           : "Efficiency",
        "ratioYlabel"      : "Ratio ",
        "ratio"            : False,
        "ratioInvert"      : False,
        "stackMCHistograms": False,
        "addMCUncertainty" : True,
        "addLuminosityText": True,
        "addCmsText"       : True,
        "errorBarsX"       : True,
        "cmsExtraText"     : "Simulation",
        "opts"             : {"ymin": opts.yMin, "ymax": opts.yMax},
        "opts2"            : {"ymin": 1e-3, "ymax": 10},
        "moveLegend"       : {"dx": -0.1 , "dy": -0.07, "dh": +dh},
        "cutBox"           : {"cutValue": 800.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True},
        "cutBoxY"          : {"cutValue": 1.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True, "mainCanvas": True, "ratioCanvas": False},
        #"xlabelsize"       : 0
        }

    # Set x-axis divisions
    n1 = 10 # primary divisions (8)
    n2 = 0 # second order divisions (5)
    n3 = 0 # third order divisions (@)
    nDivs = n1 + 100*n2 + 10000*n3
    ROOT.gStyle.SetNdivisions(nDivs, "X")
    return _kwargs


def GetSaveName(histoName):
    base = histoName.split("_")[0]
    var  = histoName.split("_")[1]
    sel  = histoName.split("_")[2]
    name = var + "_" + GetControlRegionLabel(histoName)
    return name


def SavePlot(plot, plotName, saveDir, saveFormats = [".C", ".png", ".pdf"]):
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    # Create the name under which plot will be saved
    saveName = os.path.join(saveDir, plotName.replace("/", "_"))

    # For-loop: All save formats
    for i, ext in enumerate(saveFormats):
        saveNameURL = saveName + ext
        saveNameURL = aux.convertToURL(saveNameURL, opts.url)
        Print(saveNameURL, i==0)
        plot.saveAs(saveName, formats=saveFormats)
    return


def convertHisto2TGraph(p1, binLabel, printValues=False):

    # Lists for values
    xVal   = []
    yVal   = []
    xErrL  = []
    xErrH  = []
    yErrL  = []
    yErrH  = []
    nBinsX = len(opts.datasets)
    binNum = -1 
    
    # For-loop: All histogram bins
    for d in opts.datasets:
        h = p1.histoMgr.getHisto(d).getRootHisto().Clone(d)
        mass  = float(d.split("_M")[1].split("_")[0])
        histo = GetEfficiencyHisto(h, opts.refCounter, printValues=True, hideZeros=True)
        
        # For-loop: All bins to find the correct on
        for j in range (0,  histo.GetNbinsX()+1):
            
            # Skip this counter if label is not found
            if histo.GetXaxis().GetBinLabel(j) == binLabel:
                binNum = j
            Verbose("Efficiency for %d counter \"%s\" is %s" % (j, histo.GetXaxis().GetBinLabel(j), histo.GetBinContent(j)), j==0)

        if binNum == -1:
            raise Exception("Could not find counter \"%s\"." % (binLabel) )

        # Get values
        xVal.append(mass)
        xErrL.append(0.0001)
        xErrH.append(0.0001)
        yVal.append(histo.GetBinContent(binNum))
        yErrL.append(histo.GetBinError(binNum))
        yErrH.append(histo.GetBinError(binNum))
        del histo

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
    tgraph.SetName(binLabel)

    # Construct info table (debugging)
    table  = []
    align  = "{:>7} {:>10} {:>15} {:^3} {:<15}"
    header = align.format("#", "Mass", "Efficiency", "+/-", "Error")
    hLine  = "="*50
    table.append("")
    table.append(hLine)
    table.append("{:^50}".format(binLabel))
    table.append(header)
    table.append(hLine)
    
    # For-loop: All values x-y and their errors
    for i, x in enumerate(xVal, 0):
        row = align.format(i+1, "%.0f" %  xVal[i], "%.3f" %  yVal[i], "+/-", "%.3f" %  yErrH[i])
        table.append(row)
    table.append(hLine)

    if printValues:
        for i, line in enumerate(table, 1):
            Print(line, False)
    return tgraph


def GetEfficiencyHisto(histo, refCounter, printValues=False, hideZeros=True):
    
    # Define histos here
    xMax = histo.GetXaxis().GetXmax()
    hEff = ROOT.TH1D('Eff','Eff', int(xMax), 0, xMax)
    hNum = ROOT.TH1D('Num','Num', 1, 0, 1)
    hDen = ROOT.TH1D('Den','Den', 1, 0, 1)

    # Construct info table (debugging)
    table  = []
    align  = "{:>6} {:^20} {:<50} {:>15} {:>15} {:>10} {:^3} {:<10}"
    header = align.format("Bin", "Range", "Selection", "Numerator", "Denominator", "Eff Value", "+/-", "Eff Error")
    hLine  = "="*120
    nBinsX = histo.GetNbinsX()
    table.append("{:^100}".format(histo.GetName()))
    table.append(hLine)
    table.append(header)
    table.append(hLine)

    # First determine the bin number
    binNumber = -1
    for j in range (0, nBinsX+1):

        # Skip this counter?
        binLabel = histo.GetXaxis().GetBinLabel(j)

        if binLabel == refCounter:
            binNumber = j
        else:
            continue

    if binNumber == -1:
        raise Exception("Could not find reference counter \"%s\"" % refCounter)

    if binNumber > nBinsX:
        raise Exception("Invalid bin selected (bin = %d" % (binNumber) ) 

    # First get numerator
    denValue = ROOT.Double(0.0)
    denError = ROOT.Double(0.0)
    denValue = histo.IntegralAndError(binNumber, binNumber, denError)
    ROOT.gErrorIgnoreLevel = ROOT.kFatal #kPrint = 0,  kInfo = 1000, kWarning = 2000, kError = 3000, kBreak = 4000, kSysError = 5000, kFatal = 6000   

    # For-loop: All histogram bins
    binCounter  = 0
    skipCounter = 0
    for j in range (binNumber, nBinsX+1):

        # Skip this counter?
        binLabel = histo.GetXaxis().GetBinLabel(j)
        if binLabel in opts.skipList:
            Print("Skipping counter with name \"%s\"" % (hs + binLabel + ns), skipCounter==0)
            skipCounter += 1
            continue
        
        # Declare variables
        numValue    = ROOT.Double(0.0)
        numError    = ROOT.Double(0.0)
        numValue    = histo.IntegralAndError(j, j, numError)
        effValue    = 0
        effError    = 0

        # Sanity
        if numValue < 0.0:
            raise Exception("Integral is less than zero!")
        if numError < 0.0:
            raise Exception("Error is less than zero!")
        
        # Numerator and Denominator histos
        Verbose("Evaluating efficiency for \"%s\"" % (binLabel), j==binNumber)
        hNum.SetBinContent(1, numValue)
        hNum.SetBinError(1, numError)
        hDen.SetBinContent(1, denValue)
        hDen.SetBinError(1, denError)

        # Calculate Efficiency
        teff = ROOT.TEfficiency(hNum, hDen)
        teff.SetStatisticOption(ROOT.TEfficiency.kFNormal) #statistic option is 'normal approximation'
        effValue = teff.GetEfficiency(1)
        effError = teff.GetEfficiencyErrorLow(1)

        # Bin-range or overflow bin?
        binRange = "%.1f -> %.1f" % (histo.GetXaxis().GetBinLowEdge(j), histo.GetXaxis().GetBinUpEdge(j) )
        if j >= nBinsX:
            binRange = "> %.1f"   % (histo.GetXaxis().GetBinLowEdge(j) )

        # Fill histogram
        binCounter+= 1
        hEff.SetBinContent(binCounter, effValue)
        if math.isnan(effError):
            effError = 0.0
        hEff.SetBinError(binCounter, effError)
        hEff.GetXaxis().SetBinLabel(binCounter, binLabel) #iro

        # Save information in table
        row = align.format(binCounter, binRange, binLabel, "%.1f" % numValue, "%.1f" % denValue, "%.3f" % effValue, "+/-", "%.3f" % effError)
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


def GetEventsHisto(histo, refCounter, kwargs, printValues=False, hideZeros=True):

    # Define the histo
    xMin  = 0
    xMax  = histo.GetXaxis().GetXmax()
    nBins = int(xMax)
    hEvts = ROOT.TH1D('Evts','Evts', nBins, xMin, xMax)

    # Construct info table (debugging)
    table  = []
    align  = "{:>6} {:^20} {:<35} {:>10} {:^3} {:<10}"
    header = align.format("Bin", "Range", "Selection", "Events", "+/-", "Stat.")
    hLine  = "="*120
    nBinsX = histo.GetNbinsX()
    table.append("{:^100}".format(histo.GetName()))
    table.append(hLine)
    table.append(header)
    table.append(hLine)

    # First determine the bin number of the reference counter
    binNumber = -1
    for j in range (0, nBinsX+1):

        # Skip this counter?
        binLabel = histo.GetXaxis().GetBinLabel(j)

        if binLabel == refCounter:
            binNumber = j
        else:
            continue

    if binNumber == -1:
        raise Exception("Could not find reference counter \"%s\"" % refCounter)

    if binNumber > nBinsX:
        raise Exception("Invalid bin selected (bin = %d" % (binNumber) ) 

    # For-loop: All histogram bins
    binCounter = 0
    for j in range (binNumber, nBinsX+1):

        # Skip this counter?
        binLabel = histo.GetXaxis().GetBinLabel(j)
        if binLabel in opts.skipList:
            Verbose("Skipping counter with name \"%s\"" % (binLabel), True)
            continue

        # Declare variables
        numValue    = ROOT.Double(0.0)
        numError    = ROOT.Double(0.0)
        numValue    = histo.IntegralAndError(j, j, numError)

        # Sanity
        if numValue < 0.0:
            raise Exception("Integral is less than zero!")
        if numError < 0.0:
            raise Exception("Error is less than zero!")
        
        # Bin-range or overflow bin?
        binRange = "%.1f -> %.1f" % (histo.GetXaxis().GetBinLowEdge(j), histo.GetXaxis().GetBinUpEdge(j) )
        if j >= nBinsX:
            binRange = "> %.1f"   % (histo.GetXaxis().GetBinLowEdge(j) )

        # Fill histogram
        binCounter+= 1
        hEvts.SetBinContent(binCounter, numValue)
        hEvts.SetBinError(binCounter, numError)

        # Save information in table
        row = align.format(binCounter, binRange, binLabel, "%.5f" % numValue, "+/-", "%.5f" % numError)
        table.append(row)

    # Finalise table
    table.append(hLine)

    # Print purity as function of final shape bins
    if printValues:
        for i, line in enumerate(table):
            Print(line, i==0)

    ROOT.gErrorIgnoreLevel = ROOT.kWarning #kPrint = 0,  kInfo = 1000, kWarning = 2000, kError = 3000, kBreak = 4000, kSysError = 5000, kFatal = 6000   
    return hEvts


#================================================================================================ 
# Main
#================================================================================================ 
if __name__ == "__main__":
    '''g1

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
    OPTMODE      = ""
    BATCHMODE    = True
    INTLUMI      = -1.0
    URL          = False
    SAVEDIR      = None
    SAVENAME     = None
    VERBOSE      = False
    FOLDER       = "counters/weighted/" #"counters"
    COUNTERNAME  = "counter" #"tau selection ()"
    SIGNALMASS   = [300, 700]
    REFCOUNTER   = "passed trigger" #ttree: skimCounterAll"
    SKIPLIST     = ["passed METFilter selection ()", "passed PV", "Passed tau selection and genuine ()", "#tau SF", "Fake #tau SF",  
                    "passed angular cuts Collinear ()", "b-tag SF", "passed angular cuts BackToBack ()", "Selected Events"]#, "top-tag SF"]#
    GRIDX        = False
    GRIDY        = False
    LOGX         = False
    LOGY         = False
    ERRORX       = False
    YMIN         = None
    YMAX         = None
    EFFICIENCY   = False
    COUNTERS     = None
    
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

    parser.add_option("--refCounter", dest="refCounter", type="string", default=REFCOUNTER,
                      help="Counter name to use as reference when calculating efficiencies (i.e. start point in the selections) [default: %s]" % REFCOUNTER)

    parser.add_option("--counters", dest="counters", type="string", default=COUNTERS,
                      help="Counter names to use when calculating efficiencies (Use comma-separated values) [default: %s]" % COUNTERS)

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

    parser.add_option("--errorX", dest="errorX", action="store_true", default=ERRORX,
                      help="Enable the x-axis error bar [default: %s]" % ERRORX)

    parser.add_option("--dataEra", dest="dataEra", type="string", default=DATAERA, 
                      help="Override default dataEra [default: %s]" % DATAERA)

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR, 
                      help="Directory where all pltos will be saved [default: %s]" % SAVEDIR)

    parser.add_option("--saveName", dest="saveName", type="string", default=SAVENAME, 
                      help="The Name of the histogram as it will be saved [default: %s]" % SAVENAME)

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

    parser.add_option("--counterName", dest="counterName", type="string", default = COUNTERNAME,
                      help="Counter to be plotted [default: %s]" % (COUNTERNAME) )

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
        
    # Define counters to skip (if for example empty)
    opts.skipList = SKIPLIST

    # Define save directory
    if opts.saveDir == None:
        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="", postfix="MC")
        
    # Sanity check
    allowedFolders = ["counters", "counters/weighted/"]
    if opts.folder not in allowedFolders:
        Print("Invalid folder \"%s\"! Please select one of the following:" % (opts.folder), True)
        for m in allowedFolders:
            Print(m, False)
        sys.exit()

    # Signal list
    opts.signal     = []
    opts.signalMass = SIGNALMASS
    for m in SIGNALMASS:
        #signal = "ChargedHiggs_HplusTB_HplusToTB_M_%i" % m
        signal = "CRAB_private_ChargedHiggs_HplusTB_HplusToHW_M%i_mH200_2ta_NLO" % m
        signal = "ChargedHiggs_HplusTB_HplusToHW_M%i_mH200_2ta_NLO" % m
        opts.signal.append(signal)

    # Set the counter
    opts.allCounters = ["ttree: skimCounterAll", "ttree: skimCounterPassed", "Base::AllEvents", "Base::PUReweighting", "Base::Prescale", 
                        "Base::Weighted events with top pT", "Base::Weighted events for exclusive samples", "all events", "passed trigger", 
                        "passed METFilter selection ()", "passed PV", "passed e selection (Veto)", "passed mu selection ()", 
                        "Passed tau selection ()", "Passed tau selection and genuine ()", "#tau N", "#tau OS", "#tau SF", "Fake #tau SF"
                        "#tau DM", "passed jet selection ()", "passed angular cuts Collinear ()", "passed b-jet selection ()", "b-tag SF", 
                        "passed angular cuts BackToBack ()", "passed top selection ()"]#, "top-tag SF", "Selected Events"]

    if opts.counters == None:
        # opts.counters = ["passed trigger", "passed e selection (Veto)", "passed mu selection (Veto)", "Passed tau selection (Veto)", "passed jet selection ()", "passed b-jet selection ()", "passed top selection ()"]
        opts.counters = "passed trigger,passed e selection (Veto),passed mu selection (),Passed tau selection (),#tau N,#tau OS,passed jet selection (),passed b-jet selection (),passed top selection ()"
        opts.myCounters = opts.counters.split(",")
    else:
        opts.myCounters = opts.counters.split(",")


    if opts.refCounter not in opts.allCounters:
        msg = "The reference counter \"%s\" is not supported. If you are certain it exists in the ROOT files then add it to the \"myCounters\" list and re-rerun." % opts.refCounter
        #raw_input(es + msg +ns)
        #raise Exception(es + msg + ns)

    # Set name for saving
    opts.saveName = "EffVsMass_%s" % (opts.refCounter.replace(" ", "").replace("(", "").replace(")", "").replace(":", ""))

    # Call the main function
    Print("Calculating the following efficiencies using \"%s\" as reference counter:\n\t%s" % (opts.refCounter, "\n\t".join(opts.myCounters)), True)
    main(opts)

    if not opts.batchMode:
        raw_input("=== plot_ValueVsMass.py: Press any key to quit ROOT ...")
