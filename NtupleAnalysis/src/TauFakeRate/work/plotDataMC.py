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
        if 1:
            datasetsMgr.PrintInfo()
        plots.mergeRenameReorderForDataMC(datasetsMgr, keepSourcesMC=False, analysisType="HToHW_withTop") 
        datasetsMgr.PrintInfo()

        # Set signal cross-section
        newOrder = ["Data", "DYJetsToLL", "TT", "Diboson", "SingleTop", "ttX"]

        # For-loop: All dataset objects
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
            "topSelectionMVA_",
            "config",
            "NSelectedVsRunNumber",
            ]

        # For-loop: All histograms
        for h in histoPaths:
            if "resolution" in h.lower():
                continue
            if "angular" in h.lower():
                continue
            if "backtoback" in h.lower():
                continue
            if "DeltaPhiMuonMet" in h:
                continue
            if "DeltaPhiTauMet" in h:
                continue
            myHistos.append(h)

        # For-loop: All histos
        for i, h in enumerate(myHistos, 1):
            msg   = "{:<9} {:>3} {:<1} {:<3} {:<50}".format("Histogram", "%i" % i, "/", "%s:" % (len(myHistos)), h)
            Print(ss + msg + ns, i==1)
            
            if "NprongsMatrix" in h or "etaphi" in h:
                continue

            if opts.folder == "":
                if h in skipList:
                    continue

            Verbose(h, i==1)
            PlotDataMCHistograms(datasetsMgr, h)

    Print("All plots saved under directory %s" % (ts + aux.convertToURL(opts.saveDir, opts.url) + ns), True)    
    return

def GetHistoKwargs(histoName, opts):
    h = histoName.rsplit("/")[-1]

    # Common bin settings
    _legRM  = {"dx": -10000.23, "dy": -10000.01, "dh": -0.1}
    _legNE  = {"dx": -0.08, "dy": -0.01, "dh": 0.1}
    _legNW  = {"dx": -0.52, "dy": -0.01, "dh": 0.1}
    _legSE  = {"dx": -0.23, "dy": -0.50, "dh": 0.1}
    _legSW  = {"dx": -0.52, "dy": -0.48, "dh": 0.1}
    _legCE  = {"dx": -0.08, "dy": -0.44, "dh": 0.1}
    _logY   = True
    _yLabel = "Events / %.0f "
    if "_AfterAllSelections" in h or "_AfterTopSelection" in h or "_AfterMetSelection" in h:
        _yMin = 1e-2
    else:
        _yMin = 1e-2

    if _logY:
        _yMaxF = 10
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
        #"opts2"            : {"ymin": 0.0, "ymax": 3.0},
        "opts2"            : {"ymin": 0.65, "ymax": 1.35},
        "divideByBinWidth" : False,
        "log"              : _logY,
        "moveLegend"       : _legNE,
        "moveBlindedText"  : {"dx": -0.23, "dy": +0.08, "dh": 0.0},
        }

    kwargs = copy.deepcopy(_kwargs)
    
    if "phi" in h.lower():
        kwargs["moveLegend"] = _legNW

    if "pt" in h.lower():
        kwargs["units"]      = "GeV"
        if "eSelection" in histoName:
            kwargs["rebinX"] = systematics.DataMCBinningHToHW["ElectronPt"]
        if "tauSelection" in histoName:
            kwargs["rebinX"] = systematics.DataMCBinningHToHW["TauPt"]
        kwargs["xlabel"]     = "p_{T} (%s)" % kwargs["units"]
        kwargs["ylabel"]     = _yLabel + kwargs["units"]
        kwargs["moveLegend"] = _legNE
        kwargs["divideByBinWidth"] = True

    if "eta" in h.lower():
        kwargs["rebinX"] = 2
        kwargs["xlabel"] = "#eta"
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["opts"]   = {"xmin": -2.5, "xmax": +2.5, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["moveLegend"] = _legRM #_legNE
        kwargs["divideByBinWidth"] = False

    if "tausrc" in h.lower():
        kwargs["opts"]   = {"xmax": +12, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "DeltaR" in h or "DeltaY" in h or "DR" in h:
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["opts"]   = {"xmax": +8.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "miniiso" in h.lower():
        kwargs["units"]  = ""
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["opts"] = {"xmax": +30.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["rebinX"] = systematics.DataMCBinningHToHW["MiniIso"]
        kwargs["divideByBinWidth"] = True
        if "after" in h.lower() or "passed" in h.lower():
            kwargs["opts"]   = {"xmax": +0.8, "ymin": _yMin, "ymaxfactor": _yMaxF}
            kwargs["rebinX"] = 1
            kwargs["divideByBinWidth"] = False
            
    if "reliso" in h.lower():
        kwargs["units"]  = ""
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["opts"] = {"xmax": +30.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["rebinX"] = systematics.DataMCBinningHToHW["RelIso"]
        kwargs["divideByBinWidth"] = True

    if "bdisc" in h.lower():
        kwargs["rebinX"] = 2
        kwargs["xlabel"] = "b-jet discriminant"
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["cutBox"] = {"cutValue": 0.8484, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["opts"]   = {"xmin": 0.4, "xmax": +1.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["moveLegend"] = _legSW

    if "Nvtx" in h or "IsolVtx" in h or "nvertices" in h.lower():
        kwargs["rebinX"] = systematics.DataMCBinningHToHW["Vertices"]
        kwargs["units"]  = ""
        kwargs["xlabel"] = "vertex multiplicity"
        kwargs["divideByBinWidth"] = True

    if "HT" in h:
        kwargs["rebinX"] = systematics.DataMCBinningHToHW["HT"]
        kwargs["units"]  = "GeV" #"GeV/c"
        kwargs["xlabel"] = "H_{T} (%s)" % kwargs["units"]
        kwargs["divideByBinWidth"] = True

        if "mht" in h.lower():
            kwargs["units"]      = "GeV"
            kwargs["xlabel"]     = "MHT (%s)" % kwargs["units"]
            kwargs["ylabel"]     = _yLabel + kwargs["units"]
            kwargs["moveLegend"] = _legNE
            kwargs["divideByBinWidth"] = True


    if "NTau" in h:
        kwargs["opts"]   = {"xmax": 8.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
    if "Njets" in h:
        kwargs["xlabel"] = "jet multiplicity"
    if "NBjet" in h:
        kwargs["xlabel"] = "b-jet multiplicity"
        kwargs["opts"]   = {"xmax": 8.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "JetPt" in h or "jetPt" in h:        
        if "BJetPt" in h or "BjetPt" in h:
            kwargs["rebinX"] = systematics.DataMCBinningHToHW["BJetPt"]
        else:
            kwargs["rebinX"] = systematics.DataMCBinningHToHW["JetPt"]
        kwargs["units"]  = "GeV" #"GeV/c"
        kwargs["xlabel"] = "p_{T} (%s)" % kwargs["units"]
        kwargs["ylabel"] = _yLabel + kwargs["units"]
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["divideByBinWidth"] = True

    if "counters" in opts.folder:
        ROOT.gStyle.SetLabelSize(16.0, "X")
        kwargs["moveLegend"] = {"dx": -0.08, "dy": 0.0, "dh": 0.1}
        kwargs["moveLegend"] = _legNE
        
    regex = re.compile('selectedJets.*JetPt')
    if(regex.search(h)):
        ROOT.gStyle.SetNdivisions(8, "X")
        kwargs["units"]  = "GeV"
        kwargs["xlabel"] = "p_{T} (%s)" % kwargs["units"]
        kwargs["ylabel"] = _yLabel + kwargs["units"]
        kwargs["cutBox"] = {"cutValue": 30.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}
        kwargs["rebinX"] = systematics.DataMCBinningHToHW["JetPt"]
        kwargs["divideByBinWidth"] = True
        kwargs["opts"] = {"xmin": 30.0, "xmax": +800.0}
        if "second" in h.lower():
            kwargs["opts"] = {"xmin": 30.0, "xmax": +500.0}
        if "third" in h.lower() or "fourth" in h.lower():
            kwargs["opts"] = {"xmin": 30.0, "xmax": +300.0}
        if "fifth" in h.lower() or "sixth" in h.lower():
            kwargs["opts"] = {"xmin": 30.0, "xmax": +150.0}
        
  
    if h == "btagSF":
        kwargs["xlabel"] = "b-jet SF"
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["opts"]   = {"xmin": 0.5, "xmax": 2.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["rebinX"] = 5

    if "met_" in h.lower() or h == "Met" or h == "MET":
        kwargs["units"]  = "GeV"
        kwargs["rebinX"] = systematics.DataMCBinningHToHW["Met"]
        kwargs["xlabel"] = "E_{T}^{miss} (%s)" % (kwargs["units"])
        # kwargs["opts"]   = {"xmax": 500.0, "ymaxfactor": _yMaxF}
        kwargs["divideByBinWidth"] = True        
        kwargs["ylabel"] = _yLabel + kwargs["units"]

        if "METFilter" in h:
            kwargs["rebinX"] = 1
            kwargs["xlabel"] = ""
            kwargs["ylabel"] = "Events / %.0f "
            kwargs["opts"]   = {"xmin": 0.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "deltaphi" in h.lower():
        kwargs["units"]  = "rads"
        #kwargs["rebinX"] = systematics.DataMCBinningHToHW["DeltaPhi"]
        kwargs["rebinX"] = systematics.DataMCBinningHToHW["DeltaPhiRads"]
        kwargs["xlabel"] = "#Delta#phi"
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["opts"]   = {"xmin": 0.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "deltar" in h.lower():
        kwargs["units"]  = ""
        kwargs["rebinX"] = systematics.DataMCBinningHToHW["DeltaR"]
        kwargs["xlabel"] = "#DeltaR"
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["opts"]   = {"xmin": 0.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        
    if "DecayMode" in h:
        kwargs["moveLegend"] = _legRM

    if "Nprong" in h:
        kwargs["moveLegend"] = _legNE
        kwargs["xlabel"]     = "charged particle multiplicity" # "N_{prongs}"

    if "muonN" in h:
        kwargs["xlabel"] = "#mu multiplicity"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 30.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        if "passed" in h.lower():
            kwargs["opts"]   = {"xmin": 0.0, "xmax": 8.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 1.5, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "counters" in opts.folder.lower():
        kwargs["moveLegend"] = _legRM
        kwargs["opts"]       = {"ymin": 1e-1, "ymaxfactor": 10}
        if h == "counter":
            kwargs["moveLegend"] = _legSW
            xMin =  8.0 # start from trigger (8.0)
            xMax = 27.0 # using 1 before because blinding fails there
            kwargs["opts"] = {"xmin": xMin, "xmax": xMax, "ymin": _yMin, "ymax": 5e9}
            #kwargs["opts"] = {"ymin": _yMin, "ymax": 5e9}
            kwargs["moveLegend"] = _legNE
            #kwargs["blindingRangeString"] = "0 to 100"
        elif "jet selection" in h:
            kwargs["opts"] = {"xmin": 0.0, "xmax": 7.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        elif "tau selection" in h:
            kwargs["opts"] = {"xmin": 0.0, "xmax": 11.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        else:
            pass

    if "mass" in h.lower():
        kwargs["units"] = "GeV" #"GeV/c^{2}"
        kwargs["ylabel"] = _yLabel + kwargs["units"]
        kwargs["opts"]   = {"xmax": 200.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["units"]  = "GeV"
        kwargs["xlabel"] = "m_{\\ell\\ell} (%s)" % kwargs["units"]
        kwargs["xlabel"] = None
        kwargs["ylabel"] = _yLabel + kwargs["units"]


    if kwargs["divideByBinWidth"]:
        if kwargs["units"] != "":
            kwargs["ylabel"] = "< Events / %s  > " % kwargs["units"]
        else:
            kwargs["ylabel"] = "< Events > "
    
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

def PlotDataMCHistograms(datasetsMgr, histoName):
    skipStrings = []
    if "ForDataDrivenCtrlPlots" in opts.folder:
        skipStrings = ["_Vs_", "JetEtaPhi"]

    # Skip histograms if they contain a given string
    for keyword in skipStrings:
        if keyword in histoName:
            Verbose("Skipping \"%s\" due to keyword \"%s\"" % (histoName, keyword), True)
            return
        else:
            pass

    # Get Histogram name and its kwargs
    kwargs_  = GetHistoKwargs(histoName, opts)
    saveName = histoName.rsplit("/")[-1]

    # Create the plotting object
    p = plots.DataMCPlot(datasetsMgr, histoName, saveFormats=[])

    # Overwite signal style?
    if 0:
        if opts.signalMass != 0:
            p.histoMgr.forHisto(opts.signal, styles.getSignalStyleHToTB_M(opts.signalMass))

    # Apply blinding of signal region
    if "blindingRangeString" in kwargs_:
        startBlind = float(kwargs_["blindingRangeString"].split(" to ")[1])
        endBlind   = float(kwargs_["blindingRangeString"].split(" to ")[0])
        plots.partiallyBlind(p, maxShownValue=startBlind, minShownValue=endBlind, invert=True, moveBlindedText=kwargs_["moveBlindedText"])

    # Draw and save the plot
    plots.drawPlot(p, saveName, **kwargs_) #the "**" unpacks the kwargs_ dictionary

    # Save the plots in custom list of saveFormats
    SavePlot(p, saveName, os.path.join(opts.saveDir, opts.optMode, opts.folder), opts.saveFormats)
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
    ANALYSISNAME = "TauFakeRate" #Hplus2hwAnalysisWithTop"
    SEARCHMODE   = "350to3000"
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
