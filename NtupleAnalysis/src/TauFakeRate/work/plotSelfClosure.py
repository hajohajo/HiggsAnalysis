#!/usr/bin/env python
'''
DESCRIPTIONM:
Script that plots Data/MC for all histograms under a given folder (passsed as option to the script)
Good for sanity checks for key points in the cut-flow


USAGE:
./plotSelfClosure.py -m <pseudo_mcrab_directory> [opts]


EXAMPLES:
/plotSelfClosure.py -s png --ratio --histos "NTau,NJet,NBjet,Met,Lt,DileptonPt,DileptonEta,DileptonMass" -m TauFakeRate_NewCommonPlots_Attempt1_2MuonsPt26_05Mar2020


LAST USED:
./plotSelfClosure.py -s png --ratio --histos "NTau,NJet,NBjet,Met,Lt,DileptonPt,DileptonEta,DileptonMass" -m TauFakeRate_NewCommonPlots_Attempt1_2MuonsPt26_05Mar2020

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
        plots.mergeRenameReorderForDataMC(datasetsMgr, keepSourcesMC=False, analysisType="HToHW_withTop") 
        if 0:
            datasetsMgr.PrintInfo()
        opts.intLumi = datasetsMgr.getDataset("Data").getLuminosity()

        if opts.verbose:
            datasetsMgr.PrintCrossSections()
            datasetsMgr.PrintLuminosities()
            datasetsMgr.PrintInfo()

        # Overwrite dataset order
        newOrder = ["Data", "DYJetsToLL", "TT", "Diboson", "SingleTop", "ttX"]
        datasetsMgr.selectAndReorder(newOrder)
          
        # Print dataset information
        datasetsMgr.PrintInfo()

        # Do Data-MC histograms
        folder     = opts.folder
        histoList  = datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDirectoryContent(folder)
        histoPaths = [os.path.join(folder, h) for h in histoList]
        
        for i in range(0, len(opts.histoDict["tight"]) ):
            Print("%d/%d: %s" % (i+1, len(opts.histoDict["tight"]), opts.histoDict["tight"][i]), i==0)
            PlotDataMCHistograms(datasetsMgr, [opts.histoDict["tight"][i], opts.histoDict["loose"][i]] )

    Print("All plots saved under directory %s" % (ts + aux.convertToURL(opts.saveDir, opts.url) + ns), True)    
    return

def GetHistoKwargs(histoName, opts):
    h = histoName.rsplit("/")[-1]

    # Common bin settings
    _legRM  = {"dx": -10000.23, "dy": -10000.01, "dh": -0.1}
    _legNE  = {"dx": -0.08, "dy": -0.01, "dh": -0.15}
    _legNW  = {"dx": -0.52, "dy": -0.01, "dh": -0.15}
    _legSE  = {"dx": -0.23, "dy": -0.50, "dh": -0.15}
    _legSW  = {"dx": -0.52, "dy": -0.48, "dh": -0.15}
    _legCE  = {"dx": -0.08, "dy": -0.44, "dh": -0.15}
    _logY   = True
    _yLabel = "Events / %.0f "
    _yMin   = 1e-1
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
        "ratioYlabel"      : "Data/Bkg  ",
        "ratio"            : opts.ratio, 
        "stackMCHistograms": True,
        "ratioInvert"      : False, 
        "addMCUncertainty" : True,
        "addLuminosityText": True,
        "addCmText"        : True,
        "cmsExtraText"     : "Preliminary",
        "opts"             : {"ymin": _yMin, "ymaxfactor": _yMaxF},
        "opts2"            : {"ymin": 0.65, "ymax": 1.35},
        "divideByBinWidth" : False,
        "log"              : _logY,
        "moveLegend"       : _legNE,
        "moveBlindedText"  : {"dx": -0.23, "dy": +0.08, "dh": 0.0},
        }

    kwargs = copy.deepcopy(_kwargs)
    
    if "pt" in h.lower():
        kwargs["units"]      = "GeV"
        kwargs["xlabel"]     = "p_{T} (%s)" % kwargs["units"]
        kwargs["ylabel"]     = _yLabel + kwargs["units"]
        kwargs["moveLegend"] = _legNE
        kwargs["divideByBinWidth"] = True

    if "eta_" in h.lower():
        kwargs["rebinX"] = 1
        kwargs["xlabel"] = "#eta"
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["opts"]   = {"xmin": -2.5, "xmax": +2.5, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["moveLegend"] = _legRM #_legNE
        kwargs["divideByBinWidth"] = False
        if "dilepton" in h.lower():
            kwargs["opts"]  = {"xmin": -5.0, "xmax": +5.0, "ymin": _yMin, "ymaxfactor": _yMaxF}            

    if "Lt" in h:
        kwargs["rebinX"] = systematics.DataMCBinningHToHW["LT"]
        kwargs["units"]  = "GeV"
        kwargs["divideByBinWidth"] = True

    if "HT" in h or "JT" in h:
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
    if "jetn" in h.lower() or "njet" in h.lower():
        kwargs["xlabel"] = "jet multiplicity"
        kwargs["opts"]   = {"xmax": 14.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
    if "bjetn" in h.lower() or "nbjet" in h.lower():
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

    if "met_" in h.lower() or h == "Met" or h == "MET":
        kwargs["units"]  = "GeV"
        kwargs["rebinX"] = systematics.DataMCBinningHToHW["Met"]
        kwargs["xlabel"] = "E_{T}^{miss} (%s)" % (kwargs["units"])
        kwargs["divideByBinWidth"] = True        
        kwargs["ylabel"] = _yLabel + kwargs["units"]

        if "METFilter" in h:
            kwargs["rebinX"] = 1
            kwargs["xlabel"] = ""
            kwargs["ylabel"] = "Events / %.0f "
            kwargs["opts"]   = {"xmin": 0.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "deltar" in h.lower():
        kwargs["units"]  = ""
        kwargs["rebinX"] = 1 #systematics.DataMCBinningHToHW["DeltaR"]
        kwargs["xlabel"] = "#DeltaR"
        kwargs["ylabel"] = "Events / %.2f "
        kwargs["opts"]   = {"xmax": 7.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        
    if "Nprong" in h:
        kwargs["moveLegend"] = _legNE
        kwargs["xlabel"]     = "charged particle multiplicity" # "N_{prongs}"
        kwargs["opts"]       = {"xmax": 8.0, "ymin": _yMin, "ymaxfactor": _yMaxF}

    if "MuonN" in h:
        kwargs["xlabel"] = "#mu multiplicity"
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 8.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        if "passed" in h.lower():
            kwargs["opts"]   = {"xmin": 0.0, "xmax": 8.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 2.5, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "TauN" in h:
        kwargs["opts"]   = {"xmin": 0.0, "xmax": 8.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["cutBox"] = {"cutValue": 1.5, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "mass" in h.lower():
        kwargs["units"] = "GeV" #"GeV/c^{2}"
        kwargs["ylabel"] = _yLabel + kwargs["units"]
        kwargs["opts"]   = {"xmax": 200.0, "ymin": _yMin, "ymaxfactor": _yMaxF}
        kwargs["units"]  = "GeV"
        kwargs["xlabel"] = "m_{\\ell\\ell} (%s)" % kwargs["units"]
        kwargs["xlabel"] = None
        kwargs["ylabel"] = _yLabel + kwargs["units"]

    if "TransverseMass" in h:
        kwargs["rebinX"] = systematics.DataMCBinningHToHW["Mt"]
        kwargs["units"]  = "GeV" #"GeV/c^{2}"
        kwargs["ylabel"] = _yLabel + kwargs["units"]
        kwargs["divideByBinWidth"] = True
        kwargs["opts"]   = {"ymin": _yMin, "ymaxfactor": _yMaxF}

    if kwargs["divideByBinWidth"]:
        if kwargs["units"] != "":
            kwargs["ylabel"] = "< Events / %s  > " % kwargs["units"]
        else:
            kwargs["ylabel"] = "< Events > "

    Verbose("Histagram %s will be rebinned by %s" % (h, str(kwargs["rebinX"])), True)
    #if "tauPt_" in h:
    #    kwargs["rebinX"] = None
    return kwargs
    
def PlotDataMCHistograms(datasetsMgr, hList):

    # Get Histogram name and its kwargs
    histoName = hList[0]
    saveName  = histoName.rsplit("/")[-1].rsplit("_")[0] + "_closure"
    kwargs_   = GetHistoKwargs(histoName, opts)

    # Create the plotting object
    p1 = plots.DataMCPlot(datasetsMgr, hList[0], saveFormats=[])
    p2 = plots.DataMCPlot(datasetsMgr, hList[1], saveFormats=[])

    # Draw and save the plot
    plots.drawPlot(p1, saveName + "_1", **kwargs_)
    plots.drawPlot(p2, saveName + "_2", **kwargs_)

    rhDict = {}
    rhDict["Data-1"] = p1.histoMgr.getHisto("Data").getRootHisto().Clone("Data-1")
    rhDict["Data-2"] = p2.histoMgr.getHisto("Data").getRootHisto().Clone("Data-2")

    # Create empty histogram stack list
    myStackList = []

    h1 = histograms.Histo(rhDict["Data-1"], "Data", "Tight")
    h1.setIsDataMC(isData=True, isMC=False)
    #hFakeB.setIsDataMC(isData=False, isMC=True)
    myStackList.append(h1)

    rhDict["Data-2"].Scale(0.35) #fixme:  divide inclusive histos to get factor
    styles.getFakeTauStyle().apply(rhDict["Data-2"])
    h2 = histograms.Histo(rhDict["Data-2"], "j#rightarrow #tau_{h}", "Loose")
    h2.setIsDataMC(isData=False, isMC=True)
    myStackList.append(h2)

    p = plots.DataMCPlot2( myStackList, saveFormats=[])
    p.setLuminosity(opts.intLumi)
    p.setDefaultStyles()
    plots.drawPlot(p, saveName, **kwargs_)
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
    ANALYSISNAME = "TauFakeRate_mm"
    SEARCHMODE   = "350to3000"
    DATAERA      = "Run2016"
    GRIDX        = False
    GRIDY        = False
    OPTMODE      = None
    BATCHMODE    = True
    INTLUMI      = -1.0
    SIGNALMASS   = 500
    URL          = False
    SAVEDIR      = None
    VERBOSE      = False
    #HISTOS       = "hNTau,hTauSrc,hTauSrcDM0,hTauSrcDM1,hTauSrcDM10,hNJet,hNBjet,hMet,hLt,hDileptonPt,hDileptonEta,hDileptonMass"
    HISTOS       = "NTau,NJet,NBjet,MethLt,DileptonPt,DileptonEta,DileptonMass"
    RATIO        = False
    FOLDER       = ""
    SAVEFORMATS  = "pdf,png,C"
    
    # Define the available script options
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("-m", "--mcrab", dest="mcrab", action="store", 
                      help="Path to the multicrab directory for input")

    parser.add_option("-o", "--optMode", dest="optMode", type="string", default=OPTMODE, 
                      help="The optimization mode when analysis variation is enabled  [default: %s]" % OPTMODE)

    parser.add_option("-b", "--batchMode", dest="batchMode", action="store_false", default=BATCHMODE, 
                      help="Enables batch mode (canvas creation does NOT generate a window) [default: %s]" % BATCHMODE)

    parser.add_option("-a", "--analysisName", dest="analysisName", type="string", default=ANALYSISNAME,
                      help="Override default analysisName [default: %s]" % ANALYSISNAME)

    parser.add_option("--ratio", dest="ratio", action="store_true", default=RATIO,
                      help="Enable ratio pad for Data/Bkg comparison [default: %s]" % RATIO)

    parser.add_option("--intLumi", dest="intLumi", type=float, default=INTLUMI,
                      help="Override the integrated lumi [default: %s]" % INTLUMI)

    parser.add_option("--searchMode", dest="searchMode", type="string", default=SEARCHMODE,
                      help="Override default searchMode [default: %s]" % SEARCHMODE)

    parser.add_option("--dataEra", dest="dataEra", type="string", default=DATAERA, 
                      help="Override default dataEra [default: %s]" % DATAERA)

    parser.add_option("--gridX", dest="gridX", action="store_true", default=GRIDX, 
                      help="Enable the x-axis grid lines [default: %s]" % GRIDX)

    parser.add_option("--gridY", dest="gridY", action="store_true", default=GRIDY, 
                      help="Enable the y-axis grid lines [default: %s]" % GRIDY)

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

    parser.add_option("-s", "--saveFormats", dest="formats", default = SAVEFORMATS,
                      help="Save formats for all plots [default: %s]" % SAVEFORMATS)

    parser.add_option("--histos", dest="histos", default = HISTOS,
                      help="Name of closure histograms (_LooseTau, _TightTau) to be plotted [default: %s]" % HISTOS)

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

    # Create save formats
    if "," in opts.formats:
        opts.saveFormats = opts.formats.split(",")
    else:
        opts.saveFormats = [opts.formats]
    opts.saveFormats = ["." + s for s in opts.saveFormats]

    # Defien the histograms to be plotted
    if "," in opts.histos:
        opts.histoList = opts.histos.split(",")
    else:
        opts.histoList = [opts.histos]
    opts.looseList = ["%s_LooseTau" % s for s in opts.histoList]
    opts.tightList = ["%s_TightTau" % s for s in opts.histoList]

    # Call the main function
    opts.histoDict = {}
    opts.histoDict["loose"] = opts.looseList
    opts.histoDict["tight"] = opts.tightList
    main(opts)

    if not opts.batchMode:
        raw_input("=== plotSelfClosure.py: Press any key to quit ROOT ...")
