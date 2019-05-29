#!/usr/bin/env python
'''
Description:
Script for plotting TH2 plots


USEFUL LINKS:
https://nixtricks.wordpress.com/2011/03/03/simple-loops-in-csh-ortcsh/


USAGE:
./plot_2d.py -m <pseudo_mcrab_directory> [opts]


EXAMPLES:
./plot_2d.py -m HtbKinematics_170831_085353 --url --mergeEWK -e "JetHT"


LAST USED:
./plot_2d.py --normalizeToOne --gridX --gridY --logZ -m KinematicsHToHW_190424_024614


OBSOLETE:
foreach x ( 180 200 220 250 300 350 400 500 800 1000 2000 3000 )
./plot_2d.py -m Hplus2tbAnalysis_StdSelections_TopCut100_AllSelections_NoTrgMatch_TopCut10_H2Cut0p5_170826_073257/ -i M_$x --url
end

'''

#================================================================================================ 
# Imports
#================================================================================================ 
import sys
import math
import copy
import os
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
import HiggsAnalysis.NtupleAnalysis.tools.crosssection as xsect
import HiggsAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles

#================================================================================================ 
# Shell Definitions
#================================================================================================ 
ss = ShellStyles.SuccessStyle()
ns = ShellStyles.NormalStyle()
ts = ShellStyles.NoteStyle()
ls = ShellStyles.HighlightStyle()
hs = ShellStyles.HighlightAltStyle()
es = ShellStyles.ErrorStyle()
styleList = [styles.Style(24, ROOT.kBlack)] + styles.getStyles()

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
    aux.Print(msg, printHeader)
    return

def natural_sort_key(s):
    _nsre = re.compile('([0-9]+)')
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]   

def GetLumi(datasetsMgr):
    '''
    '''
    lumi = 0.0
    for d in datasetsMgr.getAllDatasets():
        if d.isMC():
            continue
        else:
            #print "dataset = %s, lumi = %s" % (d.getName(), lumi)
            lumi += d.getLuminosity()
    #print "Luminosity = %s (pb)" % (lumi), True)
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

def GetAllDatasetsFromDir(opts):
    datasets = dataset.getDatasetsFromMulticrabDirs([opts.mcrab],
                                                    dataEra=opts.dataEra,
                                                    searchMode=opts.searchMode, 
                                                    analysisName=opts.analysisName,
                                                    optimizationMode=opts.optMode)
    return datasets
    
def getHisto(datasetsMgr, datasetName, name):

    h1 = datasetsMgr.getDataset(datasetName).getDatasetRootHisto(name)
    h1.setName("h0" + "-" + datasetName)
    return h1

def main(opts):

    # Apply TDR style
    style = tdrstyle.TDRStyle()
    style.setOptStat(False)
    style.setGridX(opts.gridX)
    style.setGridY(opts.gridY)
    style.setLogX(opts.logX)
    style.setLogY(opts.logY)
    style.setLogZ(opts.logZ)
    #style.setWide(True, 0.15)

    optModes = [""]
    if opts.optMode != None:
        optModes = [opts.optMode]
        
    # For-loop: All optimisation modes
    for opt in optModes:
        opts.optMode = opt

        # Setup & configure the dataset manager 
        datasetsMgr    = GetDatasetsFromDir(opts)
        datasetsMgr.updateNAllEventsToPUWeighted()
        # datasetsMgr.loadLuminosities() # from lumi.json
        
        # Set/Overwrite cross-sections
        for d in datasetsMgr.getAllDatasets():
            if "ChargedHiggs" in d.getName():
                datasetsMgr.getDataset(d.getName()).setCrossSection(1.0) # ATLAS 13 TeV H->tb exclusion limits
            
        if opts.verbose:
            datasetsMgr.PrintCrossSections()
            datasetsMgr.PrintLuminosities()

        # Merge histograms (see NtupleAnalysis/python/tools/plots.py) 
        plots.mergeRenameReorderForDataMC(datasetsMgr) 

        # Print dataset information
        datasetsMgr.PrintInfo()

        # Merge EWK samples
        if opts.mergeEWK:
            datasetsMgr.merge("EWK", aux.GetListOfEwkDatasets())
            plots._plotStyles["EWK"] = styles.getAltEWKStyle()

        # Do 2D histos
        histoNames  = []
        histoList   = datasetsMgr.getDataset(datasetsMgr.getAllDatasetNames()[0]).getDirectoryContent(opts.folder)
        histoPaths  = [opts.folder +"/" + h for h in histoList]

        # Axes divisions!
        ROOT.gStyle.SetNdivisions(5, "X")
        ROOT.gStyle.SetNdivisions(5, "Y")

        # For-loop: All histogram
        for index, histoName in enumerate(histoPaths, 1):

            # For-loop: All datasets
            for d in datasetsMgr.getAllDatasetNames():

                opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="", postfix=os.path.join("TH2", d))
                msg = "%s%d/%d %s" % (ss, index, len(histoPaths),  histoName.replace("/", "/"+ d + "/") + ns)
                aux.PrintFlushed(msg, index==1)
                Plot2dHistograms(datasetsMgr, d, histoName, index)

                # Avoid replacing canvas with same name warning
                #for o in gROOT.GetListOfCanvases():
                #    # print o.GetName()
                #    o.SetName(o.GetName() + "_" + d)

    print
    aux.Print("All plots saved under %s" % (hs + opts.saveDir + ns), True)
    return

def GetHistoKwargs(h, opts):
    _kwargs = {}
    _leg    = {"dx": -0.1, "dy": 0.0, "dh": -0.15}
    logY    = False
    yMin    = 0.0
    if logY:
        yMin = 0.001
        yMaxF = 10
    else:
        yMaxF = 1.0
    
    # z-axis settings
    zMin   =  0
    zMax   = None
    zLabel = "z-axis"
    if opts.normalizeToLumi:
        zLabel  = "Events"
        zMin    = 1e0
    elif opts.normalizeByCrossSection:
        zLabel  = "#sigma (pb)"
    elif opts.normalizeToOne:
        zLabel  = "Arbitrary Units"
        zMin    = 1e-5 # does not work
        zMax    = 1e+0 # does not work
    else:
        zLabel = "Unknown"

    # Cut lines/boxes
    _cutBoxX = {"cutValue": 1.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
    _cutBoxY = {"cutValue": 1.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True,
                "mainCanvas": True, "ratioCanvas": False} # box = True not working       

    _kwargs = {
        "xlabel"           : None,
        "ylabel"           : None,
        "stackMCHistograms": False,
        "addMCUncertainty" : False,
        "addLuminosityText": opts.normalizeToLumi,
        "addCmsText"       : True,
        "cmsExtraText"     : "Preliminary",
        "cmsTextPosition"  : "outframe",
        "log"              : logY,
        "moveLegend"       : _leg,
        "xtitlesize"       : 0.1,
        "ytitlesize"       : 0.1,
        "zmin"             : zMin,
        "zmax"             : zMax,
        "zlabel"           : zLabel,
        "cutBox"           : _cutBoxX,
        "cutBoxY"          : _cutBoxY,
        "rebinX"           : 1, 
        "rebinY"           : 1,
        }

    if 0:
        ROOT.gStyle.SetNdivisions(10, "X")
        ROOT.gStyle.SetNdivisions(10, "Y")
        ROOT.gStyle.SetNdivisions(10, "Z")

    if "Pt_Vs_Pt" in h:
        _kwargs["rebinX"]  =  1
        _kwargs["rebinY"]  =  1
        _kwargs["opts"]    = {"xmin": 0.0, "xmax": 300.0, "ymin": 0.0, "ymax": 300.0}

    if "Pt1_Vs_dR" in h:
        _kwargs["rebinX"]  =  1
        _kwargs["rebinY"]  =  1
        _kwargs["opts"]    = {"xmin": 0.0, "xmax": 600.0, "ymin": 0.0, "ymax": 6.0}
        _kwargs["cutBoxY"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True,
                              "mainCanvas": True, "ratioCanvas": False} # box = True not working       

    if "Pt2_Vs_dR" in h:
        _kwargs["rebinX"]  =  1
        _kwargs["rebinY"]  =  1
        _kwargs["opts"]    = {"xmin": 0.0, "xmax": 300.0, "ymin": 0.0, "ymax": 6.0}
        _kwargs["cutBoxY"] = {"cutValue": 0.4, "fillColor": 16, "box": False, "line": True, "greaterThan": True,
                              "mainCanvas": True, "ratioCanvas": False} # box = True not working       

    if "dEta_Vs_dPhi" in h:
        _kwargs["rebinX"]  =  1
        _kwargs["rebinY"]  =  1
        _kwargs["opts"]    = {"xmin": 0.0, "xmax": 5.0, "ymin": 0.0, "ymax": 3.2}
        # _kwargs["zmin"]    = 1 #fixme

    if "Eta_Vs_Eta" in h:
        _kwargs["rebinX"]  =  1
        _kwargs["rebinY"]  =  1
        _kwargs["opts"]    = {"xmin": 0.0, "xmax": 4.0, "ymin": 0.0, "ymax": 4.0}
        _kwargs["xlabel"]  = None
        _kwargs["ylabel"]  = None

    if "Phi_Vs_Phi" in h:
        _kwargs["rebinX"]  =  1
        _kwargs["rebinY"]  =  1
        _kwargs["opts"]    = {"xmin": 0.0, "xmax": 3.2, "ymin": 0.0, "ymax": 3.2}
        _kwargs["xlabel"]  = None
        _kwargs["ylabel"]  = None

    if "dRap_Vs_dPhi" in h:
        _kwargs["rebinX"]  =  1
        _kwargs["rebinY"]  =  1
        _kwargs["opts"]    = {"xmin": 0.0, "xmax": 5.0, "ymin": 0.0, "ymax": 3.2}

    if "dEta_Vs_Jet1Jet2_Mass" in h:
        _kwargs["rebinX"] =    1
        _kwargs["rebinY"] =    2
        _kwargs["opts"]   = {"xmin": 0.0, "xmax": +5.0, "ymin": 0.0, "ymax": 1000.0}

    if "dEta_Vs_Jet3Jet4_dEta" in h:
        _kwargs["rebinX"] =    1
        _kwargs["rebinY"] =    2
        _kwargs["opts"]   = {"xmin": 0.0, "xmax": +5.0, "ymin": 0.0, "ymax": 5.0}

    if "dPhi_Vs_Jet3Jet4_dPhi" in h:
        _kwargs["rebinX"] =    1
        _kwargs["rebinY"] =    1
        _kwargs["opts"]   = {"xmin": 0.0, "xmax": +3.2, "ymin": 0.0, "ymax": 3.2}

    if "dEta_Vs_Jet3Jet4_Mass" in h:
        _kwargs["rebinX"]  =    1
        _kwargs["rebinY"]  =    1
        _kwargs["opts"]    = {"xmin": 0.0, "xmax": +5.0, "ymin": 0.0, "ymax": 1000.0}
        _kwargs["cutBox"]  = {"cutValue":25.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
        _kwargs["cutBoxY"] = {"cutValue": 20.00, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
        
    return _kwargs


def GetHisto(datasetsMgr, dataset, histoName):
    h = datasetsMgr.getDataset(dataset).getDatasetRootHisto(histoName)
    return h


def Plot2dHistograms(datasetsMgr, dsetName, histoName, index):
                                    
    # Custom Filtering of datasets 
    dsetsMgr = datasetsMgr.deepCopy()
    if opts.verbose:               
        dsetsMgr.PrintInfo()       
                  
    # Get Histogram name and its kwargs
    #saveName = histoName.rsplit("/")[-1] + "_" + dsetName.split("_")[0] + dsetName.split("_")[-1]
    #saveName = dsetName + "/" + histoName.rsplit("/")[-1]
    saveName = histoName.rsplit("/")[-1]
    kwargs_  = GetHistoKwargs(saveName, opts)                                                    
    if 0:         
        print "\n%s: %s" % (histoName, kwargs_)                                                  
                  
    for i, d in enumerate(dsetsMgr.getAllDatasetNames(), 0):                                     
        if d == dsetName:          
            continue               
        else:     
            # Remove dataset from manager but do NOT close the file!                             
            dsetsMgr.remove(d, close=False)                                                      
                  
    # Sanity check
    nDatasets = len(dsetsMgr.getAllDatasets())                                                   
    if nDatasets > 1:              
        raise Exception("More than 1 datasets detected in the dataset manager! Can only support 1 dataset. Please use the -i option to choose exactly 1 dataset'")
                  
    # Get the reference histo and the list of histos to compare                                  
    datasets0 = dsetsMgr.getAllDatasets()[0].getName()                                           
    histoList = [getHisto(dsetsMgr, datasets0, histoName)]                                       
                  
    # Create the 2d plot           
    Verbose("Creating the 2d plot", True)                                                        
    if opts.normalizeToLumi:       
        p = plots.MCPlot(dsetsMgr, histoName, normalizeToLumi=opts.intLumi, saveFormats=[])      
    elif opts.normalizeByCrossSection:                                                           
        p = plots.MCPlot(dsetsMgr, histoName, normalizeByCrossSection=True, saveFormats=[], **{})
    elif opts.normalizeToOne:      
        p = plots.MCPlot(dsetsMgr, histoName, normalizeToOne=True, saveFormats=[], **{})         
    else:         
        raise Exception("One of the options --normalizeToOne, --normalizeByCrossSection, --normalizeToLumi must be enabled (set to \"True\").")
                  
    Verbose("Setting universal histo styles", True)                                              
    p.histoMgr.setHistoDrawStyleAll("COLZ")                                                      
                  
    Verbose("Customising histograms", True)                                                      
    p.histoMgr.forEachHisto(lambda h: h.getRootHisto().GetZaxis().SetTitleOffset(1.3)) #fixme    
                  
    Verbose("Setting plot styles to histograms", True)                                           
    for index, h in enumerate(p.histoMgr.getHistos()):                                           
        plots._plotStyles[p.histoMgr.getHistos()[index].getDataset().getName()].apply(p.histoMgr.getHistos()[index].getRootHisto())
        p.histoMgr.getHistos()[index].getRootHisto().SetMinimum(kwargs_["zmin"])                 
        p.histoMgr.getHistos()[index].getRootHisto().SetMaximum(kwargs_["zmax"])                 
                  
    Verbose("Drawing the plot", True)                          
    plots.drawPlot(p, saveName, **kwargs_) #the "**" unpacks the kwargs_ dictionary
                  
    # Add line?
    Verbose("Removing the legend", True)
    p.removeLegend()

                  
    Verbose("Adding text on canvas", True)
    histograms.addText(0.22, 0.89, plots._legendLabels[datasets0], 18)
                  
    Verbose("Saving the canvas", True)
    SavePlot(p, saveName, opts.saveDir, opts.saveFormats)
                  
    return                      

def HasKeys(keyList, **kwargs):
    for key in keyList:
        if key not in kwargs:
            raise Exception("Could not find the keyword \"%s\" in kwargs" % (key) )
    return 

def SavePlot(plot, saveName, saveDir, saveFormats = [".C", ".png", ".pdf"]):
    Verbose("Saving the plot in %s formats: %s" % (len(saveFormats), ", ".join(saveFormats) ) )
    
    # Check that path exists
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
        
    savePath = os.path.join(saveDir, saveName)

    # For-loop: All save formats
    for i, ext in enumerate(saveFormats, 1):
        saveNameURL = savePath + ext
        saveNameURL = saveNameURL.replace("/publicweb/a/aattikis/", "http://home.fnal.gov/~aattikis/")
        if opts.url:
            Verbose(saveNameURL, i==0)
        else:
            Verbose(savePath + ext, i==0)
        plot.saveAs(savePath, formats=saveFormats)
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
    ANALYSISNAME = "KinematicsHToHW"
    SEARCHMODE   = "80to1000"
    DATAERA      = "Run2016"
    OPTMODE      = None
    BATCHMODE    = True
    INTLUMI      = -1.0
    MERGEEWK     = False
    URL          = False
    SAVEDIR      = None
    VERBOSE      = False
    FOLDER       = "TH2"
    NORM2ONE     = False
    NORM2XSEC    = False
    NORM2LUMI    = False
    GRIDX        = False
    GRIDY        = False
    LOGX         = False
    LOGY         = False
    LOGZ         = False
    SAVEFORMATS  = "png"

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

    parser.add_option("--searchMode", dest="searchMode", type="string", default=SEARCHMODE,
                      help="Override default searchMode [default: %s]" % SEARCHMODE)

    parser.add_option("--dataEra", dest="dataEra", type="string", default=DATAERA, 
                      help="Override default dataEra [default: %s]" % DATAERA)

    parser.add_option("--mergeEWK", dest="mergeEWK", action="store_true", default=MERGEEWK, 
                      help="Merge all EWK samples into a single sample called \"EWK\" [default: %s]" % MERGEEWK)

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR, 
                      help="Directory where all pltos will be saved [default: %s]" % SAVEDIR)

    parser.add_option("--url", dest="url", action="store_true", default=URL, 
                      help="Don't print the actual save path the histogram is saved, but print the URL instead [default: %s]" % URL)
    
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=VERBOSE, 
                      help="Enables verbose mode (for debugging purposes) [default: %s]" % VERBOSE)

    parser.add_option("--folder", dest="folder", action="store", default = FOLDER,
                      help="ROOT file folder under which all histograms to be plotted are located [default: %s]" % (FOLDER) )

    parser.add_option("-i", "--includeOnlyTasks", dest="includeOnlyTasks", action="store", 
                      help="List of datasets in mcrab to include")

    parser.add_option("-e", "--excludeTasks", dest="excludeTasks", action="store", 
                      help="List of datasets in mcrab to exclude")

    parser.add_option("--normalizeToOne", dest="normalizeToOne", action="store_true", default=NORM2ONE,
                      help="Normalise plot to one [default: %s]" % NORM2ONE)
    
    parser.add_option("--normalizeByCrossSection", dest="normalizeByCrossSection", action="store_true", default=NORM2XSEC,
                      help="Normalise plot by cross-section [default: %s]" % NORM2XSEC)
    
    parser.add_option("--normalizeToLumi", dest="normalizeToLumi", action="store_true", default=NORM2LUMI,
                      help="Normalise plot to luminosity [default: %s]" % NORM2LUMI)

    parser.add_option("-s", "--saveFormats", dest="saveFormats", default = SAVEFORMATS,
                      help="Save formats for all plots [default: %s]" % SAVEFORMATS)     

    parser.add_option("--gridX", dest="gridX", action="store_true", default=GRIDX,
                      help="Enable the x-axis grid lines [default: %s]" % GRIDX)
    
    parser.add_option("--gridY", dest="gridY", action="store_true", default=GRIDY,
                      help="Enable the y-axis grid lines [default: %s]" % GRIDY)      

    parser.add_option("--logX", dest="logX", action="store_true", default=LOGX, 
                      help="Set the x-axis to log scale? [default: %s]" % LOGX)

    parser.add_option("--logY", dest="logY", action="store_true", default=LOGY, 
                      help="Set the y-axis to log scale? [default: %s]" % LOGY)

    parser.add_option("--logZ", dest="logZ", action="store_true", default=LOGZ, 
                      help="Set the z-axis to log scale? [default: %s]" % LOGZ)


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
        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="", postfix="TH2")
        Verbose("Save directory set to %s" % (ls + opts.saveDir + ns), True)

    # Call the main function
    main(opts)

    if not opts.batchMode:
        raw_input("=== plot_2d.py: Press any key to quit ROOT ...")
