#!/usr/bin/env python
'''
Description:
 
 This script creates two JSON files with the probability density functions for quark- and gluon-jets in pT bins.
 Running on QCD and/or WJetsToQQ_HT. If running on QCD remove low HT bins (at least up to HT = 200 GeV to avoid QCD mc spikes)

Usage:
./plotQGL.py -m <pseudo_mcrab> [opts]

Examples:
./plotQGL.py -m QGLAnalysis_171129_030702/ --url

Last Used:
./plotQGL.py -m QGLAnalysis_191213_055119 --url
./plotQGL.py -m QGLAnalysis_200214_150436 --url -i "QCD_HT|WJetsToQQ_HT|TTTo"


To see the results in a more user friendly environment on https://home.fnal.gov/~username/<pseudomulticrab> do:
$setenv PATH ${PATH}:$HOME/scripts/gallery
$gallery.py /publicweb/m/mkolosov/<pseudomulticrab>
'''

#================================================================================================ 
# Imports
#================================================================================================ 
import sys
import math
import copy
import json
import os
from optparse import OptionParser

import ROOT
import array

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
import HiggsAnalysis.NtupleAnalysis.tools.multicrabConsistencyCheck as consistencyCheck
import HiggsAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles

#================================================================================================ 
# Function Definition
#================================================================================================ 
def Print(msg, printHeader=False):
    fName = __file__.split("/")[-1]
    if printHeader:
        print("===", fName)
        print(msg)
    else:
        print(msg)
    return

def Verbose(msg, printHeader=True, verbose=False):
    if not opts.verbose:
        return
    aux.Print(msg, printHeader)
    return

def GetLumi(datasetsMgr):
    Verbose("Determininig Integrated Luminosity")
    
    lumi = 0.0
    for d in datasetsMgr.getAllDatasets():
        if d.isMC():
            continue
        else:
            lumi += d.getLuminosity()
    Verbose("Luminosity = %s (pb)" % (lumi), True)
    return lumi

def GetListOfQCDatasets():
    Verbose("Getting list of QCD datasets")
    QCD = []
#    QCD.append("QCD_HT50to100")
#    QCD.append("QCD_HT100to200")
#    QCD.append("QCD_HT200to300")
    QCD.append("QCD_HT300to500")
    QCD.append("QCD_HT500to700")
    QCD.append("QCD_HT700to1000")
    QCD.append("QCD_HT1000to1500")
    QCD.append("QCD_HT1500to2000")
    QCD.append("QCD_HT2000toInf")
    return QCD

def GetListOfEwkDatasets():
    Verbose("Getting list of EWK datasets")
    if "2017" in opts.dataEra:
        return ["TTToHadronic", "TTToSemiLeptonic", "TTTo2L2Nu", "WJetsToQQ_HT_800toInf_qc19_3j", "WJetsToQQ_HT600to800_qc19_3j", "WJetsToQQ_HT400to600_qc19_3j"]#, "TTWJetsToQQ"]
    elif "2016" in opts.dataEta:
        return ["WJetsToQQ_HT_600ToInf"] #"TT", "WJetsToQQ_HT_600ToInf", "DYJetsToQQHT", "SingleTop", "TTWJetsToQQ", "TTZToQQ", "Diboson", "TTTT"]
    else:
        return []
    
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

# =========================================================================================================================================
#                                                          MAIN
# =========================================================================================================================================
def main(opts):
    
    # Apply TDR style
    style = tdrstyle.TDRStyle()
    style.setGridX(False)
    style.setGridY(False)
    style.setOptStat(False)
    
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
        datasetsMgr.loadLuminosities()
        
        if 0:
            datasetsMgr.printSelections()
            
        # Set/Overwrite cross-sections
        for d in datasetsMgr.getAllDatasets():
            if "ChargedHiggs" in d.getName():
                datasetsMgr.getDataset(d.getName()).setCrossSection(1.0)
                
        if opts.verbose:
            datasetsMgr.PrintCrossSections()
            datasetsMgr.PrintLuminosities()
            
        # Filter the datasets
        datasetsMgr.remove(filter(lambda name: "Charged" in name, datasetsMgr.getAllDatasetNames()))
        if "ZJetsToQQ_HT" in datasetsMgr.getAllDatasetNames() and "DYJetsToQQ" in datasetsMgr.getAllDatasetNames():         
            # ZJets and DYJets overlap!
            Print("Cannot use both ZJetsToQQ and DYJetsToQQ due to duplicate events? Investigate. Removing ZJetsToQQ datasets for now ..", True)
            datasetsMgr.remove(filter(lambda name: "ZJetsToQQ" in name, datasetsMgr.getAllDatasetNames()))

        # Remove low QCD_HT bins
        datasetsMgr.remove(filter(lambda name: "QCD_HT50to100" in name, datasetsMgr.getAllDatasetNames()))
        datasetsMgr.remove(filter(lambda name: "QCD_HT100to200" in name, datasetsMgr.getAllDatasetNames()))
        
        datasetsMgr.PrintInfo()
            
        # Merge histograms (see NtupleAnalysis/python/tools/plots.py)
        plots.mergeRenameReorderForDataMC(datasetsMgr)
        
        # Get luminosity
        if opts.intLumi < 0:
            if "Data" in datasetsMgr.getAllDatasetNames():
                opts.intLumi = datasetsMgr.getDataset("Data").getLuminosity()
            else:
                opts.intLumi = 1.0
                
        # Remove Data
        datasetsMgr.remove(filter(lambda name: "Data" in name, datasetsMgr.getAllDatasetNames()))

        # Merge EWK samples
        #datasetsMgr.merge("EWK", GetListOfEwkDatasets())
                
        # Print post EWK-merge dataset summary
        if 1:
            datasetsMgr.PrintInfo()
            
        # Set plot Styles
        plots._plotStyles["QCD"] = styles.getQCDFillStyle()
        #plots._plotStyles["TT"] = styles.getAltEWKStyle()
        #plots._plotStyles["WJetsHT"] = styles.getEWKStyle()
        
        PtRange = ["30pt35", "35pt40", "40pt45", "45pt50", "50pt55", "55pt60", "60pt65", "65pt70", "70pt75", "75pt80", "80pt90", "90pt100", "100pt110", "110pt120",
                   "120pt140", "140pt160", "160pt180", "180pt200", "200pt250", "250pt300", "300pt350", "350pt400", "400pt450", "450pt500", "500pt550", "550pt600",
                   "600pt700", "700pt800", "800pt1000", "1000ptInf"]
        # PtRange = ["30pt40", "40pt50", "50pt65", "65pt80", "80pt100", "100pt125", "125pt160", "160pt200", "200pt250", "250pt320", "320pt400", "400pt630", "630pt800", "800ptInf"] 
        
        histos = ["UDSJetsQGL_"+ pt for pt in PtRange]
        histos.extend("GJetsQGL_"+ pt for pt in PtRange)
                
        for i, h in enumerate(histos, 1):
            msg   = "{:<9} {:>3} {:<1} {:<3} {:<50}".format("Histogram", "%i" % i, "/", "%s:" % (len(histos)), h)
            Print(ShellStyles.SuccessStyle() + msg + ShellStyles.NormalStyle(), i==1)
            PlotMC(datasetsMgr, h, opts.intLumi)
            
        DumpPDFinJSON(datasetsMgr, PtRange)
    
    Print("All plots saved under directory %s" % (ShellStyles.NoteStyle() + aux.convertToURL(opts.saveDir, opts.url) + ShellStyles.NormalStyle()), True)
    return


def DumpPDFinJSON(datasetsMgr, PtRange):
    '''
    Function to dump PDFs into JSON files
    '''
    # Loop over datasets
    for dset in datasetsMgr.getAllDatasets():
        
        # Loop over Jet type
        for JetType in ["UDS", "G"]:
            
            results = []
            jsonhistos = [JetType+"JetsQGL_"+ pt for pt in PtRange]
            
            for h in jsonhistos:
                dsetHisto = dset.getDatasetRootHisto(h)
                dsetHisto.normalizeToOne()
                
                histo = dsetHisto.getHistogram()
                
                ptbin   = h.split("_")[-1]
                minPt   = ptbin.split("pt")[0]
                maxPt   = ptbin.split("pt")[-1]
                
                if maxPt == "Inf":
                    maxPt = 9999999999.9
                    
                for k in range(1, histo.GetNbinsX()+1):
                    resultObject = {}
                    resultObject["Jet"]       = JetType
                    resultObject["QGLmin"]    = histo.GetBinLowEdge(k) 
                    resultObject["QGLmax"]    = histo.GetBinLowEdge(k)+histo.GetBinWidth(k)
                    resultObject["Ptmin"]     = minPt
                    resultObject["Ptmax"]     = maxPt
                    resultObject["prob"]      = histo.GetBinContent(k)
                    resultObject["probError"] = histo.GetBinError(k)
                    results.append(resultObject)
                    
            filename = "QGLdiscriminator_%s_%sJets_%s.json"%(dset.name, JetType, opts.dataEra)
            with open(filename, 'w') as outfile:
                json.dump(results, outfile)
                
            print("Written results to %s" % (filename))
    return

def PlotMC(datasetsMgr, histo, intLumi):
    '''
    Plot MC samples only
    '''
    kwargs = {}
    p = plots.MCPlot(datasetsMgr, histo, normalizeToOne=True, saveFormats=[], **kwargs)
    
    # Draw the histograms
    _cutBox = None
    _rebinX = 1
    _format = "%0.01f"
    _xlabel = None
    logY    = False
    _opts   = {"ymin": 1e-3, "ymaxfactor": 1.2}
    
    if "Pt" in histo:
        if "GJets" in histo:
            _xlabel = "Gluon Jets p_{T} (GeV)"
        elif "UDS" in histo:
            _xlabel = "UDS Jets p_{T} (GeV)"
            
    if "JetsN" in histo:
        if "GJets" in histo:
            _xlabel = "Gluon Jets Multiplicity"
        elif "UDS" in histo:
            _xlabel = "UDS Jets Multiplicity"
            
    if "QGL" in histo:
        _opts["xmin"] = 0.0
        _opts["xmax"] = 1.0
        _units  = ""
        _format = "%0.01f " + _units
        if "GJets" in histo:
            _xlabel = "Gluon Jets QGL"
        elif "UDS" in histo:
            _xlabel = "UDS Jets QGL"
            
    # Customise styling
    p.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetLineStyle(ROOT.kSolid))
    plots.drawPlot(p, 
                   histo,  
                   xlabel       = _xlabel,
                   ylabel       = "Arbitrary Units",# / %s" % (_format),
                   log          = logY,
                   rebinX       = _rebinX, cmsExtraText = "Simulation Preliminary 2016",
                   createLegend = {"x1": 0.20, "y1": 0.77, "x2": 0.40, "y2": 0.90},
                   opts         = _opts,
                   opts2        = {"ymin": 0.6, "ymax": 1.4},
                   cutBox       = _cutBox,
                   )
    
    if "pt" in histo:
        ptRange = histo.split("_")[1]
        ptRange = ptRange.split("pt")[0] + " < p_{T} < " + ptRange.split("pt")[-1] + " GeV"
        histograms.addText(0.55, 0.84, ptRange, size=20)
        
    # Save plot in all formats    
    saveName = histo.split("/")[-1]
    SavePlot(p, saveName, os.path.join(opts.saveDir, opts.optMode, opts.folder), [".png", ".pdf"] )
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
    ANALYSISNAME = "QGLAnalysis"
    SEARCHMODE   = "80to1000"
    DATAERA      = "Run2016"
    OPTMODE      = ""
    BATCHMODE    = True
    PRECISION    = 3
    INTLUMI      = -1.0
    SUBCOUNTERS  = False
    LATEX        = False
    URL          = False
    SAVEDIR      = None
    VERBOSE      = False
    HISTOLEVEL   = "Vital" # 'Vital' , 'Informative' , 'Debug'
    NORMALISE    = False
    FOLDER       = ""
    
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
    
    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR, 
                      help="Directory where all pltos will be saved [default: %s]" % SAVEDIR)
    
    parser.add_option("--url", dest="url", action="store_true", default=URL, 
                      help="Don't print the actual save path the histogram is saved, but print the URL instead [default: %s]" % URL)
    
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=VERBOSE, 
                      help="Enables verbose mode (for debugging purposes) [default: %s]" % VERBOSE)
    
    parser.add_option("--histoLevel", dest="histoLevel", action="store", default = HISTOLEVEL,
                      help="Histogram ambient level (default: %s)" % (HISTOLEVEL))
    
    parser.add_option("-i", "--includeOnlyTasks", dest="includeOnlyTasks", action="store", 
                      help="List of datasets in mcrab to include")
    
    parser.add_option("-e", "--excludeTasks", dest="excludeTasks", action="store", 
                      help="List of datasets in mcrab to exclude")
    
    parser.add_option("-n", "--normaliseToOne", dest="normaliseToOne", action="store_true", 
                      help="Normalise the histograms to one? [default: %s]" % (NORMALISE) )
    
    parser.add_option("--folder", dest="folder", type="string", default = FOLDER,
                      help="ROOT file folder under which all histograms to be plotted are located [default: %s]" % (FOLDER) )
    
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
        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="", postfix="")

    # Call the main function
    main(opts)

    if not opts.batchMode:
        raw_input("=== plotQGL.py: Press any key to quit ROOT ...")
