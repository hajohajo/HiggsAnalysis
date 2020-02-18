#!/usr/bin/env python
'''
DESCRIPTION:
Script for plotting TH2 histograms only.


USAGE:
./plotTBranchToTH2.py -m <pseudo-mcrab> [--options]


EXAMPLES:
./plotTBranchToTH2.py -m TopRecoTree_multFloat_AnalysisSelections --normalizeToOne --url

LAST USED:
./plotTBranchToTH2.py -m TopRecoTree_multFloat_AnalysisSelections --normalizeToOne --ref TrijetPtDR --url
./plotTBranchToTH2.py -m TopRecoTree_200217_114252_FirstNewBranches_hadronic --normalizeToOne --ref trijetPt --url --treeName treeS --logZ
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


import warnings
warnings.filterwarnings("ignore")

ROOT.gErrorIgnoreLevel = ROOT.kError

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
        Plot2DHistograms(datasetsMgr, ref, opts.treeName)
    Print("All plots saved under directory %s" % (ShellStyles.NoteStyle() + aux.convertToURL(opts.saveDir, opts.url) + ShellStyles.NormalStyle()), True)
    return


def GetHistoKwargs(h, opts):

    # Defaults
    xMin  =    0
    xMax   = 800
    yMin   =   0
    yMax   = 800
    yMaxF  =  10
    zMin   =   0
    zMax   = None
    zLabel = "z-axis"
    if opts.logX:
        xMin =  1
    if opts.logY:
        yMin =  1

    if opts.normalizeToLumi:
        zLabel  = "Events"
        zMin    = 1e0
    elif opts.normalizeByCrossSection:
        zLabel  = "#sigma (pb)"
        #zMin    = 0 #1e-3
    elif opts.normalizeToOne:
        zLabel  = "Arbitrary Units"
        zMin    = None #5e-4
        zMax    = None #1e-1
    else:
        zLabel = "Unknown"

    yLabel = "ylabel"
    xLabel = "xlabel"

    cutBox      = {"cutValue": 400.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True} #box = True works
    cutBoxY     = {"cutValue": 200.0, "fillColor": 16, "box": False, "line": True, "greaterThan": True}

    if "_vs_" not in h.lower():
        yLabel = zLabel
    else:
        yLabel = None
                   
    kwargs = {
        "stackMCHistograms": False,
        "addMCUncertainty" : False,
        "addLuminosityText": opts.normalizeToLumi,
        "addCmsText"       : True,
        "cmsExtraText"     : "  Preliminary",
        "xmin"             : xMin,
        "xmax"             : xMax,
        "ymin"             : yMin,
        "zmin"             : zMin,
        "zmax"             : zMax,
        "cutBox"           : cutBox,
        "cutBoxY"          : cutBoxY,
        "moveLegend"       : {"dx": -2.0, "dy": 0.0, "dh": 0}, #hack to remove legend (tmp)
        "xlabel"           : xLabel,
        "ylabel"           : yLabel,
        "zlabel"           : zLabel,
        }
    
    if "mass" in h.lower():
        _units  = "GeV/c^{2}"
        _format = "%0.0f " + _units
        kwargs["xlabel"] = "M (%s)" % _units

    if "pt" in h.lower():
        _units  = "GeV/c"
        _format = "%0.0f " + _units
        kwargs["xlabel"] = "p_{T} (%s)" % _units
        if "trijet" in h.lower():
            kwargs["xlabel"] = "p_{T,t} (%s)" % _units
        if "dijet" in h.lower():
            kwargs["xlabel"] = "p_{T,w} (%s)" % _units
        if "ldgjetpt" in h.lower():
            kwargs["xmax"] = 100
    if "trijetmass" in h.lower():
        _units  = "GeV/c^{2}"
        _format = "%0.0f " + _units
        kwargs["xlabel"] = "m_{top} (%s)" % _units
        kwargs["xmax"] = 805 #1005
    if "bjetmass" in h.lower():
        kwargs["xlabel"] = "m_{b-tagged jet} (%s)" % _units
        kwargs["xmax"] = 50
    if "bjetldgjet_mass" in h.lower():
        kwargs["xlabel"]= "m_{b, ldg jet} (%s)" % _units
        kwargs["xmax"] = 705
    if "bjetsubldgjet_mass" in h.lower():
        kwargs["xlabel"] = "m_{b-tagged, subldg jet} (%s)" % _units
        kwargs["xmax"] = 705
    if "jet_mass" in h.lower():
        kwargs["xmax"] = 750

    if "mult" in h.lower():
        _format = "%0.0f"        
        kwargs["xmax"] = 50
        kwargs["xlabel"] = "mult"
        if "ldg" in h.lower():
            kwargs["xlabel"] = "Leading jet mult"
        if "subldg" in h.lower():
            kwargs["xlabel"] = "Subleading jet mult"
        if "avg" in h.lower():
            kwargs["xlabel"] = "avg CvsL"

    if "cvsl" in h.lower():
        _format = "%0.2f"
        kwargs["xmax"] = 1
        kwargs["xlabel"] = "CvsL"
        if "ldg" in h.lower():
            kwargs["xlabel"] = "Leading jet CvsL"
        if "subldg" in h.lower():
             kwargs["xlabel"] = "Subleading jet CvsL"
        if "avg" in h.lower():
            kwargs["xlabel"] = "avg CvsL"
    if "axis2" in h.lower():
        _format = "%0.3f"
        kwargs["xmax"] = 0.3
        kwargs["xlabel"] = "axis2"
        if "ldg" in h.lower():
            kwargs["xlabel"] = "Leading jet axis2"
        if "subldg" in h.lower():
            kwargs["xlabel"] = "Subleading jet axis2"
        if "avg" in h.lower():
            kwargs["xlabel"] = "avg axis2"

    if "dijetmass" in h.lower():
        _units  = "GeV/c^{2}"
        _format = "%0.0f " + _units
        kwargs["xlabel"] = "m_{W} (%s)" % _units
        kwargs["xmax"] = 600

    if "ptd" in h.lower():
        _format = "%0.2f"
        kwargs["xlabel"] = "p_{T}D"
        kwargs["xmax"] = 1
        if "ldg" in h.lower():
            kwargs["xlabel"] = "Leading jet p_{T}D"
        elif "subldg" in h.lower():
            kwargs["xlabel"]= "Subleading jet p_{T}D"

    if "ptdr" in h.lower():
        kwargs["xmax"] =800
        _format = "%0.0f"
        kwargs["xlabel"] = "p_{T}#Delta R_{t}"        
        if "dijetptdr" in h.lower():
            kwargs["xlabel"] = "p_{T}#Delta R_{W}"

    if "bdisc" in h.lower():
        _format = "%0.2f"
        kwargs["xlabel"] = "CSV"
        kwargs["xmax"] = 1
        if "subldg" in h.lower():
            kwargs["xlabel"] = "Subleading jet CSV"
        elif "ldg" in h.lower():
            _xlabel = "Leading jet CSV"
        else:
            _xlabel = "b-tagged jet CSV"

    if "over" in h.lower():
         _format = "%0.2f "
         kwargs["xlabel"] = "m_{W}/m_{t}"
         kwargs["xmax"] = 1
         kwargs["xmin"] = 0

    if "eta" in h.lower():
         _format = "%0.2f "
         kwargs["xlabel"] = "|eta|"
         kwargs["xmax"] = 5
         kwargs["xmin"] = -5
        
    if "likelihood" in h.lower():
        _format = "%0.2f"
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
    else:
        pass

    if 0:
        ROOT.gStyle.SetNdivisions(8, "X")
        ROOT.gStyle.SetNdivisions(8, "Y")
        
    return kwargs
    
def CreateHistograms(inputList, ref, kwargs):
    histoMap = {}
    Print("List of branches:", True)

    for i, brName in enumerate(inputList):

        kwargs   = GetHistoKwargs(brName, opts)
        kwargs_ref   = GetHistoKwargs(ref, opts)

        Print("%s. %s" % (i, brName), False)
        histoName = "%s_vs_%s" %(brName, ref)

        xbins = 50 
        xmin = kwargs_ref["xmin"]
        xmax = kwargs_ref["xmax"]
        ybins = 50
        ymin = kwargs["xmin"]
        ymax = kwargs["xmax"]
        
        # Create histogram
        histo_vs_ref = ROOT.TH2F(histoName, histoName, xbins, xmin, xmax, ybins, ymin, ymax)

        histo_vs_ref.GetXaxis().SetTitle(kwargs_ref["xlabel"])
        histo_vs_ref.GetYaxis().SetTitle(kwargs["xlabel"])
        histo_vs_ref.GetZaxis().SetTitle(kwargs["zlabel"])

        histoMap.update({brName : histo_vs_ref})
        
    return histoMap


def GetListOfBranches(tree):
    branchList = []
    brsList = tree.GetListOfBranches()
    nBranches = len(brsList)
    for i in range(nBranches):
        if "Weight" in brsList[i].GetName():
            continue
        # Branches with old names
        if "Trijet" in brsList[i].GetName():
            continue
        branchList.append(brsList[i].GetName())        
    return branchList

def Plot2DHistograms(datasetsMgr, ref, treeName):
    
    # Get ROOT file
    f = ROOT.TFile(opts.rootFileName)
    # Get TTree 
    tree = f.Get(treeName)
    # Get number of entries
    nEntries = tree.GetEntries();

    # Define and initiale tmp variable
    var = ROOT.Double(0.0)
    
    # GetListOfBranches
    inputList = GetListOfBranches(tree)
    
    kwargs = {}
    # Create map with branch names and corresponding TH2 histograms    
    histoMap = CreateHistograms(inputList, ref, kwargs)
    # loop over entries
    for ientry in range(nEntries):    
        tree.GetEntry(ientry)        
        # loop over all branches                        
        for brName in inputList:
            # Get branch value
            # https://www.programiz.com/python-programming/methods/built-in/getattr
            xvalue = getattr(tree, ref)
            yvalue = getattr(tree, brName)
            histoMap[brName].Fill(xvalue, yvalue)
            
    for brName in inputList:
        saveName = "%s_vs_%s" % (brName, ref)
        #
        #p = plots.PlotBase(datasetRootHistos=[histoMap[brName]], saveFormats=opts.saveFormats)
        Make2D(histoMap[brName], saveName, kwargs)
    return

def AddPreliminaryText():
    # Setting up preliminary text
    tex = ROOT.TLatex(0.,0., 'Preliminary');
    tex.SetNDC();
    tex.SetX(0.27);
    tex.SetY(0.95);
    tex.SetTextFont(53);
    tex.SetTextSize(28);
    tex.SetLineWidth(2)
    return tex

def AddCMSText():
    # Settign up cms text
    texcms = ROOT.TLatex(0.,0., 'CMS');
    texcms.SetNDC();
    texcms.SetTextAlign(31);
    texcms.SetX(0.26);
    texcms.SetY(0.95);
    texcms.SetTextFont(63);
    texcms.SetLineWidth(2);
    texcms.SetTextSize(30);
    return texcms

def Make2D(h2D, saveName, kwargs):
    canvas = ROOT.TCanvas()

    canvas.SetLogx(opts.logX)
    canvas.SetLogy(opts.logY)
    canvas.SetLogz(opts.logZ)

    style = tdrstyle.TDRStyle()
    style.setOptStat(False)
    style.setGridX(False)
    style.setGridY(False)

    canvas.SetRightMargin(0.20)
    canvas.SetLeftMargin(0.16)

    #ROOT.gStyle.SetPalette(ROOT.kBlueGreenYellow)

    #h2D.GetYaxis().SetRangeUser(500, 1700)
    #h2D.GetZaxis().SetRangeUser(0.48, 1.1)

    h2D.GetZaxis().SetTitleOffset(1.2)
    h2D.GetYaxis().SetTitleOffset(1.4)
    h2D.SetMarkerColor(ROOT.kBlack)
    h2D.SetMarkerSize(1.5)
    h2D.SetMarkerStyle(8)
    
    texcms  = AddCMSText()
    texpre  = AddPreliminaryText()
    #texlumi = AddLumiText(str(round(lumi/1000.0, 1)))

    h2D.Draw("colz")

    texcms.Draw("same")
    texpre.Draw("same")

    canvas.Modified()
    canvas.Update()
        
    for f in opts.saveFormats:
        savePath = "%s%s%s" % (opts.saveDir, saveName, f)
        canvas.SaveAs(savePath)

    return

'''
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
'''
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
    REFHISTO     = "trijetPt"
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
        #opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="", postfix="TH2")
        # Get username and the initial of the username
        user    = getpass.getuser()
        initial = getpass.getuser()[0]    
        opts.saveDir = "/publicweb/%s/%s/%s/TH2/" % (initial, user, opts.mcrab)
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
