#!/usr/bin/env python
'''
Description:
This scripts creates the distributions of the input variables used for the training of the top-tagger (BDTG) 
for all the truth-matched (signal) and unmatched (background) top-candidates in the ttbar sample

Usage:
./plotSignalBackground.py -m <pseudo_mcrab> [opts]

Examples:
./plotSignalBackground.py -m <peudo_mcrab> -o "" --url --normaliseToOne

Last Used:
 ./plotSignalBackground.py -m TopRecoTree_180717_DeltaR0p3_BJetPt40_TopPtReweighting13Tev/ --normaliseToOne --url
 ./plotSignalBackground.py -m TopRecoTree_180714_DeltaR0p3_DeltaPtOverPt0p32_BJetPt40_TopPtReweighting13Tev/ --normaliseToOne --url
 ./plotSignalBackground.py -m TopRecoTree_200223_181325_hadronic --normaliseToOne --url --dataset ChargedHiggs_HplusTB_HplusToTB_M_500
'''

#================================================================================================ 
# Imports
#================================================================================================ 
import sys
import math
import copy
import os
from optparse import OptionParser
import getpass

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
import HiggsAnalysis.NtupleAnalysis.tools.multicrabConsistencyCheck as consistencyCheck
import HiggsAnalysis.NtupleAnalysis.tools.aux as aux

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

def RemoveDatasets(datasetsMgr):
    # Remove datasets
    datasetsMgr.remove(filter(lambda name: "Data" in name, datasetsMgr.getAllDatasetNames()))
    datasetsMgr.remove(filter(lambda name: "QCD-b" in name, datasetsMgr.getAllDatasetNames()))
    datasetsMgr.remove(filter(lambda name: "QCD" in name, datasetsMgr.getAllDatasetNames()))
    datasetsMgr.remove(filter(lambda name: "SingleTop" in name, datasetsMgr.getAllDatasetNames()))
    datasetsMgr.remove(filter(lambda name: "DYJetsToQQHT" in name, datasetsMgr.getAllDatasetNames()))
    datasetsMgr.remove(filter(lambda name: "TTZToQQ" in name, datasetsMgr.getAllDatasetNames()))
    datasetsMgr.remove(filter(lambda name: "TTWJetsToQQ" in name, datasetsMgr.getAllDatasetNames()))
    datasetsMgr.remove(filter(lambda name: "WJetsToQQ" in name, datasetsMgr.getAllDatasetNames()))
    datasetsMgr.remove(filter(lambda name: "Diboson" in name, datasetsMgr.getAllDatasetNames()))
    datasetsMgr.remove(filter(lambda name: "TTTT" in name, datasetsMgr.getAllDatasetNames()))
    datasetsMgr.remove(filter(lambda name: "FakeBMeasurementTrijetMass" in name, datasetsMgr.getAllDatasetNames()))
    #datasetsMgr.remove(filter(lambda name: "M_" in name and "M_" + str(opts.signalMass) not in name, datasetsMgr.getAllDatasetNames()))
    return datasetsMgr

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
    
def main(opts, signalMass):
    
    # Setup & configure the dataset manager 
    datasetsMgr = GetDatasetsFromDir(opts)
    datasetsMgr.updateNAllEventsToPUWeighted()
    #datasetsMgr.loadLuminosities() # from lumi.json
    if opts.verbose:
        datasetsMgr.PrintCrossSections()
        datasetsMgr.PrintLuminosities()

    # Set/Overwrite cross-sections
    for d in datasetsMgr.getAllDatasets():
        if "ChargedHiggs" in d.getName():
            datasetsMgr.getDataset(d.getName()).setCrossSection(1.0)
                       
    # Determine integrated Lumi before removing data
    #intLumi = datasetsMgr.getDataset("Data").getLuminosity()
    intLumi = opts.intLumi
        
    datasetsMgr = RemoveDatasets(datasetsMgr)
    
    # Merge histograms (see NtupleAnalysis/python/tools/plots.py) 
    plots.mergeRenameReorderForDataMC(datasetsMgr) 

    # Re-order datasets
    datasetOrder = []
    for d in datasetsMgr.getAllDatasets():
        datasetOrder.append(d.getName())
    for m in signalMass:
        datasetOrder.insert(0, m)
    datasetsMgr.selectAndReorder(datasetOrder)

    # Print dataset information
    datasetsMgr.PrintInfo()

    # Apply TDR style
    style = tdrstyle.TDRStyle()
    style.setOptStat(False)
    style.setGridX(opts.gridX)
    style.setGridY(opts.gridY)
    style.setLogX(opts.logX)
    style.setLogY(opts.logY)
    
    # Do the topSelection histos
    folder      = opts.folder 
    histoPaths1 = []
    if folder != "":
        histoList  = datasetsMgr.getDataset(datasetOrder[0]).getDirectoryContent(folder)
        histoPaths1 = [os.path.join(folder, h) for h in histoList]
        
    folderG     = "TrijetCandidateGenuine"
    folderFake     = "TrijetCandidateFake"

    histoListG  = datasetsMgr.getDataset(datasetOrder[0]).getDirectoryContent(folderG)
    histoListF  = datasetsMgr.getDataset(datasetOrder[0]).getDirectoryContent(folderFake)
    
    histoPathsG = [os.path.join(folderG, h) for h in histoListG]
    histoPathsF = [os.path.join(folderFake, h) for h in histoListF]

    for i in range(len(histoPathsG)):
        hG = histoPathsG[i]
        hF = histoPathsF[i]
        if "Vs" in hG: # Skip TH2D
            continue
        if "cjet" in hG.lower():
            continue
        if "all" in hG.lower():
            continue
        #if "mass" not in hG.lower():
        #    continue            
        PlotSignalBackground(datasetsMgr, hG, hF, intLumi)
        
    return

def PlotSignalBackground(datasetsMgr, hG, hF, intLumi):
    _kwargs = {}

    if opts.normaliseToOne:
        pG = plots.MCPlot(datasetsMgr, hG, normalizeToOne=True, saveFormats=[], **_kwargs) # Genuine top candidates (truth-matched)
        pF = plots.MCPlot(datasetsMgr, hF, normalizeToOne=True, saveFormats=[], **_kwargs) # Fake top candidates (unmatched)
    else:
        pG = plots.MCPlot(datasetsMgr, hG, normalizeToLumi=intLumi, saveFormats=[], **_kwargs) # Genuine top candidates (truth-matched)
        pF = plots.MCPlot(datasetsMgr, hF, normalizeToLumi=intLumi, saveFormats=[], **_kwargs) # Fake top candidates (unmatched)         
        
    myList = []
    # Customise styling
    pG.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetLineStyle(ROOT.kSolid))
    pF.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetLineStyle(ROOT.kSolid))

    # Get dataset
    dataset = datasetsMgr.getDataset(opts.dataset)

    # Get Genuine-top histogram
    h = dataset.getDatasetRootHisto(hG)
    h.normalizeToOne()
    HG = h.getHistogram()

    # Get Fake-top histogram
    h = dataset.getDatasetRootHisto(hF)
    h.normalizeToOne()
    HF = h.getHistogram()
    
    #Define Signal style
    signalStyle     = styles.StyleCompound([styles.StyleMarker(markerSize=1.2, markerColor=ROOT.kAzure-3, markerSizes=None, markerStyle=ROOT.kFullDiamond),
                                            styles.StyleLine(lineColor=ROOT.kAzure-3, lineStyle=ROOT.kSolid, lineWidth=3),
                                            styles.StyleFill(fillColor=ROOT.kAzure-3)])
    #Define Background style
    backgroundStyle = styles.StyleCompound([styles.StyleMarker(markerSize=1.2, markerColor=ROOT.kRed-4, markerSizes=None, markerStyle=ROOT.kFullDiamond),
                                            styles.StyleLine(lineColor=ROOT.kRed-4, lineStyle=ROOT.kSolid, lineWidth=3)])
    
    p = plots.ComparisonManyPlot(histograms.Histo(HF,"Fake", "p", "P"), [histograms.Histo(HG,"Genuine", "pl", "PL")], saveFormats=[])

    #Set labels
    p.histoMgr.setHistoLegendLabelMany({"Fake": "Unmatched", "Genuine": "Truth-matched"}) 
    
    if "ChargedHiggs" in opts.dataset:
        datasetName = opts.dataset
        datasetName = datasetName.replace("ChargedHiggs_HplusTB_HplusToTB_M_", "")
        datasetName = "m_{H^{#pm}} = %s GeV" % (datasetName)
    else:
        datasetName = plots._legendLabels[opts.dataset]
    
    if "hadronic" in opts.mcrab.lower():
        extra_text = "hadronic"# (%s)" % datasetName
    elif "semileptonic" in opts.mcrab.lower():
        extra_text = "semi leptonic"# (%s)" % datasetName
    else:
        extra_text = ""

    xleg = 0.595#0.665
    if (extra_text != ""):
        if (opts.dataset == "TT"):
             p.appendPlotObject(histograms.PlotText(xleg, 0.57, "%s (%s)" % (datasetName, extra_text), bold=True, size=20))
        else:
            p.appendPlotObject(histograms.PlotText(xleg, 0.57, "%s" % (datasetName), bold=True, size=20))
    else:
        p.appendPlotObject(histograms.PlotText(0.67, 0.83, plots._legendLabels[opts.dataset], bold=True, size=20))
    #Set Draw style    
    p.histoMgr.forHisto("Fake", backgroundStyle ) 
    p.histoMgr.setHistoDrawStyle("Fake", "HIST") 
    p.histoMgr.setHistoLegendStyle("Fake", "L") #F

    p.histoMgr.forHisto("Genuine", signalStyle) 
    p.histoMgr.setHistoDrawStyle("Genuine" , "HIST") 
    p.histoMgr.setHistoLegendStyle("Genuine", "L") #LP
     
    kwargs   = GetHistoKwargs(hG, opts)        
    kwargs["ylabel"] = "Arbitrary Units / %s" % str(kwargs["format"])
    #kwargs["createLegend"] = {"x1": 0.58, "y1": 0.65, "x2": 0.92, "y2": 0.82}
    kwargs["createLegend"] = {"x1": 0.58, "y1": 0.65, "x2": 0.92, "y2": 0.8}
    kwargs["opts"]   = {"ymin": 0.0, "ymaxfactor": 1.1, "xmin" : kwargs["xmin"], "xmax" : kwargs["xmax"]}

    # Save plot in all formats    
    saveName = hG.split("/")[-1]
    plots.drawPlot(p, saveName, **kwargs)
    SavePlot(p, saveName, os.path.join(opts.saveDir, opts.dataset), opts.saveFormats )

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

    if opts.normaliseToOne:
        yLabel  = "Arbitrary Units"
        yMin    = None #5e-4
        yMax    = None #1e-1
        yMaxF   = 1.1
    else:
        yLabel  = "Events"
        yMin    = 1e0
        yMaxF   = 1
        yMaxF   = 1

    xLabel = "xlabel"

    cutBox      = {"cutValue": 400.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True} #box = True works
    cutBoxY     = {"cutValue": 200.0, "fillColor": 16, "box": False, "line": False, "greaterThan": True}
                       
    _format = "%0.0f "
    
    kwargs = {
        "stackMCHistograms": False,
        "addMCUncertainty" : False,
        "addLuminosityText": (not opts.normaliseToOne),
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
        kwargs["xmax"] = 800
    if "pt" in h.lower():
        _units  = "GeV/c"
        kwargs['format'] = "%0.0f " + _units
        kwargs["xlabel"] = "p_{T} (%s)" % _units
        if "trijetpt" in h.lower():
            kwargs["xlabel"] = "p_{T,t} (%s)" % _units
        if "dijetpt" in h.lower():
            kwargs["xlabel"] = "p_{T,w} (%s)" % _units
        if "ldgjetpt" in h.lower():
            kwargs["xmax"] = 200
            kwargs["xlabel"] = "leading jet p_{T} (%s)" % _units
        if "subldgjetpt" in h.lower():
            kwargs["xmax"] = 200
            kwargs["xlabel"] = "subleading jet p_{T} (%s)" % _units
        if "bjetpt" in h.lower():
            kwargs["xlabel"] = "bjet p_{T} (%s)" % _units

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
            kwargs["xlabel"] = "avg mult"

    if "cvsl" in h.lower():
        kwargs['format'] = "%0.2f"
        kwargs["xmax"] = 1
        kwargs["xmin"] = -1
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

    if "Eta" in h:#.lower():
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
             kwargs["xmin"] = 0

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
        
    if 0:
        ROOT.gStyle.SetNdivisions(8, "X")
        ROOT.gStyle.SetNdivisions(8, "Y")
        
    return kwargs
    
def SavePlot(plot, plotName, saveDir, saveFormats = [".C", ".png", ".pdf"]):
    Verbose("Saving the plot in %s formats: %s" % (len(saveFormats), ", ".join(saveFormats) ) )

    # Check that path exists
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
    OPTMODE      = ""
    BATCHMODE    = True
    PRECISION    = 3
    SIGNALMASS   = []
    INTLUMI      = 35920
    SUBCOUNTERS  = False
    LATEX        = False
    MERGEEWK     = False
    URL          = False
    NOERROR      = True
    SAVEDIR      = None
    VERBOSE      = False
    HISTOLEVEL   = "Vital" # 'Vital' , 'Informative' , 'Debug'
    NORMALISE    = False
    FOLDER       = "" #"topSelection_" #"ForDataDrivenCtrlPlots" #"topologySelection_"
    DATASET      = "TT"
    GRIDX        = False
    GRIDY        = False
    LOGX         = False
    LOGY         = False
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

    #parser.add_option("--signalMass", dest="signalMass", type=float, default=SIGNALMASS, 
                      #help="Mass value of signal to use [default: %s]" % SIGNALMASS)

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

    parser.add_option("--dataset", dest="dataset", type="string", default = DATASET,
                      help="Name of disired dataset [default: %s]" % (DATASET) )

    parser.add_option("--gridX", dest="gridX", action="store_true", default=GRIDX, 
                      help="Enable the x-axis grid lines [default: %s]" % GRIDX)

    parser.add_option("--gridY", dest="gridY", action="store_true", default=GRIDY, 
                      help="Enable the y-axis grid lines [default: %s]" % GRIDY)

    parser.add_option("--logX", dest="logX", action="store_true", default=LOGX, 
                      help="Set x-axis to logarithm scale [default: %s]" % LOGX)

    parser.add_option("--logY", dest="logY", action="store_true", default=LOGY,
                      help="Set y-axis to logarithm scale [default: %s]" % LOGY)

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

    # Sanity check
    if opts.mergeEWK:
        Print("Merging EWK samples into a single Datasets \"EWK\"", True)

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
        opts.saveDir = "/publicweb/%s/%s/%s" % (initial, user, opts.mcrab)
    if not os.path.exists(opts.saveDir):
        os.makedirs(opts.saveDir)

    # Sanity check
    allowedMass = [180, 200, 220, 250, 300, 350, 400, 500, 800, 1000, 2000, 3000]
    signalMass = []
    for m in sorted(SIGNALMASS, reverse=True):
        signalMass.append("ChargedHiggs_HplusTB_HplusToTB_M_%.f" % m)

    # Call the main function
    main(opts, signalMass)

    if not opts.batchMode:
        raw_input("=== plotSignalBackground.py: Press any key to quit ROOT ...")
