#!/usr/bin/env python
'''
Description:

Usage:
./plotComparisons.py  -m <pseudo_mcrab> --mref <pseudo_mcrab_reference>[opts]

Examples:
./plotComparisons.py  -m TopTaggerEfficiency_180609_TopRecoBugFixed/ --mref TopTaggerEfficiency_180608_MassCut300_BeforeBugFix/ --url -v

Last Used:
./plotComparisons.py -m TopRecoTree_SemiLeptonic_EleOrMuon --mref TopRecoTree_multFloat_AnalysisSelections --url -v -i TT --normaliseToOne -s png --folder TrijetCandidate

'''

#================================================================================================ 
# Imports
#================================================================================================ 
import sys
import math
import copy
import os
import getpass
from optparse import OptionParser
import HiggsAnalysis.NtupleAnalysis.tools.aux as aux

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

def GetDatasetsFromDir(opts, mcrab):
    Verbose("Getting datasets")
    
    if (not opts.includeOnlyTasks and not opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode, 
                                                        analysisName=opts.analysisName,
                                                        optimizationMode=opts.optMode)
    elif (opts.includeOnlyTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode,
                                                        analysisName=opts.analysisName,
                                                        includeOnlyTasks=opts.includeOnlyTasks,
                                                        optimizationMode=opts.optMode)
    elif (opts.excludeTasks):
        datasets = dataset.getDatasetsFromMulticrabDirs([mcrab],
                                                        dataEra=opts.dataEra,
                                                        searchMode=opts.searchMode,
                                                        analysisName=opts.analysisName,
                                                        excludeTasks=opts.excludeTasks,
                                                        optimizationMode=opts.optMode)
    else:
        raise Exception("This should never be reached")

    datasets.updateNAllEventsToPUWeighted()
    datasets.loadLuminosities() # from lumi.json
    return datasets
    

def main(opts, signalMass):
               
    # Setup & configure the dataset manager 
    datasetsMgr = GetDatasetsFromDir(opts, opts.mcrab)
    datasetsMgrRef = GetDatasetsFromDir(opts, opts.mref)

    # Determine integrated Lumi before removing data
    #intLumi = datasetsMgr.getDataset("Data").getLuminosity()
    intLumi = opts.intLumi

    # Merge histograms (see NtupleAnalysis/python/tools/plots.py) 
    plots.mergeRenameReorderForDataMC(datasetsMgr) 
    plots.mergeRenameReorderForDataMC(datasetsMgrRef) 

    # Re-order datasets
    datasetOrder = []
    for d in datasetsMgr.getAllDatasets():
        if "M_" in d.getName():
            if d not in signalMass:
                continue
        datasetOrder.append(d.getName())
    for m in signalMass:
        datasetOrder.insert(0, m)
    datasetsMgr.selectAndReorder(datasetOrder)
    datasetsMgrRef.selectAndReorder(datasetOrder)

    # Print dataset information
    datasetsMgr.PrintInfo()

    # Apply TDR style
    style = tdrstyle.TDRStyle()
    style.setOptStat(True)
    style.setGridX(False)
    style.setGridY(False)

    # Get folder
    folderTarg     = opts.folder
    folderRef      = opts.folder

    # Get histo names
    histoListTarg = datasetsMgr.getDataset(datasetOrder[0]).getDirectoryContent(folderTarg)
    histoListRef  = datasetsMgrRef.getDataset(datasetOrder[0]).getDirectoryContent(folderRef)

    histoPathsTarg = [os.path.join(folderTarg, h) for h in histoListTarg]
    histoPathsRef  = [os.path.join(folderRef, h) for h in histoListRef]
    #histoPathsTarg = [os.path.join(folderTarg, h) for h in ["Ctrl_lepGenTop_Pt"]] #Ctrl_hadrGenTop_Pt
    #histoPathsRef  = [os.path.join(folderRef, h) for h in ["GenTop_Pt"]] #GenTop_Pt
    #histoPathsTarg = [os.path.join(folderTarg, h) for h in ["LdgJetCvsL", "SubldgJetCvsL"]] #Ctrl_hadrGenTop_Pt
    #histoPathsRef  = [os.path.join(folderRef, h) for h in ["LdgJetCvsL", "SubldgJetCvsL"]]
    for i in range(len(histoPathsRef)):            
        _continue = False
        hR = histoPathsRef[i]
        exclude = ["weighted", "alljet", "ctrl", "quarkjet", "matched"]
        for e in exclude:
            if e in hR.lower():
                _continue = True
        if _continue:
            continue
        hT = hR #histoPathsTarg[i]
        Plot_Comparisons(datasetsMgr, datasetsMgrRef, hT, hR, intLumi)
    return

def Plot_Comparisons(datasetsMgr, datasetsMgrRef, hT, hR, intLumi):
    _kwargs = {}

    #print 
    if opts.normaliseToOne:
        pT = plots.MCPlot(datasetsMgr, hT, normalizeToOne=True, saveFormats=[], **_kwargs)
        pR = plots.MCPlot(datasetsMgrRef, hR, normalizeToOne=True, saveFormats=[], **_kwargs)
    else:
        pT = plots.MCPlot(datasetsMgr, hT, normalizeToLumi=intLumi, saveFormats=[], **_kwargs)
        pR = plots.MCPlot(datasetsMgrRef, hR, normalizeToLumi=intLumi, saveFormats=[], **_kwargs)

    # Draw the histograms                                                                                                                                                                
    _cutBox = None
    _rebinX = 1
    _format = "%0.0f"
    _xlabel = None
    logY    = False
    _opts   = {"ymin": 1e-3, "ymaxfactor": 1.0}

    #print "HT", hT
    #print "HR", hR

    xMax = 200
    _units = 0
    if "pt" in hT.lower():
        _units = "GeV/c"

    _ylabel = "Events"
    _log = True
    if opts.normaliseToOne:
        _ylabel = "Arbitrary Units"
        _log = False
    
    if "mass" in hR.lower():
        _units  = "GeV/c^{2}"
        _format = "%0.0f " + _units
        _xlabel = "M (%s)" % _units

    if "trijetmass" in hR.lower():
        _units  = "GeV/c^{2}"
        _format = "%0.0f " + _units
        _xlabel = "m_{top} (%s)" % _units
        _opts["xmax"] = 805 #1005
        
    if "bjetmass" in hR.lower():
        _xlabel = "m_{b-tagged jet} (%s)" % _units
        _opts["xmax"] = 50
    if "bjetldgjet_mass" in hR.lower():
        _xlabel = "m_{b, ldg jet} (%s)" % _units
        _opts["xmax"] = 705
    if "bjetsubldgjet_mass" in hR.lower():
        _xlabel = "m_{b-tagged, subldg jet} (%s)" % _units
        _opts["xmax"] = 705
    if "jet_mass" in hR.lower():
        _opts["xmax"] = 750

    if "mult" in hR.lower():
        _format = "%0.0f"
        if "ldg" in hR.lower():
            _xlabel = "Leading jet mult"
        if "subldg" in hR.lower():
            _xlabel = "Subleading jet mult"
        if "avg" in hR.lower():
            _xlabel = "avg CvsL"

    if "cvsl" in hR.lower():
        _format = "%0.2f"
        if "ldg" in hR.lower():
            _xlabel = "Leading jet CvsL"
        if "subldg" in hR.lower():
             _xlabel = "Subleading jet CvsL"
        if "avg" in hR.lower():
            _xlabel = "avg CvsL"
    if "axis2" in hR.lower():
        _format = "%0.3f"
        if "ldg" in hR.lower():
            _xlabel = "Leading jet axis2"
        if "subldg" in hR.lower():
            _xlabel = "Subleading jet axis2"
        if "avg" in hR.lower():
            _xlabel = "avg axis2"

    if "dijetmass" in hR.lower():
        _units  = "GeV/c^{2}"
        _format = "%0.0f " + _units
        _xlabel = "m_{W} (%s)" % _units
        _opts["xmax"] = 600
        #_opts   = {"xmin": 0.0, "xmax": 605, "ymin": 1e-3, "ymaxfactor": 1.0}

    if "trijetptdr" in hR.lower():
        _opts["xmax"] =800
        _format = "%0.0f"
        _xlabel = "p_{T}#Delta R_{t}"

    if "dijetptdr" in hR.lower():
        _opts["xmax"] =800
        _format = "%0.0f"
        _xlabel = "p_{T}#Delta R_{W}"
    if "dgjetptd" in hR.lower():
        _format = "%0.2f"
        _xlabel = "Leading jet p_{T}D"
        if "subldg" in hR.lower():
            _xlabel = "Subleading jet p_{T}D"

    if "bdisc" in hR.lower():
        _format = "%0.2f"
        if "subldg" in hR.lower():
            _xlabel = "Subleading jet CSV"
            _opts   = {"ymin": 1e-3, "ymax": 0.06}
        elif "ldg" in hR.lower():
            _xlabel = "Leading jet CSV"
            _opts   = {"ymin": 1e-3, "ymax": 0.06}
        else:
            _xlabel = "b-tagged jet CSV"
            _opts   = {"ymin": 1e-3, "ymax": 0.35}

    if "over" in hR.lower():
         _format = "%0.2f "
         _xlabel = "m_{W}/m_{t}"
         _opts["xmax"] = 1
         _opts["xmin"] = 0
    if "likelihood" in hR.lower():
        _format = "%0.2f"
        if "ldg" in hR.lower():
            _xlabel = "Leading jet QGL"
        if "subldg" in hR.lower():
            _xlabel = "Subleading jet QGL"
        if "avg" in hR.lower():
            _xlabel = "avg QGL"

    else:
        pass

    if logY:
        yMaxFactor = 2.0
    else:
        yMaxFactor = 1.5

    _opts["ymaxfactor"] = yMaxFactor
    if opts.normaliseToOne:
        _opts["ymin"] = 1e-3
    else:
        _opts["ymin"] = 1e0
    
    _kwargs = {
        "xlabel"           : _xlabel,
        "ylabel"           : "%s / %s" % (_ylabel, _format),
        "ratioYlabel"      : "Ratio ",
        "ratio"            : opts.ratio,
        "ratioInvert"      : False,
        "stackMCHistograms": False,
        "addMCUncertainty" : False,
        "addLuminosityText": False,
        "addCmsText"       : True,
        "cmsExtraText"     : "Preliminary",
        "opts"             : {"ymin": 0.0, "ymaxfactor": 1.2},
        #"opts"             : {"ymin": 0.1  ,"ymaxfactor": 1.2, "xmin":0, "xmax":xMax},
        "opts2"            : {"ymin": 0.1, "ymax": 1.9},
        "log"              : _log,
        "createLegend"     : {"x1": 0.58, "y1": 0.65, "x2": 0.92, "y2": 0.82}, 
        "rebinX"           : 1,
        }
    
    myList = []
    # Customise styling
    pT.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetLineStyle(ROOT.kSolid))
    pR.histoMgr.forEachHisto(lambda h: h.getRootHisto().SetLineStyle(ROOT.kSolid))

    datasetTarg = datasetsMgr.getDataset("TT") #ChargedHiggs_HplusTB_HplusToTB_M_500
    datasetRef  = datasetsMgrRef.getDataset("TT") #ChargedHiggs_HplusTB_HplusToTB_M_500
    
    ht = datasetTarg.getDatasetRootHisto(hT)
    hr = datasetRef.getDatasetRootHisto(hR)
    if (opts.normaliseToOne):
        ht.normalizeToOne()
        hr.normalizeToOne()
    else:
        ht.normalizeToLuminosity(intLumi)
        hr.normalizeToLuminosity(intLumi)
    HT = ht.getHistogram()
    HR = hr.getHistogram()

    #print "target: integral", HT.Integral()
    #print "reference: integral", HR.Integral()

    # styles
    altSignalBDTGStyle  = styles.StyleCompound([styles.StyleMarker(markerSize=1.2, markerColor=ROOT.kAzure+9, markerSizes=None, markerStyle=ROOT.kFullDiamond),
                                                styles.StyleLine(lineColor=ROOT.kAzure+9, lineStyle=ROOT.kSolid, lineWidth=3),
                                                styles.StyleFill(fillColor=ROOT.kAzure+9, fillStyle=3001)])

    altBackgroundBDTGStyle = styles.StyleCompound([styles.StyleMarker(markerSize=1.2, markerColor=ROOT.kRed-4, markerSizes=None, markerStyle=ROOT.kFullDiamond),
                                                   styles.StyleLine(lineColor=ROOT.kRed-4, lineStyle=ROOT.kSolid, lineWidth=3),
                                                   ])

    backgroundStyle = styles.StyleCompound([styles.StyleMarker(markerSize=1.2, markerColor=ROOT.kMagenta+3, markerSizes=None, markerStyle=ROOT.kFullDiamond),
                                                   styles.StyleLine(lineColor=ROOT.kMagenta+3, lineStyle=ROOT.kSolid, lineWidth=3),
                                                   ])

    signalStyle = styles.StyleCompound([styles.StyleMarker(markerSize=1.2, markerColor=ROOT.kGreen+1, markerSizes=None, markerStyle=ROOT.kFullTriangleUp), # kTeal+2
                                        styles.StyleLine(lineColor=ROOT.kGreen+1, lineStyle=ROOT.kSolid, lineWidth=3),
                                        styles.StyleFill(fillColor=ROOT.kGreen+1, fillStyle=3001)])


    p = plots.ComparisonPlot(histograms.Histo(HT,"Target", "pl", "PL"), histograms.Histo(HR,"Reference", "pl", "PL"),) 
    p.histoMgr.setHistoLegendLabelMany({"Target": "semi-leptonic (t#bar{t})", "Reference": "hadronic (t#bar{t})"}) 

    if ("TT" in datasetTarg.getName()):    
        p.histoMgr.forHisto("Target", backgroundStyle ) 
        p.histoMgr.forHisto("Reference", signalStyle) 
    elif ("QCD" in datasetTarg.getName()):
        p.histoMgr.forHisto("Target"  , altBackgroundBDTGStyle)#styles.getABCDStyle("VR"))
        p.histoMgr.forHisto("Reference", styles.qcdFillStyle)
    elif ("Charged" in datasetTarg.getName()):
        p.histoMgr.forHisto("Target", altBackgroundBDTGStyle) 
        p.histoMgr.forHisto("Reference", signalStyle)

    p.histoMgr.setHistoDrawStyle("Target", "HIST") 
    p.histoMgr.setHistoLegendStyle("Target", "LP") #F

    p.histoMgr.setHistoDrawStyle("Reference" , "HIST") 
    p.histoMgr.setHistoLegendStyle("Reference", "F") #LP

    dName = datasetTarg.getName()
    dName = dName.replace("ChargedHiggs_HplusTB_HplusToTB_", "")

    saveName = hT.split("/")[-1] + "_" + dName
    savePath = os.path.join(opts.saveDir, hT.split("/")[0], opts.optMode)

    plots.drawPlot(p, savePath, **_kwargs)

    #leg = ROOT.TLegend(0.2, 0.8, 0.81, 0.87)
    leg = ROOT.TLegend(0.2, 0.8, 0.51, 0.87)
    leg.SetFillStyle( 0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    datasetName = datasetTarg.getName()
    datasetName = datasetName.replace("TT", "t#bar{t}")
    if "ChargedHiggs" in datasetName:
        datasetName = datasetName.replace("ChargedHiggs_HplusTB_HplusToTB_M_", "m_{H+} = ")
        datasetName = datasetName+"GeV"
    SavePlot(p, saveName, savePath)
    return
def SavePlot(plot, plotName, saveDir):#, saveFormats = [".C", ".png", ".pdf"]):
    Verbose("Saving the plot in %s formats: %s" % (len(opts.saveFormats), ", ".join(opts.saveFormats) ) )

     # Check that path exists
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    # Create the name under which plot will be saved
    saveName = os.path.join(saveDir, plotName.replace("/", "_"))
    saveName = saveName.replace("(", "_")
    saveName = saveName.replace(")", "")
    saveName = saveName.replace(" ", "")

    # For-loop: All save formats
    for i, ext in enumerate(opts.saveFormats):
        saveNameURL = saveName + ext
        saveNameURL = aux.convertToURL(saveNameURL, opts.url)
        Verbose(saveNameURL, i==0)
        plot.saveAs(saveName, formats=opts.saveFormats)
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
    SAVEDIR      = None #"/publicweb/s/skonstan/" + ANALYSISNAME
    VERBOSE      = False
    HISTOLEVEL   = "Vital" # 'Vital' , 'Informative' , 'Debug'
    NORMALISE    = False
    FOLDER       = "TrijetCandidate" #"topSelection_" #"ForDataDrivenCtrlPlots" #"topologySelection_"
    SAVEFORMATS  = "png"
    RATIO        = False
    # Define the available script options
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("-m", "--mcrab", dest="mcrab", action="store", 
                      help="Path to the multicrab directory for input")

    parser.add_option("--ref", "--mref", dest="mref", action="store", 
                      help="Path to the multicrab directory for input")

    parser.add_option("-o", "--optMode", dest="optMode", type="string", default=OPTMODE, 
                      help="The optimization mode when analysis variation is enabled  [default: %s]" % OPTMODE)

    parser.add_option("-b", "--batchMode", dest="batchMode", action="store_false", default=BATCHMODE, 
                      help="Enables batch mode (canvas creation does NOT generate a window) [default: %s]" % BATCHMODE)

    parser.add_option( "--ratio", dest="ratio", action="store_true", default=RATIO, 
                       help="Plot the ratio of the compared distributions [default: %s]" % RATIO)

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
    opts.saveFormats = [s for s in opts.saveFormats]

    for i, s in enumerate(opts.saveFormats):
        opts.saveFormats[i] = ".%s" % (s)
        print i, s, opts.saveFormats[i]

    # Sanity check
    if opts.mergeEWK:
        Print("Merging EWK samples into a single Datasets \"EWK\"", True)

    if opts.saveDir == None:
        opts.saveDir = aux.getSaveDirPath(opts.mcrab, prefix="", postfix="")

    # Sanity check
    allowedMass = [180, 200, 220, 250, 300, 350, 400, 500, 800, 1000, 2000, 3000]
    signalMass = []
    for m in sorted(SIGNALMASS, reverse=True):
        signalMass.append("ChargedHiggs_HplusTB_HplusToTB_M_%.f" % m)

    # Call the main function
    main(opts, signalMass)

    if not opts.batchMode:
        raw_input("=== plotMC_HPlusMass.py: Press any key to quit ROOT ...")

