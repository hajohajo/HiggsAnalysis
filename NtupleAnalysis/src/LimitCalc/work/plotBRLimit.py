#!/usr/bin/env python
'''
DESCRIPTION:
This is the swiss pocket knife for running Lands/Combine on large array of datacards

INSTRUCTIONS:
Produce limits using the Combine tool. The limits.json file will be the input of this
script to produced the plot

USAGE:
cd datacards_test4b/CombineResults_taujets_170913_192047
../../plotBRLimit.py [opts]


EXAMPLES:
../../../plotBRLimit.py --opaqueLegend --cutLineX 500 --gridX --gridY --yMin 1e-3 --yMax 10 --subdir StatOnly --url
../../../plotBRLimit.py --opaqueLegend --gridX --gridY --yMin 1e-1 --yMax 50 --subdir StatOnly --url
../../../plotBRLimit.py --opaqueLegend --yMin 0.1 --yMax 30 --subdir BDT0p40_Binning12_28Apr2018
/plotBRLimit_Hplus2tb.py --dir datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_190726_150229_autoMCStats/CombineResults_taujets_190726_171059/ --saveDir /publicweb/a/aattikis/Combine/ --cutLineX 300,500,600,700 --opaqueLegend--url --gridX --gridY --yMin 1e-2 --yMax 1


LAST USED:
./plotBRLimit.py --dir datacards_HToHW_EraRun2016_Search80to1000_OptNominal_limits2019_MC_mH300to700_190726_150229_autoMCStats/CombineResults_taujets_190726_171059/ --saveDir /publicweb/a/aattikis/Combine/ --opaqueLegend --yMin 1e-2 --yMax 10 --analysisType "HToHW"

'''

#================================================================================================
# Import modules
#================================================================================================
import os
import getpass
import sys
from optparse import OptionParser

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

import HiggsAnalysis.NtupleAnalysis.tools.histograms as histograms
import HiggsAnalysis.NtupleAnalysis.tools.tdrstyle as tdrstyle
import HiggsAnalysis.NtupleAnalysis.tools.plots as plots
import HiggsAnalysis.NtupleAnalysis.tools.styles as styles
import HiggsAnalysis.NtupleAnalysis.tools.aux as aux
import HiggsAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles
import HiggsAnalysis.LimitCalc.limit as limit

#================================================================================================
# Shell Types
#================================================================================================
sh_e = ShellStyles.ErrorStyle()
sh_s = ShellStyles.SuccessStyle()
sh_h = ShellStyles.HighlightStyle()
sh_a = ShellStyles.HighlightAltStyle()
sh_l = ShellStyles.AltStyle()
sh_t = ShellStyles.NoteStyle()
sh_n = ShellStyles.NormalStyle()
sh_w = ShellStyles.WarningStyle()

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


def main(opts):
    if opts.saveDir =="":
        opts.saveDir = os.getcwd()

    # Assume by default that the observed limit should be blinded
    if not opts.unblinded:
        msg  = "Working in %s mode. " % (sh_t + "Blinded" + sh_n)
        msg += "Only expected limits will be shown."
    else:
        msg  = "Working in %s mode. " % (sh_t + "Unblinded" + sh_n)
        msg += "Both expected and observed limits will be shown."
    Print(msg, True)

    Verbose("Load module for reading the BR limits from JSON file produced by Combine", True)
    limits = limit.BRLimits(opts.directory, opts.analysisType, opts.excludeMassPoints, opts.limitsFile, opts.cfgFile)

    # Enable OpenGL
    if opts.opaqueLegend:
        ROOT.gEnv.SetValue("OpenGL.CanvasPreferGL", 1)

    # Apply TDR style
    style = tdrstyle.TDRStyle()
    if not limits.isHeavyStatus:
        # Give more space for four digits on the y axis labels
        style.tdrStyle.SetPadLeftMargin(0.19)
        style.tdrStyle.SetTitleYOffset(1.6)

    # Set the paper mode
    if opts.paper:
        histograms.cmsTextMode = histograms.CMSMode.PAPER
    else:
        histograms.cmsTextMode = histograms.CMSMode.UNPUBLISHED

    Verbose("Doing the limits plots (linear scale)", True)
    doBRlimit(limits, opts.unblinded, opts, logy=False)

    Verbose("Doing the limits plots (log scale)", True)
    doBRlimit(limits, opts.unblinded, opts, logy=True)

    Verbose("Doing the limit errors", True)    
    doLimitError(limits, opts.unblinded)

    # Print the Limits?
    limits.printLimits(unblindedStatus=opts.unblinded, nDigits=opts.digits)
    
    # Save the Limits in a LaTeX table file
    limits.saveAsLatexTable(unblindedStatus=opts.unblinded, nDigits=opts.digits, savePath=os.path.join(opts.saveDir, opts.subdir) )

    # inform user of output location
    Print("Plots saved under directory %s"% (sh_s + aux.convertToURL(opts.saveDir, opts.url) + sh_n), True)
    return


def doBRlimit(limits, unblindedStatus, opts, logy=False):
    '''
    See https://twiki.cern.ch/twiki/bin/viewauth/CMS/Internal/FigGuidelines
    '''
    graphs = []
    if unblindedStatus:
        gr = limits.observedGraph()
        if gr != None:
            gr.SetPoint(gr.GetN()-1, gr.GetX()[gr.GetN()-1]-1e-10, gr.GetY()[gr.GetN()-1])
            if opts.opaqueLegend:
                graphs.append(histograms.HistoGraph(gr, "Observed", drawStyle="PL", legendStyle=None))
                excluded = gr.Clone()
                excluded.SetPoint(excluded.GetN(), excluded.GetX()[excluded.GetN()-1], 0.05)
                excluded.SetPoint(excluded.GetN(), excluded.GetX()[0], 0.05)
                limit.setExcludedStyle(excluded)
                graphs.append(histograms.HistoGraph(excluded, "Excluded", drawStyle="F", legendStyle="lpf", legendLabel="Observed"))
            else:
                graphs.append(histograms.HistoGraph(gr, "Observed", drawStyle="PL", legendStyle="lp"))

    # Add the expected lines
    graphs.extend([
            histograms.HistoGraph(limits.expectedGraph(), "Expected", drawStyle="L"),
            histograms.HistoGraph(limits.expectedBandGraph(sigma=1), "Expected1", drawStyle="F", legendStyle="f"), #fl
            histograms.HistoGraph(limits.expectedBandGraph(sigma=2), "Expected2", drawStyle="F", legendStyle="f"), #fl
            ])

    # Plot the TGraphs
    plot = plots.PlotBase(graphs, saveFormats=[])
    plot.setLuminosity(limits.getLuminosity())

    # Customise legend (NOTE: font[42]{} needed to remove bold style)
    # https://root.cern.ch/doc/master/classTAttText.html#T4)
    plot.setLegendHeader("#font[42]{95% CL upper limits}") 
    plot.histoMgr.setHistoLegendLabelMany({
            "Expected" : "#font[42]{Median expected}",
            "Expected1": "#font[42]{68% expected}",
            "Expected2": "#font[42]{95% expected}"
            })

    # Create legend
    xPos   = 0.53
    legend = getLegend(opts, xPos)
    plot.setLegend(legend)

    # Get y-min, y-max, and histogram name to be saved as
    ymin, ymax, saveName = getYMinMaxAndName(limits, "limitsBr", logy, opts)
    if opts.yMin != -1:
        ymin = opts.yMin
    if opts.yMax != -1:
        ymax = opts.yMax
    
    if len(limits.mass) == 1:
        plot.createFrame(saveName, opts={"xmin": limits.mass[0]-5.0, "xmax": limits.mass[0]+5.0, "ymin": ymin, "ymax": ymax})
    else:
        plot.createFrame(saveName, opts={"ymin": ymin, "ymax": ymax})

    # Add x-axis cut boxes
    for x in opts.cutLinesX:
        kwargs = {"greaterThan": True}
        plot.addCutBoxAndLine(cutValue=int(x), fillColor=ROOT.kRed, box=False, line=True, **kwargs)

    # Set x-axis title
    plot.frame.GetXaxis().SetTitle(limit.mHplus()) 

    # Set y-axis title
    plot.frame.GetYaxis().SetTitle(limits.getSigmaBRlimitText())

    # Enable/Disable logscale for axes
    if logy:
        plot.getPad().SetLogy(logy)
        plot.getPad().SetLogx(opts.logx)

    # Enable grids in x and y?
    plot.getPad().SetGridx(opts.gridX)
    plot.getPad().SetGridy(opts.gridY)

    # Draw the plot with standard texts
    plot.draw()    

    # Add standard CMS/Lumi/Preliminary text on canvas
    plot.addStandardTexts(addLuminosityText=True)

    # Add physics-related text on canvas
    addPhysicsText(histograms, limits, x=xPos)
    
    # Save the canvas
    plot.save()

    # Save the plots
    SavePlot(plot, saveName, os.path.join(opts.saveDir, opts.subdir), opts.saveFormats)
    return


def SavePlot(plot, plotName, saveDir, saveFormats = [".png", ".pdf", ".C"]):

    # Check that path exists
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    # Create the name under which plot will be saved
    saveName = os.path.join(saveDir, plotName.replace("/", "_"))

    # For-loop: All save formats
    for i, ext in enumerate(saveFormats, 0):
        saveNameURL = saveName + ext
        if "afs" in saveNameURL: #lxplus
            saveNameURL = saveNameURL.replace("/afs/cern.ch/user/a/attikis/public/html/", "https://cmsdoc.cern.ch/~attikis/")
        else: #lpc
            saveNameURL = saveNameURL.replace("/publicweb/a/aattikis/", "http://home.fnal.gov/~aattikis/")

        if opts.url:
            msg = saveNameURL
        else:
            msg = saveName + ext
        Verbose(msg, i==0)
        plot.saveAs(saveName, formats=saveFormats)
    return


def getLegend(opts, xLeg1=0.53):
    dy = -0.14

    # Create customised legend
    xLeg2 = 0.93
    yLeg1 = 0.70 + dy
    yLeg2 = 0.91 + dy
    if opts.unblinded:
        yLeg2 = 0.92 + dy

    # Adjust legend slightly to visually align with text
    legend = histograms.createLegend(xLeg1*.98, yLeg1, xLeg2, yLeg2) 
    legend.SetMargin(0.17)

    # Make room for the final state text
    if opts.opaqueLegend:
        legend.SetFillStyle(1001) 
    return legend


def getYMinMaxAndName(limits, name, logy, opts):
    ymin = limits.getYMin()
    ymax = limits.getYMax()

    if logy:
        name += "_logy"
        ymax *= 2
    else:
        ymin =  0.0
        ymax *= 1.2

    if opts.logx:
        name += "_logx"
    return ymin, ymax, name


def addPhysicsText(histograms, limits, x=0.45, y=0.84, size=20):
    '''
    Add physics-process text on canvas
    '''
    # Add process text on canvas
    histograms.addText(x, y+0.04, limits.getHardProcessText(), size=size)

    # Add final-state text
    histograms.addText(x, y, limits.getFinalStateText(), size=size)

    #lumi="%s" % limits.getLuminosity(), sqrts="13 TeV", addCmsText=True)#, cmsTextPosition=None, cmsExtraTextPosition=None, cmsText=Non)

    # Add BR assumption (if applicaple)
    if len(limits.getBRassumptionText()) >0:
        histograms.addText(x, y-0.05, limits.getBRassumptionText(), size=size)
    return


def doLimitError(limits, unblindedStatus):
    expRelErrors = []
    expLabels    = {}
    obsRelErrors = []
    obsLabels    = {}

    order = [0, 1, -1, 2, -2]
    expErrors = [limits.expectedErrorGraph(sigma=s) for s in order]
    if expErrors[0] != None:
        exps = [limits.expectedGraph(sigma=s) for s in order]
        expRelErrors = [(limit.divideGraph(expErrors[i], exps[i]), "ExpRelErr%d"%i) for i in xrange(len(exps))]
        expLabels = {
            "ExpRelErr0": "Expected median",
            "ExpRelErr1": "Expected +1#sigma",
            "ExpRelErr2": "Expected -1#sigma",
            "ExpRelErr3": "Expected +2#sigma",
            "ExpRelErr4": "Expected -2#sigma",
            }

    if unblindedStatus:
        obsErr = limits.observedErrorGraph()
        if obsErr != None:
            obs = limits.observedGraph()
            if obs != None:
                obsRelErrors = [(limit.divideGraph(obsErr, obs), "ObsRelErr")]
                obsLabels = {"ObsRelErr": "Observed"}


    if len(expRelErrors) == 0 and len(obsRelErrors) == 0:
        return

    # Create the plot
    plot = plots.PlotBase()
    if len(expRelErrors) > 0:
        plot.histoMgr.extendHistos([histograms.HistoGraph(x[0], x[1], drawStyle="PL", legendStyle="lp") for x in expRelErrors])
        plot.histoMgr.forEachHisto(styles.generator())
        def sty(h):
            r = h.getRootHisto()
            r.SetLineStyle(1)
            r.SetLineWidth(3)
            r.SetMarkerSize(1.4)
        plot.histoMgr.forEachHisto(sty)
        plot.histoMgr.setHistoLegendLabelMany(expLabels)
    if unblindedStatus:
        if len(obsRelErrors) > 0:
            obsRelErrors[0][0].SetMarkerSize(1.4)
            obsRelErrors[0][0].SetMarkerStyle(25)
            plot.histoMgr.insertHisto(0, histograms.HistoGraph(obsRelErrors[0][0], obsRelErrors[0][1], drawStyle="PL", legendStyle="lp"))
            plot.histoMgr.setHistoLegendLabelMany(obsLabels)

    plot.setLegend(histograms.moveLegend(histograms.createLegend(0.48, 0.75, 0.85, 0.92), dx=0.1, dy=-0.1))
 
    if len(limits.mass) == 1:
        plot.createFrame("limitsBrRelativeUncertainty", opts={"xmin": limits.mass[0]-5.0, "xmax": limits.mass[0]+5.0,  "ymin": 0, "ymaxfactor": 1.5})
    else:
        plot.createFrame("limitsBrRelativeUncertainty", opts={"ymin": 0, "ymaxfactor": 1.5})
    plot.frame.GetXaxis().SetTitle(limits.mHplus())
    plot.frame.GetYaxis().SetTitle("Uncertainty/limit")

    # Draw the plot
    plot.draw()
    plot.setLuminosity(limits.getLuminosity())
    plot.addStandardTexts(addLuminosityText=True)
    
    # Add auxiliary text
    histograms.addText(0.20, 0.88, limits.getHardProcessText() , size=20)
    histograms.addText(0.20, 0.84, limits.getFinalStateText()  , size=20)
    histograms.addText(0.20, 0.79, limits.getBRassumptionText(), size=20)
    histograms.addText(0.55, 0.88, "Toy MC relative"           , size=22)
    histograms.addText(0.55, 0.84, "statistical uncertainty"   , size=22)

    plot.save()

    # Save the plots
    SavePlot(plot, "limitError", os.path.join(opts.saveDir, opts.subdir), opts.saveFormats)

    return


if __name__ == "__main__":

    # Default Values
    BATCHMODE    = True
    CFGFILE      = "configuration.json"
    CUTLINEX     = ""
    DIGITS       = 5
    DIRECTORY    = "."
    EXCLUDE      = ""
    EXCLUDEAREA  = False
    GRIDX        = False
    GRIDY        = False
    LIMITSFILE   = "limits.json"
    LOGX         = False
    MAXY         = -1
    MINY         = -1
    PAPER        = False
    PARENTHESES  = False
    SAVEDIR      = None
    SAVEFORMATS  = "pdf,png,C"
    SUBDIR       = ""
    UNBLINDED    = False
    URL          = False
    VERBOSE      = False
    ANALYSISTYPE = "HToTauNu"

    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("-v", "--verbose", dest="verbose", default=VERBOSE, action="store_true",
                      help="Verbose mode for debugging purposes [default: %s]" % (VERBOSE))

    parser.add_option("--dir", dest="directory", default=DIRECTORY,
                      help=" Path to the multicrab task directory with the JSON files [default: %s]" % (DIRECTORY))

    parser.add_option("--limitsFile", dest="limitsFile", default=LIMITSFILE, 
                      help="JSON file containing the limits [default: %s]" % (LIMITSFILE))

    parser.add_option("--cfgFile", dest="cfgFile", default=CFGFILE, 
                      help="JSON file containing the configurations [default: %s]" % (CFGFILE))

    parser.add_option("--exclude", dest="exclude", default=EXCLUDE, 
                      help="List for mass points to exclude (comma separated WITHOUT space) [default: %s]" % (EXCLUDE))
                     
    parser.add_option("--analysisType", dest="analysisType", default=ANALYSISTYPE, 
                      help="Flag to indicate the analysis type (e.g. \"HToTauNu\", \"HToTB\", \"HToHW\") [default: %s]" % (ANALYSISTYPE) )
                     
    parser.add_option("--unblinded", dest="unblinded", default=UNBLINDED, action="store_true",
                      help="Enable unblined mode [default: %s]" % (UNBLINDED) )

    parser.add_option("--paper", dest="paper", default=PAPER, action="store_true",
                      help="Paper mode [default: %s]" % (PAPER) )

    parser.add_option("--url", dest="url", action="store_true", default=URL,
                      help="Don't print the actual save path the plots are saved, but print the URL instead [default: %s]" % URL)
    
    parser.add_option("--opaqueLegend", dest="opaqueLegend", default=EXCLUDEAREA, action="store_true",
                      help="Add excluded area as in MSSM exclusion plots [default: %s]" % (EXCLUDEAREA) )

    parser.add_option("--logx", dest="logx", action="store_true", default=LOGX, 
                      help="Plot x-axis (H+ mass) as logarithmic [default: %s]" % (LOGX) )
    
    parser.add_option("--digits", dest="digits", type="int", default=DIGITS,
                      help="Number of digits (precision) to print/save limit results [default: %s]" % (DIGITS) )

    parser.add_option("--cutLineX", dest="cutLineX", default=CUTLINEX,
                      help="List for values for x-axis lines to be drawn on canvas (comma separated WITHOUT space) [default: %s]" % (CUTLINEX))

    parser.add_option("--gridX", dest="gridX", default=GRIDX, action="store_true",
                      help="Enable the grid for the x-axis [default: %s]" % (GRIDX) )

    parser.add_option("--gridY", dest="gridY", default=GRIDY, action="store_true",
                      help="Enable the grid for the y-axis [default: %s]" % (GRIDY) )

    parser.add_option("--yMin", dest="yMin", default=MINY, type="float",
                      help="Overwrite automaticly calculated minimum value of y-axis [default: %s]" % (MINY) )

    parser.add_option("--yMax", dest="yMax", default=MAXY, type="float",
                      help="Overwrite automaticly calculated maximum value of y-axis [default: %s]" % (MAXY) )

    parser.add_option("--subdir", dest="subdir", type="string", default=SUBDIR,
                      help="Sub-directory describing additional settings used when creating the limits (e.g. no lumi) [default: %s]" % SUBDIR) 

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR,
                      help="Directory where all plots will be saved [default: %s]" % SAVEDIR)

    parser.add_option("-s", "--saveFormats", dest="saveFormats", default = SAVEFORMATS,
                      help="Save formats for all plots [default: %s]" % SAVEFORMATS)

    parser.add_option("-b", "--batchMode", dest="batchMode", action="store_false", default=BATCHMODE,
                      help="Enables batch mode (canvas creation does NOT generate a window) [default: %s]" % BATCHMODE)

    (opts, args) = parser.parse_args()


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

    opts.cutLinesX = []
    if len(opts.cutLineX) > 0:
        opts.cutLinesX = opts.cutLineX.split(",")
    if len(opts.cutLinesX) > 0:
        msg = "Adding x-axis lines \"%s\"" % (", ".join(opts.cutLinesX))
        Print(sh_t +  msg + sh_n, True)
    else:
        pass

    # Sanity check
    if not os.path.isdir(opts.directory):
        msg = "Directory %s does not exists" % (opts.directory)
        raise Exception(sh_e + msg + sh_n)
    else:
        if opts.saveDir == None:
            opts.saveDir = opts.directory

    # Create save formats
    if "," in opts.saveFormats:
        opts.saveFormats = opts.saveFormats.split(",")
    else:
        opts.saveFormats = [opts.saveFormats]
    opts.saveFormats = ["." + s for s in opts.saveFormats]
    if not opts.opaqueLegend:
        opts.saveFormats.append(".eps")

    # Call the main function
    main(opts)
    
    if not opts.batchMode:
        raw_input("=== plotBRLimit.py: Press any key to quit ROOT ...")
