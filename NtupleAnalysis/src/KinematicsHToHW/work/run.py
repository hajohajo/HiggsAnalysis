#!/usr/bin/env python
'''
INSTRUCTIONS:
The required minimum input is a multiCRAB directory with at least one dataset. If successfull
a pseudo multiCRAB with name "analysis_YYMMDD_HHMMSS/" will be created, inside which each
dataset has its own directory with the results (ROOT files with histograms). These can be later
used as input to plotting scripts to get the desired results.


PROOF:
Enable only if your analysis is CPU-limited (e.g. limit calculation) With one analyzer at
a time most probably you are I/O -limited. The limit is how much memory one process is using.


USAGE:
./run.py -m <multicrab_directory> -j <numOfCores> -i <DatasetName>
or
./run.py -m <multicrab_directory> -n 10 -e "Keyword1|Keyword2|Keyword3"


EXAMPLES:
./run.py -m /multicrab_CMSSW752_Default_07Jan2016/
./run.py -m multicrab_CMSSW752_Default_07Jan2016/ -j 16
./run.py -m multicrab_Hplus2tbAnalysis_v8014_20160818T1956 -n 1000 -e QCD
./run.py -m <multicrab-directory> -e TT_extOB
./run.py -m <multicrab_directory> -n 10 -e "QCD_bEnriched_HT300|2016|ST_"


LAST USED:
./run.py -m /uscms_data/d3/aattikis/workspace/multicrab/800k_multicrab_Hplus2hwAnalysis_v8030_20190414T2017
./run.py -m /uscms_data/d3/aattikis/workspace/multicrab/multicrab_Hplus2hwAnalysis_v8030_20181205T1455

ROOT:
The available ROOT options for the Error-Ignore-Level are (const Int_t):
        kUnset    =  -1
        kPrint    =   0
        kInfo     =   1000
        kWarning  =   2000
        kError    =   3000
        kBreak    =   4000

HistoLevel:
For the histogramAmbientLevel each DEEPER level is a SUBSET of the rest. 
For example "kDebug" will include all kDebug histos but also kInformative, kVital, kSystematics, and kNever.  
Setting histogramAmbientLevel=kSystematics will include kSystematics AND kNever.
    1. kNever = 0,
    2. kSystematics,
    3. kVital,
    4. kInformative,
    5. kDebug,
    6. kNumberOfLevels
'''

#================================================================================================
# Imports
#================================================================================================
import sys
from optparse import OptionParser
import time

from HiggsAnalysis.NtupleAnalysis.main import Process, PSet, Analyzer
from HiggsAnalysis.NtupleAnalysis.AnalysisBuilder import AnalysisBuilder


import ROOT
    
#================================================================================================
# Options
#================================================================================================
prefix      = "KinematicsHToHW"
postfix     = ""
dataEras    = ["2016"]
searchModes = ["80to1000"]

ROOT.gErrorIgnoreLevel = 0 


#================================================================================================
# Function Definition
#================================================================================================
def Verbose(msg, printHeader=False):
    if not opts.verbose:
        return

    if printHeader:
        print "=== run.py:"

    if msg !="":
        print "\t", msg
    return


def Print(msg, printHeader=True):
    if printHeader:
        print "=== run.py:"

    if msg !="":
        print "\t", msg
    return


#================================================================================================
# Setup the main function
#================================================================================================
def main():

    # Save start time (epoch seconds)
    tStart = time.time()
    Verbose("Started @ " + str(tStart), True)

    # Require at least two arguments (script-name, path to multicrab)      
    if len(sys.argv) < 2:
        Print("Not enough arguments passed to script execution. Printing docstring & EXIT.")
        print __doc__
        sys.exit(0)
    else:
        pass
        

    # ================================================================================================
    # Setup the process
    # ================================================================================================
    maxEvents = {}
    maxEvents["All"] = opts.nEvts
    process = Process(prefix, postfix, maxEvents)
                
    # ================================================================================================
    # Add the datasets (according to user options)
    # ================================================================================================
    if (opts.includeOnlyTasks):
        Verbose("Adding only dataset %s from multiCRAB directory %s" % (opts.includeOnlyTasks, opts.mcrab))
        process.addDatasetsFromMulticrab(opts.mcrab, includeOnlyTasks=opts.includeOnlyTasks)
    elif (opts.excludeTasks):
        Verbose("Adding all datasets except %s from multiCRAB directory %s" % (opts.excludeTasks, opts.mcrab))
        Verbose("If collision data are present, then vertex reweighting is done according to the chosen data era (era=2015C, 2015D, 2015) etc...")
        process.addDatasetsFromMulticrab(opts.mcrab, excludeTasks=opts.excludeTasks)
    else:
        myBlackList = []

        Print("Adding all datasets from multiCRAB directory %s" % (opts.mcrab))
        Print("If collision data are present, then vertex reweighting is done according to the chosen data era (era=2015C, 2015D, 2015) etc...")
        regex =  "|".join(myBlackList)
        if len(myBlackList)>0:
            process.addDatasetsFromMulticrab(opts.mcrab, excludeTasks=regex)
        else:
            process.addDatasetsFromMulticrab(opts.mcrab)



    # ================================================================================================
    # Overwrite Default Settings  
    # ================================================================================================
    from HiggsAnalysis.NtupleAnalysis.parameters.hplus2tbAnalysis import allSelections
    allSelections.verbose = opts.verbose
    allSelections.histogramAmbientLevel = opts.histoLevel

    allSelections.ElectronSelection.electronPtCut  = 10.0
    allSelections.ElectronSelection.electronEtaCut = 2.1

    allSelections.MuonSelection.muonPtCut  = 25.0
    allSelections.MuonSelection.muonEtaCut = 2.4
    
    allSelections.TauSelection.triggerMatchingApply = False
    allSelections.TauSelection.triggerMatchingCone  =  0.1
    allSelections.TauSelection.tauPtCut             = 30.0
    allSelections.TauSelection.tauEtaCut            =  2.1
    allSelections.TauSelection.tauLdgTrkPtCut       =  0.0
    allSelections.TauSelection.prongs               = -1
    allSelections.TauSelection.rtau                 =  0.0
    allSelections.TauSelection.againstElectronDiscr = "againstElectronTightMVA6"
    allSelections.TauSelection.againstMuonDiscr     = "againstMuonLoose3"
    allSelections.TauSelection.isolationDiscr       = "byLooseCombinedIsolationDeltaBetaCorr3Hits"

    allSelections.JetSelection.jetType                  = "Jets"
    allSelections.JetSelection.jetPtCuts                = [30.0]
    allSelections.JetSelection.jetEtaCuts               = [2.4]
    allSelections.JetSelection.numberOfJetsCutValue     = 3
    allSelections.JetSelection.numberOfJetsCutDirection = ">="
    allSelections.JetSelection.HTCutValue               = 0.0
    allSelections.JetSelection.HTCutDirection           = ">="

    allSelections.BJetSelection.triggerMatchingApply      = False
    allSelections.BJetSelection.triggerMatchingCone       = 0.1
    allSelections.BJetSelection.jetPtCuts                 = [20.0]
    allSelections.BJetSelection.jetEtaCuts                = [2.4]
    allSelections.BJetSelection.bjetDiscr                 = "pfCombinedInclusiveSecondaryVertexV2BJetTags"
    allSelections.BJetSelection.bjetDiscrWorkingPoint     = "Medium"
    allSelections.BJetSelection.numberOfBJetsCutValue     = 1
    allSelections.BJetSelection.numberOfBJetsCutDirection = ">="

    allSelections.METSelection.METCutValue                 = 0 # default: 40.0
    allSelections.METSelection.METCutDirection             = ">="
    allSelections.METSelection.METSignificanceCutValue     = -1000.0
    allSelections.METSelection.METSignificanceCutDirection = ">"
    allSelections.METSelection.METType                     = "MET_Type1"
    allSelections.METSelection.applyPhiCorrections          = False

    # ================================================================================================
    # Add Analysis Variations
    # ================================================================================================
    # selections = allSelections.clone()
    # process.addAnalyzer(prefix, Analyzer(prefix, config=selections, silent=False) ) #trigger passed from selections


    # ================================================================================================
    # Command Line Options
    # ================================================================================================ 
    # from HiggsAnalysis.NtupleAnalysis.parameters.signalAnalysisParameters import applyAnalysisCommandLineOptions
    # applyAnalysisCommandLineOptions(sys.argv, allSelections)

    
    #================================================================================================
    # Build analysis modules
    #================================================================================================
    PrintOptions(opts)
    builder = AnalysisBuilder(prefix,
                              dataEras,
                              searchModes,
                              usePUreweighting       = opts.usePUreweighting,
                              useTopPtReweighting    = opts.useTopPtReweighting,
                              doSystematicVariations = opts.doSystematics,
                              analysisType="HToTB",
                              verbose=opts.verbose)

    # Add variations (e.g. for optimisation)
    # builder.addVariation("METSelection.METCutValue", [100,120,140])
    # builder.addVariation("AngularCutsBackToBack.workingPoint", ["Loose","Medium","Tight"])
    # builder.addVariation("BJetSelection.triggerMatchingApply", [False])
    # builder.addVariation("TopSelection.ChiSqrCutValue", [5, 10, 15, 20])

    # Build the builder
    builder.build(process, allSelections)

    # ================================================================================================
    # Example of adding an analyzer whose configuration depends on dataVersion
    # ================================================================================================
    #def createAnalyzer(dataVersion):
    #a = Analyzer("ExampleAnalysis")
    #if dataVersion.isMC():
    #a.tauPtCut = 10
    #else:
    #a.tauPtCut = 20
    #return a
    #process.addAnalyzer("test2", createAnalyzer)

    
    # ================================================================================================
    # Pick events
    # ================================================================================================
    #process.addOptions(EventSaver = PSet(enabled = True,pickEvents = True))
    # ================================================================================================
    # Run the analysis
    # ================================================================================================
    # Run the analysis with PROOF? You can give proofWorkers=<N> as a parameter
    if opts.jCores:
        Print("Running process with PROOF (proofWorkes=%s)" % ( str(opts.jCores) ) )
        process.run(proof=True, proofWorkers=opts.jCores)
    else:
        Verbose("Running process")
        process.run()

    # Print total time elapsed
    tFinish = time.time()
    dt      = int(tFinish) - int(tStart)
    days    = divmod(dt,86400)      # days
    hours   = divmod(days[1],3600)  # hours
    mins    = divmod(hours[1],60)   # minutes
    secs    = mins[1]               # seconds
    Print("Total elapsed time is %s days, %s hours, %s mins, %s secs" % (days[0], hours[0], mins[0], secs), True)
    return

#================================================================================================
def PrintOptions(opts):
    if not opts.verbose:
        return
    table    = []
    msgAlign = "{:<20} {:<10} {:<10}"
    title    =  msgAlign.format("Option", "Value", "Default")
    hLine    = "="*len(title)
    table.append(hLine)
    table.append(title)
    table.append(hLine)
    #table.append( msgAlign.format("mcrab" , opts.mcrab , "") )
    table.append( msgAlign.format("jCores", opts.jCores, "") )
    table.append( msgAlign.format("includeOnlyTasks", opts.includeOnlyTasks, "") )
    table.append( msgAlign.format("excludeTasks", opts.excludeTasks, "") )
    table.append( msgAlign.format("nEvts", opts.nEvts, NEVTS) )
    table.append( msgAlign.format("verbose", opts.verbose, VERBOSE) )
    table.append( msgAlign.format("histoLevel", opts.histoLevel, HISTOLEVEL) )
    table.append( msgAlign.format("usePUreweighting", opts.usePUreweighting, PUREWEIGHT) )
    table.append( msgAlign.format("useTopPtReweighting", opts.useTopPtReweighting, TOPPTREWEIGHT) )
    table.append( msgAlign.format("doSystematics", opts.doSystematics, DOSYSTEMATICS) ) 
    table.append( hLine )

    # Print("Will run on multicrab directory %s" % (opts.mcrab), True)     
    for i, line in enumerate(table):
        Print(line, i==0)

    return


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

    # Default Values
    VERBOSE       = False
    NEVTS         = -1
    HISTOLEVEL    = "Debug" #"Informative" #"Debug"
    PUREWEIGHT    = False
    TOPPTREWEIGHT = False
    DOSYSTEMATICS = False

    parser = OptionParser(usage="Usage: %prog [options]" , add_help_option=False,conflict_handler="resolve")
    parser.add_option("-m", "--mcrab", dest="mcrab", action="store", 
                      help="Path to the multicrab directory for input")

    parser.add_option("-j", "--jCores", dest="jCores", action="store", type=int, 
                      help="Number of CPU cores (PROOF workes) to use. (default: all available)")

    parser.add_option("-i", "--includeOnlyTasks", dest="includeOnlyTasks", action="store", 
                      help="List of datasets in mcrab to include")

    parser.add_option("-e", "--excludeTasks", dest="excludeTasks", action="store", 
                      help="List of datasets in mcrab to exclude")

    parser.add_option("-n", "--nEvts", dest="nEvts", action="store", type=int, default = NEVTS,
                      help="Number of events to run on (default: %s" % (NEVTS) )

    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default = VERBOSE, 
                      help="Enable verbosity (for debugging) (default: %s)" % (VERBOSE))

    parser.add_option("-h", "--histoLevel", dest="histoLevel", action="store", default = HISTOLEVEL,
                      help="Histogram ambient level (default: %s)" % (HISTOLEVEL))

    parser.add_option("--noPU", dest="usePUreweighting", action="store_false", default = PUREWEIGHT, 
                      help="Do NOT apply Pileup re-weighting (default: %s)" % (PUREWEIGHT) )

    parser.add_option("--topPt", dest="useTopPtReweighting", action="store_true", default = TOPPTREWEIGHT, 
                      help="Do apply top-pt re-weighting (default: %s)" % (TOPPTREWEIGHT) )

    parser.add_option("--doSystematics", dest="doSystematics", action="store_true", default = DOSYSTEMATICS, 
                      help="Do systematics variations  (default: %s)" % (DOSYSTEMATICS) )

    (opts, args) = parser.parse_args()

    if opts.mcrab == None:
        raise Exception("Please provide input multicrab directory with -m")

    allowedLevels = ['Never', 'Systematics', 'Vital', 'Informative', 'Debug']
    if opts.histoLevel not in allowedLevels:
        raise Exception("Invalid ambient histogram level \"%s\"! Valid options are: %s" % (opts.histoLevel, ", ".join(allowedLevels)))

    main()