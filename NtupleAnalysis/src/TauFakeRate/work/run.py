#!/usr/bin/env python
'''
INSTRUCTIONS:
The required minimum input is a multiCRAB directory with at least one dataset. If successfull
a pseudo multiCRAB with name "analysis_YYMMDD_HHMMSS/" will be created, inside which each
dataset has its own directory with the results (ROOT files with histograms). These can be later
used as input to plotting scripts to get the desired results.


USAGE:
./run.py -m <multicrab_directory> -n 10 -e "Keyword1|Keyword2|Keyword3"


EXAMPLE:
./run.py -m <multicrab_directory> -n 10 -e "QCD_bEnriched_HT300|2016|ST_"


LAST USED:
./run.py -m  /uscms_data/d3/aattikis/workspace/multicrab/multicrab_Hplus2hwAnalysis_v8030_20190628T1421

ROOT:
The available ROOT options for the Error-Ignore-Level are (const Int_t):
        kUnset    =  -1
        kPrint    =   0
        kInfo     =   1000
        kWarning  =   2000
        kError    =   3000
        kBreak    =   4000


HISTOLEVEL:
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
prefix      = "TauFakeRate"
postfix     = ""
dataEras    = ["2016"]
searchModes = ["350to3000"] # ["80to1000"]

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
        myBlackList = ["QCD_Pt_15to20_MuEnrichedPt5",
                       "QCD_Pt_20to30_MuEnrichedPt5",
                       "QCD_Pt_30to50_MuEnrichedPt5",
                       "QCD_Pt_50to80_MuEnrichedPt5",
                       "QCD_Pt_80to120_MuEnrichedPt5",
                       "QCD_Pt_80to120_MuEnrichedPt5_ext1",
                       "QCD_Pt_120to170_MuEnrichedPt5",
                       "QCD_Pt_170to300_MuEnrichedPt5",
                       "QCD_Pt_170to300_MuEnrichedPt5_ext1",
                       "QCD_Pt_300to470_MuEnrichedPt5",
                       "QCD_Pt_300to470_MuEnrichedPt5_ext1",
                       "QCD_Pt_300to470_MuEnrichedPt5_ext2",
                       "QCD_Pt_470to600_MuEnrichedPt5",
                       "QCD_Pt_470to600_MuEnrichedPt5_ext1",
                       "QCD_Pt_470to600_MuEnrichedPt5_ext2",
                       "QCD_Pt_600to800_MuEnrichedPt5",
                       "QCD_Pt_600to800_MuEnrichedPt5_ext1",
                       "QCD_Pt_800to1000_MuEnrichedPt5",
                       "QCD_Pt_800to1000_MuEnrichedPt5_ext1",
                       "QCD_Pt_800to1000_MuEnrichedPt5_ext2",
                       "QCD_Pt_1000toInf_MuEnrichedPt5",
                       "QCD_bEnriched_HT200to300",
                       "QCD_bEnriched_HT300to500",
                       "QCD_bEnriched_HT500to700",
                       "QCD_bEnriched_HT700to1000",
                       "QCD_bEnriched_HT1000to1500",
                       "QCD_bEnriched_HT1500to2000",
                       "QCD_bEnriched_HT2000toInf",
                       "CRAB_private_ChargedHiggs_HplusTB_HplusToHW_M300_mH200_2ta_NLO",
                       "CRAB_private_ChargedHiggs_HplusTB_HplusToHW_M700_mH200_2ta_NLO",
                       #"ttHJetToGG_M125",
                       #"ttHJetToNonbb_M125_ext1",
                       #"ttHJetToTT_M125_ext4",
                       #"ttHJetTobb_M125_ext3",
                       #"TTTT",
                       #"TTWJetsToLNu_ext1",
                       #"TTWJetsToLNu_ext2",
                       #"TTZToLLNuNu_M_10_ext3",
                       #"TTZToQQ",
                       #"WJetsToLNu_HT_100To200",
                       #"WJetsToLNu_HT_100To200_ext1",
                       #"WJetsToLNu_HT_100To200_ext2",
                       #"WJetsToLNu_HT_1200To2500",
                       #"WJetsToLNu_HT_1200To2500_ext1",
                       #"WJetsToLNu_HT_200To400",
                       #"WJetsToLNu_HT_200To400_ext1",
                       #"WJetsToLNu_HT_200To400_ext2",
                       #"WJetsToLNu_HT_2500ToInf",
                       #"WJetsToLNu_HT_2500ToInf_ext1",
                       #"WJetsToLNu_HT_400To600",
                       #"WJetsToLNu_HT_400To600_ext1",
                       #"WJetsToLNu_HT_600To800",
                       #"WJetsToLNu_HT_600To800_ext1",
                       #"WJetsToLNu_HT_70To100",
                       #"WJetsToLNu_HT_800To1200",
                       #"WJetsToLNu_HT_800To1200_ext1",
                       "WWTo4Q",
                       "SingleElectron_Run2016G_03Feb2017_v1_278820_280385",
                       "SingleElectron_Run2016F_03Feb2017_v1_277932_278800",
                       "SingleElectron_Run2016E_03Feb2017_v1_276831_277420",
                       "SingleElectron_Run2016D_03Feb2017_v1_276315_276811",
                       "SingleElectron_Run2016B_03Feb2017_ver2_v2_273150_275376",
                       "SingleElectron_Run2016H_03Feb2017_ver3_v1_284036_284044",
                       "SingleElectron_Run2016H_03Feb2017_ver2_v1_281613_284035",
                       "SingleElectron_Run2016F_03Feb2017_v1_278801_278808",
                       "SingleElectron_Run2016C_03Feb2017_v1_275656_276283",
                       ]

        #if opts.doSystematics:
        #    myBlackList.append("QCD")

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
    from HiggsAnalysis.NtupleAnalysis.parameters.hplus2hwAnalysisWithTop import allSelections
    allSelections.verbose = opts.verbose
    allSelections.histogramAmbientLevel = opts.histoLevel
    allSelections.MuonSelection.applyTriggerMatching = False # cannot use it for 2 muons and single-muon trigger
    #allSelections.MuonSelection.muonPtCut = 40.0 # cannot use it for 2 muons and single-muon trigger
    #allSelections.JetSelection.numberOfJetsCutValue = 1 # [default; 0]

    # allSelections.BjetSelection.triggerMatchingApply = True
    # allSelections.BJetSelection.numberOfBJetsCutValue = 0
    # allSelections.BJetSelection.numberOfBJetsCutDirection = "=="
    
    #allSelections.TopSelectionMVA.MVACutValue           = 0.90 # [default: 0.4]
    #allSelections.FakeBTopSelectionMVA.MVACutValue      = 0.00 # [default: -1.0]
    #allSelections.FakeBMeasurement.LdgTopMVACutValue    = allSelections.TopSelectionMVA.MVACutValue
    #allSelections.FakeBMeasurement.SubldgTopMVACutValue = allSelections.TopSelectionMVA.MVACutValue

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
                              analysisType="HToHW_withTop",
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
    PUREWEIGHT    = True
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
