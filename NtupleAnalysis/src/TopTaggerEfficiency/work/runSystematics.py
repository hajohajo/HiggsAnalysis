#!/usr/bin/env python
'''
INSTRUCTIONS:
Designed for use with CONDOR

The required minimum input is a multiCRAB directory with at least one dataset. If successfull
a pseudo multiCRAB with name "analysis_YYMMDD_HHMMSS/" will be created, inside which each
dataset has its own directory with the results (ROOT files with histograms). These can be later
used as input to plotting scripts to get the desired results.


PROOF:
Enable only if your analysis is CPU-limited (e.g. limit calculation) With one analyzer at
a time most probably you are I/O -limited. The limit is how much memory one process is using.


USAGE:
./runSystematics.py -m <multicrab_directory> -j <numOfCores> -i <DatasetName>
or
./runSystematics.py -m <multicrab_directory> -n 10 -e "Keyword1|Keyword2|Keyword3"

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
prefix      = "Hplus2tbAnalysis"
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
        print "=== runSystematics.py:"

    if msg !="":
        print "\t", msg
    return


def Print(msg, printHeader=True):
    if printHeader:
        print "=== runSystematics.py:"

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
    completeList = GetDatasetCompleteList()
    whiteList    = GetDatasetWhitelist(opts)
    blackList    = GetDatasetBlackList(completeList, whiteList)
    maxEvents = {}
    for d in whiteList:
        maxEvents[d] = -1
        #maxEvents[d] = 100 #for testing
        #if  d == "ChargedHiggs_HplusTB_HplusToTB_M_650":
        #    maxEvents[d] = 4000000
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
        myBlackList = ["QCD_b",
                       "ChargedHiggs_HplusTB_HplusToTB_M_180_ext1", 
                       "ChargedHiggs_HplusTB_HplusToTB_M_200_ext1", 
                       "ChargedHiggs_HplusTB_HplusToTB_M_220_ext1", 
                       "ChargedHiggs_HplusTB_HplusToTB_M_250_ext1", 
                       "ChargedHiggs_HplusTB_HplusToTB_M_300_ext1", 
                       "ChargedHiggs_HplusTB_HplusToTB_M_350_ext1", 
                       "ChargedHiggs_HplusTB_HplusToTB_M_400_ext1", 
                       "ChargedHiggs_HplusTB_HplusToTB_M_500_ext1", 
                       "ChargedHiggs_HplusTB_HplusToTB_M_800_ext1", 
                       "ChargedHiggs_HplusTB_HplusToTB_M_1000_ext1", 
                       "ChargedHiggs_HplusTB_HplusToTB_M_2000_ext1"
                       "ChargedHiggs_HplusTB_HplusToTB_M_2500_ext1", 
                       "ChargedHiggs_HplusTB_HplusToTB_M_3000_ext1", 
                       "ChargedHiggs_HplusTB_HplusToTB_M_5000",
                       "ChargedHiggs_HplusTB_HplusToTB_M_7000",
                       "ChargedHiggs_HplusTB_HplusToTB_M_10000",
                       ]
        if opts.doSystematics:
            whiteList = GetDatasetWhitelist(opts)

        # Extend the blacklist with datasets not in the group
        myBlackList.extend(blackList)

        Verbose("Adding the following datasets from multiCRAB directory %s:\n\t%s" % (opts.mcrab, "\n\t".join(whiteList) ))
        Verbose("If collision data are present, then vertex reweighting is done according to the chosen data era (era=2015C, 2015D, 2015) etc...")
        if len(whiteList)>0:
            process.addDatasetsFromMulticrab(opts.mcrab, includeOnlyTasks="|".join(whiteList))
        elif len(myBlackList)>0:
        #if len(myBlackList)>0:
            process.addDatasetsFromMulticrab(opts.mcrab, excludeTasks="|".join(myBlackList))
        else:
            process.addDatasetsFromMulticrab(opts.mcrab)

    # ================================================================================================
    # Overwrite Default Settings  
    # ================================================================================================
    from HiggsAnalysis.NtupleAnalysis.parameters.hplus2tbAnalysis import allSelections
    allSelections.verbose = opts.verbose
    allSelections.histogramAmbientLevel = opts.histoLevel
    # allSelections.BjetSelection.triggerMatchingApply = True
    # allSelections.TopSelection.ChiSqrCutValue = 100.0
    # allSelections.BJetSelection.numberOfBJetsCutValue = 0
    # allSelections.BJetSelection.numberOfBJetsCutDirection = "=="

    print allSelections.TopSelectionMVA
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
def GetDatasetBlackList(completeList, whiteList):
    myBlacklist = []
    for d in completeList:
        if d in whiteList:
            continue
        myBlacklist.append(d)
    return myBlacklist
            
def GetDatasetCompleteList():
    myCompleteList = []
    myCompleteList.append("JetHT_Run2016B_03Feb2017_ver2_v2_273150_275376")
    myCompleteList.append("JetHT_Run2016C_03Feb2017_v1_275656_276283")
    myCompleteList.append("JetHT_Run2016D_03Feb2017_v1_276315_276811")
    myCompleteList.append("JetHT_Run2016E_03Feb2017_v1_276831_277420")
    myCompleteList.append("JetHT_Run2016F_03Feb2017_v1_277932_278800")
    myCompleteList.append("JetHT_Run2016F_03Feb2017_v1_278801_278808")
    myCompleteList.append("JetHT_Run2016G_03Feb2017_v1_278820_280385")
    myCompleteList.append("JetHT_Run2016H_03Feb2017_ver2_v1_281613_284035")
    myCompleteList.append("JetHT_Run2016H_03Feb2017_ver3_v1_284036_284044")
    #
    myCompleteList.append("ChargedHiggs_HplusTB_HplusToTB_M_500")
    myCompleteList.append("ChargedHiggs_HplusTB_HplusToTB_M_180")
    myCompleteList.append("ChargedHiggs_HplusTB_HplusToTB_M_200")
    myCompleteList.append("ChargedHiggs_HplusTB_HplusToTB_M_220")
    myCompleteList.append("ChargedHiggs_HplusTB_HplusToTB_M_250")
    myCompleteList.append("ChargedHiggs_HplusTB_HplusToTB_M_300")
    myCompleteList.append("ChargedHiggs_HplusTB_HplusToTB_M_350")
    myCompleteList.append("ChargedHiggs_HplusTB_HplusToTB_M_400")
    myCompleteList.append("ChargedHiggs_HplusTB_HplusToTB_M_500")
    myCompleteList.append("ChargedHiggs_HplusTB_HplusToTB_M_650")
    myCompleteList.append("ChargedHiggs_HplusTB_HplusToTB_M_800")
    myCompleteList.append("ChargedHiggs_HplusTB_HplusToTB_M_1000")
    myCompleteList.append("ChargedHiggs_HplusTB_HplusToTB_M_1500")
    myCompleteList.append("ChargedHiggs_HplusTB_HplusToTB_M_2000")
    myCompleteList.append("ChargedHiggs_HplusTB_HplusToTB_M_2500")
    myCompleteList.append("ChargedHiggs_HplusTB_HplusToTB_M_3000")
    myCompleteList.append("ChargedHiggs_HplusTB_HplusToTB_M_5000")
    #
    myCompleteList.append("ZZTo4Q")
    myCompleteList.append("ZJetsToQQ_HT600toInf")
    myCompleteList.append("WZ_ext1")
    myCompleteList.append("WZ")
    myCompleteList.append("WWTo4Q")
    myCompleteList.append("WJetsToQQ_HT_600ToInf")
    myCompleteList.append("TTZToQQ")
    myCompleteList.append("TTWJetsToQQ")
    myCompleteList.append("TTTT")
    myCompleteList.append("TT")
    #
    myCompleteList.append("DYJetsToQQ_HT180")
    myCompleteList.append("QCD_HT50to100")
    myCompleteList.append("QCD_HT100to200")
    myCompleteList.append("QCD_HT200to300")
    myCompleteList.append("QCD_HT200to300_ext1")
    myCompleteList.append("QCD_HT300to500")
    myCompleteList.append("QCD_HT300to500_ext1")
    myCompleteList.append("QCD_HT500to700")
    myCompleteList.append("QCD_HT500to700_ext1")
    myCompleteList.append("QCD_HT700to1000")
    myCompleteList.append("QCD_HT700to1000_ext1")
    myCompleteList.append("QCD_HT1000to1500")
    myCompleteList.append("QCD_HT1000to1500_ext1")
    myCompleteList.append("QCD_HT1500to2000")
    myCompleteList.append("QCD_HT1500to2000_ext1")
    myCompleteList.append("QCD_HT2000toInf")
    myCompleteList.append("QCD_HT2000toInf_ext1")
    myCompleteList.append("ST_s_channel_4f_InclusiveDecays")
    myCompleteList.append("ST_t_channel_antitop_4f_inclusiveDecays")
    myCompleteList.append("ST_t_channel_top_4f_inclusiveDecays")
    myCompleteList.append("ST_tW_antitop_5f_inclusiveDecays")
    myCompleteList.append("ST_tW_antitop_5f_inclusiveDecays_ext1")
    myCompleteList.append("ST_tW_top_5f_inclusiveDecays")
    myCompleteList.append("ST_tW_top_5f_inclusiveDecays_ext1")
    #
    myCompleteList.append("QCD_bEnriched_HT1000to1500")
    myCompleteList.append("QCD_bEnriched_HT1500to2000")
    myCompleteList.append("QCD_bEnriched_HT2000toInf")
    myCompleteList.append("QCD_bEnriched_HT200to300")
    myCompleteList.append("QCD_bEnriched_HT300to500")
    myCompleteList.append("QCD_bEnriched_HT500to700")
    myCompleteList.append("QCD_bEnriched_HT700to1000")
    return myCompleteList

def GetDatasetWhitelist(opts):
    myWhitelist = []
    if opts.group == "A":
        myWhitelist.append("JetHT_Run2016B_03Feb2017_ver2_v2_273150_275376")
        myWhitelist.append("JetHT_Run2016C_03Feb2017_v1_275656_276283")
        myWhitelist.append("JetHT_Run2016D_03Feb2017_v1_276315_276811")
        myWhitelist.append("JetHT_Run2016E_03Feb2017_v1_276831_277420")
        myWhitelist.append("JetHT_Run2016F_03Feb2017_v1_277932_278800")
    elif opts.group == "B":
        myWhitelist.append("JetHT_Run2016F_03Feb2017_v1_278801_278808")
        myWhitelist.append("JetHT_Run2016G_03Feb2017_v1_278820_280385")
        myWhitelist.append("JetHT_Run2016H_03Feb2017_ver2_v1_281613_284035")
        myWhitelist.append("JetHT_Run2016H_03Feb2017_ver3_v1_284036_284044")
    elif opts.group == "C":
        myWhitelist.append("ChargedHiggs_HplusTB_HplusToTB_M_180")
        myWhitelist.append("ChargedHiggs_HplusTB_HplusToTB_M_800")
    elif opts.group == "D":
        myWhitelist.append("ChargedHiggs_HplusTB_HplusToTB_M_350")
        myWhitelist.append("ChargedHiggs_HplusTB_HplusToTB_M_400")
    elif opts.group == "E":
        myWhitelist.append("ChargedHiggs_HplusTB_HplusToTB_M_650")
    elif opts.group == "F":
        myWhitelist.append("ChargedHiggs_HplusTB_HplusToTB_M_220")
        myWhitelist.append("ChargedHiggs_HplusTB_HplusToTB_M_1500")
    elif opts.group == "G":
        myWhitelist.append("ChargedHiggs_HplusTB_HplusToTB_M_250")
        myWhitelist.append("ChargedHiggs_HplusTB_HplusToTB_M_2500")
    elif opts.group == "H":
        myWhitelist.append("ZZTo4Q")
        myWhitelist.append("ZJetsToQQ_HT600toInf")
        myWhitelist.append("WZ_ext1")
        myWhitelist.append("WZ")
        myWhitelist.append("WWTo4Q")
        myWhitelist.append("WJetsToQQ_HT_600ToInf")
        #myWhitelist.append("TTZToQQ")
        #myWhitelist.append("TTWJetsToQQ")
    elif opts.group == "I":    
        myWhitelist.append("DYJetsToQQ_HT180")
        myWhitelist.append("QCD_HT50to100")
        myWhitelist.append("QCD_HT100to200")
        myWhitelist.append("QCD_HT200to300")
        myWhitelist.append("QCD_HT200to300_ext1")
        myWhitelist.append("QCD_HT300to500")
        myWhitelist.append("QCD_HT300to500_ext1")
        myWhitelist.append("QCD_HT500to700")
        myWhitelist.append("QCD_HT500to700_ext1")
        myWhitelist.append("QCD_HT700to1000")
        myWhitelist.append("QCD_HT700to1000_ext1")
        myWhitelist.append("QCD_HT1000to1500")
        myWhitelist.append("QCD_HT1000to1500_ext1")
        myWhitelist.append("QCD_HT1500to2000")
        myWhitelist.append("QCD_HT1500to2000_ext1")
        myWhitelist.append("QCD_HT2000toInf")
        myWhitelist.append("QCD_HT2000toInf_ext1")
        myWhitelist.append("ST_s_channel_4f_InclusiveDecays")
        myWhitelist.append("ST_t_channel_antitop_4f_inclusiveDecays")
        myWhitelist.append("ST_t_channel_top_4f_inclusiveDecays")
        myWhitelist.append("ST_tW_antitop_5f_inclusiveDecays")
        myWhitelist.append("ST_tW_antitop_5f_inclusiveDecays_ext1")
        myWhitelist.append("ST_tW_top_5f_inclusiveDecays")
        myWhitelist.append("ST_tW_top_5f_inclusiveDecays_ext1")
    elif opts.group == "J":
        myWhitelist.append("ChargedHiggs_HplusTB_HplusToTB_M_200")
        myWhitelist.append("ChargedHiggs_HplusTB_HplusToTB_M_2000")
    elif opts.group == "K":
        myWhitelist.append("ChargedHiggs_HplusTB_HplusToTB_M_300")
        myWhitelist.append("ChargedHiggs_HplusTB_HplusToTB_M_3000")
    elif opts.group == "L":
        myWhitelist.append("TT")
        myWhitelist.append("TTTT")
    elif opts.group == "M":
        myWhitelist.append("ChargedHiggs_HplusTB_HplusToTB_M_500")
        myWhitelist.append("ChargedHiggs_HplusTB_HplusToTB_M_5000")
    elif opts.group == "N":
        myWhitelist.append("ChargedHiggs_HplusTB_HplusToTB_M_1000")
        myWhitelist.append("ChargedHiggs_HplusTB_HplusToTB_M_10000")
    else:
        msg = "Unknown systematics submission dataset group \"%s\"%" % (opts.group)
        raise Exception(msg)
    return myWhitelist

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
    TOPPTREWEIGHT = True
    DOSYSTEMATICS = False
    GROUP         = "A"

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

    parser.add_option("--group", dest="group", default = GROUP, 
                      help="The group of datasets to run on. Capital letter from \"A\" to \"I\"  (default: %s)" % (GROUP) )

    (opts, args) = parser.parse_args()

    if opts.mcrab == None:
        raise Exception("Please provide input multicrab directory with -m")

    allowedLevels = ['Never', 'Systematics', 'Vital', 'Informative', 'Debug']
    if opts.histoLevel not in allowedLevels:
        raise Exception("Invalid ambient histogram level \"%s\"! Valid options are: %s" % (opts.histoLevel, ", ".join(allowedLevels)))

    main()
