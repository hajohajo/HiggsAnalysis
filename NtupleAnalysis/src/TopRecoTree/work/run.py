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
./run.py -m <multicrab_directory> -n 10 -e "Keyword1|Keyword2|Keyword3"


EXAMPLE:
./run.py -m /uscms_data/d3/aattikis/workspace/multicrab/multicrab_Hplus2tbAnalysis_v8030_20180508T0644/ -i "TT"
./run.py -m /uscms_data/d3/skonstan/workspace/multicrab/multicrab_TopSystBDT_v8030_20200213T0455 -i "TT" --noPU
./run.py -m /uscms_data/d3/skonstan/workspace/multicrab/multicrab_TopSystBDT_v8030_20200213T0455 --noPU -i "TT|M_250|M_500|M_800"

LAST USED:
./run.py -m /uscms_data/d3/aattikis/workspace/multicrab/multicrab_Hplus2tbAnalysis_v8030_20180508T0644/ -i "TT" --noPU (hadronic selections)
./run.py -m /uscms_data/d3/skonstan/workspace/multicrab/multicrab_TopSystBDT_v8030_20200213T0455 -i "TT" --noPU (semi leptonic selections)

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

from HiggsAnalysis.NtupleAnalysis.main import Process, PSet, Analyzer
from HiggsAnalysis.NtupleAnalysis.AnalysisBuilder import AnalysisBuilder

import copy
import ROOT
    
#================================================================================================
# Options
#================================================================================================
prefix      = "TopRecoTree"
postfix     = ""
dataEras    = ["2016"]
searchModes = ["80to1000"]

ROOT.gErrorIgnoreLevel = 0 

# Define colors
ns = "\033[0;0m"
ts = "\033[1;34m"

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
        Print("Adding only dataset %s from multiCRAB directory %s" % (opts.includeOnlyTasks, opts.mcrab))
        process.addDatasetsFromMulticrab(opts.mcrab, includeOnlyTasks=opts.includeOnlyTasks)
    elif (opts.excludeTasks):
        Print("Adding all datasets except %s from multiCRAB directory %s" % (opts.excludeTasks, opts.mcrab))
        Print("If collision data are present, then vertex reweighting is done according to the chosen data era (era=2015C, 2015D, 2015) etc...")
        process.addDatasetsFromMulticrab(opts.mcrab, excludeTasks=opts.excludeTasks)
    else:
        myBlackList = [
            "ChargedHiggs_HplusTB_HplusToTB_M_1000",
            "hargedHiggs_HplusTB_HplusToTB_M_10000",
            "hargedHiggs_HplusTB_HplusToTB_M_1000_ext1",
            "hargedHiggs_HplusTB_HplusToTB_M_1500",
            "hargedHiggs_HplusTB_HplusToTB_M_1500_ext1",
            "hargedHiggs_HplusTB_HplusToTB_M_180",
            "hargedHiggs_HplusTB_HplusToTB_M_180_ext1",
            "hargedHiggs_HplusTB_HplusToTB_M_200",
            "hargedHiggs_HplusTB_HplusToTB_M_2000",
            "ChargedHiggs_HplusTB_HplusToTB_M_2000_ext1",
            "ChargedHiggs_HplusTB_HplusToTB_M_200_ext1",
            "ChargedHiggs_HplusTB_HplusToTB_M_220",
            "ChargedHiggs_HplusTB_HplusToTB_M_220_ext1",
            "ChargedHiggs_HplusTB_HplusToTB_M_250",
            "ChargedHiggs_HplusTB_HplusToTB_M_2500",
            "ChargedHiggs_HplusTB_HplusToTB_M_2500_ext1",
            "ChargedHiggs_HplusTB_HplusToTB_M_250_ext1",
            "ChargedHiggs_HplusTB_HplusToTB_M_300",
            "ChargedHiggs_HplusTB_HplusToTB_M_3000",
            "ChargedHiggs_HplusTB_HplusToTB_M_3000_ext1",
            "ChargedHiggs_HplusTB_HplusToTB_M_300_ext1",
            "ChargedHiggs_HplusTB_HplusToTB_M_350",
            "ChargedHiggs_HplusTB_HplusToTB_M_350_ext1",
            "ChargedHiggs_HplusTB_HplusToTB_M_400",
            "ChargedHiggs_HplusTB_HplusToTB_M_400_ext1",
            "ChargedHiggs_HplusTB_HplusToTB_M_500",
            "ChargedHiggs_HplusTB_HplusToTB_M_5000",
            "ChargedHiggs_HplusTB_HplusToTB_M_500_ext1",
            "ChargedHiggs_HplusTB_HplusToTB_M_650",
            "ChargedHiggs_HplusTB_HplusToTB_M_7000",
            "ChargedHiggs_HplusTB_HplusToTB_M_800",
            "ChargedHiggs_HplusTB_HplusToTB_M_800_ext1",
            "DYJetsToQQ_HT180",
            "JetHT_Run2016B_03Feb2017_ver2_v2_273150_275376",
            "JetHT_Run2016C_03Feb2017_v1_275656_276283",
            "JetHT_Run2016D_03Feb2017_v1_276315_276811",
            "JetHT_Run2016E_03Feb2017_v1_276831_277420",
            "JetHT_Run2016F_03Feb2017_v1_277932_278800",
            "JetHT_Run2016F_03Feb2017_v1_278801_278808",
            "JetHT_Run2016G_03Feb2017_v1_278820_280385",
            "JetHT_Run2016H_03Feb2017_ver2_v1_281613_284035",
            "JetHT_Run2016H_03Feb2017_ver3_v1_284036_284044",
            "QCD_HT1000to1500",
            "QCD_HT1000to1500_ext1",
            "QCD_HT100to200",
            "QCD_HT1500to2000",
            "QCD_HT1500to2000_ext1",
            "QCD_HT2000toInf",
            "QCD_HT2000toInf_ext1",
            "QCD_HT200to300",
            "QCD_HT200to300_ext1",
            "QCD_HT300to500",
            "QCD_HT300to500_ext1",
            "QCD_HT500to700",
            "QCD_HT500to700_ext1",
            "QCD_HT50to100",
            "QCD_HT700to1000",
            "QCD_HT700to1000_ext1",
            "ST_s_channel_4f_InclusiveDecays",
            "ST_tW_antitop_5f_inclusiveDecays",
            "ST_tW_antitop_5f_inclusiveDecays_ext1",
            "ST_tW_top_5f_inclusiveDecays",
            "ST_tW_top_5f_inclusiveDecays_ext1",
            "ST_t_channel_antitop_4f_inclusiveDecays",
            "ST_t_channel_top_4f_inclusiveDecays",
            #"TT",
            "TTTT",
            "TTWJetsToQQ",
            "TTZToQQ",
            "WJetsToQQ_HT_600ToInf",
            "WWTo4Q",
            "WZ",
            "WZ_ext1",
            "ZJetsToQQ_HT600toInf",
            "ZZTo4Q",
            "QCD_b"]

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
    from HiggsAnalysis.NtupleAnalysis.main import PSet
    RecoTopMVASelection = PSet(
        SelectionsType           = "hadronic", # Default value. This changes according to the mcrab name. Options ["hadronic", "semiLeptonic"]
        LepBJetDRCutValue        = "9999.99",
        LepBJetDRCutDirection    = "<=",
        MiniIsoCutValue          = "0.1",
        MiniIsoCutDirection      = "<=",
        METCutValue              = "-9999.99",
        METCutDirection          = ">=",
        DeltaPtOverPtCutValue    = "0.32",
        DeltaRCutValue           = "0.3"
        )

    # Hadronic selections
    from HiggsAnalysis.NtupleAnalysis.parameters.hplus2tbAnalysis import allSelections as _hadronicSelections

    _hadronicSelections.verbose = opts.verbose
    _hadronicSelections.histogramAmbientLevel = opts.histoLevel
    _hadronicSelections.RecoTopMVASelection = RecoTopMVASelection

    _hadronicSelections.JetSelection.numberOfJetsCutValue  = 7 # [default: 7]
    _hadronicSelections.BJetSelection.numberOfBJetsCutValue = 3 # [default: 3]
    
    hadronicSelections = copy.deepcopy(_hadronicSelections)

    # Semi leptonic selections
    from HiggsAnalysis.NtupleAnalysis.parameters.jetTriggers import allSelections as _semiLeptonicSelections
    
    _semiLeptonicSelections.verbose = opts.verbose
    _semiLeptonicSelections.histogramAmbientLevel = opts.histoLevel
    _semiLeptonicSelections.RecoTopMVASelection = RecoTopMVASelection
    # Muon
    _semiLeptonicSelections.MuonSelection.muonPtCut = 26
    # Electron
    _semiLeptonicSelections.ElectronSelection.electronPtCut = 29
    # Jets
    _semiLeptonicSelections.JetSelection.numberOfJetsCutValue = 5
    _semiLeptonicSelections.JetSelection.jetPtCuts = [40] #[40.0, 40.0, 40.0, 30.0]
    
    # Trigger
    _semiLeptonicSelections.Trigger.triggerOR = ["HLT_IsoMu24", "HLT_Ele27_eta2p1_WPTight_Gsf"]
    
    # Bjets
    _semiLeptonicSelections.BJetSelection.jetPtCuts = [40.0, 40.0]
    _semiLeptonicSelections.BJetSelection.numberOfBJetsCutValue = 2
    
    semiLeptonicSelections = copy.deepcopy(_semiLeptonicSelections)

    # Define selections type based on the mcrab name
    if "Hplus2tbAnalysis" in opts.mcrab:
        RecoTopMVASelection.SelectionsType = "hadronic"
    elif "TopSyst" in opts.mcrab:
        RecoTopMVASelection.SelectionsType = "semiLeptonic"
        
    # Get allSelections 
    if (RecoTopMVASelection.SelectionsType == "hadronic"):
        allSelections = copy.deepcopy(hadronicSelections)
    else:
        allSelections = copy.deepcopy(semiLeptonicSelections)

    Print("Selections Type: %s %s %s" %(ts, RecoTopMVASelection.SelectionsType, ns), True)

    #Print ("Hadronic Selections", True)
    #print hadronicSelections.Trigger
    #Print ("SemiLeptonic Selections", True)
    #print semiLeptonicSelections.Trigger
    #Print ("All Selections", True)
    #print allSelections.Trigger

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
                              )

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
        Print("Running process")
        process.run()


#================================================================================================
def PrintOptions(opts):
    '''
    '''
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
