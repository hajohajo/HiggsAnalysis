#! /usr/bin/env python
'''
DESCRIPTION:
This is a script for generating datacards to be used by combine (or LimitOMatic) to 
produce exclusion limits.


INSTRUCTIONS:
To run this datacard generator you first need to run your signal analysis (and 
the data-driven backdround estimation) to produce a pseudo-multicrab directory with the
results. The go to HiggsAnalysis/NtupleAnalysis/src/LimitCalc/work and create a new directory
(e.g. "exampleDirectory"). If not already the case, Rename the pseudo-multicrab directory so
that it starts with "SignalAnalysis_" and move it under the "exampleDirectory". 
Set the input paremters (datasets, signal rates, backgound rates, systematics, etc..) for
datacard generation, defined in your template datacard python file 
(e.g. LimitCalc/work/dcardHplustb2017Datacard). Run the datacard generator 
by providing a minimum of two arguments:
1) The datacard python file with all the datacard definitions & settings
2) The pseudo-multicrab directory with all the results


USAGE:
./dcardGenerator.py -x <datacard> --dir <pseudo-multicrab-dir> [opts]


EXAMPLES:
./dcardGenerator_v2.py -x dcardDefault_h2tb_2016.py --dir limits2016/ --analysisType HToHW
./dcardGenerator_v3.py --dir limits2019/ --analysisType HToHW --datacard datacard_HToHW.py


LAST USED:
./dcardGenerator_v3.py --analysisType HToHW --datacard datacard_HToHW.py --dir limits2019/ --dir limits2016/


'''

#================================================================================================
# Import modules
#================================================================================================
import os
import sys
import imp
from optparse import OptionParser
import gc
import cPickle
import time
import ROOT
ROOT.gROOT.SetBatch(True) # no flashing canvases
ROOT.PyConfig.IgnoreCommandLineOptions = True
import tarfile
import cProfile

import HiggsAnalysis.LimitCalc.MulticrabPathFinder as PathFinder
import HiggsAnalysis.NtupleAnalysis.tools.analysisModuleSelector as analysisModuleSelector
import HiggsAnalysis.LimitCalc.DataCardGenerator as DataCard
import HiggsAnalysis.NtupleAnalysis.tools.dataset as dataset
import HiggsAnalysis.NtupleAnalysis.tools.aux as aux
import HiggsAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles
import HiggsAnalysis.NtupleAnalysis.tools.multicrab as multicrab


#===============================================================================================
# Shell Types
#================================================================================================
sh_e = ShellStyles.ErrorStyle()
sh_s = ShellStyles.SuccessStyle()
sh_h = ShellStyles.HighlightStyle()
sh_a = ShellStyles.HighlightAltStyle()
sh_t = ShellStyles.NoteStyle()
sh_n = ShellStyles.NormalStyle()
sh_w = ShellStyles.WarningStyle()


#================================================================================================
# Function definition
#================================================================================================
def PrintOptions(moduleSelector):
    nEras       = len(moduleSelector.getSelectedEras())
    nModes      = len(moduleSelector.getSelectedSearchModes())
    nOptModes   = len(moduleSelector.getSelectedOptimizationModes())
    nDatacards  = nEras*nModes*nOptModes

    msg  = "Producing %d set(s) of datacards" % (nDatacards)
    msg += " [%d era(s) x %d search mode(s) x %d optimization mode(s)]" % (nEras, nModes, nOptModes)
    Verbose(sh_a + msg + sh_n)
    return


def PrintDsetInfo(dsetCreatorList, verbose=False):
    if not verbose:
        return

    # For-loop: All dataset creators
    for i, dc in enumerate(dsetCreatorList, 1):        
        Verbose("Dataset creator with label %s (directory %s) has these datasets:" % (sh_t + dc.getLabel() + sh_n, sh_a + dc.getBaseDirectory() + sh_n), True)
        for dset in dc.getDatasetNames():
            Verbose(dset, False)
    return


def SetModuleSelectorSources(moduleSelector, primaryDsetCreators, secondaryDsetCreators, opts):

    # Definitions
    mySources = []

    # Sanity check
    if len(primaryDsetCreators) > 1:
        msg = "Currently only support a maximum of 1 primary dataset creators. Found %d and will ignore all extra ones." % (len(primaryDsetCreators))
        Print(sh_e + msg + sh_n, True)
        raw_input("Press any key to continue ...")

    # For loop: All Primary dataset creators in list
    for i, p in enumerate(primaryDsetCreators, 1):
        if i > 1:            
            break
        moduleSelector.setPrimarySource(p.getLabel(), p)
        mySources.append(p.getLabel())

    # For loop: All Secondary dataset creators in list
    for s in secondaryDsetCreators:
        moduleSelector.addOtherSource(s.getLabel(), s)
        mySources.append(s.getLabel())

    Verbose("Added %i sources to the module selector with labels: %s" % (len(mySources), ", ".join(mySources)), True)
    moduleSelector.doSelect(opts, printSelections=opts.verbose)
    moduleSelector.closeFiles()
    return


def CheckOptions(config):
    '''
    Check various options defined in the config (=opts.datacard) and warns user if some flags
    are disabled
    '''
    msgs = []
    if not config.OptionIncludeSystematics:
        msg = sh_t + "Skipping of systematic uncertainties and will only use statistical uncertainties"  + sh_n + " (flag \"OptionIncludeSystematics\" in the datacard file)"
        msgs.append(msg)

    if not config.OptionDoControlPlots:
        msg = sh_t + "Skipping of data-driven control plot generation" + sh_n + " (flag \"OptionDoControlPlots\" in the datacard file)"
        msgs.append(msg)

    if not config.BlindAnalysis:
        msg = sh_e + "Unblinding analysis results!" + sh_n + " (flag \"BlindAnalysis\" in the datacard file)"
        msgs.append(msg)

    if len(msgs) < 1:
        return

    # Print all warnings!
    for i, m in enumerate(msgs, 1):
        Print(m + sh_n, i==1)
    return

def getSignalDsetCreator(multicrabPaths, label, mcrabInfoOutput):

    multicrabPath = multicrabPaths.getSignalPath()
    if multicrabPath == "":
        msg = "Signal analysis multicrab directory not found!"
        raise Exception(sh_e + msg + sh_n)
    
    # Get the dataset creator
    signalDsetCreator = getDsetCreator(label, multicrabPath, mcrabInfoOutput)
    msg = "The signal with label %s will be estimated from multicrab directory %s" % (sh_t + label + sh_n, sh_a + multicrabPath + sh_n)
    Print(msg, True)
    return signalDsetCreator


def getBkgDsetCreator(multicrabPaths, bkgLabel, fakesFromData, mcrabInfoOutput):

    bkgDsetCreator = None
    multicrabPath  = None
    myPaths = {}
    myPaths["EWK Genuine-b"]   = multicrabPaths.getGenuineBPath()
    myPaths["Fake-b"]          = multicrabPaths.getFakeBPath()
    myPaths["EWK Genuine-tau"] = multicrabPaths.getEWKPath()
    myPaths["Fake-tau"]        = multicrabPaths.getQCDInvertedPath()
    myPaths["MC"]              = multicrabPaths.getSignalPath()


    print myPaths
    print bkgLabel

    if bkgLabel in myPaths.keys():
        multicrabPath = myPaths[bkgLabel]
    else:
        msg = "Unsupported background type %s. Currently supported backbround types are: \"%s" % (bkgLabel, "\", \"".join(myPaths.keys()) + "\"")
        raise Exception(sh_e + msg + sh_n)

    # Sanity check
    if multicrabPath == "":
        raise Exception(sh_e + "Could not find the path for the background dataset with label %s" + sh_n % (sh_n + bkgLabel + sh_n) )

    # Get the dataset creator
    msg = "The background with label %s will be estimated from multicrab directory %s" % (sh_t + bkgLabel + sh_n, sh_a + multicrabPath + sh_n)
    Verbose(msg, True)
    bkgDsetCreator = getDsetCreator(bkgLabel, multicrabPath, mcrabInfoOutput, True)

    mcrabInfoOutput.append("- " + msg)
    Print(msg)
    return bkgDsetCreator

def GetDatasetLabels(config, fakesFromData, opts):
    #iro Print(sh_w + "GetDatasetLabels() is too specific. Must be made more generic" + sh_n, True)
    signalLabels = ["Signal Analysis"]
    bkgLabels    = []

    if fakesFromData:
        if opts.analysisType == "HToTB":
            #bkg1Label = "EWK Genuine-b"
            #bkg2Label = "Fake-b"
            bkgLabels.append("EWK Genuine-b")
            bkgLabels.append("Fake-b")
        elif opts.analysisType in ["HToTauNu", "HToHW"]:
            #bkg1Label = "EWK Genuine-tau"
            #bkg2Label = "Fake-tau"
            bkgLabels.append("EWK Genuine-tau")
            bkgLabels.append("Fake-tau")
        else:
            pass            
    else:
        bkgLabels.append("MC")
        #bkg1Label = "EWK MC"
        #bkg2Label = "QCD MC"
    return signalLabels, bkgLabels
    #return signalLabel, bkg1Label, bkg2Label


def memoryDump():
    dump = open("memory_pickle.txt", 'w')
    for obj in gc.get_objects():
        i = id(obj)
        size = sys.getsizeof(obj, 0)
        #    referrers = [id(o) for o in gc.get_referrers(obj) if hasattr(o, '__class__')]
        referents = [id(o) for o in gc.get_referents(obj) if hasattr(o, '__class__')]
        if hasattr(obj, '__class__'):
            cls = str(obj.__class__)
            cPickle.dump({'id': i, 'class': cls, 'size': size, 'referents': referents}, dump)
    return


def getDsetCreator(label, mcrabPath, mcrabInfoOutput, enabledStatus=True):
    if enabledStatus:
        if mcrabPath == "":
            msg = "The dataset with label %s is not present!" % (sh_t + label + sh_n)
            mcrabInfoOutput.append(msg)
            Print(msg, True)
        else:
            msg = sh_n + "The dataset with label %s was found in directory %s" % (sh_t + label + sh_n, sh_t + mcrabPath + sh_n)
            mcrabInfoOutput.append(msg)
            Verbose(sh_t + msg + sh_n, True)
            return dataset.readFromMulticrabCfg(directory=mcrabPath)
    else:
        msg = "The dataset with label %s will not be considered" % (sh_t + label + sh_n)
        mcrabInfoOutput.append(msg)
        Print(msg, True)
    return None


def Verbose(msg, printHeader=False):
    '''
    Calls Print() only if verbose options is set to true
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

def CreateTarball(myOutputDirectories, opts):
    '''
    Creats a tarball with entire output
    making easier to transfer results other machines
    '''
    if not opts.tarball:
        return

    # Create filename with a time stamp
    myTimestamp = time.strftime("%y%m%d_%H%M%S", time.gmtime(time.time()))
    # myFilename  = "datacards_archive_%s.tgz" % (myTimestamp) # Adopt directory name instead!
    myFilename  = myOutputDirectories[0] + ".tgz"
    fTarball    = tarfile.open(myFilename, mode="w:gz")
    
    # For-loop: All output dirs
    for d in myOutputDirectories:
        fTarball.add(d)

    # Close the tarball
    fTarball.close()

    # Inform user
    msg = "Created archive of results directories to "
    Print(msg + SuccessStyle() + myFilename + sh_n, False)
    return

def main(opts, moduleSelector, multipleDirs):

    # Fixme: Is this used/needed?
    if 0:
        gc.set_debug(gc.DEBUG_LEAK | gc.DEBUG_STATS)
        gc.set_debug(gc.DEBUG_STATS)
        ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)
        gc.set_debug(gc.DEBUG_STATS)

    # Catch any errors in the input datacard (options added to avoid double printouts from the template card; if there are any!)
    msg = "Checking the syntax of the datacard file %s by compiling it." % (opts.datacard) 

    Verbose(sh_a + msg + sh_n, True)
#    os.system("python -m py_compile \"%s\" >& tmp_validate.txt" % opts.datacard)
#    if os.stat("tmp_validate.txt").st_size != 0:
#        msg = "Datacard \"%s\" is invalid." % (opts.datacard)
#        Print(sh_e + msg + sh_n, True)
#        os.system("cat tmp_validate.txt")
#    else:
#        msg = "Datacard %s is valid" % (sh_h + opts.datacard + sh_n)
#        Print(msg, True)
#        os.system("rm -f tmp_validate.txt")
        

    # Load the datacard
    Verbose("Loading datacard \"%s\"." % (opts.datacard), True)
    config = aux.load_module(opts.datacard)

    # Replace source directory if necessary
    if multipleDirs != None:
        config.Path = multipleDirs
    else:
        try:
            Print("Will produce datacards using results in directory %s (defined by variable \"Path\" in input datacard)." % (sh_h + config.Path+ sh_n), True)
        except AttributeError:
            msg = "The imported file %s has no attribute \"Path\". You must either define this parameter in the imported module or use the --dir option. Exit" % (config.__file__)
            raise AttributeError(msg, True) # raise AttributeError(sh_e + msg + sh_n, True)

    # Obtain dataset creators (also check multicrab directory existence)
    Verbose("Checking input multicrab directory presence", True)
    multicrabPaths  = PathFinder.MulticrabPathFinder(config.Path, opts.analysisType, opts.verbose)
    Verbose("Input directory is \"%s\"" % (config.Path), True)
    if opts.verbose:
        multicrabPaths.PrintInfo()

    # fixme - is this list needed?
    mcrabInfoOutput = []
    mcrabInfoOutput.append("Input directories:")

    # Determine dataset labels.
    if opts.analysisType == "HToTB":
        try:
            fakesFromData = (config.OptionFakeBMeasurementSource == "DataDriven")
        except AttributeError:
            msg = "The imported file \"%s\" has no attribute \"OptionFakeBMeasurementSource\". Please define this parameter in the imported module." % (config.__file__)
            raise AttributeError(msg, True)
    elif opts.analysisType in ["HToTauNu", "HToHW"]:
        try:
#            fakesFromData = (config.OptionGenuineTauBackgroundSource == "DataDriven")
#            fakesFromData = (config.OptionFakeTauMeasurementSource == "DataDriven")
	     fakesFromData = True
        except AttributeError:
            msg = "The imported file \"%s\" has no attribute \"OptionGenuineTauBackgroundSource\". Please define this parameter in the imported module." % (config.__file__)
            raise AttributeError(msg, True)


    # Define lists for signal/bkg dataset creators
    signalDsetCreators = []
    bkgDsetCreators    = []

    # Get the signal and background labels (analysis-specific)
    signalLabels, bkgLabels = GetDatasetLabels(config, fakesFromData, opts)

    Verbose("Creating signal datasets", True)
    sigDset = getSignalDsetCreator(multicrabPaths, signalLabels[0], mcrabInfoOutput)
    sigDset.setLabel(signalLabels[0])
    signalDsetCreators.append(sigDset)

    Verbose("Creating the background datasets", True)
    # For-loop: All background labels
    for i, bkg in enumerate(bkgLabels, 1):
        Verbose("%d) bkg = %s" % (i, bkg), i==1)
        bkgDset = getBkgDsetCreator(multicrabPaths, bkg, fakesFromData, mcrabInfoOutput)
        bkgDset.setLabel(bkg)
        bkgDsetCreators.append(bkgDset)

    # Print the datasets contained in each of the dataset creators
    PrintDsetInfo(signalDsetCreators + bkgDsetCreators, True)
    
    # Check options that are affecting the validity of the results
    CheckOptions(config)

    # Find list of available eras, search modes, and optimization modes common for all multicrab directories
    Verbose("Find list of available eras, search modes, and optimization modes common for all multicrab directories", True)
    SetModuleSelectorSources(moduleSelector, signalDsetCreators, bkgDsetCreators, opts)
    
    # Summarise the consequences of the user choises
    PrintOptions(moduleSelector)

    # Produce the datacards
    startTime = time.time()
    myOutputDirectories = []
    nDatacards = 0
    nEras      = len(moduleSelector.getSelectedEras())
    nModes     = len(moduleSelector.getSelectedSearchModes())
    nOpts      = len(moduleSelector.getSelectedOptimizationModes())

    # For-loop: Data eras
    for i, era in enumerate(moduleSelector.getSelectedEras(), 1):        
        msg  = "{:<9} {:>3} {:<1} {:<3} {:<50}".format("Era", "%i" % i, "/", "%s:" % (nEras), era)
        Print(sh_a + msg + sh_n, i==1)

        # For-loop: Search modes
        for j, searchMode in enumerate(moduleSelector.getSelectedSearchModes(), 1):
            msg = "{:<9} {:>3} {:<1} {:<3} {:<50}".format("Mode", "%i" % j, "/", "%s:" % (nModes), searchMode)
            Print(sh_a + msg + sh_n, False)

            # For-loop: Optimization modes            
            for k, optimizationMode in enumerate(moduleSelector.getSelectedOptimizationModes(), 1):
                nDatacards +=1

                msg = "{:<9} {:>3} {:<1} {:<3} {:<50}".format("Opt", "%i" % k, "/", "%s:" % (nOpts), optimizationMode)
                Print(sh_a + msg + sh_n, False)
                
                Verbose("Create the datacard generator & check config file contents", True)
                dcgen = DataCard.DataCardGenerator(opts, config)
                dcgen.setDsetMgrCreators(signalDsetCreators, bkgDsetCreators)

                Verbose("Do the heavy stuff", True)
                myDir = dcgen.doDatacard(era, searchMode, optimizationMode, mcrabInfoOutput)
                myOutputDirectories.append(myDir)

    # Timing calculations
    endTime = time.time()
    totTime = (endTime-startTime)
    avgTime = (totTime)/float(nDatacards)
    Verbose("Running took on average %.1f s per datacard (elapsed time = %.1f s,  datacards = %d) " % (avgTime, totTime, nDatacards) )
    
    # Generate plots for systematics
    if opts.systAnalysis:
        for d in myOutputDirectories:
            Print("Generating systematics plots for %s" (d) )
            os.chdir(d)
            os.system("../plotShapes.py")
            os.chdir("..")

    # Inform user
    for d in myOutputDirectories:
        msg = "Created results dir %s" % (sh_s + d + sh_n)
        Print(msg)

    # Optionally, create a tarball with all the results
    CreateTarball(myOutputDirectories, opts)
    
    if 0:
        gc.collect()
        ROOT.SetMemoryPolicy( ROOT.kMemoryHeuristics)
        memoryDump()
    return
    

if __name__ == "__main__":

    # Default Values
    HELP          = False
    COMBINE       = True
    SYSTANALYSIS  = False
    SHAPESENSITIV = False
    SHOWCARD      = False
    DATASETDEBUG  = False
    CONFIGDEBUG   = False
    MININGDEBUG   = False
    QCDDEBUG      = False
    SHAPEDEBUG    = False
    CTRLPLOTDEBUG = False
    VERBOSE       = False    
    TARBALL       = False
    # New
    ANALYSISTYPE  = "HToTauNu"
    DATACARD      = None
    PROFILER      = False
    DIRECTORIES   = None


    # Object for selecting data eras, search modes, and optimization modes
    myModuleSelector = analysisModuleSelector.AnalysisModuleSelector() 

    parser = OptionParser(usage="Usage: %prog [options]", add_help_option=False, conflict_handler="resolve")

    parser.add_option("-h", "--help", dest="help", action="store_true", default=HELP, 
                      help="Show this help message and exit [default: %s]" % (HELP) )

    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=VERBOSE,
                      help="Print more information [default: %s]" % (VERBOSE) )

    parser.add_option("--datacard", dest="datacard", action="store", 
                      help="Name (incl. path) of the datacard to be used as an input [default: %s]" % (DATACARD) )

    myModuleSelector.addParserOptions(parser)

    parser.add_option("--combine", dest="combine", action="store_true", default=COMBINE,
                      help="Generate datacards for Combine [default=%s]" % (COMBINE) )

    parser.add_option("--tarball", dest="tarball", action="store_true", default=TARBALL,
                      help="In addition to a dir with all the output, create also a tarball for easy transfer of results to other machines [default=%s]" % (TARBALL) )

    parser.add_option("--dir", dest="directories", action="append", 
                      help="Each time this argument is called it appends a new directory name to a directory list. The script will produce datacards for each dir in the list [default=%s]" % (DIRECTORIES) )

    parser.add_option("--systAnalysis", dest="systAnalysis", action="store_true", default=SYSTANALYSIS, 
                      help="Runs the macro for generating systematic uncertainties plots [default: %s]" % (SYSTANALYSIS) )

    parser.add_option("--testShapeSensitivity", dest="testShapeSensitivity", action="store_true", default=SHAPESENSITIV,
                      help="Creates datacards for varying each shape nuisance up and down by 1 sigma [default: %s]" % (SHAPESENSITIV) )

    parser.add_option("--showcard", dest="showDatacard", action="store_true", default=SHOWCARD, 
                      help="Print datacards also to screen [default: %s]" % (SHOWCARD) )

    parser.add_option("--debugDatasets", dest="debugDatasets", action="store_true", default=DATASETDEBUG, 
                      help="Enable debugging print for datasetMgr contents [default: %s]" % (DATASETDEBUG) )

    parser.add_option("--debugConfig", dest="debugConfig", action="store_true", default=CONFIGDEBUG,
                      help="Enable debugging print for config parsing [default: %s]" % (CONFIGDEBUG) )
    
    parser.add_option("--debugMining", dest="debugMining", action="store_true", default=MININGDEBUG, 
                      help="Enable debugging print for data mining [default: %s]" % (MININGDEBUG) )
    
    parser.add_option("--debugQCD", dest="debugQCD", action="store_true", default=QCDDEBUG, 
                      help="Enable debugging print for QCD measurement [default: %s]" % (QCDDEBUG) )

    parser.add_option("--debugControlPlots", dest="debugControlPlots", action="store_true", default=CTRLPLOTDEBUG,
                      help="Enable debugging print for data-driven control plots [default: %s]" % (CTRLPLOTDEBUG) )

    parser.add_option("--profiler", dest="profiler", action="store_true", default=PROFILER,
                      help="Enable running of profiler to provide set of statistics that describe how often and for how long various parts of the program executed [default: %s]" % (PROFILER) )

    parser.add_option("--analysisType", dest="analysisType", default=ANALYSISTYPE,
                      help="Flag to indicate the analysis typeh2tb analysis (e.g. \"HToTauNu\", \"HToTB\", \"HToHW\") [default: %s]" % (ANALYSISTYPE) )
    
    (opts, args) = parser.parse_args()


    # Print a help message and exit
    if opts.help:
        parser.print_help()
        sys.exit()

    # Sanity check (analysis type)
    myAnalyses = ["HToTauNu", "HToTB", "HToHW"]
    if opts.analysisType not in myAnalyses:
        msg = "Invalid analysis type \"%s\". Please selected one of the following: \"%s" % (opts.analysisType, "\", \"".join(myAnalyses) + "\"")
        raise Exception(sh_e + msg + sh_n)
                      
    # Sanity check (datacard)
    if opts.datacard == None:
        msg = "No datacard file provided. Will use analysis-specific datacard."
        Print(sh_h + msg + sh_n, True)
        opts.datacard = "datacard_%s.py" % (opts.analysisType)
    if os.path.isfile(opts.datacard):
        msg = "Datacard is \"%s\"." % (opts.datacard)
        Verbose(sh_a + msg + sh_n, True)
    else:
        msg = "File not found. Please check that the datacard file \"%s\" exists." % (opts.datacard)
        Print(sh_e + msg + sh_n, True)
        sys.exit()

    # Execute the program with 
    if opts.directories == None:
        Verbose("Input directory will be read from the datacard (fron variable \"Path\"). ", True)
        if opts.profiler:
            cProfile.run("main(opts, myModuleSelector, multipleDirs=None)")
            sys.exit()
        else:
            main(opts, myModuleSelector, multipleDirs=None)
            sys.exit()
    else:
        Verbose("Input directory provided as argument in script execution (Will overwrite the directory (if) defined in the input datacard as \"Path\").", True)

    # Get a list containing the names of all contents in the current working directory
    myDirs = os.listdir(".")
    nDirs = len(opts.directories)
        
    Verbose("Creating datacards for %i directories" % (nDirs) )
    # For-loop: All directories
    for i, d in enumerate(opts.directories, 1):
        
        # Remove "/" from end of dir
        d = aux.rchop(d, "/")

        # Print progress
        msg   = "{:<9} {:>3} {:<1} {:<3} {:<50}".format("Directory", "%i" % i, "/", "%s:" % (nDirs), d)
        Print(sh_t + msg + sh_n, i==1)

        # Sanity check (directory is located under current working directory)
        if d not in myDirs:
            msg = "Directory not found. Please check that directory \"%s\" exists under current working directory." % (d)
            raise Exception(sh_e + msg + sh_n)
       
        if opts.profiler:
            cProfile.run("main(opts, myModuleSelector, multipleDirs=d)")
        else:
            main(opts, myModuleSelector, multipleDirs=d)            
