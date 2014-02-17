## \package CommonLimitTools
# Python interface common for LandS and Combine and running them with multicrab
#
# The interface for casual user is provided by the functions
# generateMultiCrab() (for LEP-CLs and LHC-CLs) and
# produceLHCAsymptotic (for LHC-CLs asymptotic).
#
# The multicrab configuration generation saves various parameters to
# taskdir/configuration.json, to be used in by landsMergeHistograms.py
# script. The script uses tools from this module, which write the
# limit results to taskdir/limits.json. I preferred simple text format
# over ROOT files due to the ability to read/modify the result files
# easily. Since the amount of information in the result file is
# relatively small, the performance penalty should be negligible.

import os
import re
import sys
import glob
import stat
import json
import time
import random
import shutil
import subprocess

import multicrab
import multicrabWorkflows
import git
import aux

## Deduces from directory listing the mass point list
def obtainMassPoints(pattern):
    myPattern = pattern%" "
    mySplit = myPattern.split(" ")
    prefix = mySplit[0]
    suffix = mySplit[1]
    myList = os.listdir(".")
    myMasses = []
    for item in myList:
        if prefix in item:
            myMasses.append(item.replace(prefix,"").replace(suffix,""))
    myMasses.sort()
    return myMasses

## Reads luminosity from a datacard and returns it
def readLuminosityFromDatacard(myPath, filename):
    lumi_re = re.compile("luminosity=\s*(?P<lumi>\d+\.\d+)")
    fname = os.path.join(myPath, filename)
    f = open(fname)
    myLuminosity = None
    for line in f:
        match = lumi_re.search(line)
        if match:
            #self.lumi = str(1000*float(match.group("lumi"))) # 1/fb -> 1/pb
            # Nowadays the luminosity is in 1/pb
            myLuminosity = match.group("lumi")
            f.close()
            return myLuminosity
    raise Exception("Did not find luminosity information from '%s'" % fname)

## Returns true if mass list contains only heavy H+
def isHeavyHiggs(massList):
    for m in massList:
        if int(m) < 175:
            return False
    return True

## Returns true if mass list contains only light H+
def isLightHiggs(massList):
    for m in massList:
        if int(m) > 175:
            return False
    return True

## Create OptionParser, and add common LandS options to OptionParser object
#
# \param lepDefault     Boolean for the default value of --lep switch (if None, switch is not added)
# \param lhcDefault     Boolean for the default value of --lhc switch (if None, switch is not added)
# \param lhcasyDefault  Boolean for the default value of --lhcasy switch (if None, switch is not added)
#
# \return optparse.OptionParser object
def createOptionParser(lepDefault=None, lhcDefault=None, lhcasyDefault=None):
    parser = OptionParser(usage="Usage: %prog [options] [datacard directories]\nDatacard directories can be given either as the last arguments, or with -d.")

    # Switches for different CLs flavours, the interpretation of these
    # is in the generate* scripts (i.e. in the caller responsibility)
    if lepDefault != None:
        parser.add_option("--lep", dest="lepType", default=lepDefault, action="store_true",
                          help="Use hybrid LEP-CLs (default %s)" % str(lepDefault))
    if lhcDefault != None:
        parser.add_option("--lhc", dest="lhcType", default=lhcDefault, action="store_true",
                          help="Use hybrid LHC-CLs (default %s)" % str(lhcDefault))
    if lhcasyDefault != None:
        parser.add_option("--lhcasy", dest="lhcTypeAsymptotic", default=lhcasyDefault, action="store_true",
                          help="Use asymptotic LHC-CLs (default %s)" % str(lhcasyDefault))

    # Datacard directories
    parser.add_option("-d", "--dir", dest="dirs", type="string", action="append", default=[],
                      help="Datacard directories to create the LandS MultiCrab tasks into (default: use the working directory")
    parser.add_option("--create", dest="multicrabCreate", action="store_true", default=False,
                      help="Run 'multicrab -create' for each multicrab task directory")

    return parser

## Parse OptionParser object
#
# \param parser   optparse.OptionParser object
#
# \return Options object
def parseOptionParser(parser):
    (opts, args) = parser.parse_args()
    opts.dirs.extend(args)
    if len(opts.dirs) == 0:
        opts.dirs = ["."]
    return opts

## Class to hold the limit results
class Result:
    ## Constructor
    def __init__(self, mass):
        self.mass                = mass
        self.observed            = None
        self.expected            = None

    ## Check if the result is empty, i.e. no limits has been assigned
    def empty(self):
        return self.observed == None and self.expected == None

## Collection of Result objects
class ResultContainer:
    ## Constructor
    #
    # \param path  Path to the multicrab directory (where configuration.json exists)
    def __init__(self, path):
        self.path = path

        # Read task configuration json file
        configFile = os.path.join(path, "configuration.json")
        f = open(configFile)
        self.config = json.load(f)
        f.close()

        # Read the luminosity, use tau+jets one if that's available, if not, use the first one
        datacards = self.config["datacards"]
        taujetsDc = None
        for dc in datacards:
            if "hplushadronic" in dc:
                taujetsDc = dc
                break
        if taujetsDc != None:
            self._readLuminosityTaujets(taujetsDc % self.config["masspoints"][0])
        else:
            self._readLuminosityLeptonic(self.config["datacards"][0] % self.config["masspoints"][0])

        self.results = []

    ## Append a result object to the list
    #
    # \param obj   Result object
    def append(self, obj):
        self.results.append(obj)

    ## Read luminosity from a datacard assuming it follows the tau+jets convention
    #
    # \param filename  Name of the datacard file inside the multicrab directory
    def _readLuminosityTaujets(self, filename):
        self.lumi = readLuminosityFromDatacard(self.path, filename)

    ## Read luminosity from a datacard assuming it follows the leptonic convention
    #
    # \param filename  Name of the datacard file inside the multicrab directory
    #
    # \todo This needs to be updated, it get's the scale wrong
    def _readLuminosityLeptonic(self, filename):
        scale_re = re.compile("lumi scale (?P<scale>\S+)")
        lumi_re = re.compile("lumi=(?P<lumi>\S+)")
        scale = None
        lumi = None

        fname = os.path.join(self.path, filename)
        f = open(fname)
        for line in f:
            match = scale_re.search(line)
            if match:
                scale = float(match.group("scale"))
                continue
            match = lumi_re.search(line)
            if match:
                lumi = float(match.group("lumi"))
                break
        f.close()
        if lumi == None:
            raise Exception("Did not find luminosity information from '%s'" % fname)
        if scale != None:
            lumi *= scale
        self.lumi = str(lumi)

    ## Get the integrated luminosity as a string in 1/pb
    def getLuminosity(self):
        return self.lumi

    ## Print the limits
    def print2(self,unblindedStatus=False):
        print
        print "                  Expected"
        print "Mass  Observed    Median       -2sigma     -1sigma     +1sigma     +2sigma"
        format = "%3s:  %-9s   %-10s   %-10s  %-10s  %-10s  %-10s"
        massIndex = [(int(self.results[i].mass), i) for i in range(len(self.results))]
        massIndex.sort()
        for mass, index in massIndex:
            result = self.results[index]
            if result.empty():
                continue
            if unblindedStatus:
                print format % (result.mass, result.observed, result.expected, result.expectedMinus2Sigma, result.expectedMinus1Sigma, result.expectedPlus1Sigma, result.expectedPlus2Sigma)
            else:
                print format % (result.mass, "BLINDED", result.expected, result.expectedMinus2Sigma, result.expectedMinus1Sigma, result.expectedPlus1Sigma, result.expectedPlus2Sigma)
        print

    ## Store the results in a limits.json file
    #
    # \param data   Dictionary of additional data to be stored
    def saveJson(self, data={}, unblindedStatus=False):
        output = {}
        output.update(data)
        output.update({
                "luminosity": self.getLuminosity(),
                "masspoints": {}
                })

        massIndex = [(int(self.results[i].mass), i) for i in range(len(self.results))]
        massIndex.sort()
        for mass, index in massIndex:
            result = self.results[index]
            if result.empty():
                continue

            if unblindedStatus:
                output["masspoints"][result.mass] = {
                    "mass": result.mass,
                    "observed": result.observed,
                    "expected": {
                        "-2sigma": result.expectedMinus2Sigma,
                        "-1sigma": result.expectedMinus1Sigma,
                        "median": result.expected,
                        "+1sigma": result.expectedPlus1Sigma,
                        "+2sigma": result.expectedPlus2Sigma,
                        }
                    }
            else:
                output["masspoints"][result.mass] = {
                    "mass": result.mass,
                    "observed": 0,
                    "expected": {
                        "-2sigma": result.expectedMinus2Sigma,
                        "-1sigma": result.expectedMinus1Sigma,
                        "median": result.expected,
                        "+1sigma": result.expectedPlus1Sigma,
                        "+2sigma": result.expectedPlus2Sigma,
                        }
                    }

            if unblindedStatus:
                if hasattr(result, "observedError"):
                    output["masspoints"][result.mass]["observed_error"] = result.observedError
            if hasattr(result, "expectedError"):
                output["masspoints"][result.mass]["expected"].update({
                        "-2sigma_error": result.expectedMinus2SigmaError,
                        "-1sigma_error": result.expectedMinus1SigmaError,
                        "median_error": result.expectedError,
                        "+1sigma_error": result.expectedPlus1SigmaError,
                        "+2sigma_error": result.expectedPlus2SigmaError,
                        })


        fname = os.path.join(self.path, "limits.json")
        f = open(fname, "wb")
        json.dump(output, f, sort_keys=True, indent=2)
        f.close()
        return fname


## Base class to generate (LEP-CLs, LHC-CLs) multicrab configuration, or run (LHC-CLs asymptotic)
#
# The class is not intended to be used directly by casual user, but
# from generateMultiCrab() and produceLHCAsymptotic()
class LimitMultiCrabBase:
    ## Constructor
    #
    # \param directory          Datacard directory
    # \param massPoints         List of mass points to calculate the limit for
    # \param datacardPatterns   List of datacard patterns to include in the
    #                           limit calculation
    # \param rootfilePatterns   List of shape ROOT file patterns to include
    #                           in the limit calculation
    # \param clsType            Object defining the CLs flavour (should be either
    #                           LEPType, or LHCType).
    def __init__(self, directory, massPoints, datacardPatterns, rootfilePatterns, clsType):
        self.datacardDirectory = directory
        self.massPoints = massPoints
        self.datacardPatterns = datacardPatterns
        self.rootfilePatterns = rootfilePatterns
        self.clsType = clsType.clone()
        self.jobsCreated = False
        self.datacards = {}
        self.rootfiles = {}
        self.scripts   = []
        self.configuration = {}

        if not os.path.isdir(directory):
            raise Exception("Datacard directory '%s' does not exist" % directory)

        if clsConfig != None:
            self.configuration["clsConfig"] = clsConfig

        for mass in self.massPoints:
            for dc in datacardPatterns:
                fname = os.path.join(self.datacardDirectory, dc % mass)
                if not os.path.isfile(fname):
                    raise Exception("Datacard file '%s' does not exist!" % fname)

                aux.addToDictList(self.datacards, mass, fname)

            for rf in rootfilePatterns:
                fname = os.path.join(self.datacardDirectory, rf % mass)
                if not os.path.isfile(fname):
                    raise Exception("ROOT file (for shapes) '%s' does not exist!" % fname)

                aux.addToDictList(self.rootfiles, mass, fname)

    ## Create the multicrab task directory
    #
    # \param postfix   Additional string to be included in the directory name
    def _createMultiCrabDir(self, prefix, postfix):
        if len(postfix) > 0:
            prefix += "_"+postfix
        self.dirname = multicrab.createTaskDir(prefix=prefix, path=self.datacardDirectory)
        self.clsType.setDirectory(self.dirname)

    ## Copy input files for LandS (datacards, rootfiles) to the multicrab directory
    def copyInputFiles(self):
        for d in [self.datacards, self.rootfiles]:
            for mass, files in d.iteritems():
                for f in files:
                    shutil.copy(f, self.dirname)
        shutil.copy(self.exe, self.dirname)

    ## Write  shell scripts to the multicrab directory
    def writeScripts(self):
        for mass, datacardFiles in self.datacards.iteritems():
            self.clsType.createScripts(mass, datacardFiles)

    ## Write crab.cfg to the multicrab directory
    # \param crabScheduler      CRAB scheduler to use
    # \param crabOptions        Dictionary for specifying additional CRAB
    #                           options. The keys correspond to the
    #                           sections in crab.cfg. The values are lists
    #                           containing lines to be appended to the
    #                           section.
    # \param outputFile         list of output files
    def writeCrabCfg(self, crabScheduler, crabOptions, outputFile):
        filename = os.path.join(self.dirname, "crab.cfg")
        fOUT = open(filename,'w')
        fOUT.write("[CRAB]\n")
        fOUT.write("jobtype                 = cmssw\n")
        fOUT.write("scheduler               = %s\n" % crabScheduler)
        fOUT.write("use_server              = 0\n")
        if "CRAB" in crabOptions:
            for line in crabOptions["CRAB"]:
                fOUT.write(line+"\n")
        fOUT.write("\n")

        fOUT.write("[CMSSW]\n")
        fOUT.write("datasetpath             = none\n")
        fOUT.write("pset                    = none\n")
        fOUT.write("number_of_jobs          = 1\n")
        fOUT.write("output_file             = %s\n"%",".join(outputFile))
        if "CMSSW" in crabOptions:
            for line in crabOptions["CMSSW"]:
                fOUT.write(line+"\n")
        fOUT.write("\n")

        fOUT.write("[USER]\n")
        fOUT.write("return_data             = 1\n")
        fOUT.write("copy_data               = 0\n")
        if "USER" in crabOptions:
            for line in crabOptions["USER"]:
                fOUT.write(line+"\n")
        fOUT.write("\n")

        if "GRID" in crabOptions:
            fOUT.write("[GRID]\n")
            for line in crabOptions["GRID"]:
                fOUT.write(line+"\n")
            fOUT.write("\n")

        fOUT.close()

    ## Write multicrab.cfg to the multicrab directory
    #
    # \param numberOfJobs   ValuePerMass object holding the information
    #                       of number of crab jobs (per mass point)
    def writeMultiCrabCfg(self, numberOfJobs):
        filename = os.path.join(self.dirname, "multicrab.cfg")
        fOUT = open(filename,'w')
        fOUT.write("[COMMON]\n")
        fOUT.write("CRAB.use_server              = 0\n")
        fOUT.write("CMSSW.datasetpath            = none\n")
        fOUT.write("CMSSW.pset                   = none\n")
        fOUT.write("\n")

        exe = self.exe.split("/")[-1]
        for mass in self.massPoints:
            inputFiles = [exe]+self.datacards[mass]
            if len(self.rootfiles) > 0:
                inputFiles += self.rootfiles[mass]
            self.clsType.writeMultiCrabConfig(fOUT, mass, inputFiles, numberOfJobs.getValue(mass))
            fOUT.write("\n\n")

        f = open(os.path.join(self.dirname, "configuration.json"), "wb")
        json.dump(self.configuration, f, sort_keys=True, indent=2)
        f.close()

    ## Run 'multicrab create' inside the multicrab directory
    def createMultiCrab(self):
        print "Creating multicrab task %s" % self.dirname
        print
        print "############################################################"
        print
        pwd = os.getcwd()
        os.chdir(self.dirname)
        ret = subprocess.call(["multicrab", "-create"])
        if ret != 0:
            raise Exception("'multicrab -create' failed with exit code %d in directory '%s'" % (ret, self.dirname))
        os.chdir(pwd)
        print
        print "############################################################"
        print

        self.jobsCreated = True

    def printInstruction(self):
        if self.jobsCreated:
            print "Multicrab cfg and jobs created. Type"
            print "cd %s && multicrab -submit" % self.dirname
        else:
            print "Multicrab cfg created. Type"
            print "cd %s && multicrab -create" % self.dirname


def generateMultiCrab(opts,
                      massPoints=defaultMassPoints,
                      datacardPatterns=defaultDatacardPatterns,
                      rootfilePatterns=defaultRootfilePatterns,
                      clsType=None,
                      numberOfJobs=None,
                      crabScheduler="arc",
                      crabOptions={},
                      postfix=""
                      ):
def produceLHCAsymptotic(opts,
                         massPoints=defaultMassPoints,
                         datacardPatterns=defaultDatacardPatterns,
                         rootfilePatterns=defaultRootfilePatterns,
                         clsType = None,
                         postfix=""
                         ):

## Helper class to manage mass-specific configuration values
#class LEPType:
#class LHCType:
#class LHCTypeAsymptotic:
#class ParseLandsOutput:
#def parseLandsMLOutput(outputFileName):
#class LandSInstaller:

