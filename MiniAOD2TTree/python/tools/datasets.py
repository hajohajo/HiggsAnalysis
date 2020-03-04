'''
DESCRIPTION:
All the DAS datasets used for the  HToTauNu, HToTB, HToHW, etc.. analyses.


USAGE:
Can be loaded and used when running on MiniAOD when submitting CRAB jobs


VALIDATE:
To validate a given datasets enable the "validate" boolean in the construction (disabled by default)

'''
#================================================================================================ 
# Import modules
#================================================================================================ 
import os
import subprocess
import json 
from subprocess import Popen, PIPE
import shlex
import HiggsAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles


#================================================================================================ 
# Variable definition
#================================================================================================ 
ss  = ShellStyles.SuccessStyle()
ns  = ShellStyles.NormalStyle()
ts  = ShellStyles.NoteStyle()
hs  = ShellStyles.HighlightAltStyle()
ls  = ShellStyles.HighlightStyle()
es  = ShellStyles.ErrorStyle()
cs  = ShellStyles.CaptionStyle()
cys = ShellStyles.CyanStyle()


#================================================================================================ 
# Class Definition
#================================================================================================ 
class DatasetGroup:
    def __init__(self, analysis):
        self.verbose   = False
        self.analysis  = analysis
        self.GroupDict = {}
        self.CreateGroups()
        return

    def SetVerbose(verbose):
        self.verbose = verbose
        return

    def Verbose(self, msg, printHeader=False):
        '''
        Simple print function. If verbose option is enabled prints, otherwise does nothing.
        '''
        if not self.verbose:
            return
        self.Print(msg, printHeader)
        return

    def Print(self, msg, printHeader=True):
        '''
        Simple print function. If verbose option is enabled prints, otherwise does nothing.
        '''
        fName = __file__.split("/")[-1]
        cName = self.__class__.__name__
        name  = fName + ": " + cName
        if printHeader:
                print "=== ", name
        print "\t", msg
        return

    def CreateGroups(self):
        '''
        Create dataset grouping in a dictionary for easy access.
        '''

        analyses = ["SignalAnalysis", "Hplus2tbAnalysis", "Hplus2hwAnalysis", "Hplus2hwWithTopAnalysis", "HToHWTrgEff",
                    "JetTriggers", "TopSystBDT", "TauLeg", "METLeg", "L1Study"]
        
        analysisTypes = ["HToTauNu", "HToTB", "HToHW", "HToHW_withTop"]
        analyses.extend(analysisTypes)
        if self.analysis not in analyses:
            raise Exception("Unknown analysis \"%s\". Please select one of the following: \"%s" % (self.analysis, "\", \"".join(analyses) + "\".") )

        # Add the default analysis Types as keys
        for a in analysisTypes:
            self.GroupDict[a] = dsetGroups_2016[a]

        self.GroupDict["SignalAnalysis"]   = dsetGroups_2016["HToTauNu"]
        self.GroupDict["Hplus2tbAnalysis"] = dsetGroups_2016["HToTB"]
        self.GroupDict["Hplus2hwAnalysis"] = dsetGroups_2016["HToHW"]
        self.GroupDict["Hplus2hwWithTopAnalysis"] = dsetGroups_2016["HToHW_withTop"]
        self.GroupDict["TauLeg"]           = dsetGroups_2016["TauLeg"]
        self.GroupDict["METLeg"]           = dsetGroups_2016["METLeg"]
        self.GroupDict["L1Study"]          = dsetGroups_2016["L1Study"]
        self.GroupDict["JetTriggers"]      = dsetGroups_2016["JetTriggers"]
        self.GroupDict["TopSystBDT"]       = dsetGroups_2016["TopTagSys"]
        self.GroupDict["HToHWTrgEff"]      = dsetGroups_2016["HToHWTrgEff"]
        return

    def GetDatasetNames(self):
        '''
        Return the dataset list according to the analysis name. 
        Uses pre-defined dictionary mapping: analysis->dataset list
        '''
        dsetNames = [d.getDataset() for d in self.GroupDict[self.analysis]]
        return dsetNames

    def GetDatasets(self):
        '''
        Return the dataset list according to the analysis name. 
        Uses pre-defined dictionary mapping: analysis->dataset list
        '''
        datasets = [d for d in self.GroupDict[self.analysis]]
        return datasets

    def PrintDatasets(self, printHeader=False):
        '''
        Print all datasets for given analysis
        '''
        if printHeader:
            self.Print("The datasets for analysis %s are:" % (hs + self.analysis + ns), True)
            
        for i, d in enumerate(self.GetDatasetNames(), 1):
            self.Print("%d) %s" % (i, d), False)
        return

class Dataset:
    def __init__(self, PDName, campaign, pString, PU, globalTag, ext, version, dataTier, dbs, dataVersion, lumiMask=None, validate=False, verbose=False):
        self.verbose       = verbose
        self.PDName        = PDName     # Primary Dataset Name
        self.Campaign      = campaign   # Run<Year><Era> for data
        self.ProcessString = pString    # Process string + Prefix to Global Tag (not all campaign have it)
        self.Pileup        = PU         
        self.GlobalTag     = globalTag
        self.Extension     = ext        # 0 for no extension
        self.Version       = version
        self.DataTier      = dataTier
        self.Instance      = "prod/%s" % (dbs)
        self.DBS           = dbs
        self.DataVersion   = dataVersion        
        self.IsData_       = self.getIsData()
        self.IsMC_         = self.getIsMC()
        self.Dataset       = self.getDataset()
        self.CheckedProxy  = False        
        self.Files         = []
        self.Sites         = []
        self.Summary       = []
        self.NBlocks       = 0
        self.NEvents       = 0
        self.NLumis        = 0
        self.NumBlock      = 0
        self.NumEvent      = 0
        self.NumFile       = 0
        self.NumLumi       = 0
        self.URL           = "https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod/%s&input=dataset=%s" % (self.DBS, self.getDataset())
        self.LumiMask      = self.getLumiMask_(lumiMask)
        self.RequestName   = "" # handled in multicrab.py
        if validate:
            self.validate()
        return

    def Verbose(self, msg, printHeader=False):
        if not self.verbose:
            return
        if printHeader:
            print "=== checkCondor.py:"
        if msg !="":
            print "\t", msg
        return

    def GetFName(self):
        fName = __file__.split("/")[-1]
        fName = fName.replace(".pyc", ".py")
        return fName

    def Print(self, msg, printHeader=True):
        fName = self.GetFName()
        if printHeader:
            print "=== ", fName
        if msg !="":
            print "\t", msg
        return

    def PrintFlushed(self, msg, printHeader=True):
        '''
        Useful when printing progress in a loop
        '''
        msg = "\r\t" + msg
        if printHeader:
            print "=== ", GetFName()
        sys.stdout.write(msg)
        sys.stdout.flush()
        return     

    def getPrimaryDatasetName(self):
        '''
        Convention:

        PROCESS_RANGETYPE-RANGELOWToRANGEHIGH_FILTER_TUNE_COMMENT_COMENERGY-GENERATO

        Links:
        https://monte-carlo-production-tools.gitbook.io/project/mccontact/rules-for-dataset-names
        '''
        return self.PDName

    def getRequestName(self):
        '''
        For usage when creating CRAB task.
        
        This is handled later with the multicrab.py script (for now)
        '''
        return self.RequestName

    def getDataset(self):
        '''
        data:
        /<dataset_name>/<AcquisitionEra>-<flow process string>_<request process string>_<request global tag>_<ext N if N>1>-v<M>/<data tier> 

        MC:
        /<primary_dataset_name>/<campaign>-<request process string_PU>_<request global tag>_<ext N if N>1>-v<version>/<data tier> 

        Links: 
        https://twiki.cern.ch/twiki/bin/view/CMS/PdmVMcMGlossary#requests_output_dataset
        https://twiki.cern.ch/twiki/bin/view/CMS/ProductionDataSetNames
        '''

        if not isinstance(self.PDName, str):        
            msg = "Please provide a string for the primary dataset name"
            raise Exception(es + msg + ns)        
        if not isinstance(self.Campaign, str):
            msg = "Please provide a string for the dataset campaign"
            raise Exception(es + msg + ns)        

        # Construct the name
        dataset = "/%s/%s-" % (self.PDName, self.Campaign)
        
        # Add ProcessString (if any)
        if isinstance(self.ProcessString, str):
            if self.ProcessString != "":
                dataset += "%s" % (self.ProcessString)

        if isinstance(self.Pileup, str):
            if self.Pileup != "":
                if len(self.ProcessString) > 0:
                    dataset += "_%s" % (self.Pileup)
                else:
                    dataset += "%s" % (self.Pileup)

        if isinstance(self.GlobalTag, str):
            dataset += "_%s" % (self.GlobalTag)

        if self.Extension > 0:
            dataset += "_ext%d" % (self.Extension)

        if isinstance(self.Version, int):
            dataset += "-v%d" % (self.Version)
        else:
            #msg = "Please provide an integer for the dataset version"
            #raise Exception(es + msg + ns)        
            pass

        if isinstance(self.DataTier, str):
            dataset += "/%s" % (self.DataTier)
        else:
            msg = "Please provide a valid string for the dataset TIER (e.g. \"MINIAOD\", \"MINIAODSIM\")."
            raise Exception(es + msg + ns)        

        # Sometimes it can happen that a string is empty resulting in "--" character. Simple fix below
        dataset = dataset.replace("--", "-")
        return dataset

    def getCampaign(self):
        return self.Campaign

    def getProcessString(self):
        '''
        A process string is used to customize 
        a dataset name in cases when datasets
        with identical PD name exist
        '''
        return self.ProcessString

    def getPileup(self):
        return self.Pileup
    
    def getGlobalTag(self):
       return self.GlobalTag

    def getExtension(self):
        '''        
        Extension number has to be changed if
        request is a statistics extension of 
        existing dataset.
        '''
        return self.Extension
    
    def getVersion(self):
        return self.Version

    def getDataTier(self):
        return self.DataTier
    
    def getIsData(self):
        if "data" in self.DataVersion:
            return True
        return False

    def getDataVersion(self):
        return self.DataVersion

    def getIsMC(self):
        return not self.getIsData()

    def IsData(self):
        return self.IsData_

    def IsMC(self):
        return self.IsMC_

    def getPrimaryName(self):
	return self.PDName

    def getName(self):
	return ls + self.Dataset + ns

    def getDBS(self):
        return self.DBS

    def checkVomsProxy(self):
        if self.CheckedProxy:
            return
        else:
           self.CheckedProxy = True

        # Setup the shell command
        cmds = []
        cmds.append("voms-proxy-info")
        cmds.append("--timeleft")

        self.Verbose(" ".join(cmds), True)
        process = Popen(cmds, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = process.communicate()
        if len(err) > 0:
            raise Exception(es + err + ns)
        else:
            msg = "Found %svalid%s VOMS proxy (%s seconds left)" % (ss, ns, output.replace("\n",""))
            self.Verbose(msg, True)
        return

    def getSummary(self):
        '''
        https://github.com/dmwm/dasgoclient
        '''
        self.checkVomsProxy()

        # Setup the shell command ("--json" option to return results in JSON data-format)
        if self.IsMC():
            cmd = "dasgoclient --query='summary dataset=%s instance=%s' --json" % (self.Dataset, self.Instance)
        else:
            cmd = "dasgoclient --query='summary dataset=%s' --json" % (self.Dataset)
        self.Verbose(cmd, True)
        
        # Provide an array of arguments (not a string of arguments) to subprocess
        cmds = shlex.split(cmd)
        process = Popen(cmds, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = process.communicate()
        if len(err) > 0:
            raise Exception(es + err + ns)
        else:
            parsed = json.loads(output)
            if len(parsed) == 0:
                msg  = "Invalid dataset \"%s\". " % (self.getDataset())
                msg += "Please double-check the primary dataset name, campaign, pString, globalTag, ext, version, dataTier, and dbs parameters!\n"
                msg += hs + cmd
                #msg += self.getURL()
                raise Exception(es + msg + ns)

            self.Verbose(json.dumps(parsed[0], indent=0, sort_keys=True), True)
            self.Verbose(json.dumps(parsed[0]["das"]["expire"], indent=0, sort_keys=True), True)

            self.NBlocks  = json.dumps(parsed[0]["summary"][0]["nblocks"], indent=0, sort_keys=True)
            self.NEvents  = json.dumps(parsed[0]["summary"][0]["nevents"], indent=0, sort_keys=True)
            self.NLumis   = json.dumps(parsed[0]["summary"][0]["nlumis"], indent=0, sort_keys=True)
            self.NumBlock = json.dumps(parsed[0]["summary"][0]["num_block"], indent=0, sort_keys=True)
            self.NumEvent = json.dumps(parsed[0]["summary"][0]["num_event"], indent=0, sort_keys=True)
            self.NumFile  = json.dumps(parsed[0]["summary"][0]["num_file"], indent=0, sort_keys=True)
            self.NumLumi  = json.dumps(parsed[0]["summary"][0]["num_lumi"], indent=0, sort_keys=True)
            self.Summary  = output.split("\n")        
        return self.Summary

    def getSites(self):
        '''
        https://github.com/dmwm/dasgoclient
        '''
        self.checkVomsProxy()
        if len(self.Sites) > 0:
            return self.Sites

        # Setup the shell command
        cmd = "dasgoclient --query='site dataset=%s instance=%s'" % (self.Dataset, self.Instance)
        self.Verbose(cmd, True)
        
        # Provide an array of arguments (not a string of arguments) to subprocess
        cmds = shlex.split(cmd)
        process = Popen(cmds, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = process.communicate()
        if len(err) > 0:
            raise Exception(es + err + ns)
        else:
            self.Sites = output.split("\n")
        return self.Sites
   
    def getNumOfSites(self):
        return len(self.getSites())

    def getFiles(self):
        '''
        https://github.com/dmwm/dasgoclient
        '''
        self.checkVomsProxy()
        if len(self.Files) > 0:
            return self.Files

        # Setup the shell command
        cmd = "dasgoclient --query='file dataset=%s instance=%s'" % (self.Dataset, self.Instance)
        self.Verbose(cmd, True)
        
        # Provide an array of arguments (not a string of arguments) to subprocess
        files = []
        cmds  = shlex.split(cmd)
        process = Popen(cmds, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = process.communicate()
        if len(err) > 0:
            raise Exception(es + err + ns)
        else:
            files = output.split("\n")

        # Remove empty strings ("")
        self.Files = [f for f in files if f]

        # Sanity check
        #for f in self.Files:
        #    print "\"%s\"" % (f)
        
        return self.Files

    def getNumOfFiles(self):
        return len(self.getFiles())
        
    def getNumOfEvents(self):
        if self.NEvents == 0:
            self.getSummary()
        return self.NEvents

    def getNumOfLumis(self):
        if self.NLumis == 0:
            self.getSummary()
        return self.NLumis
            
    def printFiles(self, printHeader=True):
        if len(self.Files) == 0:
            self.getFiles()
        for i, f in enumerate(self.Files, 1):
            self.Print(f, i==1 and printHeader==True)
        return

    def printSites(self, printHeader=True):
        if len(self.Sites) == 0:
            self.getSites()
        for i, f in enumerate(self.Sites, 1):
            self.Print(f, i==1 and printHeader==True)
        return

    def printSummary(self, printHeader=True):
        if len(self.Summary) == 0:
            self.getSummary()
        for i, line in enumerate(self.Summary, 1):
            self.Verbose(line, i==1)
            
        self.Print("NBlocks  = %s" % self.NBlocks, printHeader)
        self.Print("NEvents  = %s" % self.NEvents, False)
        self.Print("NLumis   = %s" % self.NLumis, False)
        self.Print("NumBlock = %s" % self.NumBlock, False)
        self.Print("NumEvent = %s" % self.NumEvent, False)
        self.Print("NumFile  = %s" % self.NumFile, False)
        self.Print("NumLumi  = %s" % self.NumLumi, False)
        return

    def getURL(self):
        return hs + self.URL + ns

    def printURL(self, printHeader=True):
        self.Print(self.getURL(), printHeader)

    def printName(self, printHeader=False):
        self.Print(self.getName(), printHeader)
        
    def validate(self):
        self.Verbose("Validating dataset %s" % (self.getName()), True)
        self.printName(True)
        self.printSummary(True)
        #self.printFiles(True)
        #self.printInfo()
        self.printSites(True)
        self.printURL(True)
        self.Print("", False)
        return
        
    def getLumiMask_(self, lumiMask):
        '''
        For GOLDEN JSON files see: https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2016Analysis

        lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt" # ICHEP dataset 271036-276811
        lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt" # without L1T certification?
        lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt" #PromptReco
        lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
        '''
        if lumiMask == None:
            return lumiMask

        if self.IsMC() == True:
            if lumiMask != None:
                msg = "The dataset \"%s\" is of type MC but the lumi-mask file is not None (=\"%s\")" % (self.getName(), lumiMask)
                raise Exception(es + msg + ns)
            else:
                return lumiMask
        
        lumiMaskFile = os.path.join(os.environ['CMSSW_BASE'],"src/HiggsAnalysis/MiniAOD2TTree/data", lumiMask)
        if not os.path.exists(lumiMaskFile):
            msg = "The lumi-mask file \"%s\" for dataset %s does not exist!" % (lumiMaskFile, self.getName())
            raise Exception(es + msg + ns)
        else:
            if not os.path.isfile(lumiMaskFile):
                msg = "The path \"%s\" used as lumi-mask file for dataset %s is not a file!" % (lumiMaskFile, self.getName() + es)
                self.Print(es + msg + ns, True)
            else:
                pass
        return lumiMaskFile
    
    def getLumiMaskFileBase(self):
        if self.LumiMask == None:
            return "N/A"
        else:
            return os.path.basename(self.LumiMask)

    def getLumiMask(self):
        return self.getLumiMaskFileBase()

    def getLumiMaskFile(self):
        if self.LumiMask == None:
            return "N/A"
        else:
            return self.LumiMask

    def getRunRange(self):
        '''
        https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGoodLumiSectionsJSONFile
        https://raw.githubusercontent.com/attikis/cmssw/l1t-integration-CMSSW_10_3_1/FWCore/PythonUtilities/python/LumiList.py
        '''
        if self.getIsMC():
            msg = "Cannot get run-range for dataset %s (self.getIsMC() = %s)" % (self.getName() + es, self.getIsMC())
            raise Exception(es + msg + ns)

        import LumiList as LumiList
        l = LumiList.LumiList(filename = self.getLumiMaskFile())

        # Extract the information
        runs = l.getRuns()
        runRange = "%s-%s" % (runs[0], runs[-1])
        return runRange

    def printInfo(self):
        table = []
        if self.getIsMC():
            table = self.getInfoMC()
        else:
            table = self.getInfoData()

        for i, line in enumerate(table, 1):
            self.Print(line, i==1)
        return

    def getInfoMC(self):
        msgAlign = "{:^30} {:^10} {:^25} {:^40} {:^5} {:^5} {:^15} {:^10} {:^15} {:^10} {:^10} {:^10}"
        header   = msgAlign.format("Campaign", "p-String", "PU", "Global Tag", "Ext", "Ver", "DataTier", "DBS", "DataVersion", "#Events", "#Files", "#Sites")
        hLine    =  "="*200
        table    = []
        table.append("{:^200}".format(self.getName()))
        table.append(hLine)
        table.append(header)
        table.append(hLine)
        table.append(msgAlign.format(self.getCampaign(), self.getProcessString(), self.getPileup(), self.getGlobalTag(), 
                                     self.getExtension(), self.getVersion(), self.getDataTier(), self.getDBS(), 
                                     self.getDataVersion(), self.getNumOfEvents(), self.getNumOfLumis(),
                                     self.getNumOfFiles(), self.getNumOfSites()))
        return table

    def getInfoData(self):
        msgAlign = "{:^30} {:^15} {:^15} {:^10} {:^10} {:^5} {:^5} {:^15} {:^15} {:^10} {:^10} {:^10}"
        header   = msgAlign.format("<Run><Year><Era>", "Run-Range", "p-String", "Global Tag", "Ext", "Ver", "DataTier", "DBS", "DataVersion", "#Events", "#Lumis", "#Sites")
        hLine    =  "="*190
        table    = []
        table.append("{:^190}".format(self.getName()))
        table.append(hLine)
        table.append(header)
        table.append(hLine)
        table.append(msgAlign.format(self.getCampaign(), self.getRunRange(), self.getProcessString(), self.getGlobalTag(), 
                                     self.getExtension(), self.getVersion(), self.getDataTier(), self.getDBS(),
                                     self.getDataVersion(), self.getNumOfEvents(), self.getNumOfLumis(),
                                     self.getNumOfFiles(), self.getNumOfSites()))
        return table

#================================================================================================ 
# Definitions
#================================================================================================ 
datasets_2016 = {}
datasets_2017 = {}
datasets_2018 = {}

dsetGroups_2016 = {}
dsetGroups_2017 = {}
dsetGroups_2018 = {}

#================================================================================================ 
# [MC] RunIISummer16MiniAODv2 campaign
#================================================================================================ 
campaign = "RunIISummer16MiniAODv2"
pString  = ""
PU       = "PUMoriond17"
gTag     = "80X_mcRun2_asymptotic_2016_TrancheIV_v6"
dataTier = "MINIAODSIM"
dbs      = "global"
dataVer  = "80Xmc"

datasets_2016["HplusToTB"] = []
# For-loop: Initial samples (1.5M / sample except M650 which had 10M)
for mass in [180, 200, 220, 250, 300, 350, 400, 500, 650, 800, 1000, 1500, 2000, 2500, 3000, 5000, 7000, 10000]:
    pdName = "ChargedHiggs_HplusTB_HplusToTB_M-%d_13TeV_amcatnlo_pythia8" % (mass)
    dset   = Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc")
    datasets_2016["HplusToTB"].append(dset)

# For-loop: Additional signal samples (10M / sample)
for mass in [180, 200, 220, 250, 300, 350, 400, 500, 800, 1000, 1500, 2000, 2500, 3000]:
    pdName = "ChargedHiggs_HplusTB_HplusToTB_M-%d_13TeV_amcatnlo_pythia8" % (mass)
    ext    = 1
    ver    = 1
    if mass == 250 or mass == 500:
        ver = 2
    dset   = Dataset(pdName, campaign, pString, PU, gTag, ext, ver, "MINIAODSIM", "global", "80Xmc")
    datasets_2016["HplusToTB"].append(dset)

datasets_2016["HplusToTauNu"] = []
# Heavy
for mass in [180, 200, 220, 250, 300, 400, 500, 750, 1000, 2000, 3000]:
    pdName = "ChargedHiggs_HplusTB_HplusToTauNu_M-%s_13TeV_amcatnlo_pythia8" % (mass)
    ext    = 0
    ver    = 1
    if mass == 500:
        ext = 1
    if mass > 500:
        pdName = "ChargedHiggs_HplusTB_HplusToTauNu_HeavyMass_M-%d_13TeV_amcatnlo_pythia8" % (mass)
    dset = Dataset(pdName, campaign, pString, PU, gTag, ext, ver, "MINIAODSIM", "global", "80Xmc")
    datasets_2016["HplusToTauNu"].append(dset)

# Intermediate (No Neutrals)
for mass in [145, 150, 155, 160, 165, 170, 175, 180, 190, 200]:
    pdName = "ChargedHiggs_HplusTB_HplusToTauNu_IntermediateMassNoNeutral_M-%d_13TeV_amcatnlo_pythia8" % (mass)
    dset   = Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc")
    datasets_2016["HplusToTauNu"].append(dset)

# Intermediate (With Neutrals)
versions = [3, 5, 5, 5, 6, 3, 5, 5]
for i, mass in enumerate( [145, 150, 155, 165, 170, 175, 180, 200], 0):
    pdName = "ChargedHiggs_HplusTB_HplusToTauNu_IntermediateMassWithNeutral_M-%d_13TeV_amcatnlo_pythia8" % (mass)
    ver    = versions[i]
    dset   = Dataset(pdName, campaign, pString, PU, gTag, 0, ver, "MINIAODSIM", "global", "80Xmc")
    datasets_2016["HplusToTauNu"].append(dset)

# Light
for mass in [80, 90, 100, 120, 140, 150, 155, 160]:
    pdName = "ChargedHiggs_TTToHplusBWB_HplusToTauNu_M-%d_13TeV_amcatnlo_pythia8" % (mass)
    dset   = Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc")
    datasets_2016["HplusToTauNu"].append(dset)

# Private
dset = Dataset("HplusToTauNu_M_200_TuneCUETP8M1_tauola_13TeV_pythia8", "amarini-GEN-SIM-71", "6edf4b210aa48b81088c0de44a7af6f5", None, None, None, None, "USER", "phys03", "80Xmc", None, False)
datasets_2016["HplusToTauNu"].append(dset)
for mass in [140, 160, 180, 220]:
    pdName = "ChargedHToTauNu_M%d_13TeV_pythia8" % (mass)
    dset   = Dataset(pdName, "amarini-PUMoriond17-MINIAODSIM-v1", "28028af67189b3de7224b79195bd0e1d",  None, None, None, None, "USER", "phys03", "80Xmc")
    datasets_2016["HplusToTauNu"].append(dset)

for mass in [140, 150, 155, 160]:
    pdName = "ChargedHToTauNu_lowmass_M%d_13TeV_pythia8" % (mass)
    dset   = Dataset(pdName, "amarini-PUMoriond17-MINIAODSIM-v1", "28028af67189b3de7224b79195bd0e1d",  None, None, None, None, "USER", "phys03", "80Xmc")
    datasets_2016["HplusToTauNu"].append(dset)


datasets_2016["HplusToHW"] = []
#datasets_2016["HplusToHW"].append( Dataset("CRAB_PrivateMC", "mlotti-Hplus2hw_2ta_PAT_m350_f", "71a58a62b6d71fe0f402eda209b9b80f", None, None, None, None, "USER", "phys03", "80Xmc", None, False) )
#datasets_2016["HplusToHW"].append( Dataset("CRAB_PrivateMC", "mlotti-Hplus2hw_2ta_PAT_mhp350_mh150", "71a58a62b6d71fe0f402eda209b9b80f", None, None, None, None, "USER", "phys03", "80Xmc", None, False) )
#datasets_2016["HplusToHW"].append( Dataset("CRAB_PrivateMC", "mlotti-CRAB3_Hplus_PAT", "71a58a62b6d71fe0f402eda209b9b80f", None, None, None, None, "USER", "phys03", "80Xmc", None, False) )
#datasets_2016["HplusToHW"].append( Dataset("CRAB_PrivateMC", "mlotti-Hplus2hw_4l_PAT_m350_f", "71a58a62b6d71fe0f402eda209b9b80f", None, None, None, None, "USER", "phys03", "80Xmc", None, False) )
#datasets_2016["HplusToHW"].append( Dataset("CRAB_PrivateMC", "mlotti-Hplus2hw_4l_PAT_m1500_f", "71a58a62b6d71fe0f402eda209b9b80f", None, None, None, None, "USER", "phys03", "80Xmc", None, False) )
datasets_2016["HplusToHW"].append( Dataset("CRAB_private_ChargedHiggs_HplusTB_HplusToHW_M300_mH200_2ta_NLO", "gkole-Hplus2hw_2ta_PAT_m300_mh200", "5514f561afe549de22a98231d58a7b06", None, None, None, None, "USER", "phys03", "80Xmc", None, False) )
datasets_2016["HplusToHW"].append( Dataset("CRAB_private_ChargedHiggs_HplusTB_HplusToHW_M700_mH200_2ta_NLO", "gkole-Hplus2hw_2ta_PAT_m700_mh200", "5514f561afe549de22a98231d58a7b06", None, None, None, None, "USER", "phys03", "80Xmc", None, False) )
#datasets_2016["HplusToHW"].append( Dataset("CRAB_private_ChargedHiggs_HplusTB_HplusToHW_M300_mH200_2ta_NLO", "mlotti-Hplus2hw_2ta_w2all_NLO_mhp300_mh200_PAT", "b023a05347da3c2a313d219710ff54e9", None, None, None, None, "USER", "phys03", "80Xmc", None, False) )
#datasets_2016["HplusToHW"].append( Dataset("CRAB_private_ChargedHiggs_HplusTB_HplusToHW_M700_mH200_2ta_NLO", "mlotti-Hplus2hw_2ta_w2all_NLO_mhp700_mh200_PAT", "b023a05347da3c2a313d219710ff54e9", None, None, None, None, "USER", "phys03", "80Xmc", None, False) )


# Single Top
datasets_2016["SingleTop"] = []
pdName = "ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1"
dset   = Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False)
datasets_2016["SingleTop"].append(dset)

pdName = "ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1"
dset   = Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc")
datasets_2016["SingleTop"].append(dset)

pdName = "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4"
dset   = Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc")
datasets_2016["SingleTop"].append(dset)

pdName = "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1"
dset   = Dataset(pdName, campaign, pString, PU, gTag, 1, 1, "MINIAODSIM", "global", "80Xmc")
datasets_2016["SingleTop"].append(dset)

pdName = "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1"
dset   = Dataset(pdName, campaign, pString, PU, gTag, 1, 1, "MINIAODSIM", "global", "80Xmc")
datasets_2016["SingleTop"].append(dset)

pdName = "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4"
dset   = Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc")
datasets_2016["SingleTop"].append(dset)

pdName = "ST_s-channel_4f_InclusiveDecays_13TeV-amcatnlo-pythia8"
dset   = Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc")
datasets_2016["SingleTop"].append(dset)

# TT
datasets_2016["TT"] = []
pdName = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8"
dset   = Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False)
datasets_2016["TT"].append(dset)

# TT with variations 
datasets_2016["TopTagSys"] = []
pdName = "TT_TuneCUETP8M2T4_13TeV-powheg-fsrdown-pythia8"
for i in [0, 1, 2]:
    dset = Dataset(pdName, campaign, pString, PU, gTag, i, 1, "MINIAODSIM", "global", "80Xmc", None, False)
    datasets_2016["TopTagSys"].append(dset)

pdName = "TT_TuneCUETP8M2T4_13TeV-powheg-fsrup-pythia8"
for i in [0, 1, 2]:
    dset = Dataset(pdName, campaign, pString, PU, gTag, i, 1, "MINIAODSIM", "global", "80Xmc", None, False)
    datasets_2016["TopTagSys"].append(dset)

pdName = "TT_TuneCUETP8M2T4_13TeV-powheg-isrdown-pythia8"
for i in [0, 1, 2]:
    dset = Dataset(pdName, campaign, pString, PU, gTag, i, 1, "MINIAODSIM", "global", "80Xmc", None, False)
    datasets_2016["TopTagSys"].append(dset)

pdName = "TT_TuneCUETP8M2T4_13TeV-powheg-isrup-pythia8"
for i in [1, 2]:
    dset   = Dataset(pdName, campaign, pString, PU, gTag, i, 1, "MINIAODSIM", "global", "80Xmc", None, False)
    datasets_2016["TopTagSys"].append(dset)

pdName = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8-evtgen"
dset   = Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False)
datasets_2016["TopTagSys"].append(dset)

pdName = "TT_TuneCUETP8M2T4_erdON_13TeV-powheg-pythia8"
dset   = Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False)
datasets_2016["TopTagSys"].append(dset)

for mt in [1665, 1695, 1715, 1735, 1755, 1785]:
    pdName = "TT_TuneCUETP8M2T4_mtop%d_13TeV-powheg-pythia8" % (mt)
    if mt == 1695:
        dset = Dataset(pdName, campaign, pString, PU, gTag, 1, 1, "MINIAODSIM", "global", "80Xmc", None, False)
        dset = Dataset(pdName, campaign, pString, PU, gTag, 2, 1, "MINIAODSIM", "global", "80Xmc", None, False)
    elif mt == 1755:
        dset = Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False)
        dset = Dataset(pdName, campaign, pString, PU, gTag, 1, 1, "MINIAODSIM", "global", "80Xmc", None, False)
        dset = Dataset(pdName, campaign, pString, PU, gTag, 2, 1, "MINIAODSIM", "global", "80Xmc", None, False)
    else:
        dset = Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False)
    datasets_2016["TopTagSys"].append(dset)

pdName = "TT_TuneEE5C_13TeV-powheg-herwigpp"
for i in [0, 2, 3]:
    dset = Dataset(pdName, campaign, pString, PU, gTag, i, 1, "MINIAODSIM", "global", "80Xmc", None, False)
    datasets_2016["TopTagSys"].append(dset)

pdName = "TT_TuneCUETP8M2T4_erdON_13TeV-powheg-pythia8"
dset = Dataset(pdName, campaign, pString, PU, gTag, 1, 1, "MINIAODSIM", "global", "80Xmc", None, False)
datasets_2016["TopTagSys"].append(dset)

pdName = "TT_hdampDOWN_TuneCUETP8M2T4_13TeV-powheg-pythia8"
for i in [0, 1]:
    dset = Dataset(pdName, campaign, pString, PU, gTag, i, 1, "MINIAODSIM", "global", "80Xmc", None, False)
    datasets_2016["TopTagSys"].append(dset)

pdName = "TT_hdampUP_TuneCUETP8M2T4_13TeV-powheg-pythia8"
for i in [0, 1]:
    dset = Dataset(pdName, campaign, pString, PU, gTag, i, 1, "MINIAODSIM", "global", "80Xmc", None, False)
    datasets_2016["TopTagSys"].append(dset)


# TTJets (-ve weights)
datasets_2016["TTJets_HT"] = []
for htRange in ["600to800", "800to1200", "1200to2500", "2500toInf"]:
    pdName = "TTJets_HT-%s_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" % (htRange)
    datasets_2016["TTJets_HT"].append( Dataset(pdName, campaign, pString, PU, gTag, 1, 1, "MINIAODSIM", "global", "80Xmc", None, False) )


# TTTT 
datasets_2016["TTTT"] = []
datasets_2016["TTTT"].append( Dataset("TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8"  , campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["TTTT"].append( Dataset("TTTT_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )


# TTW
datasets_2016["TTWJetsToLNu"] = []
pdName = "TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8"
datasets_2016["TTWJetsToLNu"].append( Dataset(pdName, campaign, pString, PU, gTag, 1, 3, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["TTWJetsToLNu"].append( Dataset(pdName, campaign, pString, PU, gTag, 2, 1, "MINIAODSIM", "global", "80Xmc", None, False) )

datasets_2016["TTWJetsToQQ"] = []
pdName = "TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8"                                      
datasets_2016["TTWJetsToQQ"].append( Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )


# TTZ
datasets_2016["TTZToLLNuNu"] = []
pdName = "TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8"
datasets_2016["TTZToLLNuNu"].append(Dataset(pdName, campaign, pString, PU, gTag, 3, 1, "MINIAODSIM", "global", "80Xmc", None, False) )

datasets_2016["TTZToQQ"] = []
pdName = "TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"
datasets_2016["TTZToQQ"].append( Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )


# ttH
datasets_2016["ttHJet"] = []
datasets_2016["ttHJet"].append( Dataset("ttHJetTobb_M125_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 3, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet"].append( Dataset("ttHJetToTT_M125_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 4, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet"].append( Dataset("ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix", campaign, pString, PU, gTag, 1, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet"].append( Dataset("ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_v2", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )

datasets_2016["ttHJet_Syst"] = []
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToTT_M130_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToTT_M120_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M95_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M90_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M85_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M80_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M75_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M70_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M65_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M60_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M130_13TeV_amcatnloFXFX_madspin_pythia8_v2", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M130_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M127_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M126_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_UpPS", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_DownPS", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_CUETP8M1Up", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_CUETP8M1Down", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M124_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M123_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M120_13TeV_amcatnloFXFX_madspin_pythia8_v2", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M120_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M110_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M105_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["ttHJet_Syst"].append( Dataset("ttHJetToGG_M100_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )


# Zprime
datasets_2016["Zprime"] = []
for mass in ["500", "1000", "3000"]:
    pdName = "ZprimeToTauTau_M-%s_TuneCUETP8M1_13TeV-pythia8-tauola" % (mass)
    dset   = Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) 
    datasets_2016["Zprime"].append(dset)


# QCD
datasets_2016["QCD"] = []
for ptRange in ["15to30", "30to50", "50to80", "80to120", "120to170",
                "170to300", "300to470", "470to600", "600to800", "800to1000", 
                 "1000to1400", "1400to1800", "1800to2400", "2400to3200", "3200toInf"]:
    
    pdName = "QCD_Pt_%s_TuneCUETP8M1_13TeV_pythia8" % (ptRange)
    ver = 1
    if ptRange == "3200toInf":
        ver = 3
    dset = Dataset(pdName, campaign, pString, PU, gTag, 0, ver, "MINIAODSIM", "global", "80Xmc", None, False)
    datasets_2016["QCD"].append(dset)

    if ptRange == "80to120":
        dset = Dataset(pdName, campaign, pString, PU, gTag, 2, 1, "MINIAODSIM", "global", "80Xmc", None, False)
        datasets_2016["QCD"].append(dset)
    if ptRange in ["120to170", "170to300", "300to470", "600to800", "800to1000", "1400to1800", "1800to2400", "2400to3200"]:
        dset = Dataset(pdName, campaign, pString, PU, gTag, 1, 1, "MINIAODSIM", "global", "80Xmc", None, False)
        datasets_2016["QCD"].append(dset)

# QCD (mu-enriched)
datasets_2016["QCD_MuEnriched"] = []
for ptRange in ["15to20", "20to30", "30to50", "50to80", "80to120",
                "120to170", "170to300", "300to470", "470to600", 
                "600to800", "800to1000", "1000toInf"]:
    
    pdName = "QCD_Pt-%s_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8" % (ptRange)
    ver = 1
    dset = Dataset(pdName, campaign, pString, PU, gTag, 0, ver, "MINIAODSIM", "global", "80Xmc", None, False)
    datasets_2016["QCD_MuEnriched"].append(dset)

    if ptRange in ["80to120", "170to300", "300to470", "470to600", "600to800", "800to1000"]:
        if ptRange == "80to120":
            ver = 3
        dset = Dataset(pdName, campaign, pString, PU, gTag, 1, ver, "MINIAODSIM", "global", "80Xmc", None, False)
        datasets_2016["QCD_MuEnriched"].append(dset)

    if ptRange in ["300to470", "470to600", "800to1000"]:
        dset = Dataset(pdName, campaign, pString, PU, gTag, 2, ver, "MINIAODSIM", "global", "80Xmc", None, False)
        datasets_2016["QCD_MuEnriched"].append(dset)

# QCD (b-enriched)
datasets_2016["QCD_bEnriched"] = []
for htRange in ["200to300", "300to500", "500to700", "700to1000", "1000to1500", "1500to2000", "2000toInf"]:
    
    pdName = "QCD_bEnriched_HT%s_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" % (htRange)
    ver = 1
    dset = Dataset(pdName, campaign, pString, PU, gTag, 0, ver, "MINIAODSIM", "global", "80Xmc", None, False)
    datasets_2016["QCD_bEnriched"].append(dset)

# QCD (HT binned)
datasets_2016["QCD_HT"] = []
for htRange in ["50to100", "100to200", "200to300", "300to500", "500to700", "700to1000", "1000to1500", "1500to2000", "2000toInf"]:
    
    pdName = "QCD_HT%s_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" % (htRange)
    dset = Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False)
    datasets_2016["QCD_HT"].append(dset)

    if htRange in ["200to300", "300to500", "500to700", "700to1000", "1000to1500", "1500to2000", "2000toInf"]:
        if htRange == "500to700":
            ver = 2
        else:
            ver = 1
        dset = Dataset(pdName, campaign, pString, PU, gTag, 1, ver, "MINIAODSIM", "global", "80Xmc", None, False)
        datasets_2016["QCD_HT"].append(dset)
    dset   = None


# Drell-Yan
datasets_2016["DYJetsToLL"] = []
pdName = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
dset = Dataset(pdName, campaign, pString, PU, gTag, 1, 2, "MINIAODSIM", "global", "80Xmc", None, False)
datasets_2016["DYJetsToLL"].append(dset)

pdName = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
dset = Dataset(pdName, campaign, pString, PU, gTag, 2, 1, "MINIAODSIM", "global", "80Xmc", None, False)
datasets_2016["DYJetsToLL"].append(dset)

datasets_2016["DYJetsToLL_HT"] = []
for htRange in ["70to100", "100to200", "200to400", "400to600", "600to800", "800to1200", "1200to2500", "2500toInf"]:
    pdName = "DYJetsToLL_M-50_HT-%s_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" % (htRange)
    ver = 1
    if htRange == "600to800":
        ver = 2
    dset = Dataset(pdName, campaign, pString, PU, gTag, 0, ver, "MINIAODSIM", "global", "80Xmc", None, False)
    datasets_2016["DYJetsToLL_HT"].append(dset)
    if htRange in ["100to200", "200to400", "400to600"]:
        dset = Dataset(pdName, campaign, pString, PU, gTag, 1, 1, "MINIAODSIM", "global", "80Xmc", None, False)
        datasets_2016["DYJetsToLL_HT"].append(dset)

datasets_2016["DYJetsToQQ"] = []
pdName = "DYJetsToQQ_HT180_13TeV-madgraphMLM-pythia8"
dset   = Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) 
datasets_2016["DYJetsToQQ"].append(dset)


# Z+jets
datasets_2016["ZJetsToQQ"] = []
dset = Dataset("ZJetsToQQ_HT600toInf_13TeV-madgraph", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False)
datasets_2016["ZJetsToQQ"].append(dset)


# W+jets
datasets_2016["WJetsToLNu"] = []
datasets_2016["WJetsToLNu"].append(Dataset("WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False))
datasets_2016["WJetsToLNu"].append(Dataset("WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", campaign, pString, PU, gTag, 2, 1, "MINIAODSIM", "global", "80Xmc", None, False))

datasets_2016["WJetsToLNu_HT"] = []
for htRange in ["70To100", "100To200", "200To400", "400To600", "600To800", "800To1200", "1200To2500", "2500ToInf"]:
    pdName = "WJetsToLNu_HT-%s_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" % (htRange)
    dset = Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False)

    datasets_2016["WJetsToLNu_HT"].append(dset)
    if htRange in ["100To200", "200To400", "400To600", "600To800", "800To1200", "1200To2500", "2500ToInf"]:
        dset = Dataset(pdName, campaign, pString, PU, gTag, 1, 1, "MINIAODSIM", "global", "80Xmc", None, False)
        datasets_2016["WJetsToLNu_HT"].append(dset)
    if htRange in ["100To200", "200To400"]:
        dset = Dataset(pdName, campaign, pString, PU, gTag, 2, 1, "MINIAODSIM", "global", "80Xmc", None, False)
        datasets_2016["WJetsToLNu_HT"].append(dset)


# WJetsToQQ
datasets_2016["WJetsToQQ"] = []
pdName = "WJetsToQQ_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
datasets_2016["WJetsToQQ"].append( Dataset(pdName, campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )


# Diboson
datasets_2016["Diboson"] = []
datasets_2016["Diboson"].append( Dataset("WZ_TuneCUETP8M1_13TeV-pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["Diboson"].append( Dataset("WZ_TuneCUETP8M1_13TeV-pythia8", campaign, pString, PU, gTag, 1, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["Diboson"].append( Dataset("ZZ_TuneCUETP8M1_13TeV-pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["Diboson"].append( Dataset("ZZ_TuneCUETP8M1_13TeV-pythia8", campaign, pString, PU, gTag, 1, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["Diboson"].append( Dataset("WW_TuneCUETP8M1_13TeV-pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["Diboson"].append( Dataset("WW_TuneCUETP8M1_13TeV-pythia8", campaign, pString, PU, gTag, 1, 1, "MINIAODSIM", "global", "80Xmc", None, False) )

datasets_2016["DibosonTo4Q"] = []
datasets_2016["DibosonTo4Q"].append( Dataset("WWTo4Q_13TeV-powheg", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )
datasets_2016["DibosonTo4Q"].append( Dataset("ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )

datasets_2016["DibosonTo2L2Nu"] = []
datasets_2016["DibosonTo2L2Nu"].append (Dataset("WWTo2L2Nu_13TeV-powheg", campaign, pString, PU, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )

datasets_2016["DibosonLNuQQ"] = []
datasets_2016["DibosonLNuQQ"].append (Dataset("WWToLNuQQ_13TeV-powheg", campaign, pString, PU, gTag, 1, 1, "MINIAODSIM", "global", "80Xmc", None, False) )


# Neutrino
datasets_2016["SingleNeutrino"] = []
datasets_2016["SingleNeutrino"].append( Dataset("SingleNeutrino", "RunIIHighPUTrainsMiniAODv2", "HighPUTrains", None, gTag, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False, False) )


#================================================================================================ 
# MC, Trigger Development
#================================================================================================ 
datasets_2016["HplusToTauNu_TrgDev"] = []
cTrgDev  = "PhaseIFall16MiniAOD"
pTrgDev  = ""
puTrgDev = "FlatPU28to62HcalNZSRAW_PhaseIFall16"
gTrgDev  = "81X_upgrade2017_realistic_v26"
for mass in [80, 160]:
    pdName = "ChargedHiggs_TTToHplusBWB_HplusToTauNu_M-%d_13TeV_amcatnlo_pythia8" % (mass) 
    datasets_2016["HplusToTauNu_TrgDev"].append( Dataset(pdName, cTrgDev, pTrgDev, puTrgDev, gTrgDev, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )

for mass in [200, 500, 1000]:
    pdName = "ChargedHiggs_HplusTB_HplusToTauNu_M-%d_13TeV_amcatnlo_pythia8" % (mass) 
    datasets_2016["HplusToTauNu_TrgDev"].append( Dataset(pdName, cTrgDev, pTrgDev, puTrgDev, gTrgDev, 0, 1, "MINIAODSIM", "global", "80Xmc", None, False) )


#================================================================================================ 
# [Data] 03Feb2017
#================================================================================================ 
# Tau Primary Dataset
datasets_2016["Tau-03Feb2017"] = []
datasets_2016["Tau-03Feb2017"].append( Dataset("Tau", "Run2016B-03Feb2017_ver2",  None, None, None, 0, 2, "MINIAOD", "global", "80Xdata", "Cert_273150-275376_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016B.txt", False) )
datasets_2016["Tau-03Feb2017"].append( Dataset("Tau", "Run2016C-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_275656-276283_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016C.txt", False) )
datasets_2016["Tau-03Feb2017"].append( Dataset("Tau", "Run2016D-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_276315-276437_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016D_MET80.txt", False) )
datasets_2016["Tau-03Feb2017"].append( Dataset("Tau", "Run2016D-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_276453-276811_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016D_noMET80.txt", False) )
datasets_2016["Tau-03Feb2017"].append( Dataset("Tau", "Run2016E-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_276831-277420_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016E.txt", False) )
datasets_2016["Tau-03Feb2017"].append( Dataset("Tau", "Run2016F-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_277932-278800_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016F_HIP.txt", False) )
datasets_2016["Tau-03Feb2017"].append( Dataset("Tau", "Run2016F-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_278801-278808_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016F_HIPfixed.txt", False) )
datasets_2016["Tau-03Feb2017"].append( Dataset("Tau", "Run2016G-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_278820-280385_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016G.txt", False) )
datasets_2016["Tau-03Feb2017"].append( Dataset("Tau", "Run2016H-03Feb2017_ver2",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_281613-284035_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016H.txt", False) )
datasets_2016["Tau-03Feb2017"].append( Dataset("Tau", "Run2016H-03Feb2017_ver3",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_284036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016H.txt", False) )

datasets_2016["Tau-18Apr2017"] = []
datasets_2016["Tau-18Apr2017"].append( Dataset("Tau", "Run2016B-18Apr2017_ver2",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata18Apr", "Cert_273150-275376_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016B.txt", False) )
datasets_2016["Tau-18Apr2017"].append( Dataset("Tau", "Run2016C-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata18Apr", "Cert_275656-276283_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016C.txt", False) )
datasets_2016["Tau-18Apr2017"].append( Dataset("Tau", "Run2016D-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata18Apr", "Cert_276315-276437_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016D_MET80.txt", False) )
datasets_2016["Tau-18Apr2017"].append( Dataset("Tau", "Run2016D-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata18Apr", "Cert_276453-276811_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016D_noMET80.txt", False) )
datasets_2016["Tau-18Apr2017"].append( Dataset("Tau", "Run2016E-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata18Apr", "Cert_276831-277420_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016E.txt", False) )
datasets_2016["Tau-18Apr2017"].append( Dataset("Tau", "Run2016F-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata18Apr", "Cert_277932-278800_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016F_HIP.txt", False) )
datasets_2016["Tau-18Apr2017"].append( Dataset("Tau", "Run2016F-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata18Apr", "Cert_278801-278808_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016F_HIPfixed.txt", False) )
datasets_2016["Tau-18Apr2017"].append( Dataset("Tau", "Run2016G-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata18Apr", "Cert_278820-280385_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016G.txt", False) )
datasets_2016["Tau-18Apr2017"].append( Dataset("Tau", "Run2016H-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata18Apr", "Cert_281613-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016H.txt", False) )

datasets_2016["Tau-07Aug2017"] = []
datasets_2016["Tau-07Aug2017"].append( Dataset("Tau", "Run2016B-07Aug17_ver1",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata07Aug", "Cert_273150-275376_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016B.txt", False) )
datasets_2016["Tau-07Aug2017"].append( Dataset("Tau", "Run2016C-07Aug17",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata07Aug", "Cert_275656-276283_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016C.txt", False) )
datasets_2016["Tau-07Aug2017"].append( Dataset("Tau", "Run2016E-07Aug17",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata07Aug", "Cert_276831-277420_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016E.txt", False) )


# JetHT Primary Dataset
datasets_2016["JetHT-03Feb2017"] = []
datasets_2016["JetHT-03Feb2017"].append( Dataset("JetHT", "Run2016B-03Feb2017_ver2",  None, None, None, 0, 2, "MINIAOD", "global", "80Xdata", "Cert_273150-275376_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016B.txt", False) )
datasets_2016["JetHT-03Feb2017"].append( Dataset("JetHT", "Run2016C-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_275656-276283_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016C.txt", False) )
datasets_2016["JetHT-03Feb2017"].append( Dataset("JetHT", "Run2016D-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_276315-276811_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016D.txt", False) )
datasets_2016["JetHT-03Feb2017"].append( Dataset("JetHT", "Run2016E-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_276831-277420_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016E.txt", False) )
datasets_2016["JetHT-03Feb2017"].append( Dataset("JetHT", "Run2016F-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_277932-278800_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016F_HIP.txt", False) )
datasets_2016["JetHT-03Feb2017"].append( Dataset("JetHT", "Run2016F-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_278801-278808_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016F_HIPfixed.txt", False) )
datasets_2016["JetHT-03Feb2017"].append( Dataset("JetHT", "Run2016G-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_278820-280385_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016G.txt", False) )
datasets_2016["JetHT-03Feb2017"].append( Dataset("JetHT", "Run2016H-03Feb2017_ver2",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_281613-284035_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016H.txt", False) )
datasets_2016["JetHT-03Feb2017"].append( Dataset("JetHT", "Run2016H-03Feb2017_ver3",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_284036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016H.txt", False) )

datasets_2016["JetHT-18Apr2017"] = []
datasets_2016["JetHT-18Apr2017"].append( Dataset("JetHT", "Run2016B-18Apr2017_ver2",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_273150-275376_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016B.txt", False) )
datasets_2016["JetHT-18Apr2017"].append( Dataset("JetHT", "Run2016C-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_275656-276283_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016C.txt", False) )
datasets_2016["JetHT-18Apr2017"].append( Dataset("JetHT", "Run2016D-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_276315-276811_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016D.txt", False) )
datasets_2016["JetHT-18Apr2017"].append( Dataset("JetHT", "Run2016E-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_276831-277420_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016E.txt", False) )
datasets_2016["JetHT-18Apr2017"].append( Dataset("JetHT", "Run2016F-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_277932-278800_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016F_HIP.txt", False) )
datasets_2016["JetHT-18Apr2017"].append( Dataset("JetHT", "Run2016F-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_278801-278808_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016F_HIPfixed.txt", False) )
datasets_2016["JetHT-18Apr2017"].append( Dataset("JetHT", "Run2016G-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_278820-280385_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016G.txt", False) )
datasets_2016["JetHT-18Apr2017"].append( Dataset("JetHT", "Run2016H-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_281613-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016H.txt", False) )

datasets_2016["SingleMuon-03Feb2017"] = []
datasets_2016["SingleMuon-03Feb2017"].append( Dataset("SingleMuon", "Run2016B-03Feb2017_ver2",  None, None, None, 0, 2, "MINIAOD", "global", "80Xdata", "Cert_273150-275376_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016B.txt", False) )
datasets_2016["SingleMuon-03Feb2017"].append( Dataset("SingleMuon", "Run2016C-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_275656-276283_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016C.txt", False) )
datasets_2016["SingleMuon-03Feb2017"].append( Dataset("SingleMuon", "Run2016D-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_276315-276811_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016D.txt", False) )
datasets_2016["SingleMuon-03Feb2017"].append( Dataset("SingleMuon", "Run2016D-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_276315-276811_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016D.txt", False) )
datasets_2016["SingleMuon-03Feb2017"].append( Dataset("SingleMuon", "Run2016E-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_276831-277420_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016E.txt", False) )
datasets_2016["SingleMuon-03Feb2017"].append( Dataset("SingleMuon", "Run2016F-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_277932-278800_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016F_HIP.txt", False) )
datasets_2016["SingleMuon-03Feb2017"].append( Dataset("SingleMuon", "Run2016F-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_278801-278808_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016F_HIPfixed.txt", False) )
datasets_2016["SingleMuon-03Feb2017"].append( Dataset("SingleMuon", "Run2016G-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_278820-280385_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016G.txt", False) )
datasets_2016["SingleMuon-03Feb2017"].append( Dataset("SingleMuon", "Run2016H-03Feb2017_ver2",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_281613-284035_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016H.txt", False) )
datasets_2016["SingleMuon-03Feb2017"].append( Dataset("SingleMuon", "Run2016H-03Feb2017_ver3",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_284036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016H.txt", False) )


datasets_2016["SingleMuon-18Apr2017"] = []
datasets_2016["SingleMuon-18Apr2017"].append( Dataset("SingleMuon", "Run2016B-18Apr2017_ver2",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_273150-275376_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016B.txt", False) )
datasets_2016["SingleMuon-18Apr2017"].append( Dataset("SingleMuon", "Run2016C-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_275656-276283_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016C.txt", False) )
datasets_2016["SingleMuon-18Apr2017"].append( Dataset("SingleMuon", "Run2016D-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_276315-276811_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016D.txt", False) )
datasets_2016["SingleMuon-18Apr2017"].append( Dataset("SingleMuon", "Run2016E-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_276831-277420_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016E.txt", False) )
datasets_2016["SingleMuon-18Apr2017"].append( Dataset("SingleMuon", "Run2016F-18Apr2017",  None, None, None, 0, 2, "MINIAOD", "global", "80Xdata", "Cert_277932-278800_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016F_HIP.txt", False) )
datasets_2016["SingleMuon-18Apr2017"].append( Dataset("SingleMuon", "Run2016F-18Apr2017",  None, None, None, 0, 2, "MINIAOD", "global", "80Xdata", "Cert_278801-278808_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016F_HIPfixed.txt", False) )
datasets_2016["SingleMuon-18Apr2017"].append( Dataset("SingleMuon", "Run2016G-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_278820-280385_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016G.txt", False) )
datasets_2016["SingleMuon-18Apr2017"].append( Dataset("SingleMuon", "Run2016H-18Apr2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_281613-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016H.txt", False) )


datasets_2016["SingleElectron-03Feb2017"] = []
datasets_2016["SingleElectron-03Feb2017"].append( Dataset("SingleElectron", "Run2016B-03Feb2017_ver2",  None, None, None, 0, 2, "MINIAOD", "global", "80Xdata", "Cert_273150-275376_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016B.txt", False) )
datasets_2016["SingleElectron-03Feb2017"].append( Dataset("SingleElectron", "Run2016C-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_275656-276283_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016C.txt", False) )
datasets_2016["SingleElectron-03Feb2017"].append( Dataset("SingleElectron", "Run2016D-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_276315-276811_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016D.txt", False) )
datasets_2016["SingleElectron-03Feb2017"].append( Dataset("SingleElectron", "Run2016E-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_276831-277420_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016E.txt", False) )
datasets_2016["SingleElectron-03Feb2017"].append( Dataset("SingleElectron", "Run2016F-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_277932-278800_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016F_HIP.txt", False) )
datasets_2016["SingleElectron-03Feb2017"].append( Dataset("SingleElectron", "Run2016F-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_278801-278808_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016F_HIPfixed.txt", False) )
datasets_2016["SingleElectron-03Feb2017"].append( Dataset("SingleElectron", "Run2016G-03Feb2017",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_278820-280385_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016G.txt", False) )
datasets_2016["SingleElectron-03Feb2017"].append( Dataset("SingleElectron", "Run2016H-03Feb2017_ver2",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_281613-284035_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016H.txt", False) )
datasets_2016["SingleElectron-03Feb2017"].append( Dataset("SingleElectron", "Run2016H-03Feb2017_ver3",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", "Cert_284036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_Run2016H.txt", False ) )

datasets_2016["ZeroBias-PromptReco"] = []
datasets_2016["ZeroBias-PromptReco"].append( Dataset("ZeroBias", "Run2016G-PromptReco",  None, None, None, 0, 1, "MINIAOD", "global", "80Xdata", None, False) )
datasets_2016["ZeroBias-PromptReco"].append( Dataset("ZeroBias", "Run2016H-PromptReco",  None, None, None, 0, 2, "MINIAOD", "global", "80Xdata2016H", None, False) )
datasets_2016["ZeroBias-PromptReco"].append( Dataset("ZeroBias", "Run2016H-PromptReco",  None, None, None, 0, 3, "MINIAOD", "global", "80Xdata2016H", None, False) )

#================================================================================================ 
# Dataset Grouping
#================================================================================================ 
dsetGroups_2016["TauLeg"] = []
dsetGroups_2016["TauLeg"].extend(datasets_2016["SingleMuon-03Feb2017"])
dsetGroups_2016["TauLeg"].extend(datasets_2016["DYJetsToLL"])
dsetGroups_2016["TauLeg"].extend(datasets_2016["DYJetsToLL_HT"])
dsetGroups_2016["TauLeg"].extend(datasets_2016["Zprime"])
dsetGroups_2016["TauLeg"].extend(datasets_2016["TT"])
dsetGroups_2016["TauLeg"].extend(datasets_2016["SingleTop"])
dsetGroups_2016["TauLeg"].extend(datasets_2016["WJetsToLNu"] )
dsetGroups_2016["TauLeg"].extend(datasets_2016["QCD_MuEnriched"])
# dsetGroups_2016["TauLeg"].extend(datasetsH125)


dsetGroups_2016["METLeg"] = []
dsetGroups_2016["METLeg"].extend(datasets_2016["Tau-03Feb2017"])
dsetGroups_2016["METLeg"].extend(datasets_2016["DYJetsToLL"])
dsetGroups_2016["METLeg"].extend(datasets_2016["DYJetsToLL_HT"])
dsetGroups_2016["METLeg"].extend(datasets_2016["TT"])
dsetGroups_2016["METLeg"].extend(datasets_2016["SingleTop"])
dsetGroups_2016["METLeg"].extend(datasets_2016["WJetsToLNu"] )
dsetGroups_2016["METLeg"].extend(datasets_2016["QCD"])

dsetGroups_2016["L1Study"] = []
dsetGroups_2016["L1Study"].extend(datasets_2016["ZeroBias-PromptReco"])
dsetGroups_2016["L1Study"].extend(datasets_2016["SingleNeutrino"])
dsetGroups_2016["L1Study"].extend(datasets_2016["QCD"])

dsetGroups_2016["HToTauNu"] = []
dsetGroups_2016["HToTauNu"].extend(datasets_2016["Tau-03Feb2017"])
#dsetGroups_2016["HToTauNu"].extend(datasets_2016["Tau-18Apr2017"])
dsetGroups_2016["HToTauNu"].extend(datasets_2016["TT"])
dsetGroups_2016["HToTauNu"].extend(datasets_2016["SingleTop"])
dsetGroups_2016["HToTauNu"].extend(datasets_2016["WJetsToLNu"] )  
dsetGroups_2016["HToTauNu"].extend(datasets_2016["Diboson"])
dsetGroups_2016["HToTauNu"].extend(datasets_2016["QCD_HT"])
dsetGroups_2016["HToTauNu"].extend(datasets_2016["HplusToTauNu"])
#dsetGroups_2016["HToTauNu"].extend(datasets_2016["HplusToTauNu_TrgDev"])

dsetGroups_2016["HToTB"] = []
dsetGroups_2016["HToTB"].extend(datasets_2016["JetHT-03Feb2017"])
dsetGroups_2016["HToTB"].extend(datasets_2016["HplusToTB"])
dsetGroups_2016["HToTB"].extend(datasets_2016["TT"])
dsetGroups_2016["HToTB"].extend(datasets_2016["SingleTop"])
dsetGroups_2016["HToTB"].extend(datasets_2016["DYJetsToQQ"])
dsetGroups_2016["HToTB"].extend(datasets_2016["WJetsToQQ"])
dsetGroups_2016["HToTB"].extend(datasets_2016["ZJetsToQQ"])
dsetGroups_2016["HToTB"].extend(datasets_2016["DibosonTo4Q"])
dsetGroups_2016["HToTB"].extend(datasets_2016["QCD_HT"])
dsetGroups_2016["HToTB"].extend(datasets_2016["QCD_bEnriched"])
dsetGroups_2016["HToTB"].extend(datasets_2016["TTWJetsToQQ"])
dsetGroups_2016["HToTB"].extend(datasets_2016["TTTT"]) 
dsetGroups_2016["HToTB"].extend(datasets_2016["TTZToQQ"])

dsetGroups_2016["HToHW"] = []
dsetGroups_2016["HToHW"].extend(datasets_2016["HplusToHW"])
dsetGroups_2016["HToHW"].extend(datasets_2016["SingleMuon-03Feb2017"])
dsetGroups_2016["HToHW"].extend(datasets_2016["SingleElectron-03Feb2017"])
dsetGroups_2016["HToHW"].extend(datasets_2016["TT"])
dsetGroups_2016["HToHW"].extend(datasets_2016["SingleTop"])
#dsetGroups_2016["HToHW"].extend(datasets_2016["DYJetsToLL"])
dsetGroups_2016["HToHW"].extend(datasets_2016["DYJetsToLL_HT"])
dsetGroups_2016["HToHW"].extend(datasets_2016["WJetsToLNu"] ) #for doing the HT binned 0To70
dsetGroups_2016["HToHW"].extend(datasets_2016["WJetsToLNu_HT"] )
#dsetGroups_2016["HToHW"].extend(datasets_2016["ZJetsToQQ"])
dsetGroups_2016["HToHW"].extend(datasets_2016["Diboson"])
#dsetGroups_2016["HToHW"].extend(datasets_2016["QCD"])
dsetGroups_2016["HToHW"].extend(datasets_2016["QCD_MuEnriched"])
#dsetGroups_2016["HToHW"].extend(datasetsQCDbEnriched) 
dsetGroups_2016["HToHW"].extend(datasets_2016["TTWJetsToLNu"])
dsetGroups_2016["HToHW"].extend(datasets_2016["TTTT"])
dsetGroups_2016["HToHW"].extend(datasets_2016["TTZToLLNuNu"])
dsetGroups_2016["HToHW"].extend(datasets_2016["TTZToQQ"])
dsetGroups_2016["HToHW"].extend(datasets_2016["ttHJet"])

dsetGroups_2016["HToHWTrgEff"] = []
dsetGroups_2016["HToHWTrgEff"].extend(datasets_2016["HplusToHW"])

dsetGroups_2016["HToHW_withTop"] = []
dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["HplusToHW"])
dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["SingleMuon-03Feb2017"])
dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["SingleElectron-03Feb2017"])
dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["TT"])
dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["SingleTop"])
dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["DYJetsToLL"])
dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["DYJetsToLL_HT"])
dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["WJetsToLNu"] ) #for doing the HT binned 0To70
dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["WJetsToLNu_HT"] )
dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["Diboson"])
#dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["QCD"])
dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["QCD_MuEnriched"])
#dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["QCD_bEnriched"])
dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["TTWJetsToLNu"])
dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["TTTT"])
dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["TTZToLLNuNu"])
dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["TTZToQQ"])
dsetGroups_2016["HToHW_withTop"].extend(datasets_2016["ttHJet"])

dsetGroups_2016["JetTriggers"] = []
dsetGroups_2016["JetTriggers"].extend(datasets_2016["SingleMuon-03Feb2017"])
dsetGroups_2016["JetTriggers"].extend(datasets_2016["TT"])
dsetGroups_2016["JetTriggers"].extend(datasets_2016["QCD_MuEnriched"])

dsetGroups_2016["TopTagSys"] = []
dsetGroups_2016["TopTagSys"].extend(datasets_2016["JetHT-03Feb2017"])
dsetGroups_2016["TopTagSys"].extend(datasets_2016["TT"])
dsetGroups_2016["TopTagSys"].extend(datasets_2016["TopTagSys"])


#================================================================================================ 
# Testing (python datasets.py)
#================================================================================================ 
#DatasetGroup("SignalAnalysis").GetDatasetList()
#DatasetGroup("HToTauNu").GetDatasetList()

#datasets_2016["TT"][0].printSummary()
#datasets_2016["TT"][0].printInfo()

#for i, d in enumerate(datasets_2016["JetHT-03Feb2017"], 0):
#    print d.printInfo()

# datasets_2016["JetHT-03Feb2017"][0].printInfo()
# datasets_2016["JetHT-03Feb2017"][1].printInfo()
# datasets_2016["JetHT-03Feb2017"][2].printInfo()
# datasets_2016["JetHT-03Feb2017"][3].printInfo()
# datasets_2016["JetHT-03Feb2017"][4].printInfo()
