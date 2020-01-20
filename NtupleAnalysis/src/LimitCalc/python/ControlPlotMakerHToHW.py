'''
## \package ControlPlotMaker

Classes for making control plots
'''
#================================================================================================ 
# Imports
#================================================================================================ 
from HiggsAnalysis.LimitCalc.DatacardColumn import DatacardColumn
import HiggsAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles
from HiggsAnalysis.NtupleAnalysis.tools.dataset import Count,RootHistoWithUncertainties
import HiggsAnalysis.NtupleAnalysis.tools.aux as aux

import HiggsAnalysis.NtupleAnalysis.tools.histograms as histograms
import HiggsAnalysis.NtupleAnalysis.tools.plots as plots
import HiggsAnalysis.NtupleAnalysis.tools.counter as counter
import HiggsAnalysis.NtupleAnalysis.tools.tdrstyle as tdrstyle
import HiggsAnalysis.NtupleAnalysis.tools.styles as styles

from math import pow,sqrt,log10
import os
import sys
import ROOT

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
# Definitions
#================================================================================================ 
_legendLabelQCDMC = "QCD (MC)" 
_legendLabelEWKMC = "EWK (MC)" 
_legendLabelQCDdata = "Mis-ID. #tau_{h} (data)"

drawPlot = plots.PlotDrawer(ratio             = True, 
                            ratioYlabel       = "Data / Bkg  ",
                            ratioCreateLegend = True,
                            ratioType         = "errorScale",
                            ratioErrorOptions = {"numeratorStatSyst": False},
                            stackMCHistograms = True, 
                            addMCUncertainty  = True, 
                            addLuminosityText = True,
                            cmsTextPosition   = "outframe")

drawPlotPaper = plots.PlotDrawer(ratio             = True, 
                                 ratioYlabel       = "Data / Bkg  ",
                                 ratioCreateLegend = True,
                                 ratioType         = "errorScalePaper",
                                 ratioErrorOptions = {"numeratorStatSyst": False},
                                 stackMCHistograms = True, 
                                 addMCUncertainty  = True, 
                                 addLuminosityText = True,
                                 cmsTextPosition   = "outframe")

drawPlot2D = plots.PlotDrawer(opts2={"ymin": 0.5, "ymax": 1.5},
                              ratio=False, 
                              # ratioYlabel="Data/Bkg.", 
                              # ratioCreateLegend=True,
                              # ratioType="errorScale", 
                              # ratioErrorOptions={"numeratorStatSyst": False},
                              # stackMCHistograms=True, addMCUncertainty=True, 
                              addLuminosityText=True,
                              cmsTextPosition="outframe")


#================================================================================================ 
# Definitions
#================================================================================================ 
#class ControlPlotMakerHToHW:
class ControlPlotMaker:
    def __init__(self, opts, config, dirname, luminosity, observation, datasetGroups, verbose=False):
        self._validateDatacard(config)
        self._config  = config
        self._verbose = verbose
        self._opts    = opts
        self._dirname = dirname
        self._luminosity    = luminosity
        self._observation   = observation
        self._datasetGroups = datasetGroups
        
        # Define label options
        myStyle = tdrstyle.TDRStyle()
        myStyle.setOptStat(False)
        if self._config.OptionPaper:
            plots._legendLabels["MCStatError"] = "Stat unc"
            plots._legendLabels["MCStatSystError"] = "Stat #oplus syst unc"
            plots._legendLabels["BackgroundStatError"] = "Stat #oplus syst unc"
            plots._legendLabels["BackgroundStatSystError"] = "Stat #oplus syst unc"
        else:
            plots._legendLabels["MCStatError"] = "Bkg. stat."
            plots._legendLabels["MCStatSystError"] = "Bkg. stat.#oplussyst."
            plots._legendLabels["BackgroundStatError"] = "Bkg. stat. unc"
            plots._legendLabels["BackgroundStatSystError"] = "Bkg. stat.#oplussyst. unc."
            
        # Make control plots
        self.Verbose(sh_h+ "Generating control plots" + sh_n, True)

        # Definitions
        massPoints = []
        massPoints.extend(self._config.MassPoints)
        if self._config.OptionDoWithoutSignal:
            massPoints.append(-1) # for plotting with no signal
        nMasses = len(massPoints)
        nPlots  = len(self._config.ControlPlots)
        counter = 0

        # For-loop: All mass points
        for m in massPoints:

            # Initialize flow plot            
            selectionFlow = SelectionFlowPlotMaker(self._opts, self._config, m)

            # For-loop: All control plots
            for i in range(0, nPlots):
                counter += 1                

                # Skip if control plot does not exist
                if observation.getControlPlotByIndex(i) == None:
                    continue
                
                # Get the control plot
                myCtrlPlot = self._config.ControlPlots[i]

                # The case m < 0 is for plotting hitograms without any signal
                if m > 0:
                    saveName = "%s/%s_M%d_%s" % (self._dirname, self._config.OptionPlotNamePrefix, m, myCtrlPlot.title) #xenios - iro
                    msg      = "Control Plot %d/%d: mH=%s GeV"  % (counter, nMasses*nPlots, str(m))
                else:
                    saveName = "%s/%s_%s" % (self._dirname, self._config.OptionPlotNamePrefix, myCtrlPlot.title)
                    msg      = "Control Plot %d/%d: no signal"  % (counter, nMasses*nPlots)

                # Initialize histograms
                hData     = None
                hSignal   = None
                hQCDdata    = None
                hQCDMC    = None
                stackList = []

                # Inform user of progress
                self.PrintFlushed(sh_l + msg + sh_n, counter==1)
                if counter == len(massPoints)*nPlots:
                    print

                # For-loop: All dataset columns (to find histograms)
                for i, c in enumerate(self._datasetGroups, 1):
                    msg = "DatasetGroup %d/%d) %s" % (i, len(self._datasetGroups), c.getLabel())
                    if 0:
                        self.PrintFlushed(sh_l + msg + sh_n, False)

                    msg  = "Dataset is %s for plot %s. " % (sh_h + c.getLabel() + sh_n, sh_a + myCtrlPlot.histoName + sh_n)
                    msg += "The title is \"%s\" (Serves as dictionary key to get binning from systematics.py)" % (myCtrlPlot.title)
                    self.Verbose(msg, True)

                    # Skip plot (NOTE! DatasetGroup is signal then all other mass points should be disabled. This concerns bValidMass)
                    bValidMass   = (m < 0 or c.isActiveForMass(m, self._config)) 
                    bEmptyCol    = (c.typeIsEmptyColumn())
                    #bInvalidPlot = (c.getControlPlotByIndex(i) == None) # broken (for unknown reason)
                    bInvalidPlot = (c.getControlPlotByTitle(myCtrlPlot.histoName, myKey="shape") == None) 
                    bSkipPlot    = (not bValidMass or bEmptyCol or bInvalidPlot)

                    if bSkipPlot:
                        self.Verbose("Skipping plot %s (%s): bValidMass=%s, bEmptyCol=%s, bInvalidPlot=%s" % (sh_h + c.getLabel() + sh_n, myCtrlPlot.histoName, bValidMass, bEmptyCol, bInvalidPlot), False)
                        continue

                    # Print available control plots?
                    if 0:
                        c.PrintControlPlotTitles(myKey="shape")

                    # Clone histo
                    h = c.getControlPlotByTitle(myCtrlPlot.histoName, myKey="shape").Clone()
                    #h = c.getControlPlotByIndex(i)["shape"].Clone() # broken (for unknown reason)
                    
                    if bSkipPlot:
                        self.Print("Skipping plot %s (bValidMass=%s, bEmptyCol=%s, bInvalidPlot=%s)" % (sh_h + c.getLabel() + sh_n, bValidMass, bEmptyCol, bInvalidPlot), False)
                        continue

                    msg = "Histogram %s is of type %s" % (sh_h + c.getLabel() + sh_n, c.type())
                    if c.typeIsSignal():
                        msg += ". Signal"
                        self.Verbose(msg, True)

                        if self._config.OptionBr != 1.0:
                            msg = "Scaling %s histograms with BR=%.2f" % (sh_h + c.getLabel() + sh_n, self._config.OptionBr)
                            self.Print(msg, False)
                        h.Scale(self._config.OptionBr)

                        if hSignal == None:
                            hSignal = h.Clone()
                        else:
                            hSignal.Add(h)
                    elif c.typeIsQCDinverted():
                        msg += ". QCDdata"
                        self.Verbose(msg, True)
                        if hQCDdata == None:
                            hQCDdata = h.Clone()
                        else:
                            hQCDdata.Add(h)
                    elif c.typeIsQCDMC():
                        msg += ". QCDMC"
                        self.Verbose(msg, True)
                        if hQCDMC== None:
                            hQCDMC = h.Clone()
                        else:
                            hQCDMC.Add(h)
                    elif c.typeIsEWKMC() or c.typeIsGenuineB():
                        msg += ". EWKMC"
                        self.Verbose(msg, False)
                        myHisto = histograms.Histo(h, c._datasetMgrColumn)
                        myHisto.setIsDataMC(isData=False, isMC=True)
                        stackList.append(myHisto)
                    else:
                        msg = "Histogram %s has unexpected type" % (c.getLabel())
                        raise Exception(sh_e + msg + sh_n)
                    
                # If stackList is empty skip the next steps
                if len(stackList) < 1:
                    continue

                # Stack all the histograms
                if hQCDdata != None:
                    self.Verbose("Stacking QCDdata", True)
                    myHisto = histograms.Histo(hQCDdata, "QCDdata", legendLabel=_legendLabelQCDdata)
                    myHisto.setIsDataMC(isData=False, isMC=True)
                    stackList.insert(0, myHisto)
                elif hQCDMC != None:
                    self.Verbose("Stacking QCDMC", True)
                    myHisto = histograms.Histo(hQCDMC,"QCDMC",legendLabel=_legendLabelQCDMC)
                    myHisto.setIsDataMC(isData=False, isMC=True)
                    stackList.insert(0, myHisto)
                else:
                    self.Verbose("The stack list has %d datasets" % len(stackList), True)
                    for i, s in enumerate(stackList, 1):
                        self.Verbose(s.getName(), False)                        

                    msg = "This should never be reached (type=%s)" % (c.type())
                    #raise Exception(sh_e + msg + sh_n)

                #hData = observation.getControlPlotByIndex(i)["shape"].Clone() # broken (for unknown reason)
                hData = observation.getControlPlotByTitle(myCtrlPlot.histoName, myKey="shape").Clone()
                hDataUnblinded = hData.Clone()

                # Apply blinding & Get blinding string
                myBlindingString = self._applyBlinding(myCtrlPlot, stackList, hData, hSignal)

                # Data
                myDataHisto = histograms.Histo(hData,"Data")
                myDataHisto.setIsDataMC(isData=True, isMC=False)
                stackList.insert(0, myDataHisto)
                
                # Add signal
                if m > 0:
                    try:
                        signalLabel = self._config.SignalName % (m)
                    except AttributeError:
                        msg  = "The imported file %s has no attribute \"SignalName\". "
                        msg += "Make sure you define this parameter in the imported module (e.g. SignalName = ChargedHiggs_HplusTB_HplusToHW_M\%s_mH200_2ta_NLO)"
                        raise AttributeError(msg, True)
                    
                    self.Verbose("Signal label is %s" % (sh_h + signalLabel + sh_n), True)
                    myHisto = histograms.Histo(hSignal, signalLabel)
                    myHisto.setIsDataMC(isData=False, isMC=True)
                    stackList.insert(1, myHisto)
            
                self.Verbose("Add data to selection flow plot", True)
                selectionFlow.addColumn(myCtrlPlot.flowPlotCaption, hDataUnblinded, stackList[1:])

                # Make plot
                stackPlot = None
                myParams  = myCtrlPlot.details.copy()
                self.Verbose("Creating DataMCPlot2 object using stackList of length %d: %s" % ( len(stackList), ", ".join([s.getName() for s in stackList])), True)
                stackPlot = plots.DataMCPlot2(stackList)

                self.Verbose("Setting stack plot values to: luminosity = %s, energy = %s" % ( self._luminosity, self._config.OptionSqrtS), True)
                stackPlot.setLuminosity(self._luminosity)
                stackPlot.setEnergy("%d" % self._config.OptionSqrtS)
                stackPlot.setDefaultStyles(paperStyle=self._config.OptionPaper)
                try:
                    stackPlot.setLegendHeader(self._config.PlotLegendHeader)
                except AttributeError:
                    msg  = "The imported file %s has no attribute \"PlotLegendHeader\". "
                    msg += "Make sure you define this parameter in the imported module (e.g. LegendHeader = Associated Production)"
                    raise AttributeError(msg, True)
                    
                # Tweak paramaters
                if not "unit" in myParams.keys():
                    myParams["unit"] = ""

                if myParams["unit"] != "":
                    myParams["xlabel"] = "%s (%s)" % (myParams["xlabel"], myParams["unit"])

                # Apply various settings to my parameters
                self._setBlingingString(myBlindingString, myParams)
                self._setYlabelWidthSuffix(hData, myParams)
                self._setLegendPosition(myParams)
                self._setRatioLegendPosition(myParams)

                # Remove non-dientified keywords
                del myParams["unit"]

                # Ratio axis
                if not "opts2" in myParams.keys():
                    myParams["opts2"] = {"ymin": 0.3, "ymax": 1.7}

                # Make sure BR is indicated if anyting else but BR=1.0
                if m > 0 and self._config.OptionBr != 1.0:
                    stackPlot.histoMgr.setHistoLegendLabelMany({
                            #mySignalLabel : "H^{+} m_{H^{+}}=%d GeV (x %s)" % (m, self._config.OptionBr)
                            signalLabel : "m_{H^{+}}=%d GeV (x %s)" % (m, self._config.OptionBr)
                            })

                # Do plotting
                if self._config.OptionPaper:
                    drawPlotPaper(stackPlot, saveName, **myParams)
                else:
                    drawPlot(stackPlot, saveName, **myParams)
                    
                
            # Do selection flow plot
            selectionFlow.makePlot(self._dirname, m, len(self._config.ControlPlots), self._luminosity)
        return
    
    def GetFName(self):
        fName = __file__.split("/")[-1]
        fName = fName.replace(".pyc", ".py")
        return fName

    def _validateDatacard(self, config):
        if config.ControlPlots == None:
            return

        if config.OptionSqrtS == None:
            raise Exception(ShellStyles.ErrorLabel()+"Please set the parameter OptionSqrtS = <integer_value_in_TeV> in the config file!"+sh_n)
        return

    def _setBlingingString(self, myBlindingString, myParams):
        if myBlindingString != None:
            if myParams["unit"] != "" and myParams["unit"][0] == "^":
                myParams["blindingRangeString"] = "%s%s" % (myBlindingString, myParams["unit"])
            else:
                myParams["blindingRangeString"] = "%s %s" % (myBlindingString, myParams["unit"])
        return

    def _setYlabelWidthSuffix(self, histo, myParams):

        ylabelBinInfo = True
        if "ylabelBinInfo" in myParams:
            ylabelBinInfo = myParams["ylabelBinInfo"]
            del myParams["ylabelBinInfo"]

        if ylabelBinInfo:
            minBinWidth, maxBinWidth = self._getMinMaxBinWidth(histo)

            widthSuffix = ""
            minBinWidthString = "%d" % minBinWidth
            maxBinWidthString = "%d" % maxBinWidth
            
            if minBinWidth < 1.0:
                myFormat = "%%.%df" % (abs(int(log10(minBinWidth)))+1)
                minBinWidthString = myFormat % minBinWidth
            if maxBinWidth < 1.0:
                myFormat = "%%.%df" % (abs(int(log10(maxBinWidth)))+1)
                maxBinWidthString = myFormat % maxBinWidth

            widthSuffix = "%s-%s" % (minBinWidthString, maxBinWidthString)
            if abs(minBinWidth-maxBinWidth) < 0.001:
                widthSuffix = "%s" % (minBinWidthString)
            if (myParams["unit"] == "" and widthSuffix == "1"):
                return
            elif "< Events / " in myParams["ylabel"]: # e.g. < Events / GeV >
                return                
            else:
                myParams["ylabel"] = "%s / %s %s" % (myParams["ylabel"], widthSuffix, myParams["unit"])
            return

    def _getMinMaxBinWidth(self, histo):

        minWidth = 10000.0
        maxWidth = 0.0
        
        # For-loop: All bins
        for j in range(1, histo.getRootHisto().GetNbinsX()+1):
            w = histo.getRootHisto().GetBinWidth(j)
            if w < minWidth:
                minWidth = w
            if w > maxWidth:
                maxWidth = w
        return minWidth, maxWidth

    def _setRatioLegendPosition(self, myParams):
        if "ratioLegendPosition" in myParams.keys():
            dy = 0.03
            if self._config.OptionPaper:
                dh = -0.12
            else:
                dh = 0.0
            self.Verbose("Setting the legend (ratio) position to \"%s\"" % ( myParams["ratioLegendPosition"] ) )
            if myParams["ratioLegendPosition"] == "left":
                myParams["ratioMoveLegend"] = {"dx": -0.51, "dy": dy, "dh": dh}
            elif myParams["ratioLegendPosition"] == "right":
                myParams["ratioMoveLegend"] = {"dx": -0.01, "dy": dy, "dh": dh}
            elif myParams["ratioLegendPosition"] == "SE":
                myParams["ratioMoveLegend"] = {"dx": -0.08, "dy": -0.33, "dh": dh}
            else:
                raise Exception("Unknown value for option ratioLegendPosition: %s!", myParams["ratioLegendPosition"])
            del myParams["ratioLegendPosition"]
        else:
            if not "ratioMoveLegend" in myParams:
                myParams["ratioMoveLegend"] = {"dx": -0.51, "dy": 0.03} # default: left
        return

    def _setLegendPosition(self, myParams):
        if "legendPosition" in myParams.keys():
            self.Verbose("Setting the legend position to \"%s\"" % ( myParams["legendPosition"] ) )
            dh = 0.14
            if self._config.OptionBr == 1.0:
                dx = -0.10
            else:
                dx = -0.12
            if self._config.OptionPaper:
                dx = -0.20
                dh = +0.05
            if myParams["legendPosition"] == "NE":
                myParams["moveLegend"] = {"dx": dx, "dy": -0.02, "dh": dh}
            elif myParams["legendPosition"] == "SE":
                myParams["moveLegend"] = {"dx": dx, "dy": -0.40, "dh": dh}
            elif myParams["legendPosition"] == "SW":
                myParams["moveLegend"] = {"dx": -0.53, "dy": -0.40, "dh": dh}
            elif myParams["legendPosition"] == "NW":
                myParams["moveLegend"] = {"dx": -0.53, "dy": -0.02, "dh": dh}
            elif myParams["legendPosition"] == "RM":
                myParams["moveLegend"] = {"dx": +10.0, "dy": +10.0, "dh": -100.0}
            else:
                raise Exception("Unknown value for option legendPosition: %s!", myParams["legendPosition"])
            del myParams["legendPosition"]
        elif not "moveLegend" in myParams:
            myParams["moveLegend"] = {"dx": -0.10, "dy": -0.02} # default: NE
        return

    def Verbose(self, msg, printHeader=True):
        '''
        Calls Print() only if verbose options is set to true
        '''                                                                                                                                                                         
        if not self._verbose:
            return
        self.Print(msg, printHeader)
        return

    def Print(self, msg, printHeader=True):
        '''
        Simple print function. If verbose option is enabled prints, otherwise does nothing
        '''
        fName = __file__.split("/")[-1]
        fName = fName.replace(".pyc", ".py")
        if printHeader:
            print "=== ", fName
        print "\t", msg
        return

    def PrintFlushed(self, msg, printHeader=True):
        '''
        Useful when printing progress in a loop
        '''
        msg = "\r\t" + msg
        if printHeader:
            print "=== ", self.GetFName()
        sys.stdout.write(msg)
        sys.stdout.flush()
        return

    def _applyCustomBlinding(self, myObject, blindedRange = []):
        '''
        Blind observed (i.e. data) histogram bins for the given 
        range

        Returns blinding string 
        '''
        myMin = None
        myMax = None
        myHisto = myObject.getRootHisto()
        
        # For-loop: All histo bins
        for i in range (1, myHisto.GetNbinsX()+1):
            myUpEdge  = myHisto.GetXaxis().GetBinUpEdge(i)
            myLowEdge = myHisto.GetXaxis().GetBinLowEdge(i)

            # Blind if any edge of the current bin is inside the blinded range or if bin spans over the blinded range
            if ((myLowEdge >= blindedRange[0] and myLowEdge <= blindedRange[1]) or
                (myUpEdge >= blindedRange[0] and myUpEdge <= blindedRange[1]) or 
                (myLowEdge <= blindedRange[0] and myUpEdge >= blindedRange[1])):
                if myMin == None or myLowEdge < myMin:
                    myMin = myLowEdge
                if myMax == None or myUpEdge > myMax:
                    myMax = myUpEdge
                myHisto.SetBinContent(i, -1.0)
                myHisto.SetBinError(i, 0.0)

        if myMin == None:
            return None
        myMinFormat = "%"+"d"
        myMaxFormat = "%"+"d"
        if abs(myMin) < 1.0 and abs(myMin) > 0.00000001:
            myMinFormat = "%%.%df"%(abs(int(log10(myMin)))+1)
        if abs(myMax) < 1.0  and abs(myMax) > 0.00000001:
            myMaxFormat = "%%.%df"%(abs(int(log10(myMax)))+1)
        s = myMinFormat%myMin+"-"+myMaxFormat%myMax
        return s

    def _applyBlinding(self, myCtrlPlot, stackList, hData, hSignal):
        myBlindingString = None
        if self._config.BlindAnalysis:
            if len(myCtrlPlot.blindedRange) > 0:                        
                myBlindingString = self._applyCustomBlinding(hData, myCtrlPlot.blindedRange)

        if self._config.OptionBlindThreshold != None:
            # For-loop: All histogram bins
            for k in xrange(1, hData.GetNbinsX()+1):
                myExpValue = 0.0
                for item in stackList:
                    myExpValue += item.getRootHisto().GetBinContent(k)
                # Fixme: Loop over ALL mass points and do this check! Not just the current one
                if hSignal.getRootHisto().GetBinContent(k) >= myExpValue * self._config.OptionBlindThreshold:
                    hData.getRootHisto().SetBinContent(k, -1.0)
                    hData.getRootHisto().SetBinError(k, 0.0)
        return myBlindingString
                        
class SignalAreaEvaluator:
    def __init__(self):
        self._output = ""

    def addEntry(self,mass,title,evaluationRange,hSignal,hQCD,hEmbedded,hEWKfake):
        # Obtain event counts
        mySignal = self._evaluate(evaluationRange,hSignal)
        myQCD = self._evaluate(evaluationRange,hQCD)
        myEmbedded = self._evaluate(evaluationRange,hEmbedded)
        myEWKfake = self._evaluate(evaluationRange,hEWKfake)
        # Produce output
        myOutput = "%s, mass=%d, range=%d-%d\n"%(title,mass,evaluationRange[0],evaluationRange[1])
        myOutput += "  signal: %f +- %f\n"%(mySignal.value(),mySignal.uncertainty())
        myOutput += "  QCD: %f +- %f\n"%(myQCD.value(),myQCD.uncertainty())
        myOutput += "  EWKtau: %f +- %f\n"%(myEmbedded.value(),myEmbedded.uncertainty())
        myOutput += "  EWKfake: %f +- %f\n"%(myEWKfake.value(),myEWKfake.uncertainty())
        myExpected = Count(0.0, 0.0)
        myExpected.add(myQCD)
        myExpected.add(myEmbedded)
        myExpected.add(myEWKfake)
        myOutput += "  Total expected: %f +- %f\n"%(myExpected.value(),myExpected.uncertainty())
        mySignal.divide(myExpected)
        myOutput += "  signal/expected: %f +- %f\n"%(mySignal.value(),mySignal.uncertainty())
        myOutput += "\n"
        self._output += myOutput

    def save(self,dirname):
        myFilename = dirname+"/signalAreaEvaluation.txt"
        myFile = open(myFilename, "w")
        myFile.write(self._output)
        myFile.close()
        print sh_h + "Signal area evaluation written to: " + sh_n + myFilename
        self._output = ""

    def _evaluate(self,evaluationRange,h):
        myResult = 0.0
        myError = 0.0
        for i in range(1,h.GetNbinsX()+1):
            if (h.GetXaxis().GetBinLowEdge(i) >= evaluationRange[0] and
                h.GetXaxis().GetBinUpEdge(i) <= evaluationRange[1]):
                myResult += h.GetBinContent(i)
                myError += pow(h.GetBinError(i),2)
        return Count(myResult,sqrt(myError))

class SelectionFlowPlotMaker:
    def __init__(self, opts, config, mass, verbose=False):
        self._opts = opts
        self._config = config
        self._mass = mass
        self._verbose = verbose
        self._createEmptyHisto(mass)
        self._expectedList = []
        self._expectedLabelList = []
        self._expectedListSystUp = []
        self._expectedListSystDown = []
        self._data = None
        self._myCurrentColumn = 1
        self._pickStatus = False
        self._pickLabel = ""
        return

    def Verbose(self, msg, printHeader=True):
        '''
        Calls Print() only if verbose options is set to true
        '''                                                                                                                                                                         
        if not self._verbose:
            return
        Print(msg, printHeader)
        return

    def GetFName(self):
        fName = __file__.split("/")[-1]
        fName = fName.replace(".pyc", ".py")
        return fName

    def Print(self, msg, printHeader=True):
        '''
        Simple print function. If verbose option is enabled prints, otherwise does nothing
        '''
        if printHeader:
            print "=== ", self.GetFName()
        print "\t", msg
        return

    def _createEmptyHisto(self, mass):
        myBinList = []

        # For-loop: All control plots (1 bin for each control plot)
        for i, c in enumerate(self._config.ControlPlots, 1):
            if (c.flowPlotCaption != "" and c.flowPlotCaption != "final"):
                self.Verbose("Selection flow. Ignoring plot \"%s\"" % (c.title))
                continue
            else:
                self.Verbose("Appending selection flow bin #%d (%s)" % (i, c.title), i==1)
                myBinList.append(c.flowPlotCaption)

        # Make an empty frame
        nBins = len(myBinList)
        myPlotName = "SelectionFlow_%d" % mass
        self._hFrame = ROOT.TH1F(myPlotName, myPlotName, nBins, 0, nBins)

        # For-loop: All (selection flow) bins
        for i, b in enumerate(myBinList, 1):
            self.Verbose("Setting label for bin #%d to \"%s\"" % (i, b), i==1)
            self._hFrame.GetXaxis().SetBinLabel(i, b)
        return

    def delete(self):
        for h in self._expectedList:
            h.Delete()
        self._expectedList = None
        for h in self._expectedListSystUp:
            h.Delete()
        self._expectedListSystUp = None
        for h in self._expectedListSystDown:
            h.Delete()
        self._expectedListSystDown = None
        self._expectedLabelList = None
        self._data.Delete()

    def addColumn(self,label,data,expectedList):
        '''
        System to pick the correct input for correct label
        '''
        if label == "":
            return

        # Create histograms if necessary
        if self._data == None:
            self._createHistograms(data, expectedList)
            return

        # Add expected
        for i in range(0, len(expectedList)):
            myRate = expectedList[i].getRootHistoWithUncertainties().getRate()
            self._expectedList[i].SetBinContent(self._myCurrentColumn, myRate)
            self._expectedList[i].SetBinError(self._myCurrentColumn, expectedList[i].getRootHistoWithUncertainties().getRateStatUncertainty())
            uncertUp = 0.0
            uncertDown = 0.0
            if myRate > 0.0:
                (uncertUp,uncertDown) = expectedList[i].getRootHistoWithUncertainties().getRateSystUncertainty()
                self._expectedListSystUp[i].SetBinContent(self._myCurrentColumn, uncertUp/myRate)
                self._expectedListSystDown[i].SetBinContent(self._myCurrentColumn, -uncertDown/myRate)
            if self._opts.debugControlPlots:
                s = "debugControlPlots:,"
                s += "After "+self._pickLabel
                s += ","+expectedList[i].getName()
                s += ",%f"%myRate
                s += ",+-,%f,(stat.)"%self._expectedList[i].GetBinError(self._myCurrentColumn)
                s += ",+,%f"%uncertUp
                s += ",-,%f,(syst.)"%uncertDown
                print s

        # Add data
        if data != None:
            self._data.SetBinContent(self._myCurrentColumn, data.getRate())
            self._data.SetBinError(self._myCurrentColumn, data.getRateStatUncertainty())
        else:
            self._data.SetBinContent(self._myCurrentColumn, -1)
        self._myCurrentColumn += 1
        # Refresh pick status
        self._pickLabel = label
        return

    def _createHistograms(self, data, expectedList):
        for e in expectedList:
            self._expectedList.append(aux.Clone(self._hFrame))
            self._expectedList[len(self._expectedList)-1].Reset()
            self._expectedListSystUp.append(aux.Clone(self._hFrame))
            self._expectedListSystUp[len(self._expectedListSystUp)-1].Reset()
            self._expectedListSystDown.append(aux.Clone(self._hFrame))
            self._expectedListSystDown[len(self._expectedListSystDown)-1].Reset()
            self._expectedLabelList.append(e.name)
        self._data = aux.Clone(self._hFrame)
        self._data.Reset()
        return

    def makePlot(self, dirname, m, index, luminosity):
        if self._data == None:
            self.Verbose("Did not find any data. Skipping creation of selection flow histogram", True)
            return

        self.Print("creating data-driven selection flow histogram!", True)
        stackList = []
        # For-loop: All expected histos
        for i in range(0,len(self._expectedList)):
            myRHWU = RootHistoWithUncertainties(self._expectedList[i])
            myRHWU.addShapeUncertaintyRelative("syst", th1Plus=self._expectedListSystUp[i], th1Minus=self._expectedListSystDown[i])
            myRHWU.makeFlowBinsVisible()
            if self._expectedLabelList[i] == "QCD":
                myHisto = histograms.Histo(myRHWU, self._expectedLabelList[i], legendLabel=_legendLabelQCD)
            elif self._expectedLabelList[i] == "QCDdata":
                myHisto = histograms.Histo(myRHWU, self._expectedLabelList[i], legendLabel=_legendLabelQCDdata)
            elif self._expectedLabelList[i] == "EWKfakes":
                myHisto = histograms.Histo(myRHWU, self._expectedLabelList[i], legendLabel=_legendLabelEWKFakes)
            else:
                myHisto = histograms.Histo(myRHWU, self._expectedLabelList[i])
            myHisto.setIsDataMC(isData=False, isMC=True)
            stackList.append(myHisto)

        # Data
        myRHWU = RootHistoWithUncertainties(self._data)
        myRHWU.makeFlowBinsVisible()
        myHisto = histograms.Histo(myRHWU, "Data")
        myHisto.setIsDataMC(isData=True, isMC=False)
        stackList.insert(0, myHisto)

        # Make plot
        stackPlot = plots.DataMCPlot2(stackList)
        stackPlot.setLuminosity(luminosity)
        stackPlot.setEnergy("%d"%self._config.OptionSqrtS)
        stackPlot.setDefaultStyles()

        myParams = {}
        myParams["ylabel"] = "Events"
        myParams["log"] = True
        myParams["cmsTextPosition"] = "right"
        myParams["opts"] = {"ymin": 0.9}
        myParams["opts2"] = {"ymin": 0.3, "ymax":1.7}
        myParams["moveLegend"] = {"dx": -0.53, "dy": -0.52, "dh":0.05} # for data-driven
        myParams["ratioMoveLegend"] = {"dx": -0.51, "dy": 0.03}
        drawPlot(stackPlot, "%s/%s_M%d_%02d_SelectionFlow"%(dirname, self._config.OptionPlotNamePrefix, m, index), **myParams) #xenios
        return
