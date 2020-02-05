#!/usr/bin/env python
'''
DESCRIPTION:                                                                                                                                                                         
This script is used to develop a mass-decorrelated top-quark tagger using a Neural Network architecture.                                                                             
The model identifies hadronically decaying top-quarks using as input a list of discriminating variables of                                                                           
truth-matched (signal) or unmatched (background) trijet combinations.                                                                                                                
For the mass decorrelated tagger we develop an Adversarial Neural Network that consists of two models:                                                                               
 1. Classifer: The nominal model that takes as input the set of input variables and gives as output                                                                                  
               a value from 0 to 1 that characterizes an object as signal or background.                                                                                             
               Loss function: L_c                                                                                                                                                    
 2. Regressor: This model is introduced in order to prevent the classifier from learning the top-quark mass.                                                                         
               It takes as input the output of the classifier and its output is a prediction of the top-quark mass.                                                                  
               Loss function: L_r                                                                                                                                                    
Combining the two models, we get a new Neural Network that minimizes a joint Loss functions: L = L_c - lambda*L_r                                                                    
The factor lambda controls the mass independence and the performance of the classifier.                                                                                              
During a training epoch, both classifier and regressor are trained and after each iteration, the models use                                                                          
the updated weights.                                                                                                                                                                 
During the training, we expect an increase of the loss value of the classifier and the regressor, and at the same                                                                    
time the minimization of the joint loss function.                                                                                                                                    

EXAMPLE: 
./adversarial.py --activ_clf relu,relu,relu,sigmoid --neurons_clf 32,32,32,1 --epochs 500 --batchSize 32
./adversarial.py --activ_clf relu,relu,relu,sigmoid --neurons_clf 32,32,32,1 --epochs 5 --batchSize 50000 --entrystop 20000
./adversarial.py --activ_clf relu,relu,relu,sigmoid --neurons_clf 32,32,32,1 --activ_reg tanh,relu,relu,relu --neurons_reg 32,32,32,1 --epochs 2000 --entrystop 500000 -s pdf 
./adversarial.py --activ_clf relu,relu,relu,sigmoid --neurons_clf 32,32,32,1 --activ_reg tanh,relu,relu,relu --neurons_reg 32,32,32,1 --epochs 2000 --entrystop 500000 -s pdf --standardise Robust
LAST USED:                                                                                                                                                                           
./adversarial.py --activ_clf relu,relu,relu,sigmoid --neurons_clf 32,32,32,1 --activ_reg tanh,relu,relu,relu --neurons_reg 32,32,32,1 --epochs 2000 --entrystop 500000 -s pdf
'''
#================================================================================================ 
# Imports
#================================================================================================ 
# print "=== Importing KERAS"
import keras
from keras import backend as backend
import uproot
import numpy
from numpy import digitize
import pandas
import ROOT
import array
import os
import psutil
import math
import json
import random as rn

import tensorflow as tf
from keras.models import Sequential, Model
from keras.layers import Dense, Activation, Flatten, Input
from sklearn.utils.class_weight import compute_sample_weight
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.externals import joblib

import subprocess
import plot
import func
import tdrstyle
from jsonWriter import JsonWriter 

import sys
import time
from datetime import datetime 
from optparse import OptionParser
import getpass
import socket

#================================================================================================
# Variable definition
#================================================================================================
# https://en.wikipedia.org/wiki/ANSI_escape_code#Colors  
ss = "\033[92m"
ns = "\033[0;0m"
ts = "\033[1;34m"
hs = "\033[0;35m"   
ls = "\033[0;33m"
es = "\033[1;31m"
cs = "\033[0;44m\033[1;37m"

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' #Disable AVX/FMA Warning

#Run in single CPU: this ensures reproducible results!                                                                                                                               
config = tf.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1,
                        allow_soft_placement=True, device_count = {'CPU': 1})
session = tf.Session(config=config)
backend.set_session(session)

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

def PrintFlushed(msg, printHeader=True):
    '''
    Useful when printing progress in a loop
    '''
    msg = "\r\t" + msg
    ERASE_LINE = '\x1b[2K'
    if printHeader:
        print "=== aux.py"
    sys.stdout.write(ERASE_LINE)
    sys.stdout.write(msg)
    sys.stdout.flush()
    return

def Verbose(msg, printHeader=True, verbose=False):
    if not opts.verbose:
        return
    Print(msg, printHeader)
    return

def PrintNetworkSummary(opts):
    table    = []
    msgAlign = "{:^10} {:^10} {:>12} {:>10}"
    title    =  msgAlign.format("Layer #", "Neurons", "Activation", "Type")
    hLine    = "="*len(title)
    table.append(hLine)
    table.append(title)
    table.append(hLine)
    table.append("=== Classifier")
    for i, n in enumerate(opts.neurons_clf, 0): 
        layerType = "unknown"
        if i == 0:
            layerType = "input"
        elif i+1 == len(opts.neurons_clf):
            layerType = "output"
        else:
            layerType = "hidden"
        table.append( msgAlign.format(i+1, opts.neurons_clf[i], opts.activ_clf[i], layerType) )
    table.append("=== Regressor")
    for i, n in enumerate(opts.neurons_reg, 0): 
        layerType = "unknown"
        if i == 0:
            layerType = "input"
        elif i+1 == len(opts.neurons_reg):
            layerType = "output"
        else:
            layerType = "hidden"
        table.append( msgAlign.format(i+1, opts.neurons_reg[i], opts.activ_reg[i], layerType) )
    table.append("")

    Print("Will construct a DNN with the following architecture", True)
    for r in table:
        Print(r, False)
    return

def GetDataFramesRowsColumns(df_sig, df_bkg):
    nEntries_sig, nColumns_sig = GetDataFrameRowsColumns(df_sig)
    nEntries_bkg, nColumns_bkg = GetDataFrameRowsColumns(df_bkg)

    if nColumns_sig != nColumns_bkg:
        msg = "The number of columns for signal (%d variables) does not match the correspondigng number for background (%d variables)" % (nColumns_sig, nColumns_bkg)
        raise Exception(es + msg + ns)
    else:
        msg  = "The signal has a total %s%d entries%s for each variable." % (ts, nEntries_sig, ns)
        msg += "The background has %s %d entries%s for each variable. " % (hs, nEntries_bkg, ns)
        msg += "(%d variables in total)" % (nColumns_sig)
        Verbose(msg, True)

    return nEntries_sig, nColumns_sig, nEntries_bkg, nColumns_bkg

def GetDataFrameRowsColumns(df):
    '''
    https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html
    '''

    # Get the index (row labels) of the DataFrame object (df).  
    rows = df.index.values
    Verbose("Printing indices (i.e. rows) of DataFrame:\n%s" % (rows), True)

    # Get the data types in the DataFrame
    dtypes = df.dtypes
    Verbose("Printing dtypes of DataFrame:\n%s" % (dtypes), True)

    # Get the columns (labels) of the DataFrame.
    columns = df.columns.values 
    Verbose("Printing column labels of DataFrame:\n%s" % (columns), True)

    nRows    = rows.size    #len(rows)     # NOTE: This is the number of Entries for each TBranch (variable) in the TTree
    nColumns = columns.size #len(columns)  # NOTE: This is the number of TBranches (i.e. variables) in the TTree
    Verbose("DataFrame has %s%d rows%s and %s%d columns%s" % (ts, nRows, ns, hs, nColumns, ns), True)
    return nRows, nColumns

def GetModelWeightsAndBiases(inputList, neuronsList):
    '''
    https://keras.io/models/about-keras-models/

    returns a dictionary containing the configuration of the model. 
    The model can be reinstantiated from its config via:
    model = Model.from_config(config)
    '''
    nParamsT  = 0
    nBiasT    = 0 
    nWeightsT = 0
    nInputs   = len(inputList)

    for i, n in enumerate(neuronsList, 0):
        nParams = 0
        nBias   = neuronsList[i]
        if i == 0:
            nWeights = nInputs * neuronsList[i]
        else:
            nWeights = neuronsList[i-1] * neuronsList[i]
        nParams += nBias + nWeights

        nParamsT  += nParams
        nWeightsT += nWeights
        nBiasT += nBias
    return nParamsT, nWeightsT, nBiasT
        
def GetModelParams(model):
    '''
    https://stackoverflow.com/questions/45046525/how-can-i-get-the-number-of-trainable-parameters-of-a-model-in-keras
    '''
    trainable_count = int(
        numpy.sum([backend.count_params(p) for p in set(model.trainable_weights)]))

    non_trainable_count = int(
        numpy.sum([backend.count_params(p) for p in set(model.non_trainable_weights)]))
    
    total_count = trainable_count + non_trainable_count
    return total_count, trainable_count, non_trainable_count

def GetModelConfiguration(myModel, verbose):
    '''
    NOTE: For some older releases of Keras the
    Model.get_config() method returns a list instead of a dictionary.
    This is fixed in newer releases. See details here:
    https://github.com/keras-team/keras/pull/11133
    '''
    
    config = myModel.get_config()
    for i, c in enumerate(config, 1):
        Verbose( "%d) %s " % (i, c), True)
    return config


def GetKwargs(var, standardise=False):

    kwargs  = {
        "normalizeToOne": True,
        "xTitle" : "DNN output",
        "yTitle" : "a.u.",
        "xMin"   :  0.0,
        "xMax"   : +1.0,
        "nBins"  : 200,
        "log"    : True,
        }

    if "output" in var.lower() or var.lower() == "output":
        kwargs["normalizeToOne"] = True
        kwargs["nBins"]  = 50
        kwargs["xMin"]   =  0.0
        kwargs["xMax"]   =  1.0
        kwargs["yMin"]   =  1e-5
        kwargs["yMax"]   =  2.0
        kwargs["xTitle"] = "DNN output"
        kwargs["yTitle"] = "a.u."
        kwargs["log"]    = True
        if var.lower() == "output":
            kwargs["legHeader"]  = "all data"
        if var.lower() == "outputpred":
            kwargs["legHeader"]  = "all data"
            kwargs["legEntries"] = ["train", "test"]
        if var.lower() == "outputtrain":
            kwargs["legHeader"]  = "training data"
        if var.lower() == "outputtest":
            kwargs["legHeader"]  = "test data"
        return kwargs

    if "efficiency" in var:
        kwargs["normalizeToOne"] = False
        kwargs["xMin"]   = 0.0
        kwargs["xMax"]   = 1.0
        kwargs["nBins"]  = 200
        kwargs["xTitle"] = "DNN output"
        kwargs["yTitle"] = "efficiency" # (\varepsilon)"
        kwargs["yMin"]   = 0.0
        kwargs["log"]    = False

    if "roc" in var:
        kwargs["normalizeToOne"] = False
        kwargs["xMin"]   = 0.0
        kwargs["xMax"]   = 1.0
        kwargs["nBins"]  = 200
        kwargs["xTitle"] = "signal efficiency"
        kwargs["yTitle"] = "background efficiency"
        kwargs["yMin"]   = 1e-4
        kwargs["yMax"]   = 10
        kwargs["log"]    = True

    if var == "significance":
        kwargs["normalizeToOne"] = False
        kwargs["xMin"]   = 0.0
        kwargs["xMax"]   = 1.0
        kwargs["nBins"]  = 200
        kwargs["xTitle"] = "DNN output"
        kwargs["yTitle"] = "significance"
        kwargs["yMin"]   = 0.0
        kwargs["log"]    = False

    if var == "loss":
        kwargs["normalizeToOne"] = False
        kwargs["xMin"]   =  0.
        kwargs["xMax"]   = opts.epochs
        kwargs["nBins"]  = opts.epochs
        kwargs["xTitle"] = "epoch"
        kwargs["yTitle"] = "loss"
        kwargs["yMin"]   = 0.0
        kwargs["log"]    = False

    if var == "acc":
        kwargs["normalizeToOne"] = False
        kwargs["xMin"]   =  0.0
        kwargs["xMax"]   = opts.epochs
        kwargs["nBins"]  = opts.epochs
        kwargs["xTitle"] = "epoch"
        kwargs["yTitle"] = "accuracy"
        kwargs["yMin"]   = 0.0
        kwargs["log"]    = False

    if var == "TrijetMass":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   = 900 #900.0
        kwargs["nBins"]  = 180 #450
        kwargs["xTitle"] = "m_{t} [GeV]"
        kwargs["yMin"]   = 1e-3
        kwargs["yMax"]   = 1.0

    if standardise:
        kwargs["xMin"]  =  -5.0
        kwargs["xMax"]  =  +5.0
        kwargs["nBins"] = 500 #1000

    return kwargs

def GetTime(tStart):
    tFinish = time.time()
    dt      = int(tFinish) - int(tStart)
    days    = divmod(dt,86400)      # days
    hours   = divmod(days[1],3600)  # hours
    mins    = divmod(hours[1],60)   # minutes
    secs    = mins[1]               # seconds
    return days, hours, mins, secs

def SaveModelParameters(myModel, opts):
    nParams = myModel.count_params()
    total_count, trainable_count, non_trainable_count  = GetModelParams(myModel)
    nParams, nWeights, nBias = GetModelWeightsAndBiases(opts.inputList, opts.neurons_clf)
    opts.modelParams = total_count
    opts.modelParamsTrain = trainable_count
    opts.modelParamsNonTrainable = non_trainable_count
    opts.modelWeights = nWeights
    opts.modelBiases  = nBias

    ind  = "{:<6} {:<30}"
    msg  = "The model has a total of %s%d parameters%s (neuron weights and biases):" % (hs, opts.modelParams, ns)
    msg += "\n\t" + ind.format("%d" % opts.modelWeights, "Weights")
    msg += "\n\t" + ind.format("%d" % opts.modelBiases, "Biases")
    msg += "\n\t" + ind.format("%d" % opts.modelParamsTrain, "Trainable")
    msg += "\n\t" + ind.format("%d" % opts.modelParamsNonTrainable, "Non-trainable")
    Print(msg, True)
    return

def PrintDataset(myDset):
    Verbose("%sPrinting 1 instance of the Numpy representation of the DataFrame:%s" % (ns, es), True)

    # For-loop: All TBranch entries
    for vList in myDset:
        # For-loop: All input variable values for given entry
        for v in vList:
            Verbose(v, False)
        break
    Verbose("%s" % (ns), True)
    return

def writeCfgFile(opts):

    # Write to json file
    jsonWr = JsonWriter(saveDir=opts.saveDir, verbose=opts.verbose)
    jsonWr.addParameter("keras version", opts.keras)
    jsonWr.addParameter("host name", opts.hostname)
    jsonWr.addParameter("python version", opts.python)
    jsonWr.addParameter("model", "adversarial")
    jsonWr.addParameter("model parameters (total)", opts.modelParams)
    jsonWr.addParameter("model parameters (trainable)", opts.modelParamsTrain)
    jsonWr.addParameter("model parameters (non-trainable)", opts.modelParamsNonTrainable)
    jsonWr.addParameter("model weights", opts.modelWeights)
    jsonWr.addParameter("model biases", opts.modelBiases)
    jsonWr.addParameter("standardised datasets", opts.standardise)
    jsonWr.addParameter("rndSeed", opts.rndSeed)
    jsonWr.addParameter("layers", len(opts.neurons_clf))
    jsonWr.addParameter("hidden layers", len(opts.neurons_clf)-2)
    jsonWr.addParameter("activation functions", [a for a in opts.activ_clf])
    jsonWr.addParameter("neurons", [n for n in opts.neurons_clf])
    jsonWr.addParameter("loss function", opts.loss_clf)
    jsonWr.addParameter("optimizer", opts.optimizer)
    jsonWr.addParameter("epochs", opts.epochs)
    jsonWr.addParameter("batch size", opts.batchSize)
    jsonWr.addParameter("train sample", opts.trainSample)
    jsonWr.addParameter("test sample" , opts.testSample)
    jsonWr.addParameter("ROOT file" , opts.rootFileName)
    jsonWr.addParameter("input variables", len(opts.inputList))
    jsonWr.addParameter("elapsed time", opts.elapsedTime)
    for i,b in enumerate(opts.inputList, 1):
        jsonWr.addParameter("var%d"% i, b)
    jsonWr.write(opts.cfgJSON)
    return

def writeGitFile(opts):

    # Write to json file
    path  = os.path.join(opts.saveDir, "gitBranch.txt")
    gFile = open(path, 'w')
    gFile.write(opts.gitBranch)

    path  = os.path.join(opts.saveDir, "gitStatus.txt")
    gFile = open(path, 'w')
    gFile.write(opts.gitStatus)

    path  = os.path.join(opts.saveDir, "gitDiff.txt")
    gFile = open(path, 'w')
    gFile.write(opts.gitDiff)
    return

def PlotLossFunction(TrainLosses, ValLosses, nepochs, saveDir):
    '''                                                                                                                                                                              
    Plot the loss function vs the number of the epoch                                                                                                                                
    '''
    def ApplyStyle(h, ymin, ymax, color, lstyle=ROOT.kSolid):
        h.SetLineColor(color)
        h.SetLineWidth(2)
        h.SetLineStyle(lstyle)
        if (lstyle != ROOT.kSolid):
            h.SetLineWidth(3)
        h.SetMaximum(ymax)
        h.SetMinimum(ymin)
        h.GetXaxis().SetTitle("# epoch")
        h.GetYaxis().SetTitle("Loss")

    ROOT.gStyle.SetOptStat(0)
    canvas = ROOT.TCanvas()
    canvas.cd()
    leg=plot.CreateLegend(0.6, 0.65, 0.9, 0.88)

    xarr = []
    for i in range(1, nepochs+1):
        xarr.append(i)

    if (opts.lam > 0):
        hLf_t  = ROOT.TGraph(len(xarr), array.array('d', xarr), array.array('d', TrainLosses["L_f"]))
        hLr_t  = ROOT.TGraph(len(xarr), array.array('d', xarr), array.array('d', TrainLosses["L_r"]))
        hLfr_t = ROOT.TGraph(len(xarr), array.array('d', xarr), array.array('d', TrainLosses["L_f - L_r"]))

        hLf_v  = ROOT.TGraph(len(xarr), array.array('d', xarr), array.array('d', ValLosses["L_f"]))
        hLr_v  = ROOT.TGraph(len(xarr), array.array('d', xarr), array.array('d', ValLosses["L_r"]))
        hLfr_v = ROOT.TGraph(len(xarr), array.array('d', xarr), array.array('d', ValLosses["L_f - L_r"]))

        # Find minimum and maximum
        ymin = 999.99
        ymax = 0
        for h in [hLf_t, hLr_t, hLfr_t, hLf_v, hLr_v, hLfr_v]:
            ymin = min(ymin, h.GetHistogram().GetMinimum())
            ymax = max(ymax, h.GetHistogram().GetMaximum())
        ymax = ymax + 1
        ymin = ymin - 1
        
        ApplyStyle(hLf_t, ymin, ymax, ROOT.kRed, ROOT.kDashed)
        ApplyStyle(hLr_t, ymin, ymax, ROOT.kBlue, ROOT.kDashed)
        ApplyStyle(hLfr_t, ymin, ymax, ROOT.kGreen+1, ROOT.kDashed)
        ApplyStyle(hLf_v, ymin, ymax, ROOT.kRed)
        ApplyStyle(hLr_v, ymin, ymax, ROOT.kBlue)
        ApplyStyle(hLfr_v, ymin, ymax, ROOT.kGreen+1)

        for i, h in enumerate([hLf_t, hLr_t, hLfr_t, hLf_v, hLr_v, hLfr_v], 0):
            if (i == 0):
                h.Draw("la")
            else:
                h.Draw("l same")

        leg.AddEntry(hLf_t, "L_{c} (train)","l")
        leg.AddEntry(hLr_t,"L_{r} (train)","l")
        leg.AddEntry(hLfr_t,"L_{c} - L_{r} (train)","l")
        leg.AddEntry(hLf_v, "L_{c} (test)","l")
        leg.AddEntry(hLr_v,"L_{r} (test)","l")
        leg.AddEntry(hLfr_v,"L_{c} - L_{r} (test)","l")

    else:
        h_TrainLoss = ROOT.TGraph(len(xarr), array.array('d', xarr), array.array('d', TrainLosses["loss"]))
        h_ValLoss   = ROOT.TGraph(len(xarr), array.array('d', xarr), array.array('d', ValLosses["val_loss"]))

        ymax = max(h_TrainLoss.GetHistogram().GetMaximum(), h_ValLoss.GetHistogram().GetMaximum()) + 1
        ymin = min(h_TrainLoss.GetHistogram().GetMinimum(), h_ValLoss.GetHistogram().GetMinimum()) - 1
        ApplyStyle(h_TrainLoss, ymin, ymax,  ROOT.kGreen +1, ROOT.kDashed)
        ApplyStyle(h_ValLoss, ymin, ymax,  ROOT.kGreen +1)
                
        h_TrainLoss.Draw("la")
        h_ValLoss.Draw("l same")
        leg.AddEntry(h_TrainLoss, "L_{c} (train)","l")
        leg.AddEntry(h_ValLoss, "L_{c} (test)","l")

    leg.Draw()
    saveName = "Loss_lam%s_epochs%s" %(opts.lam, len(xarr))
    plot.SavePlot(canvas, saveDir, saveName)
    canvas.Close()
    return

def getClassifier(model_clf, nInputs):
    # Add classifier layers           
    Print("Classifier network", True)
    for iLayer, n in enumerate(opts.neurons_clf, 0):
        layer = "layer#%d" % (int(iLayer)+1)
        if iLayer == len(opts.neurons_clf)-1:
            layer += " (output Layer)" # Layers of nodes between the input and output layers. There may be one or more of these layers.
        else:
            layer += " (hidden layer)" # A layer of nodes that produce the output variables.
            
        Print("Adding %s, with %s%d neurons%s and activation function %s" % (hs + layer + ns, ls, n, ns, ls + opts.activ_clf[iLayer] + ns), False)

        # First layer demands input_dim                                                                                                                                              
        if (iLayer == 0):
            model_clf.add(Dense(opts.neurons_clf[iLayer], input_dim =  nInputs))
            model_clf.add(Activation(opts.activ_clf[iLayer]))
        else:
            model_clf.add(Dense(opts.neurons_clf[iLayer]))
            model_clf.add(Activation(opts.activ_clf[iLayer]))
    return model_clf

def getRegressor(model_clf, model_reg):
    # Add classifier                                                                                                                                                                 
    model_reg.add(model_clf)

    # Add regressor layers                                                                                                                                                           
    for iLayer, n in enumerate(opts.neurons_reg, 0):
        model_reg.add(Dense(opts.neurons_reg[iLayer]))
        model_reg.add(Activation(opts.activ_reg[iLayer]))

    return model_reg

def make_trainable(model, isTrainable):
    '''
    Make layers trainiable / not trainable
    For not trainable layers, the weights are not updated 
    '''
    model.trainable = isTrainable
    for l in model.layers:
        l.trainable = isTrainable

def getLoss_clf(c):
    def loss_C(y_true, y_pred):
        return c * backend.binary_crossentropy(y_true, y_pred)
    return loss_C

def getLoss_reg(c):
    def loss_R(z_true, z_pred):
        if opts.loss_reg == "msle":
            loss =  c * keras.losses.mean_squared_logarithmic_error(z_true, z_pred) #c: controls the mass dependence
        elif opts.loss_reg == "mse":
            loss =  c * keras.losses.mean_squared_error(z_true, z_pred) #c: controls the mass dependence
        return loss
    return loss_R

def SaveTheModel(model_clf, last_epoch):
    modelName = "model_trained_%sEpochs.h5" % (last_epoch)
    model_clf.save(os.path.join(opts.saveDir,  modelName) )
    
    # Serialize model to JSON (contains arcitecture of model)
    model_json = model_clf.to_json() # model_clf.to_yaml()
    with open(opts.saveDir + "/model_architecture_%sEpochs.json" % (last_epoch), "w") as json_file:
        json_file.write(model_json)
    # Serialize weights to HDF5
    model_clf.save_weights(os.path.join(opts.saveDir, 'model_weights_%sEpochs.h5' % (last_epoch)), overwrite=True)
    model_clf.save(os.path.join(opts.saveDir, modelName))
        
    # Write weights and architecture in txt file
    modelFilename = os.path.join(opts.saveDir, "model_%sEpochs.txt" % (last_epoch))
    func.WriteModel(model_clf, model_json, opts.inputList, modelFilename, verbose=False)

def main(opts): 

    # Save start time (epoch seconds)
    tStart = time.time()
    Verbose("Started @ " + str(tStart), True)

    # Do not display canvases & disable screen output info
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")

    # Setup the style
    style = tdrstyle.TDRStyle() 
    style.setOptStat(False) 
    style.setGridX(opts.gridX)
    style.setGridY(opts.gridY)

    stdPlots = opts.standardise

    # Open the ROOT file
    Verbose("Opening ROOT file %s" %  (opts.rootFileName), True)
    ROOT.TFile.Open(opts.rootFileName)

    # Setting the seed for numpy-generated random numbers
    numpy.random.seed(opts.rndSeed)
    # Setting the seed for python random numbers
    rn.seed(opts.rndSeed)

    # Determine the number of threads running
    p = psutil.Process(os.getpid())
    Verbose("Number of threads is %s" % (ts + str(p.num_threads()) + ns), True)
    if 0:
        sess = tf.Session()
        Print("Number of threads after starting TF session is %s" % (ts + str(p.num_threads()) + ns), True)

    # Setting tensorflow random seed
    tf.set_random_seed(opts.rndSeed)

    # Open the signal and background TTrees with uproot (uproot allows one to read ROOT data, in python, without using ROOT)
    Verbose("Opening the signal and background TTrees with uproot using ROOT file %s" % (ts + opts.rootFileName + ns), True)
    signal     = uproot.open(opts.rootFileName)["treeS"]
    background = uproot.open(opts.rootFileName)["treeB"]    

    # Construct signal and background dataframes using a list of TBranches (a Dataframe is a two dimensional structure representing data in python)
    # https://uproot.readthedocs.io/en/latest/ttree-handling.html#id7
    nInputs   = len(opts.inputList)
    Verbose("Constucting dataframes for signal and background with %d input variables:\n\t%s%s%s" % (nInputs, ss, "\n\t".join(opts.inputList), ns), True)
    df_sig = signal.pandas.df(opts.inputList, entrystop=opts.entrystop) # call an array-fetching method to fill a Pandas DataFrame.
    df_bkg = background.pandas.df(opts.inputList, entrystop=opts.entrystop)
    # For testing:
    if 0:
        Print("Printing the signal DataFrame: %s" % df_sig, True)
        print
        Print("Printing the background DataFrame: %s" % df_bkg, True)
        print
        Print(df_sig['TrijetMass'], True)
        Print(df_sig['TrijetMass'].values, True)

    # Get the number of entries and columns of the DataFrames. Perform sanity checks
    nEntries_sig, nColumns_sig, nEntries_bkg, nColumns_bkg = GetDataFramesRowsColumns(df_sig, df_bkg)
    # Entries will be limited to the number of entries available from the signal and background samples
    nEntries = min(nEntries_sig, nEntries_bkg)

    # Apply rule-of-thumb to prevent over-fitting!
    #checkNeuronsPerLayer(nEntries, opts)

    # Assigning new column in our DataFrames DataFrames with label "signal". It takes value of "1" for signal and "0" for background.
    # Total number of columns is now increased by 1 (19 + 1 = 20)
    Verbose("Assigning a new column to our signal and bacbkround DataFrame labelled \"signal\" that will designate whether the entry is signal or not. Once this is done we can merge the DataFrames into a single one", True)
    df_sig = df_sig.assign(signal=1)
    df_bkg = df_bkg.assign(signal=0)
    Verbose("Printing tabular data for signal:\n%s%s%s"     % (ss, df_sig,ns), True)
    Verbose("Printing tabular data for background:\n%s%s%s" % (hs, df_bkg,ns), True)

    # Create a DataFrame list
    df_all  = pandas.concat( [df_sig, df_bkg] )
    Verbose("Printing the combined tabular data (signal first, background appended after signal):\n%s%s%s" % (ts, df_all, ns), True)
    if opts.standardise:
        msg  = "Standardising dataset features with the %sScaler%s" % (ls + opts.standardise, ns)
        Print(msg, True)    
        # WARNING! Don't cheat - fit only on training https://scikit-learn.org/stable/modules/neural_networks_supervised.html)
        scaler_sig, df_sig = func.GetStandardisedDataFrame(df_sig, opts.inputList, scalerType=opts.standardise) # fixme - iro - should fit/transform only the TEST data
        scaler_bkg, df_bkg = func.GetStandardisedDataFrame(df_bkg, opts.inputList, scalerType=opts.standardise)
        scaler_all, df_all = func.GetStandardisedDataFrame(df_all, opts.inputList, scalerType=opts.standardise)

    # Get a Numpy representation of the DataFrames for signal and background datasets (again, and AFTER assigning signal and background)
    Verbose("Getting a numpy representation of the DataFrames for signal and background datasets", True)
    dset_all = df_all.values

    # For future use
    jsonWr = JsonWriter(saveDir=opts.saveDir, verbose=opts.verbose)    

    # Define keras models
    # Classifier (characterizes events as signal or background)                                                                                                                        
    # Inputs: List of input variables                                                                                                                                                  
    # Output: classifier                                                                                                                                                               
    # Regressor: Predicts the top-quark mass from the classifiers output  

    Verbose("Creating the Keras models", True)
    model_clf = Sequential()    
    model_clf = getClassifier(model_clf, nInputs)
    model_reg = Sequential()
    model_reg = getRegressor(model_clf, model_reg)

    # Print a summary representation of your model
    if opts.verbose:
        Print("Printing model summary:", True)
        Print("Classifier", False)
        model_clf.summary()
        Print("Regressor", False)
        model_reg.summary()

    # Get the number of parameters of the model
    SaveModelParameters(model_clf, opts)

    # Get a dictionary containing the configuration of the model. 
    model_cfg = GetModelConfiguration(model_clf, opts.verbose)
        
    # Split data into input (X) and output (Y), using an equal number for signal and background. 
    # This is by creating the dataset with double the number of rows as is available for signal. Remember
    # that the background entries were appended to those of the signal.
    X      = dset_all[:2*nEntries, 0:nInputs]           # rows: 0 -> 2*Entries, columns: 0 -> 19 (all variables but not the "signal" column)
    Y      = dset_all[:2*nEntries, nInputs:]            # rows: 0 -> 2*Entries, column : 20  [i.e. everything after column 19 which is only the "signal" column (0 or 1)]
    Z      = dset_all[:2*nEntries, 7]                   #fixme! 7th row corresponce to top mass!
    X_sig  = dset_all[:nEntries, 0:nInputs]             # contains all 19 variables (total of "nEntries" values per variable) - from "treeS"
    X_bkg  = dset_all[nEntries:2*nEntries, 0:nInputs]   # contains all 19 variables (total of "nEntries" values per variable) - from "treeB"
    Print("Signal dataset has %s%d%s rows. Background dataset has %s%d%s rows" % (ss, len(X_sig), ns, es, len(X_bkg), ns), True)

    # Split the datasets (X= 19 inputs, Y=output variable, Z=target - trijet mass). test_size 0.5 means half for training half for testing. Shuffle the entry order?
    X_train, X_test, Y_train, Y_test, Z_train, Z_test = train_test_split(X, Y, Z, test_size=0.5, random_state = opts.rndSeed, shuffle=True)
    
    # Save size of test & training samples for future reference
    opts.testSample  = len(X_test)
    opts.trainSample = len(X_train)
    # batch size equal to the training batch size (See https://machinelearningmastery.com/use-different-batch-sizes-training-predicting-python-keras/)
    if opts.batchSize == None:
        opts.batchSize = len(X_train)/2

    # Early stop. Stop training when a monitored quantity has stopped improving.
    # Show patience of "50" epochs with a change in the loss function smaller than "min_delta" before stopping procedure
    # https://stackoverflow.com/questions/43906048/keras-early-stopping
    #earlystop = keras.callbacks.EarlyStopping(monitor='val_loss', min_delta=0.0, patience=10, verbose=1, mode='auto')
    earlystop = keras.callbacks.EarlyStopping(monitor='val_loss', min_delta=0.0001, patience=50)
    callbacks = [earlystop]

    # [Loss function is used to understand how well the network is working (compare predicted label with actual label via some function)]
    # Optimizer function is related to a function used to optimise the weights
    Print("Compiling the classifier with the loss function %s and optimizer %s " % (ls + opts.loss_clf + ns, ls + opts.optimizer + ns), True)
    model_clf.compile(loss=opts.loss_clf, optimizer=opts.optimizer, metrics=['accuracy'])
    
    # Define combined models
    inputs   = Input(shape=(X_train.shape[1],))
    optimizer_comb = opts.optimizer
    optimizer_adv  = opts.optimizer

    # Customise the optimiser settings
    if 0: #fixme! Add lr option 
        if opts.optimizer == "adam":
            optimizer_comb = keras.optimizers.Adam(learning_rate=0.001, beta_1=0.9, beta_2=0.999, amsgrad=False)
            optimizer_adv  = keras.optimizers.Adam(learning_rate=0.001, beta_1=0.9, beta_2=0.999, amsgrad=False)


    # Combined network: updates the classifier's weights                                                                                                                             
    # Input: list of discriminating variables                                                                                                                                        
    # Output: Classifier's output, top-quark mass prediction                                                                                                                          
    model_comb = Model(inputs=[inputs], outputs=[model_clf(inputs), model_reg(inputs)])

    # Regression layers are not trainable !                                                                                                                                          
    make_trainable(model_reg, False)
    make_trainable(model_clf, True)
    # Compile the combined model                                                                                                                                                     
    # L_comb = L_clf - lambda*L_reg                                                                                                                                                  
    model_comb.compile(loss=[getLoss_clf(c=1.0), getLoss_reg(c=-opts.lam)], optimizer=optimizer_comb)
    
    # Adversary network: predicts top-quark mass                                                                                                                                     
    model_adv = Model(inputs=[inputs], outputs=[model_reg(inputs)])

    # Classification layers are not trainable !                                                                                                                                      
    make_trainable(model_reg, True)
    make_trainable(model_clf, False)
    # Compile the model                                                                                                                                                              
    # L_adv = L_reg                                                                                                                                                                  
    model_adv.compile(loss=[getLoss_reg(c=1.0)], optimizer=optimizer_adv)

    # Pretrain the Classifier 
    make_trainable(model_reg, False)
    make_trainable(model_clf, True)

    Print("=== Classifier: Fit", True)
    model_clf.fit(X_train, Y_train, validation_data=(X_test, Y_test), epochs=10, batch_size=opts.batchSize, verbose = opts.debug)
    
    # Pretrain the adversarial network
    # Here the output is the top-quark mass (target) 
    make_trainable(model_reg, True)
    make_trainable(model_clf, False)

    Print("=== Adversarial Network: Fit", True)
    model_adv.fit(X_train, Z_train, validation_data=(X_test, Z_test), epochs=10, batch_size=opts.batchSize, verbose=opts.debug)

    Verbose("Number of threads before fitting model is %s" % (ts + str(p.num_threads()) + ns), True)    
    ###################################################
    # Train the model
    ###################################################
    
    # Dictionary with loss functions                                                                                                                                                 
    TrainLosses = {"L_f": [], "L_r": [], "L_f - L_r": []}
    ValLosses   = {"L_f": [], "L_r": [], "L_f - L_r": []}
    
    Print("Train the combined model")
    if (opts.lam > 0):
        for i in range(opts.epochs):
            # Evaluate the loss function after each epoch                            
            l_train = model_comb.evaluate(X_train, [Y_train, Z_train], verbose=opts.debug)                                                                                                
            l_val  = model_comb.evaluate(X_test, [Y_test, Z_test], verbose=opts.debug)

            # Store the value of the loss function after each epoch                                                                                                                  
            # Print("Metrics names: %s" %(model_comb.metrics_names), True) 
            # Metrics names: ['loss', 'sequential_1_loss', 'sequential_2_loss']
            TrainLosses["L_f - L_r"].append(l_train[0][None][0])
            TrainLosses["L_f"].append(l_train[1][None][0])
            TrainLosses["L_r"].append(-l_train[2][None][0])

            ValLosses["L_f - L_r"].append(l_val[0][None][0])
            ValLosses["L_f"].append(l_val[1][None][0])
            ValLosses["L_r"].append(-l_val[2][None][0])

            if (opts.epochs < 10) or (i % (opts.epochs/10) == 0):
                lf  = l_val[1][None][0]
                lr  = -l_val[2][None][0]
                lfr = l_val[0][None][0]
                Print("=== Epoch: %s / %s. Losses: Lf = %.3f, Lr = %.3f, (L_f - L_r) = %.3f" % (i, opts.epochs, lf, lr, lfr), False)
                # Save the model every "opts.epochs/10" epochs
                if i > 0:#(opts.epochs >= 100):                    
                    SaveTheModel(model_clf, i)
                    PlotLossFunction(TrainLosses, ValLosses, i, opts.saveDir)

            # Fit Classifier (with updated weights) to minimize joint loss function                                                                                                  
            make_trainable(model_reg, False)
            make_trainable(model_clf, True)
            indices = numpy.random.permutation(len(X_train))[:opts.batchSize]

            # Train and test on a batch                                                                                                                                              
            model_comb.train_on_batch(X_train[indices], [Y_train[indices], Z_train[indices]])
            model_comb.test_on_batch(X_test[indices], [Y_test[indices], Z_test[indices]])

            # Fit Regressor (with updated weights) to minimize Lr                                                                                                                    
            make_trainable(model_reg, True)
            make_trainable(model_clf, False)
            model_adv.fit(X_train, Z_train, validation_data=(X_test, Z_test), batch_size=opts.batchSize, epochs=1, callbacks=callbacks, verbose=opts.debug)
    else: #lambda == 0
        make_trainable(model_clf, True)
        hist = model_clf.fit(X_train, Y_train, validation_data=(X_test, Y_test), batch_size=opts.batchSize, epochs=opts.epochs, callbacks=callbacks, verbose=opts.debug)

    # Last epoch
    last_epoch = earlystop.stopped_epoch
    if (last_epoch == 0):
        last_epoch = opts.epochs
    epochList     = range(last_epoch)

    # Plot the input variables
    if opts.plotInputs:
        std_ = (opts.standardise != None)
        Verbose("Plotting all %d input variables for signal and bacgkround" % (len(opts.inputList)), True)
        for i, var in enumerate(opts.inputList, 0):
            
            # Get the lists
            sigList   = X[0:nEntries, i:i+1]          # first nEntries is signal. i:i+1 get the column (variable) of inteerest
            bkgList   = X[nEntries:2*nEntries, i:i+1] # after the first nEntries there is an equal number of background entries
            trainList = X_train[0:nEntries, i:i+1]
            testList  = X_test[0:nEntries, i:i+1]

            # Make the plots
            func.PlotInputs(sigList  , bkgList , var, "%s/%s" % (opts.saveDir, "sigVbkg")   , opts.saveFormats, pType="sigVbkg"   , standardise=std_)
            func.PlotInputs(trainList, testList, var, "%s/%s" % (opts.saveDir, "trainVtest"), opts.saveFormats, pType="trainVtest", standardise=std_)

    # Save the model
    SaveTheModel(model_clf, last_epoch)

    #https://keras.io/visualization/
    #https://machinelearningmastery.com/display-deep-learning-model-training-history-in-keras/
    if "fnal" not in opts.hostname:
        from keras.utils import plot_model
        modelFilename = "model.png"
        print molelFilename
        Print("plot a graph of the model and save it to the file %s" % (hs + os.path.basename(modelFilename) + ns), True)    
        plot_model(model_clf, to_file=modelFilename)

    # Produce method score (i.e. predict output value for given input dataset). Computation is done in batches.
    # https://stackoverflow.com/questions/49288199/batch-size-in-model-fit-and-model-predict
    Print("Generating output predictions (numpy arrays) for the input samples", True)
    pred_train  = model_clf.predict(X_train, batch_size=None, verbose=1, steps=None) # DNN output for training data (for both signal & bkg)
    pred_test   = model_clf.predict(X_test , batch_size=None, verbose=1, steps=None) # DNN output for test data (for both signal & bkg)
    pred_signal = model_clf.predict(X_sig  , batch_size=None, verbose=1, steps=None) # DNN output for signal only (all data)
    pred_bkg    = model_clf.predict(X_bkg  , batch_size=None, verbose=1, steps=None) # DNN output for data only (all data)

    # Join a sequence of arrays (X and Y) along an existing axis (1). In other words, add the ouput variable (Y) to the input variables (X)
    XY_train = numpy.concatenate((X_train, Y_train), axis=1)
    XY_test  = numpy.concatenate((X_test , Y_test ), axis=1)

    # Pick events with output = 1
    Verbose("Select events/samples which have an output variable Y (last column) equal to 1 (i.e. prediction is combatible with signal)", True)
    X_train_S = XY_train[XY_train[:,nInputs] == 1]; X_train_S = X_train_S[:,0:nInputs]
    X_test_S  = XY_test[XY_test[:,nInputs] == 1];   X_test_S  = X_test_S[:,0:nInputs]

    Verbose("Select events/samples which have an output variable Y (last column) equal to 0 (i.e. prediction is NOT combatible with signal)", False)
    X_train_B = XY_train[XY_train[:,nInputs] == 0]; X_train_B = X_train_B[:,0:nInputs]
    X_test_B  = XY_test[XY_test[:,nInputs] == 0];   X_test_B  = X_test_B[:,0:nInputs]
    
    # Produce method score for signal (training and test) and background (training and test)
    pred_train_S =  model_clf.predict(X_train_S, batch_size=None, verbose=1, steps=None)
    pred_train_B =  model_clf.predict(X_train_B, batch_size=None, verbose=1, steps=None)
    pred_test_S  =  model_clf.predict(X_test_S , batch_size=None, verbose=1, steps=None)
    pred_test_B  =  model_clf.predict(X_test_B , batch_size=None, verbose=1, steps=None)

    # Inform user of early stop
    stopEpoch = earlystop.stopped_epoch
    if stopEpoch != 0 and stopEpoch < opts.epochs:
        msg = "Early stop occured after %d epochs!" % (stopEpoch)
        opts.epochs = stopEpoch
        Print(cs + msg + ns, True)

    # Create json file
    writeGitFile(opts)

    # Plot selected output and save to JSON file for future use
    func.PlotAndWriteJSON(pred_signal , pred_bkg    , opts.saveDir, "Output"       , jsonWr, opts.saveFormats, **GetKwargs("Output"     )) # DNN score (predicted): Signal Vs Bkg (all data)
    func.PlotAndWriteJSON(pred_train  , pred_test   , opts.saveDir, "OutputPred"   , jsonWr, opts.saveFormats, **GetKwargs("OutputPred" )) # DNN score (sig+bkg)  : Train Vs Predict
    func.PlotAndWriteJSON(pred_train_S, pred_train_B, opts.saveDir, "OutputTrain"  , jsonWr, opts.saveFormats, **GetKwargs("OutputTrain")) # DNN score (training) : Sig Vs Bkg
    func.PlotAndWriteJSON(pred_test_S , pred_test_B , opts.saveDir, "OutputTest"   , jsonWr, opts.saveFormats, **GetKwargs("OutputTest" )) # DNN score (predicted): Sig Vs Bkg  

    rocDict = {}
    WPs     = [float(x)/float(40) for x in range(0, 40, 1)] 

    # For-loop: All branches to be plotted for various DNN score cuts (v. slow, especially for large number of "entrystop")
    for i, var in enumerate(opts.inputList, 0):
        if var != "TrijetMass":
            continue
        else:
            PrintFlushed("Variable %s (%d/%d)" % (var, i, len(opts.inputList)), True) #i==0
        
        Verbose("Plotting variable %s (DNN score = 0.0)" % (var), True)
        func.PlotAndWriteJSON(dset_all[0:nEntries, i:i+1], dset_all[nEntries:2*nEntries, i:i+1], opts.saveDir, var , jsonWr, opts.saveFormats, **GetKwargs(var, stdPlots))
        
        # Definitions
        xMin     = 140.0
        xMax     = 200.0
        xVals    = []
        xErrs    = []
        yVals    = []
        yErrs    = []
        print
        # For-loop: All Working Points (WPs)
        for wp in WPs:
            PrintFlushed("Plotting variable %s (DNN score = %.2f)" % (var, wp), wp==0.0)
            hList = func.PlotAndWriteJSON_DNNscore(pred_signal, pred_bkg, wp, dset_all[0:nEntries, i:i+1], dset_all[nEntries:2*nEntries, i:i+1], opts.saveDir, var , jsonWr, opts.saveFormats, **GetKwargs(var, stdPlots))
            rocDict[wp] = hList

            sig_all  = 0
            sig_pass = 0
            bkg_all  = 0
            bkg_pass = 0

            for j, x in enumerate(dset_all[0:nEntries, i:i+1], 0):
                sig_all += 1
                if pred_signal[j] < wp:
                    continue
                if x < xMin:
                    continue
                if x > xMax:
                    continue
                #print "x_s = %.2f, DNN = %.2f" % (x, pred_signal[j])
                sig_pass+=1

            for k, x in enumerate(dset_all[nEntries:2*nEntries, i:i+1], 0):
                bkg_all += 1
                if pred_bkg[k] < wp:
                    continue
                if x < xMin:
                    continue
                if x > xMax:
                    continue
                #print "x_b = %.2f, DNN = %.2f" % (x, pred_bkg[k])
                bkg_pass+=1
               
            xValue = 0.0
            yValue = 0.0
            if sig_all > 0.0:
                xValue = float(sig_pass)/float(sig_all)
            if bkg_all > 0.0:
                yValue = float(bkg_pass)/float(bkg_all)
            xVals.append(xValue)
            xErrs.append(0.001)
            yVals.append(yValue)
            yErrs.append(0.001)
            #print "%.2f WP) all_s = %.1f, pass_s = %.1f, x = %.2f | all_b = %.2f, pass_b = %.2f, y = %.2f" % (wp, sig_all, sig_pass, xValue, bkg_all, bkg_pass, yValue)

        gr = plot.GetGraph(xVals, yVals, xErrs, xErrs, yErrs, yErrs)
        gDict = {"graph" : [gr], "name" : [os.path.basename(opts.saveDir)]}
        rocName = "ROC_%s_GE%.0f_LE%.0f" % (var, xMin, xMax)
        #func.PlotROC(gDict, opts.saveDir, rocName, opts.saveFormats)
        func.PlotTGraph(xVals, xErrs, yVals, yErrs, opts.saveDir, rocName, jsonWr, opts.saveFormats, **GetKwargs("roc") )
        print

    # Plot overtraining test
    htrain_s, htest_s, htrain_b, htest_b = func.PlotOvertrainingTest(pred_train_S, pred_test_S, pred_train_B, pred_test_B, opts.saveDir, "OvertrainingTest", opts.saveFormats)

    # Plot summary plot (efficiency & singificance)
    func.PlotEfficiency(htest_s, htest_b, opts.saveDir, "Summary", opts.saveFormats)

    # Write efficiencies (signal and bkg)
    xVals_S, xErrs_S, effVals_S, effErrs_S  = func.GetEfficiency(htest_s)
    xVals_B, xErrs_B, effVals_B, effErrs_B  = func.GetEfficiency(htest_b)
    func.PlotTGraph(xVals_S, xErrs_S, effVals_S, effErrs_S, opts.saveDir, "EfficiencySig", jsonWr, opts.saveFormats, **GetKwargs("efficiency") )
    func.PlotTGraph(xVals_B, xErrs_B, effVals_B, effErrs_B, opts.saveDir, "EfficiencyBkg", jsonWr, opts.saveFormats, **GetKwargs("efficiency") )

    xVals, xErrs, sig_def, sig_alt = func.GetSignificance(htest_s, htest_b)
    func.PlotTGraph(xVals, xErrs, sig_def, effErrs_B, opts.saveDir, "SignificanceA", jsonWr, opts.saveFormats, **GetKwargs("significance") )
    func.PlotTGraph(xVals, xErrs, sig_alt, effErrs_B, opts.saveDir, "SignificanceB", jsonWr, opts.saveFormats, **GetKwargs("significance") )

    # Plot Loss functions
    xErr = [0.0 for i in range(0, opts.epochs)]
    yErr = [0.0 for i in range(0, opts.epochs)]
    
    # === Loss function vs # epoch                                                                                                                                                   
    if (opts.lam > 0):
        func.PlotTGraph(epochList, xErr, TrainLosses["L_f"]      , yErr , opts.saveDir, "TrainLossClf"    , jsonWr, opts.saveFormats, **GetKwargs("loss") )
        func.PlotTGraph(epochList, xErr, ValLosses["L_f"]        , yErr , opts.saveDir, "ValLossClf"      , jsonWr, opts.saveFormats, **GetKwargs("loss") )
        func.PlotTGraph(epochList, xErr, TrainLosses["L_r"]      , yErr , opts.saveDir, "TrainLossReg"    , jsonWr, opts.saveFormats, **GetKwargs("loss") )
        func.PlotTGraph(epochList, xErr, ValLosses["L_r"]        , yErr , opts.saveDir, "ValLossReg"      , jsonWr, opts.saveFormats, **GetKwargs("loss") )
        func.PlotTGraph(epochList, xErr, TrainLosses["L_f - L_r"], yErr , opts.saveDir, "TrainLossComb"   , jsonWr, opts.saveFormats, **GetKwargs("loss") )
        func.PlotTGraph(epochList, xErr, ValLosses["L_f - L_r"]  , yErr , opts.saveDir, "ValLossComb"     , jsonWr, opts.saveFormats, **GetKwargs("loss") )
        PlotLossFunction(TrainLosses, ValLosses, last_epoch, opts.saveDir)
    else:
        func.PlotTGraph(epochList, xErr, hist.history['loss'][0:last_epoch]     , yErr , opts.saveDir, "TrainLoss"    , jsonWr, opts.saveFormats, **GetKwargs("loss") )
        func.PlotTGraph(epochList, xErr, hist.history['val_loss'][0:last_epoch] , yErr , opts.saveDir, "ValLoss"      , jsonWr, opts.saveFormats, **GetKwargs("loss") )
        PlotLossFunction(hist.history, hist.history, last_epoch, opts.saveDir)
    
    # Plot ROC curve
    gSig  = func.GetROC(htest_s, htest_b)
    gDict = {"graph" : [gSig], "name" : [os.path.basename(opts.saveDir)]}
    func.PlotROC(gDict, opts.saveDir, "ROC", opts.saveFormats)

    # Write the resultsJSON file!
    jsonWr.write(opts.resultsJSON)
    
    # Print total time elapsed
    days, hours, mins, secs = GetTime(tStart)
    dt = "%s days, %s hours, %s mins, %s secs" % (days[0], hours[0], mins[0], secs)
    Print("Elapsed time: %s" % (hs + dt + ns), True)
    opts.elapsedTime = dt
    writeCfgFile(opts)
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
    #ROOTFILENAME = "/uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-TT_19var.root"               #  92 000
    ROOTFILENAME = "/uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-TT_19var_6Jets_2BJets.root"  # 460 000
    #ROOTFILENAME = "/uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-TT_19var_5Jets_1BJets.root"  # 875 000
    NOTBATCHMODE = False
    SAVEDIR      = None
    SAVEFORMATS  = "png,pdf"
    URL          = False
    STANDARDISE  = None
    ENTRYSTOP    = None
    VERBOSE      = False
    DEBUG        = False
    RNDSEED      = 1234
    EPOCHS       = 100
    BATCHSIZE    = None  # See: http://stats.stackexchange.com/questions/153531/ddg#153535
    LAMBDA       = 10
    ACTIV_CLF    = 'relu,relu,relu,sigmoid'
    ACTIV_REG    = 'tanh,relu'
    NEURONS_CLF  = '32,32,32,1'
    NEURONS_REG  = '16,1'
    LOSS_REG     = 'msle'
    LOSS_CLF     = 'binary_crossentropy'
    LOG          = False
    GRIDX        = False
    GRIDY        = False
    OPTIMIZER    = 'adam'
    CFGJSON      = "config.json"
    RESULTSJSON  = "results.json"
    PLOTINPUTS   = True
    INPUTVARS    = "TrijetPtDR,TrijetDijetPtDR,TrijetBjetMass,TrijetLdgJetBDisc,TrijetSubldgJetBDisc,TrijetBJetLdgJetMass,TrijetBJetSubldgJetMass,TrijetMass,TrijetDijetMass,TrijetBJetBDisc,TrijetSoftDrop_n2,TrijetLdgJetCvsL,TrijetSubldgJetCvsL,TrijetLdgJetPtD,TrijetSubldgJetPtD,TrijetLdgJetAxis2,TrijetSubldgJetAxis2,TrijetLdgJetMult,TrijetSubldgJetMult"


    # Define the available script options
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("--notBatchMode", dest="notBatchMode", action="store_true", default=NOTBATCHMODE, 
                      help="Disable batch mode (opening of X window) [default: %s]" % NOTBATCHMODE)

    parser.add_option("--standardise", dest="standardise", default=STANDARDISE,
                      help="Standardizing a dataset involves rescaling the distribution of INPUT values so that the mean of observed values is 0 and the standard deviation is 1 (e.g. StandardScaler) [default: %s]" % STANDARDISE)

    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=VERBOSE, 
                      help="Enable verbose mode (for debugging purposes mostly) [default: %s]" % VERBOSE)

    parser.add_option("--debug", dest="debug", action="store_true", default=DEBUG, help="Enable debugging mode [default: %s]" % DEBUG)

    parser.add_option("--plotInputs", dest="plotInputs", action="store_true", default=PLOTINPUTS, 
                      help="Enable plotting of input variables [default: %s]" % PLOTINPUTS)

    parser.add_option("--rootFileName", dest="rootFileName", type="string", default=ROOTFILENAME, 
                      help="Input ROOT file containing the signal and backbground TTrees with the various TBrances *variables) [default: %s]" % ROOTFILENAME)

    parser.add_option("--resultsJSON", dest="resultsJSON", default=RESULTSJSON,
                      help="JSON file containing the results [default: %s]" % (RESULTSJSON))

    parser.add_option("--cfgJSON", dest="cfgJSON", default=CFGJSON,
                      help="JSON file containing the configurations [default: %s]" % (CFGJSON))

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR,
                      help="Directory where all pltos will be saved [default: %s]" % SAVEDIR)

    parser.add_option("--url", dest="url", action="store_true", default=URL, 
                      help="Don't print the actual save path the histogram is saved, but print the URL instead [default: %s]" % URL)

    parser.add_option("-s", "--saveFormats", dest="saveFormats", default = SAVEFORMATS,
                      help="Save formats for all plots [default: %s]" % SAVEFORMATS)
    
    parser.add_option("--rndSeed", dest="rndSeed", type=int, default=RNDSEED, 
                      help="Value of random seed (integer) [default: %s]" % RNDSEED)
    
    parser.add_option("--epochs", dest="epochs", type=int, default=EPOCHS, 
                      help="Number of \"epochs\" to be used (how many times you go through your training set) [default: %s]" % EPOCHS)

    parser.add_option("--batchSize", dest="batchSize", type=int, default=BATCHSIZE,
                      help="The \"batch size\" to be used (= a number of samples processed before the model is updated). Batch size impacts learning significantly; typically networks train faster with mini-batches. However, the batch size has a direct impact on the variance of gradients (the larger the batch the better the appoximation and the larger the memory usage). [default: %s]" % BATCHSIZE)

    parser.add_option("--lam",dest="lam",  default=LAMBDA,  type = int,   help="Lambda (Default: %s)" % LAMBDA)
    
    parser.add_option("--loss_reg", dest="loss_reg",  default=LOSS_REG,   help="Regressor loss function (Default: %s)" % LOSS_REG)

    parser.add_option("--loss_clf", dest="loss_clf",  default=LOSS_CLF,   help="Classifier loss function (Default: %s)" % LOSS_CLF)

    parser.add_option("--activ_clf",dest="activ_clf", type="string", default=ACTIV_CLF,
                      help="List of activation functions of the classifier (comma-separated) [default: %s]" % ACTIV_CLF)

    parser.add_option("--activ_reg", dest="activ_reg", type="string", default=ACTIV_REG,
                      help="List of activation functions of the regressor (comma-separated) [default: %s]" % ACTIV_REG)

    parser.add_option("--neurons_clf", dest="neurons_clf", type="string", default=NEURONS_CLF,
                      help="List of neurons to use for each classification layer (comma-separated integers)  [default: %s]" % NEURONS_CLF)

    parser.add_option("--neurons_reg", dest="neurons_reg", type="string", default=NEURONS_REG,
                      help="List of neurons to use for each regression layer (comma-separated integers)  [default: %s]" % NEURONS_REG)

    parser.add_option("--inputVariables", dest="inputVariables", type="string", default=INPUTVARS,
                      help="List of input variables (TBranches) to use for the classifier (comma-separated integers)  [default: %s]" % INPUTVARS)
    
    parser.add_option("--entrystop", dest="entrystop", type="int", default=ENTRYSTOP,
                      help="Entry at which the (TBranch of TTree) reading stops. If ROOT file itoo big it may result to a \"memory error\", depending on the resources available in the machine used.  If \"None\", stop at the end of the branch. [default: %s]" % ENTRYSTOP)

    parser.add_option("--optimizer", dest="optimizer", type="string", default=OPTIMIZER,
                      help="Name of optimizer function; one of the two parameters required to compile a model: [default: %s]" % OPTIMIZER)

    parser.add_option("--log", dest="log", action="store_true", default=LOG,
                      help="Boolean for redirect sys.stdout to a log file [default: %s]" % LOG)

    parser.add_option("--gridX", dest="gridX", action="store_true", default=GRIDX,
                      help="Enable x-axis grid [default: %s]" % GRIDX)

    parser.add_option("--gridY", dest="gridY", action="store_true", default=GRIDY,
                      help="Enable y-axis grid [default: %s]" % GRIDY)

    (opts, parseArgs) = parser.parse_args()
    
    # Require at least two arguments (script-name, ...)
    if len(sys.argv) < 1:
        parser.print_help()
        sys.exit(1)
    else:
        pass
        
    # Create save formats
    if "," in opts.saveFormats:
        opts.saveFormats = opts.saveFormats.split(",")
    else:
        opts.saveFormats = [opts.saveFormats]
    opts.saveFormats = [s for s in opts.saveFormats]
    Verbose("Will save all output in %d format(s): %s" % (len(opts.saveFormats), ss + ", ".join(opts.saveFormats) + ns), True)

    # Create specification lists
    print "=== Classifier"
    if "," in opts.activ_clf:
        opts.activ_clf = opts.activ_clf.split(",")
    else:
        opts.activ_clf = [opts.activ_clf]
    Print("Activation = %s" % (opts.activ_clf), False)
    if "," in opts.neurons_clf:
        opts.neurons_clf = list(map(int, opts.neurons_clf.split(",")) )
    else:
        opts.neurons_clf = list(map(int, [opts.neurons_clf]))
    Print("Neurons = %s" % (opts.neurons_clf), False)
    
    print 
    print "=== Regressor"
    if "," in opts.activ_reg:
        opts.activ_reg = opts.activ_reg.split(",")
    else:
        opts.activ_reg = [opts.activ_reg]
    Print("Activation = %s" % (opts.activ_reg), False)
    if "," in opts.neurons_reg:
        opts.neurons_reg = list(map(int, opts.neurons_reg.split(",")) )
    else:
        opts.neurons_reg = list(map(int, [opts.neurons_reg]))
    Print("Neurons = %s" % (opts.neurons_reg), False)

        
    # Sanity checks (One activation function for each layer)
    if len(opts.neurons_clf) != len(opts.activ_clf):
        msg = "== Classifier: The list of neurons (size=%d) is not the same size as the list of activation functions (=%d)" % (len(opts.neurons_clf), len(opts.activ_clf))
        raise Exception(es + msg + ns)  
    if len(opts.neurons_reg) != len(opts.activ_reg):
        msg = "== Regressor: The list of neurons (size=%d) is not the same size as the list of activation functions (=%d)" % (len(opts.neurons_reg), len(opts.activ_reg))
        raise Exception(es + msg + ns)  

    # Sanity check (Last layer)
    if opts.neurons_clf[-1] != 1:
        msg = "Classification: The number of neurons for the last layer should be equal to 1 (=%d instead)" % (opts.neurons_clf[-1])
        raise Exception(es + msg + ns)   
    if opts.neurons_reg[-1] != 1:
        msg = "Regression: The number of neurons for the last layer should be equal to 1 (=%d instead)" % (opts.neurons_reg[-1])
        raise Exception(es + msg + ns)   

        #Print(es + msg + ns, True) 
    if opts.activ_clf[-1] != "sigmoid":
        msg = "Classification: The activation function for the last layer should be set to \"sigmoid\"  (=%s instead)" % (opts.activ_clf[-1])        
        Print(es + msg + ns, True)  #Print(es + msg + ns, True)

    # Define dir/logfile names
    specs = "%dLayers" % (len(opts.neurons_clf))
    for i,n in enumerate(opts.neurons_clf, 0):
        specs+= "_%s%s" % (opts.neurons_clf[i], opts.activ_clf[i])
    specs+= "_%sEpochs_%sBatchSize" % (opts.epochs, opts.batchSize)

    # Input list of discriminatin variables (TBranches)
    opts.inputList = opts.inputVariables.split(",")
    if opts.inputList < 1:
        raise Exception("At least one input variable needed to create the DNN. Only %d provided" % (len(opts.inputList)) )

    # Get the current date and time
    now    = datetime.now()
    nDay   = now.strftime("%d")
    nMonth = now.strftime("%h")
    nYear  = now.strftime("%Y")
    nTime  = now.strftime("%Hh%Mm%Ss") # w/ seconds
    nDate  = "%s%s%s_%s" % (nDay, nMonth, nYear, nTime)
    if opts.entrystop != None:
        sName  = "Adversarial_%s_%s" % (specs, str(opts.entrystop) + "Entrystop")
    else:
        sName  = "Adversarial_%s" % (specs)
    if opts.standardise != None:
        scalerTypes = ["Standard", "Robust", "MinMax"]
        if not opts.standardise in scalerTypes:
            msg = "Unsupported scaler type \"%s\". Please select one of the following (case-sensitive): %s" % (opts.standardise, ", ".join(scalerTypes) )
            raise Exception(es + msg + ns)
        #sName += "_Standardised"
        sName += "_%sScaler" % (opts.standardise)

    # Add the number of input variables
    sName += "_%dInputs" % (len(opts.inputList))
    # Add the lamda value
    sName += "_Lamda%s" % (opts.lam)
    # Add the time-stamp last
    sName += "_%s" % (nDate)
    
    # Determine path for saving plots
    if opts.saveDir == None:
        usrName = getpass.getuser()
        usrInit = usrName[0]
        myDir   = ""
        if "lxplus" in socket.gethostname():
            myDir = "/afs/cern.ch/user/%s/%s/public/html/" % (usrInit, usrName)
        else:
            myDir = os.path.join(os.getcwd())
        opts.saveDir = os.path.join(myDir, sName)
    
    # Create dir if it does not exist
    if not os.path.exists(opts.saveDir):
        os.mkdir(opts.saveDir)    

    # Create logfile
    opts.logFile = "stdout.log" #sName + ".log"
    opts.logPath = os.path.join(opts.saveDir, opts.logFile)
    bak_stdout   = sys.stdout
    log_file     = None
    if opts.log:
        # Keep a copy of the original "stdout"
        log_file   = open(opts.logPath, 'w')
        sys.stdout = log_file
        Print("Log file %s created" % (ls + opts.logFile + ns), True)

    # Inform user of network stup
    PrintNetworkSummary(opts)
    Print("A total of %s%d input variables%s will be used:\n\t%s%s%s" % (ls, len(opts.inputList),ns, ls, "\n\t".join(opts.inputList), ns), True)

    # See https://keras.io/activations/
    actList = ["elu", "softmax", "selu", "softplus", "softsign", "PReLU", "LeakyReLU",
               "relu", "tanh", "sigmoid", "hard_sigmoid", "exponential", "linear"] # Loukas used "relu" for resolved top tagger
    # Sanity checks
    for a in opts.activ_clf:
        if a not in actList:
            msg = "Unsupported activation function %s. Please select on of the following:%s\n\t%s" % (a, ss, "\n\t".join(actList))
            raise Exception(es + msg + ns)    
    for a in opts.activ_reg:
        if a not in actList:
            msg = "Unsupported activation function %s. Please select on of the following:%s\n\t%s" % (a, ss, "\n\t".join(actList))
            raise Exception(es + msg + ns)    

    # See https://keras.io/losses/
    lossList = ["binary_crossentropy", "is_categorical_crossentropy",
                "mean_squared_error", "mean_absolute_error", "mean_absolute_percentage_error", "mean_squared_logarithmic_error", "squared_hinge",
                "hinge", "categorical_hinge", "logcosh", "huber_loss", "categorical_crossentropy", "sparse_categorical_crossentropy", 
                "kullback_leibler_divergenc", "poisson", "cosine_proximity"]
    bLossList = ["binary_crossentropy", "is_categorical_crossentropy"]
    rLossList = ["mse", "mae", "mape", "msle"] #"mean_squared_error", "mean_absolute_error", "mean_absolute_percentage_error", "mean_squared_logarithmic_error"
    # Sanity checks
    if opts.loss_clf not in lossList:
        msg = "Unsupported loss function %s. Please select on of the following:%s\n\t%s" % (opts.loss_clf, ss, "\n\t".join(lossList))
        raise Exception(es + msg + ns)    
    elif opts.loss_clf not in bLossList:
        msg = "Binary output currently only supports the following loss fucntions: %s" % ", ".join(bLossList)
        raise Exception(es + msg + ns)
    else:
        pass

    if opts.loss_reg not in lossList + rLossList:
        msg = "Unsupported loss function %s. Please select on of the following:%s\n\t%s" % (opts.loss_reg, ss, "\n\t".join(lossList))
        raise Exception(es + msg + ns)    

    # See https://keras.io/optimizers/. Also https://www.dlology.com/blog/quick-notes-on-how-to-choose-optimizer-in-keras/
    optList = ["sgd", "rmsprop", "adagrad", "adadelta", "adam", "adamax", "nadam"]
    # optList = ["SGD", "RMSprop", "Adagrad", "Adadelta", "Adam", "Adamax", "Nadam"]
    if opts.optimizer not in optList:
        msg = "Unsupported loss function %s. Please select on of the following:%s\n\t%s" % (opts.optimizer, ss, "\n\t".join(optList))
        raise Exception(es + msg + ns)    

    # Determine and check the number of layers
    #checkNumberOfLayers(opts)

    # Get some basic information
    opts.keras      = keras.__version__
    opts.tensorflow = tf.__version__
    opts.hostname   = socket.gethostname()
    opts.python     = "%d.%d.%d" % (sys.version_info[0], sys.version_info[1], sys.version_info[2])
    opts.gitBranch  = subprocess.check_output(["git", "branch", "-a"])
    opts.gitStatus  = subprocess.check_output(["git", "status"])
    opts.gitDiff    = subprocess.check_output(["git", "diff"])

    # Call the main function
    Print("Hostname is %s" % (hs + opts.hostname + ns), True)
    Print("Using Keras %s" % (hs + opts.keras + ns), False)
    Print("Using Tensorflow %s" % (hs + opts.tensorflow + ns), False)

    main(opts)

    Print("Directory %s created" % (ls + opts.saveDir + ns), True)

    # Restore "stdout" to its original state and close the log file
    if opts.log:
        sys.stdout =  bak_stdout
        log_file.close()

    if opts.notBatchMode:
        raw_input("=== adversarial.py: Press any key to quit ROOT ...")
