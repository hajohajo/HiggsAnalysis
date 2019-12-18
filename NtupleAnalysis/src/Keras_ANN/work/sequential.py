#!/usr/bin/env python
'''
DESCRIPTION:
Keras is a modular, powerful and intuitive open-source Deep Learning library built on Theano and TensorFlow. 
Thanks to its minimalist, user-friendly interface, it has become one of the most popular packages to build, 
train, and test neural networks. This script provides an interface for conducting deep learning studies in 
python with Keras. The step included are:
1. Load Data.
2. Define Keras Model.
3. Compile Keras Model.
4. Fit Keras Model.
5. Evaluate Keras Model.
6. Tie It All Together.
7. Make Predictions

ENV SETUP:
cd ~/scratch0/CMSSW_10_2_11/src/ && cmsenv && setenv SCRAM_ARCH slc7_amd64_gcc700 && cd /afs/cern.ch/user/a/attikis/workspace/Keras_ANN/


OPTIMISATION:
In two words; trial & error. 
The number of layers and the number of nodes in Keras (aka "hyperparameters") should be tuned on a validation set 
(split from your original data into train/validation/test). Tuning just means trying different combinations of parameters
 and keep the one with the lowest loss value or better accuracy on the validation set, depending on the problem. 
There are two basic methods:
1) Grid search: For each parameter, decide a range and steps into that range, like 8 to 64 neurons, in powers of two 
(e.g. 8, 16, 32, 64), and try each combination of the parameters. 
2) Random search: Do the same but just define a range for each parameter and try a random set of parameters, drawn from 
an uniform distribution over each range.
Unfortunately there is no other way to tune such parameters. Just keep in mind that a deep model provides a hierarchy of 
layers that build up increasing levels of abstraction from the space of the input variables to the output variables.

About layers having different number of neurons, that could come from the tuning process, or you can also see it as 
dimensionality reduction, like a compressed version of the previous layer.


ACTIVATION FUNCTIONS:
Activation functions are really important for a Artificial Neural Network to learn and make sense of something really 
complicated and Non-linear complex functional mappings between the inputs and response variable. They introduce non-linear
properties to our Network.Their main purpose is to convert a input signal of a node in a A-NN to an output signal. 
That output signal now is used as a input in the next layer in the stack. Specifically in A-NN we do the sum of products of
inputs (x_{i}) and their corresponding Weights(w_{i}) and apply a Activation function f(x_{i}) to it to get the output of that 
layer and feed it as an input to the next layer. The question arises that why can't we do it without activating the input signal?
If we do not apply a Activation function then the output signal would simply be a simple linear function.A linear function is 
just a polynomial of one degree. Now, a linear equation is easy to solve but they are limited in their complexity and have less 
power to learn complex functional mappings from data. A Neural Network without Activation function would simply be a Linear 
regression Model, which has limited power and does not perform good most of the times. We want our Neural Network to not just 
learn and compute a linear function but something more complicated than that. Also without activation function our Neural network
 would not be able to learn and model other complicated kinds of data such as images, videos , audio , speech etc. That is why
 we use Artificial Neural network techniques such as Deep learning to make sense of something complicated ,high dimensional, 
non-linear -big datasets, where the model has lots and lots of hidden layers in between and has a very complicated architecture
 which helps us to make sense and extract knowledge form such complicated big datasets.


USAGE:
./sequential.py [opts]


EXAMPLES:
./sequential.py --activation relu,relu,sigmoid --neurons 36,19,1 --epochs 10000 --batchSize 50000 -s pdf --log
./sequential.py --activation relu,relu,sigmoid --neurons 19,190,1 --epochs 50 --batchSize 500 -s pdf --saveDir ~/public/html/Test
./sequential.py --activation relu,relu,sigmoid --neurons 50,50,1 --epochs 200 --batchSize 2000 -s pdf --entrystop 600000 --inputVariables "TrijetPtDR,TrijetDijetPtDR,TrijetBjetMass,TrijetLdgJetBDisc,TrijetSubldgJetBDisc,TrijetBJetLdgJetMass,TrijetBJetSubldgJetMass,TrijetMass,TrijetDijetMass,TrijetBJetBDisc,TrijetSoftDrop_n2,TrijetLdgJetCvsL,TrijetSubldgJetCvsL,TrijetLdgJetPtD,TrijetSubldgJetPtD,TrijetLdgJetAxis2,TrijetSubldgJetAxis2,TrijetLdgJetMult,TrijetSubldgJetMult"
./sequential.py --activation relu,relu,sigmoid --neurons 50,50,1 --epochs 200 --batchSize 50000 -s pdf --entrystop 600000 --standardise --scaleBack


LAST USED:
./sequential.py --activation relu,relu,sigmoid --neurons 50,50,1 --epochs 5 --batchSize 20000 -s pdf --entrystop 20000 --standardise Robust


GITHUB:
https://github.com/skonstantinou/Keras_ANN/
https://github.com/attikis/Keras_ANN


URL:
https://stats.stackexchange.com/questions/181/how-to-choose-the-number-of-hidden-layers-and-nodes-in-a-feedforward-neural-netw <----- Really good
https://towardsdatascience.com/activation-functions-neural-networks-1cbd9f8d91d6
https://theffork.com/activation-functions-in-neural-networks/?source=post_page-----4dfb9c7ce9c9----------------------
http://www.faqs.org/faqs/ai-faq/neural-nets/part2/ <----- Really good
https://github.com/Kulbear/deep-learning-nano-foundation/wiki/ReLU-and-Softmax-Activation-Functions
https://ml-cheatsheet.readthedocs.io/en/latest/loss_functions.html
https://gombru.github.io/2018/05/23/cross_entropy_loss/
https://stackoverflow.com/questions/36950394/how-to-decide-the-size-of-layers-in-keras-dense-method
https://machinelearningmastery.com/how-to-configure-the-number-of-layers-and-nodes-in-a-neural-network/
https://machinelearningmastery.com/tutorial-first-neural-network-python-keras/
https://keras.io/activations/
https://keras.io/getting-started/
https://keras.io/getting-started/faq/#what-does-sample-batch-epoch-mean
https://indico.cern.ch/event/749214/attachments/1754601/2844298/CERN_TH_Seminar_Nov2018.pdf#search=keras%20AND%20StartDate%3E%3D2018%2D10%2D01
https://www.indico.shef.ac.uk/event/11/contributions/338/attachments/281/319/rootTutorialWeek5_markhod_2018.pdf
https://indico.cern.ch/event/683620/timetable/#20190512.detailed
https://cds.cern.ch/record/2157570

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
from keras.models import Sequential
from keras.layers import Dense, Activation, Flatten
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

def PrintXYTestTrain(x_tr, x_te, y_tr, y_te, w_tr, w_te, nEntries, verbose=False):
    if not verbose:
        return
    table = []
    align = "{:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}"
    table.append("="*110)
    title = align.format("index", "X_train", "X_test", "Y_train", "Y_test", "W_train", "W_test")
    table.append(title)
    table.append("="*110)
    # For-loop: All entries
    for i,x in enumerate(x_tr, 0):
        msg = align.format("%d" % i, "%.2f" %  x_tr[i], "%.2f" %  x_te[i], "%.2f" %  y_tr[i], "%.2f" %  y_te[i], "%.2f" % w_tr[i], "%.2f" % w_te[i])
        if i < (nEntries/2):
            table.append(ss + msg + ns)
        else:
            table.append(es + msg + ns)
    table.append("="*100)
    for i,r in enumerate(table, 0):
        Print(r, i==0)
    return

def PrintNetworkSummary(opts):
    table    = []
    msgAlign = "{:^10} {:^10} {:>12} {:>10}"
    title    =  msgAlign.format("Layer #", "Neurons", "Activation", "Type")
    hLine    = "="*len(title)
    table.append(hLine)
    table.append(title)
    table.append(hLine)
    for i, n in enumerate(opts.neurons, 0): 
        layerType = "unknown"
        if i == 0:
            layerType = "input"
        elif i+1 == len(opts.neurons):
            layerType = "output"
        else:
            layerType = "hidden"
        table.append( msgAlign.format(i+1, opts.neurons[i], opts.activation[i], layerType) )
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

def eval(myModel, x_test, y_test, opts):
    '''
    https://keras.io/models/model/#evaluate

    NOTE: computation is done in batches
    '''
    metrics = myModel.evaluate(
        x=x_test,
        y=y_test,
        batch_size=opts.batchSize,
        verbose=1,
        sample_weight=None
    )
    print "*"*100
    print ('Loss: {:.3f}, Accuracy: {:.3f}'.format(metrics[0], metrics[1]))
    print "myModel.metrics_names = ", myModel.metrics_names
    print "len(metrics) = ", len(metrics)
    print 
    print "*"*100
    return

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

    if var == "TrijetPtDR":
        kwargs["xMin"]   =  0.0
        kwargs["xMax"]   = 800.0
        kwargs["nBins"]  =  80
        kwargs["xTitle"] = "p_{T}#DeltaR_{t}"
        kwargs["yMin"]   = 1e-3        
        #kwargs["yTitle"] = "a.u. / %0.0f"

    if var == "TrijetDijetPtDR":
        kwargs["xMin"]   =  0.0
        kwargs["xMax"]   = 800.0
        kwargs["nBins"]  = 160
        kwargs["xTitle"] = "p_{T}#DeltaR_{W}"
        kwargs["yMin"]   = 1e-3
        #kwargs["yTitle"] = "a.u. / %0.f"

    if var == "TrijetBjetMass":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   = 100.0
        kwargs["nBins"]  = 200
        kwargs["xTitle"] = "m_{b} [GeV]"
        #kwargs["yTitle"] = "a.u. / %0.f GeV"

    if var == "TrijetLdgJetBDisc":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   =   1.0
        kwargs["nBins"]  = 100
        kwargs["xTitle"] = "b-tag discr."
        #kwargs["yTitle"] = "a.u. / %0.2f"

    if var == "TrijetSubldgJetBDisc":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   =   1.0
        kwargs["nBins"]  = 200
        kwargs["xTitle"] = "b-tag discr."
        #kwargs["yTitle"] = "a.u. / %0.2f"

    if var == "TrijetBJetLdgJetMass":
        kwargs["xMin"]   =  0.0
        kwargs["xMax"]   = 600.0
        kwargs["nBins"]  = 300
        kwargs["xTitle"] = "m_{b,j_{1}} [GeV]"
        kwargs["yMin"]   = 1e-3

    if var == "TrijetBJetSubldgJetMass":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   = 600.0
        kwargs["nBins"]  = 300
        kwargs["xTitle"] = "m_{b,j_{2}} [GeV]"
        kwargs["yMin"]   = 1e-3
        #kwargs["yTitle"] = "a.u. / %0.0f GeV"

    if var == "TrijetMass":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   = 900 #900.0
        kwargs["nBins"]  = 180 #450
        kwargs["xTitle"] = "m_{t} [GeV]"
        kwargs["yMin"]   = 1e-3
        kwargs["yMax"]   = 1.0
        #kwargs["yTitle"] = "a.u. / %0.0f GeV"

    if var == "TrijetDijetMass":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   = 500.0
        kwargs["nBins"]  = 250
        kwargs["xTitle"] = "m_{W} [GeV]"
        kwargs["yMin"]   = 1e-3
        #kwargs["yTitle"] = "a.u. / %0.0f GeV"

    if var == "TrijetBJetBDisc":
        kwargs["xMin"]   =  0.0
        kwargs["xMax"]   =  1.0
        kwargs["nBins"]  = 200
        kwargs["xTitle"] = "b-tag discr."
        #kwargs["yTitle"] = "a.u. / %0.2f"

    if var == "TrijetSoftDrop_n2":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   =   2.0
        kwargs["nBins"]  = 400
        kwargs["xTitle"] = "soft-drop n_{2}"
        kwargs["yMin"]   = 1e-3
        #kwargs["yTitle"] = "a.u. / %0.2f"

    if var == "TrijetLdgJetCvsL":
        kwargs["xMin"]   =  -1.0
        kwargs["xMax"]   =  +1.0
        kwargs["nBins"]  = 400
        kwargs["xTitle"] = "CvsL discr."
        kwargs["yMin"]   = 1e-3
        kwargs["yMax"]   = 1e-1
        #kwargs["yTitle"] = "a.u. / %0.2f"

    if var == "TrijetSubldgJetCvsL":
        kwargs["xMin"]   =  -1.0
        kwargs["xMax"]   =  +1.0
        kwargs["nBins"]  = 400
        kwargs["xTitle"] = "CvsL discr."
        kwargs["yMin"]   = 1e-3
        kwargs["yMax"]   = 1e-1
        #kwargs["yTitle"] = "a.u. / %0.2f"

    if var == "TrijetLdgJetPtD":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   =   1.0
        kwargs["nBins"]  = 200
        kwargs["xTitle"] = "p_{T}D"
        #kwargs["yTitle"] = "a.u. / %0.2f"

    if var == "TrijetSubldgJetPtD":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   =   1.0
        kwargs["nBins"]  = 200
        kwargs["xTitle"] = "p_{T}D"
        #kwargs["yTitle"] = "a.u. / %0.2f"

    if var == "TrijetLdgJetAxis2":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   =   0.16
        kwargs["nBins"]  = 320
        kwargs["xTitle"] = "axis2"
        #kwargs["yTitle"] = "a.u. / %0.3f"

    if var == "TrijetSubldgJetAxis2":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   =   0.16
        kwargs["nBins"]  = 320
        kwargs["xTitle"] = "axis2"
        #kwargs["yTitle"] = "a.u. / %0.3f"

    if var == "TrijetLdgJetMult":
        kwargs["xMin"]  =  0.0
        kwargs["xMax"]   = 40.0
        kwargs["nBins"]  = 80
        kwargs["xTitle"] = "constituent multiplicity"
        #kwargs["yTitle"] = "a.u. / %0.0f"

    if var == "TrijetSubldgJetMult":
        kwargs["xMin"]   =  0.0
        kwargs["xMax"]   = 40.0
        kwargs["nBins"]  = 80
        kwargs["xTitle"] = "constituent multiplicity"
        #kwargs["yTitle@"] = "a.u. / %0.0f"

    if standardise:
        kwargs["xMin"]  =  -5.0
        kwargs["xMax"]  =  +5.0
        kwargs["nBins"] = 500 #1000

    return kwargs


def CrossEntropy(yHat, y):
    ''' 
    For future class implementation

    https://ml-cheatsheet.readthedocs.io/en/latest/loss_functions.html
    ''' 
    if y == 1:
      return -log(yHat)
    else:
      return -log(1 - yHat)

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
    nParams, nWeights, nBias = GetModelWeightsAndBiases(opts.inputList, opts.neurons)
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

def checkNumberOfLayers(opts):

    opts.nLayers       = len(opts.neurons)+1
    opts.nHiddenLayers = len(opts.neurons)-1

    if opts.nHiddenLayers < 1:
        msg = "The NN has %d hidden layers. It must have at least 1 hidden layer" % (opts.nHiddenLayers)
        raise Exception(es + msg + ns)    
    else:
        msg = "The NN has %s%d hidden layers%s (in addition to input and output layers)." % (hs, opts.nHiddenLayers, ns)
        Verbose(msg, True) 
        msg = "[NOTE: One hidden layer is sufficient for the large majority of problems. In general, the optimal size of the hidden layer is usually between the size of the input and size of the output layers.]"
        Verbose (msg, False)         
    return

def checkNeuronsPerLayer(nsignal, opts):
    '''
    https://stats.stackexchange.com/questions/181/how-to-choose-the-number-of-hidden-layers-and-nodes-in-a-feedforward-neural-netw
    '''
    # For-loop: All neurons
    for i, o in enumerate(opts.neurons, 1):
        # Skip first (input layer) and last (output) layers (i.e. do only hidden layers)
        if i==0 or i>=len(opts.neurons)-1:
            continue
        # Determine in/out neurons and compare to min-max range allowed
        nN   = opts.neurons[i]
        nIn  = opts.neurons[0]
        nOut = opts.neurons[-1]
        dof  = (2*nsignal) * (nIn + nOut)     # 2*nsignal = number of samples in training data set
        nMax = (2*nsignal) / (2*(nIn + nOut))
        nMin = (2*nsignal) / (10*(nIn + nOut))
        if nN > nIn or nN < nOut:
            msg = "The number of OUT neurons for layer #%d (%d neurons) is not within the recommended range(%d <= N <= %d). This may result in over-fitting" % (i, nN, nIn, nOut)
            Print(hs + msg + ns, True)
    return

def writeCfgFile(opts):

    # Write to json file
    jsonWr = JsonWriter(saveDir=opts.saveDir, verbose=opts.verbose)
    jsonWr.addParameter("keras version", opts.keras)
    jsonWr.addParameter("host name", opts.hostname)
    jsonWr.addParameter("python version", opts.python)
    jsonWr.addParameter("model", "sequential")
    jsonWr.addParameter("model parameters (total)", opts.modelParams)
    jsonWr.addParameter("model parameters (trainable)", opts.modelParamsTrain)
    jsonWr.addParameter("model parameters (non-trainable)", opts.modelParamsNonTrainable)
    jsonWr.addParameter("model weights", opts.modelWeights)
    jsonWr.addParameter("model biases", opts.modelBiases)
    jsonWr.addParameter("standardised datasets", opts.standardise)
    jsonWr.addParameter("decorrelate",  opts.decorrelate)
    jsonWr.addParameter("rndSeed", opts.rndSeed)
    jsonWr.addParameter("layers", len(opts.neurons))
    jsonWr.addParameter("hidden layers", len(opts.neurons)-2)
    jsonWr.addParameter("activation functions", [a for a in opts.activation])
    jsonWr.addParameter("neurons", [n for n in opts.neurons])
    jsonWr.addParameter("loss function", opts.lossFunction)
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

    if opts.scaleBack:
        stdPlots = False
    else:
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
        Print("Printing the signal DaaFrame: %s" % df_sig, True)
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
    checkNeuronsPerLayer(nEntries, opts)

    # Optional numpy array of weights for the training & testing samples, used for weighting the loss function (during training only). 
    # Can be used to decorrelate a variable before training; alternative approach to adversarial neural network (ANN) for combatting mass sculpting
    if opts.decorrelate != None:        
        msg = "Applying sample reweighting so that signal and background distributions of variable \"%s\" match"  % opts.decorrelate
        Print(cs + msg + ns, True)
        weights = func.doSampleReweighing(df_sig, df_bkg, "TrijetMass", _verbose=opts.verbose, **GetKwargs("TrijetMass"))

    # Assigning new column in our DataFrames DataFrames with label "signal". It takes value of "1" for signal and "0" for background.
    # Total number of columns is now increased by 1 (19 + 1 = 20)
    Verbose("Assigning a new column to our signal and bacbkround DataFrame labelled \"signal\" that will designate whether the entry is signal or not. Once this is done we can merge the DataFrames into a single oneOB", True)
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

    # Define keras model as a linear stack of layers. Add layers one at a time until we are happy with our network architecture.
    Verbose("Creating the sequential Keras model", True)
    myModel = Sequential()

    # The best network structure is found through a process of trial and error experimentation. Generally, you need a network large enough to capture the structure of the problem.
    # The Dense function defines each layer - how many neurons and mathematical function to use.
    for iLayer, n in enumerate(opts.neurons, 0):
        layer = "layer#%d" % (int(iLayer)+1)
        if iLayer == len(opts.neurons)-1:
            layer += " (output Layer)" # Layers of nodes between the input and output layers. There may be one or more of these layers.
        else:            
            layer += " (hidden layer)" # A layer of nodes that produce the output variables.
            
        Print("Adding %s, with %s%d neurons%s and activation function %s" % (hs + layer + ns, ls, n, ns, ls + opts.activation[iLayer] + ns), iLayer==0)
        if iLayer == 0:
            # Only first layer demands input_dim. For the rest it is implied.
            myModel.add( Dense(opts.neurons[iLayer], input_dim = nInputs) ) #, weights = [np.zeros([692, 50]), np.zeros(50)] OR bias_initializer='zeros',  or bias_initializer=initializers.Constant(0.1)
            myModel.add( Activation(opts.activation[iLayer]) )
            # myModel.add( Dense(opts.neurons[iLayer], input_dim = nInputs, activation = opts.activation[iLayer]) ) # her majerty Soti requested to break this into 2 lines
        else:
            if 0: #opts.neurons[iLayer] < nInputs:
                msg = "The number of  neurons (=%d) is less than the number of input variables (=%d). Please set the number of neurons to be at least the number of inputs!" % (opts.neurons[iLayer], nInputs)
                raise Exception(es + msg + ns)  
            myModel.add( Dense(opts.neurons[iLayer]))
            myModel.add( Activation(opts.activation[iLayer]) )
            # myModel.add( Dense(opts.neurons[iLayer], activation = opts.activation[iLayer]) ) 

    # Print a summary representation of your model
    if opts.verbose:
        Print("Printing model summary:", True)
        myModel.summary()
    
    # Get the number of parameters of the model
    SaveModelParameters(myModel, opts)

    # Get a dictionary containing the configuration of the model. 
    model_cfg = GetModelConfiguration(myModel, opts.verbose)
        
    # Split data into input (X) and output (Y), using an equal number for signal and background. 
    # This is by creating the dataset with double the number of rows as is available for signal. Remember
    # that the background entries were appended to those of the signal.
    X     = dset_all[:2*nEntries, 0:nInputs]           # rows: 0 -> 2*Entries, columns: 0 -> 19 (all variables but not the "signal" column)
    Y     = dset_all[:2*nEntries, nInputs:]            # rows: 0 -> 2*Entries, column : 20  [i.e. everything after column 19 which is only the "signal" column (0 or 1)]
    X_sig = dset_all[0:nEntries, 0:nInputs]            # contains all 19 variables (total of "nEntries" values per variable) - from "treeS"
    X_bkg = dset_all[nEntries:2*nEntries, 0:nInputs]   # contains all 19 variables (total of "nEntries" values per variable) - from "treeB"
    Y_sig = dset_all[:nEntries, nInputs:]              # contains 1's (total of "nEntries" values)   - NOT USED
    Y_bkg = dset_all[:nEntries, nInputs:]              # contains 0's (total ODof "nEntries" values) - N0T USED
    Print("Signal dataset has %s%d%s rows. Background dataset has %s%d%s rows" % (ss, len(X_sig), ns, es, len(X_bkg), ns), True)

    # Split the datasets (X= 19 inputs, Y=output variable). test_size 0.5 means half for training half for testing. Shuffle the entry order?
    if opts.decorrelate != None:
        X_train, X_test, Y_train, Y_test, W_train, W_test = train_test_split(X, Y, weights["all"], test_size=0.5, random_state = opts.rndSeed, shuffle=True)
        x_tr = X_train[0:nEntries, 7].tolist()
        x_te = X_test[0:nEntries, 7].tolist() # iro - fixme: index ambiguous! 
        y_tr = Y_train[0:nEntries, 0].tolist()
        y_te = Y_test[0:nEntries, 0].tolist()
        w_tr = W_train.tolist()
        w_te = W_test.tolist()
        PrintXYTestTrain(x_tr, x_te, y_tr, y_te, w_tr, w_te, nEntries, opts.verbose)
    else:
        X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.5, random_state = opts.rndSeed, shuffle=True)

    # Save size of test & training samples for future reference
    opts.testSample  = len(X_test)
    opts.trainSample = len(X_train)
    # batch size equal to the training batch size (See https://machinelearningmastery.com/use-different-batch-sizes-training-predicting-python-keras/)
    if opts.batchSize == None:
        opts.batchSize = len(X_train)/2


    # Early stop? Stop training when a monitored quantity has stopped improving.
    # Show patience of "50" epochs with a change in the loss function smaller than "min_delta" before stopping procedure
    # https://stackoverflow.com/questions/43906048/keras-early-stopping
    #earlystop = keras.callbacks.EarlyStopping(monitor='val_loss', min_delta=0.0, patience=10, verbose=1, mode='auto')
    earlystop = keras.callbacks.EarlyStopping(monitor='val_loss', min_delta=0.0001, patience=50)
    callbacks = [earlystop]


    # [Loss function is used to understand how well the network is working (compare predicted label with actual label via some function)]
    # Optimizer function is related to a function used to optimise the weights
    Print("Compiling the model with the loss function %s and optimizer %s " % (ls + opts.lossFunction + ns, ls + opts.optimizer + ns), True)
    myModel.compile(loss=opts.lossFunction, optimizer=opts.optimizer, metrics=['accuracy'])
    

    # Customise the optimiser settings?
    if 0: #opts.optimizer == "adam":  # does not work (required more recent version of Keras?)
        # Default parameters follow those provided in the original paper.
        # learning_rate: float >= 0. Learning rate.
        # beta_1: float, 0 < beta < 1. Generally close to 1.
        # beta_2: float, 0 < beta < 1. Generally close to 1.
        # amsgrad: boolean. Whether to apply the AMSGrad variant of this algorithm from the paper "On the Convergence of Adam and Beyond".
        keras.optimizers.Adam(learning_rate=0.001, beta_1=0.9, beta_2=0.999, amsgrad=False)


    # Fit the model with our data
    # (An "epoch" is an arbitrary cutoff, generally defined as "one iteration of training on the whole dataset", 
    # used to separate training into distinct phases, which is useful for logging and periodic evaluation.)
    sampleWeights  = None
    validationData = (X_test, Y_test)
    if opts.decorrelate != None:
        # Optional array of the same length as X_train, containing weights to apply to the model's loss for each sample.
        sampleWeights  = W_train
        validationData = (X_test, Y_test, W_test) # sample weights for validation data
    else:
        pass
    Verbose("Number of threads before fitting model is %s" % (ts + str(p.num_threads()) + ns), True)
    seqModel = myModel.fit(X_train,
                           Y_train,
                           validation_data=validationData,
                           epochs     = opts.epochs,    # a full pass over all of your training data
                           batch_size = opts.batchSize, # a set of N samples (https://stats.stackexchange.com/questions/153531/what-is-batch-size-in-neural-network)
                           shuffle    = False,
                           verbose    = 1,              # 0=silent, 1=progress, 2=mention the number of epoch
                           callbacks  = callbacks,
                           sample_weight=sampleWeights
                           )
    Verbose("Number of threads after fitting model is %s" % (ts + str(p.num_threads()) + ns), True)
        
    # Retrieve  the training / validation loss / accuracy at each epoch
    Verbose("The available history objects of the model are: %s" % (", ".join(seqModel.history.keys())), True)
    # For plotting: https://machinelearningmastery.com/display-deep-learning-model-training-history-in-keras/
    trainLossList = seqModel.history['loss']
    trainAccList  = seqModel.history['acc']
    valLossList   = seqModel.history['val_loss']
    valAccList    = seqModel.history['val_acc']
    epochList     = range(opts.epochs)

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
            if opts.decorrelate != None:
                func.PlotInputs(sigList  , bkgList , var, "%s/%s" % (opts.saveDir, "sigVbkg_weighted")   , opts.saveFormats, pType="sigVbkg"   , standardise=std_, w1=weights["sig"], w2=weights["bkg"])
                func.PlotInputs(trainList, testList, var, "%s/%s" % (opts.saveDir, "trainVtest_weighted"), opts.saveFormats, pType="trainVtest", standardise=std_, w1=W_train, w2=W_test)

    # Write the model
    modelName = "model_trained.h5"
    myModel.save(os.path.join(opts.saveDir,  modelName) )
        
    # Serialize model to JSON (contains arcitecture of model)
    model_json = myModel.to_json() # myModel.to_yaml()
    with open(opts.saveDir + "/model_architecture.json", "w") as json_file:
        json_file.write(model_json)
    # Serialize weights to HDF5
    myModel.save_weights(os.path.join(opts.saveDir, 'model_weights.h5'), overwrite=True)
    myModel.save(os.path.join(opts.saveDir, modelName))
        
    # Write weights and architecture in txt file
    modelFilename = os.path.join(opts.saveDir, "model.txt")
    Print("Writing the model (weights and architecture) in the file %s" % (hs + os.path.basename(modelFilename) + ns), True)
    func.WriteModel(myModel, model_json, opts.inputList, modelFilename, verbose=False)

    #https://keras.io/visualization/
    #https://machinelearningmastery.com/display-deep-learning-model-training-history-in-keras/
    molelFilename = modelFilename.replace(".txt", ".png")
    Print("plot a graph of the model and save it to the file %s" % (hs + os.path.basename(modelFilename) + ns), True)
    if "fnal" not in opts.hostname:
        from keras.utils import plot_model
        plot_model(myModel, to_file=modelFilename)


    # Produce method score (i.e. predict output value for given input dataset). Computation is done in batches.
    # https://stackoverflow.com/questions/49288199/batch-size-in-model-fit-and-model-predict
    Print("Generating output predictions (numpy arrays) for the input samples", True)
    pred_train  = myModel.predict(X_train, batch_size=None, verbose=1, steps=None) # DNN output for training data (for both signal & bkg)
    pred_test   = myModel.predict(X_test , batch_size=None, verbose=1, steps=None) # DNN output for test data (for both signal & bkg)
    pred_signal = myModel.predict(X_sig  , batch_size=None, verbose=1, steps=None) # DNN output for signal only (all data)
    pred_bkg    = myModel.predict(X_bkg  , batch_size=None, verbose=1, steps=None) # DNN output for data only (all data)

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
    pred_train_S =  myModel.predict(X_train_S, batch_size=None, verbose=1, steps=None)
    pred_train_B =  myModel.predict(X_train_B, batch_size=None, verbose=1, steps=None)
    pred_test_S  =  myModel.predict(X_test_S , batch_size=None, verbose=1, steps=None)
    pred_test_B  =  myModel.predict(X_test_B , batch_size=None, verbose=1, steps=None)

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

    # Scale back the data to the original representation
    if opts.scaleBack:
        msg = "Performing inverse-transform (%s) to all variables to get their original representation" % (hs + "scaleBack" + ns)
        Print(msg, True)
        df_sig   = func.GetOriginalDataFrame(scaler_sig, df_sig, opts.inputList)
        df_bkg   = func.GetOriginalDataFrame(scaler_bkg, df_bkg, opts.inputList)
        df_all   = func.GetOriginalDataFrame(scaler_all, df_all, opts.inputList)
        dset_all = df_all.values

    rocDict = {}
    WPs     = [float(x)/float(20) for x in range(0, 20, 1)]

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

    # Plot some metrics
    xErr = [0.0 for i in range(0, opts.epochs)]
    yErr = [0.0 for i in range(0, opts.epochs)]
    func.PlotTGraph(epochList[0:opts.epochs], xErr, trainLossList[0:opts.epochs-1], yErr , opts.saveDir, "TrainLoss"    , jsonWr, opts.saveFormats, **GetKwargs("loss") )
    func.PlotTGraph(epochList[0:opts.epochs], xErr, trainAccList[0:opts.epochs-1] , yErr , opts.saveDir, "TrainAccuracy", jsonWr, opts.saveFormats, **GetKwargs("acc") )
    func.PlotTGraph(epochList[0:opts.epochs], xErr, valLossList[0:opts.epochs-1]  , yErr , opts.saveDir, "ValLoss"      , jsonWr, opts.saveFormats, **GetKwargs("loss") )
    func.PlotTGraph(epochList[0:opts.epochs], xErr, valAccList[0:opts.epochs-1]   , yErr , opts.saveDir, "ValAccuracy"  , jsonWr, opts.saveFormats, **GetKwargs("acc") )
        
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
    #ROOTFILENAME = "/uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-TT_19var_6Jets_2BJets.root"  # 460 000
    ROOTFILENAME = "/uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-TT_19var_5Jets_1BJets.root"  # 875 000
    NOTBATCHMODE = False
    SAVEDIR      = None
    SAVEFORMATS  = "png"
    URL          = False
    STANDARDISE  = None
    SCALEBACK    = False
    DECORRELATE  = None
    ENTRYSTOP    = None
    VERBOSE      = False
    RNDSEED      = 1234
    EPOCHS       = 100
    BATCHSIZE    = None  # See: http://stats.stackexchange.com/questions/153531/ddg#153535
    ACTIVATION   = "relu" # "relu" or PReLU" or "LeakyReLU"
    NEURONS      = "36,19,1"
    LOG          = False
    GRIDX        = False
    GRIDY        = False
    LOSSFUNCTION = 'binary_crossentropy'
    OPTIMIZER    = 'adam'
    CFGJSON      = "config.json"
    RESULTSJSON  = "results.json"
    PLOTINPUTS   = True
    INPUTVARS    = "TrijetPtDR,TrijetDijetPtDR,TrijetBjetMass,TrijetLdgJetBDisc,TrijetSubldgJetBDisc,TrijetBJetLdgJetMass,TrijetBJetSubldgJetMass,TrijetMass,TrijetDijetMass,TrijetBJetBDisc,TrijetSoftDrop_n2,TrijetLdgJetCvsL,TrijetSubldgJetCvsL,TrijetLdgJetPtD,TrijetSubldgJetPtD,TrijetLdgJetAxis2,TrijetSubldgJetAxis2,TrijetLdgJetMult,TrijetSubldgJetMult"


    # Define the available script options
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("--notBatchMode", dest="notBatchMode", action="store_true", default=NOTBATCHMODE, 
                      help="Disable batch mode (opening of X window) [default: %s]" % NOTBATCHMODE)

    parser.add_option("--decorrelate", dest="decorrelate", default=DECORRELATE,
                      help="Calculate weights to decorrelate a variable from the training. This is done by reweighting the branches so that signal and background have similar mass distributions [default: %s]" % DECORRELATE)

    parser.add_option("--standardise", dest="standardise", default=STANDARDISE,
                      help="Standardizing a dataset involves rescaling the distribution of INPUT values so that the mean of observed values is 0 and the standard deviation is 1 (e.g. StandardScaler) [default: %s]" % STANDARDISE)

    parser.add_option("--scaleBack", dest="scaleBack", action="store_true", default=SCALEBACK,
                      help="Scale back the data to the original representation (before the standardisation). i.e. Performing inverse-transform to all variables to get their original representation. [default: %s]" % SCALEBACK)

    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=VERBOSE, 
                      help="Enable verbose mode (for debugging purposes mostly) [default: %s]" % VERBOSE)

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

    parser.add_option("--activation", dest="activation", type="string", default=ACTIVATION,
                      help="Type of transfer function that will be used to map the output of one layer to another [default: %s]" % ACTIVATION)

    parser.add_option("--neurons", dest="neurons", type="string", default=NEURONS,
                      help="List of neurons to use for each sequential layer (comma-separated integers)  [default: %s]" % NEURONS)

    parser.add_option("--inputVariables", dest="inputVariables", type="string", default=INPUTVARS,
                      help="List of input variables (TBranches) to use for the classifier (comma-separated integers)  [default: %s]" % INPUTVARS)
    
    parser.add_option("--entrystop", dest="entrystop", type="int", default=ENTRYSTOP,
                      help="Entry at which the (TBranch of TTree) reading stops. If ROOT file itoo big it may result to a \"memory error\", depending on the resources available in the machine used.  If \"None\", stop at the end of the branch. [default: %s]" % ENTRYSTOP)

    parser.add_option("--lossFunction", dest="lossFunction", type="string", default=LOSSFUNCTION,
                      help="One of the two parameters required to compile a model. The weights will take on values such that the loss function is minimized [default: %s]" % LOSSFUNCTION)

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
    #opts.saveFormats = ["." + s for s in opts.saveFormats]
    opts.saveFormats = [s for s in opts.saveFormats]
    Verbose("Will save all output in %d format(s): %s" % (len(opts.saveFormats), ss + ", ".join(opts.saveFormats) + ns), True)

    # Create specification lists
    if "," in opts.activation:
        opts.activation = opts.activation.split(",")
    else:
        opts.activation = [opts.activation]
    Verbose("Activation = %s" % (opts.activation), True)
    if "," in opts.neurons:
        opts.neurons = list(map(int, opts.neurons.split(",")) )
    else:
        opts.neurons = list(map(int, [opts.neurons]))
    Verbose("Neurons = %s" % (opts.neurons), True)
        
    # Sanity checks (One activation function for each layer)
    if len(opts.neurons) != len(opts.activation):
        msg = "The list of neurons (size=%d) is not the same size as the list of activation functions (=%d)" % (len(opts.neurons), len(opts.activation))
        raise Exception(es + msg + ns)  
    # Sanity check (Last layer)
    if opts.neurons[-1] != 1:
        msg = "The number of neurons for the last layer should be equal to 1 (=%d instead)" % (opts.neurons[-1])
        raise Exception(es + msg + ns)   
        #Print(es + msg + ns, True) 
    if opts.activation[-1] != "sigmoid":
        msg = "The activation function for the last layer should be set to \"sigmoid\"  (=%s instead)" % (opts.activation[-1])        
        #raise Exception(es + msg + ns)  #Print(es + msg + ns, True)
        Print(es + msg + ns, True)  #Print(es + msg + ns, True)

    # Define dir/logfile names
    specs = "%dLayers" % (len(opts.neurons))
    for i,n in enumerate(opts.neurons, 0):
        specs+= "_%s%s" % (opts.neurons[i], opts.activation[i])
    specs+= "_%sEpochs_%sBatchSize" % (opts.epochs, opts.batchSize)

    # Input list of discriminatin variables (TBranches)
    opts.inputList = opts.inputVariables.split(",")
    if opts.inputList < 1:
        raise Exception("At least one input variable needed to create the DNN. Only %d provided" % (len(opts.inputList)) )
    if opts.decorrelate != None:
        if opts.decorrelate not in opts.inputList:
            msg = "Cannot apply sample reweighting. The input variable \"%s\" is not in the inputList." % (opts.decorrelate)
            raise Exception(es + msg + ns)

    # Get the current date and time
    now    = datetime.now()
    nDay   = now.strftime("%d")
    nMonth = now.strftime("%h")
    nYear  = now.strftime("%Y")
    nTime  = now.strftime("%Hh%Mm%Ss") # w/ seconds
    nDate  = "%s%s%s_%s" % (nDay, nMonth, nYear, nTime)
    if opts.entrystop != None:
        sName  = "Keras_%s_%s" % (specs, str(opts.entrystop) + "Entrystop")
    if opts.decorrelate != None:
        sName += "_%sDecorrelated" % opts.decorrelate
    if opts.standardise != None:
        scalerTypes = ["Standard", "Robust", "MinMax"]
        if not opts.standardise in scalerTypes:
            msg = "Unsupported scaler type \"%s\". Please select one of the following (case-sensitive): %s" % (opts.standardise, ", ".join(scalerTypes) )
            raise Exception(es + msg + ns)
        #sName += "_Standardised"
        sName += "_%sScaler" % (opts.standardise)
        if opts.scaleBack:
            sName += "_ScaleBack"

    # Add the number of input variables
    sName += "_%dInputs" % (len(opts.inputList))
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
    for a in opts.activation:
        if a not in actList:
            msg = "Unsupported activation function %s. Please select on of the following:%s\n\t%s" % (opts.activation, ss, "\n\t".join(actList))
            raise Exception(es + msg + ns)    
    
    # See https://keras.io/losses/
    lossList = ["binary_crossentropy", "is_categorical_crossentropy",
                "mean_squared_error", "mean_absolute_error", "mean_absolute_percentage_error", "mean_squared_logarithmic_error", "squared_hinge",
                "hinge", "categorical_hinge", "logcosh", "huber_loss", "categorical_crossentropy", "sparse_categorical_crossentropy", 
                "kullback_leibler_divergenc", "poisson", "cosine_proximity"]
    bLossList = ["binary_crossentropy", "is_categorical_crossentropy"]
    # Sanity checks
    if opts.lossFunction not in lossList:
        msg = "Unsupported loss function %s. Please select on of the following:%s\n\t%s" % (opts.lossFunction, ss, "\n\t".join(lossList))
        raise Exception(es + msg + ns)    
    elif opts.lossFunction not in bLossList:
        msg = "Binary output currently only supports the following loss fucntions: %s" % ", ".join(bLossList)
        raise Exception(es + msg + ns)
    else:
        pass

    # See https://keras.io/optimizers/. Also https://www.dlology.com/blog/quick-notes-on-how-to-choose-optimizer-in-keras/
    optList = ["sgd", "rmsprop", "adagrad", "adadelta", "adam", "adamax", "nadam"]
    # optList = ["SGD", "RMSprop", "Adagrad", "Adadelta", "Adam", "Adamax", "Nadam"]
    if opts.optimizer not in optList:
        msg = "Unsupported loss function %s. Please select on of the following:%s\n\t%s" % (opts.optimizer, ss, "\n\t".join(optList))
        raise Exception(es + msg + ns)    

    # Determine and check the number of layers
    checkNumberOfLayers(opts)

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

    # Sanity check
    if not opts.standardise:
        opts.scaleBack = False

    main(opts)

    Print("Directory %s created" % (ls + opts.saveDir + ns), True)

    # Restore "stdout" to its original state and close the log file
    if opts.log:
        sys.stdout =  bak_stdout
        log_file.close()

    if opts.notBatchMode:
        raw_input("=== sequential.py: Press any key to quit ROOT ...")
