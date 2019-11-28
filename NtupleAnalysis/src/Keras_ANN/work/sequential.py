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
./sequential.py --test --activation relu -s png,pdf,root,C
./sequential.py --test --activation relu,relu,sigmoid --neurons 36,19,1 -s pdf
./sequential.py --activation relu,relu,relu,relu,sigmoid --neurons 36,35,19,19,1 --epochs 10000 --batchSize 50000 -s pdf
./sequential.py --activation relu,relu,sigmoid --neurons 36,19,1 --epochs 100 --batchSize 500 -s pdf --saveDir Test
./sequential.py --activation relu,relu,sigmoid --neurons 36,19,1 --epochs 10000 --batchSize 50000 -s pdf --log
./sequential.py --activation relu,relu,sigmoid --neurons 19,190,1 --epochs 50 --batchSize 500 -s pdf --log --saveDir ~/public/html/Test
./sequential.py --activation relu,relu,sigmoid --neurons 19,190,1 --epochs 50 --batchSize 500 -s pdf --saveDir ~/public/html/Test
./sequential.py --activation relu,sigmoid --neurons 100,1 --epochs 100 --batchSize 2000 -s pdf --entrystop 600000
./sequential.py --activation relu,sigmoid --neurons 100,1 --epochs 100 --batchSize 2000 -s pdf --entrystop 600000 --decorrelate  
./sequential.py --activation relu,sigmoid --neurons 100,1 --epochs 100 --batchSize 2000 -s pdf --entrystop 600000 --standardise
./sequential.py --activation relu,sigmoid --neurons 100,1 --epochs 100 --batchSize 2000 -s pdf --entrystop 600000 --decorrelate --standardise


LAST USED:
./sequential.py --activation relu,relu,sigmoid --neurons 50,50,1 --epochs 200 --batchSize 2000 -s pdf --entrystop 600000
./sequential.py --activation relu,relu,sigmoid --neurons 50,50,1 --epochs 200 --batchSize 2000 -s pdf --entrystop 600000 --decorrelate  
./sequential.py --activation relu,relu,sigmoid --neurons 50,50,1 --epochs 200 --batchSize 2000 -s pdf --entrystop 600000 --standardise
./sequential.py --activation relu,relu,sigmoid --neurons 50,50,1 --epochs 200 --batchSize 2000 -s pdf --entrystop 600000 --decorrelate --standardise


GITHUB:
https://github.com/skonstantinou/Keras_ANN/
https://github.com/attikis/Keras_ANN


URL:
https://stats.stackexchange.com/questions/181/how-to-choose-the-number-of-hidden-layers-and-nodes-in-a-feedforward-neural-netw <----- Really good
https://towardsdatascience.com/activation-functions-neural-networks-1cbd9f8d91d6
https://theffork.com/activation-functions-in-neural-networks/?source=post_page-----4dfb9c7ce9c9----------------------
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
ts = "\033[0;35m"   
hs = "\033[1;34m"
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

def split_list(a_list, firstHalf=True):
    half = len(a_list)//2
    if firstHalf:
        return a_list[:half]
    else:
        return a_list[half:]

def GetModelWeightsAndBiases(inputList, neuronsList):
    '''
    Customly-defined function to estimate the model
    parameters (total, weights, biases)
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

def eval(model, x_test, y_test, opts):
    '''
    https://keras.io/models/model/#evaluate

    NOTE: computation is done in batches
    '''
    metrics = model.evaluate(
        x=x_test,
        y=y_test,
        batch_size=opts.batchSize,
        verbose=1,
        sample_weight=None
    )
    print "*"*100
    print ('Loss: {:.3f}, Accuracy: {:.3f}'.format(metrics[0], metrics[1]))
    print "model.metrics_names = ", model.metrics_names
    print "len(metrics) = ", len(metrics)
    print 
    print "*"*100
    return

def GetKwargs(var):

    kwargs  = {
        "normalizeToOne": True,
        "xTitle" : "DNN Output",
        "yTitle" : "a.u.",
        "xMin"   : 0.0,
        "xMax"   : 1.0,
        "nBins"  : 50,
        "log"    : True,
        }

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
        kwargs["nBins"]  =  80
        kwargs["xTitle"] = "p_{T}#DeltaR_{W}"
        kwargs["yMin"]   = 1e-3
        #kwargs["yTitle"] = "a.u. / %0.f"

    if var == "TrijetBjetMass":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   = 100.0
        kwargs["nBins"]  = 100
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
        kwargs["nBins"]  = 100
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
        kwargs["xMin"]   =    0.0
        kwargs["xMax"]   = 1000.0
        kwargs["nBins"]  =  200
        kwargs["xTitle"] = "m_{t} [GeV]"
        kwargs["yMin"]   = 1e-3
        #kwargs["yTitle"] = "a.u. / %0.0f GeV"

    if var == "TrijetDijetMass":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   = 500.0
        kwargs["nBins"]  = 100
        kwargs["xTitle"] = "m_{W} [GeV]"
        kwargs["yMin"]   = 1e-3
        #kwargs["yTitle"] = "a.u. / %0.0f GeV"

    if var == "TrijetBJetBDisc":
        kwargs["xMin"]   =  0.0
        kwargs["xMax"]   =  1.0
        kwargs["nBins"]  = 100
        kwargs["xTitle"] = "b-tag discr."
        #kwargs["yTitle"] = "a.u. / %0.2f"

    if var == "TrijetSoftDrop_n2":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   =   2.0
        kwargs["nBins"]  = 200
        kwargs["xTitle"] = "soft-drop n_{2}"
        kwargs["yMin"]   = 1e-3
        #kwargs["yTitle"] = "a.u. / %0.2f"

    if var == "TrijetLdgJetCvsL":
        kwargs["xMin"]   =  -1.0
        kwargs["xMax"]   =  +1.0
        kwargs["nBins"]  = 200*2
        kwargs["xTitle"] = "CvsL discr."
        kwargs["yMin"]   = 1e-3
        kwargs["yMax"]   = 1e-1
        #kwargs["yTitle"] = "a.u. / %0.2f"

    if var == "TrijetSubldgJetCvsL":
        kwargs["xMin"]   =  -1.0
        kwargs["xMax"]   =  +1.0
        kwargs["nBins"]  = 200*2
        kwargs["xTitle"] = "CvsL discr."
        kwargs["yMin"]   = 1e-3
        kwargs["yMax"]   = 1e-1
        #kwargs["yTitle"] = "a.u. / %0.2f"

    if var == "TrijetLdgJetPtD":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   =   1.0
        kwargs["nBins"]  = 100
        kwargs["xTitle"] = "p_{T}D"
        #kwargs["yTitle"] = "a.u. / %0.2f"

    if var == "TrijetSubldgJetPtD":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   =   1.0
        kwargs["nBins"]  = 100
        kwargs["xTitle"] = "p_{T}D"
        #kwargs["yTitle"] = "a.u. / %0.2f"

    if var == "TrijetLdgJetAxis2":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   =   0.16
        kwargs["nBins"]  = 160
        kwargs["xTitle"] = "axis2"
        #kwargs["yTitle"] = "a.u. / %0.3f"

    if var == "TrijetSubldgJetAxis2":
        kwargs["xMin"]   =   0.0
        kwargs["xMax"]   =   0.16
        kwargs["nBins"]  = 160
        kwargs["xTitle"] = "axis2"
        #kwargs["yTitle"] = "a.u. / %0.3f"

    if var == "TrijetLdgJetMult":
        kwargs["xMin"]  =  0.0
        kwargs["xMax"]   = 40.0
        kwargs["nBins"]  = 40
        kwargs["xTitle"] = "constituent multiplicity"
        #kwargs["yTitle"] = "a.u. / %0.0f"

    if var == "TrijetSubldgJetMult":
        kwargs["xMin"]   =  0.0
        kwargs["xMax"]   = 40.0
        kwargs["nBins"]  = 40
        kwargs["xTitle"] = "constituent multiplicity"
        #kwargs["yTitle@"] = "a.u. / %0.0f"

    kwargs["nBins"]  *= 2
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

def printDataset(myDset):
    Verbose("%sPrinting 1 instance of the Numpy representation of the DataFrame%s:" % (ns, es), True)

    # For-loop: All TBranch entries
    for vList in myDset:
        # For-loop: All input variable values for given entry
        for v in vList:
            Verbose(v, False)
        break
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
        #nIn  = opts.neurons[-1] #?
        #nOut = opts.neurons[i]  #?
        nN  = opts.neurons[i]
        nIn  = opts.neurons[0]
        nOut = opts.neurons[-1]
        dof  = (2*nsignal) * (nIn + nOut)     # 2*nsignal = number of samples in training data set
        nMax = (2*nsignal) / (2*(nIn + nOut))
        nMin = (2*nsignal) / (10*(nIn + nOut))
        #if nOut > nMax or nOut < nMin:
        #    msg = "The number of OUT neurons for layer #%d (IN=%d, OUT=%d) is not within the recommended range(%d-%d). This may result in over-fitting" % (i, nIn, nOut, nMin, nMax)
        #    Print(hs + msg + ns, True)
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
    jsonWr.addParameter("decorrelate mass",  opts.decorrelate)
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

    # Open the ROOT file
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
    Print("Opening the signal and background TTrees with uproot using ROOT file %s" % (ts + opts.rootFileName + ns), True)
    signal     = uproot.open(opts.rootFileName)["treeS"]
    background = uproot.open(opts.rootFileName)["treeB"]

    # Input list of discriminatin variables (TBranches)
    inputList = []
    inputList.append("TrijetPtDR")
    inputList.append("TrijetDijetPtDR")
    inputList.append("TrijetBjetMass")
    inputList.append("TrijetLdgJetBDisc")
    inputList.append("TrijetSubldgJetBDisc")
    inputList.append("TrijetBJetLdgJetMass")
    inputList.append("TrijetBJetSubldgJetMass")
    inputList.append("TrijetMass")
    inputList.append("TrijetDijetMass")
    inputList.append("TrijetBJetBDisc")
    inputList.append("TrijetSoftDrop_n2")
    inputList.append("TrijetLdgJetCvsL")
    inputList.append("TrijetSubldgJetCvsL")
    inputList.append("TrijetLdgJetPtD")
    inputList.append("TrijetSubldgJetPtD")
    inputList.append("TrijetLdgJetAxis2")
    inputList.append("TrijetSubldgJetAxis2")
    inputList.append("TrijetLdgJetMult")
    inputList.append("TrijetSubldgJetMult")
    opts.inputList = inputList
    nInputs = len(inputList)
    
    # Construct signal and background dataframes using a list of TBranches (a Dataframe is a two dimensional structure representing data in python)
    # https://uproot.readthedocs.io/en/latest/ttree-handling.html#id7
    Print("Constucting dataframes for signal and background with %d input variables:\n\t%s%s%s" % (nInputs, ss, "\n\t".join(inputList), ns), True)
    df_signal = signal.pandas.df(inputList, entrystop=opts.entrystop) # call an array-fetching method to fill a Pandas DataFrame.
    df_bkg    = background.pandas.df(inputList, entrystop=opts.entrystop)

    # Get the index (row labels) of the DataFrame.
    nsignal = len(df_signal.index)    
    nbkg    = len(df_bkg.index)
    Verbose("Signal has %s%d%s row labels. Background has %s%d%s row labels" % (ss, nsignal, ns, es, nbkg, ns), True)

    # Apply rule-of-thumb to prevent over-fitting
    checkNeuronsPerLayer(nsignal, opts)
    
    # Sanity check
    columns = list(df_signal.columns.values)
    Verbose("The signal columns are :\n\t%s%s%s" % (ss, "\n\t".join(columns), ns), True)    

    # Sanity check:
    if 0:
        Print(df_signal['TrijetMass'], True)
        Print(df_signal['TrijetMass'].values, True)
    
    # Optional Numpy array of weights for the training samples, used for weighting the loss function (during training only). 
    sampleWeights = None
    if opts.decorrelate:
        msg = "De-correlating the mass by applying weights"
        Print(ts + msg + ns, True)
        # Return evenly spaced numbers over a specified interval.
        minPt       =    0.0
        maxPt       = 1000.0
        massBinning = numpy.linspace(minPt, maxPt, 500)
        massRatio   = [s/b for s, b in zip(df_signal['TrijetMass'].values, df_bkg['TrijetMass'].values)] # unused. keep for reference
        # Use Mass as the variable to be decorrelated
        if 1:
            #massEvents  = df_signal['TrijetMass'].values + df_bkg['TrijetMass'].values # WRONG! It adds the two arrays! does not extend! Use concatenate instead!
            massEvents  = df_bkg['TrijetMass'].values
        else:
            massEvents  = split_list(df_bkg['TrijetMass'].values)
            
        # Digitize will return numbers from 1 to len(bins) depending on which bin the event belongs to. However it won't handle
        # the situation if value is over the maximum bin edge, so you'll want to clip your values (or adjust the binning) accordingly
        # In other words; get the bin indices (convert values to bin indices according to the massBinning variable)
        digitizedMass = digitize( numpy.clip(massEvents, minPt, maxPt-1.0), bins=massBinning, right=False ) # binIndices
        
        # These are the weights you can give to keras.fit() function in parameter "sample_weight". The idea is to reweight the braches to get 
        # both signal and bkg to have similar mass (decouple mass from learning) 
        if 1:
            sampleWeights  = compute_sample_weight('balanced', digitizedMass) # Flatten-out both signal and bkg mass
        else:
            signalWeights = numpy.ones(len( split_list(df_signal['TrijetMass'].values)))
            bkgWeights    = compute_sample_weight('balanced', digitizedMass) # Flatten-out both signal and bkg mass
            sampleWeights = numpy.concatenate( (signalWeights, bkgWeights), axis=0)

        # Flatten out only the bkg (and not the signal)
        #sampleWeights[isSignal==0]  = compute_sample_weight('balanced', digitizedMass[isSignal==0])

    # Get a Numpy representation of the DataFrames for signal and background datasets
    Verbose("Getting a numpy representation of the DataFrames for signal and background datasets", True)
    dset_signal = df_signal.values
    dset_bkg    = df_bkg.values
    # Alternative fix in case of "memory error"
    if 0: 
        df_bkg   = df_bkg.iloc[0:nsignal] 
        dset_bkg = df_bkg.values

    # Print datasets?
    printDataset(dset_signal)
    printDataset(dset_bkg)

    # Construct the pandas DataFrames (2D size-mutable tabular data structure with labeled axes i.e. rows and columns)
    Verbose("Constructing pandas DataFrames for signal and background", True)
    ds_signal = pandas.DataFrame(data=dset_signal,columns=inputList)
    ds_bkg    = pandas.DataFrame(data=dset_bkg,columns=inputList)

    # Construct pandas DataFrames (2D size-mutable tabular data structure with labeled axes i.e. rows and columns)
    Verbose("Constructing pandas DataFrames", True)
    df_signal = df_signal.assign(signal=1)
    df_bkg    = df_bkg.assign(signal=0)
    Verbose("Printing tabular data for signal:\n%s%s%s" % (ss, ds_signal,ns), True)
    Verbose("Printing tabular data for background:\n%s%s%s" % (ss, ds_bkg,ns), True)
    
    # Create dataframe lists
    df_list = [df_signal, df_bkg]
    df_all  = pandas.concat(df_list)

    # Get a Numpy representation of the DataFrames for signal and background datasets (again, and AFTER assigning signal and background)
    dset_signal = df_signal.values
    dset_bkg    = df_bkg.values
    dset_all    = df_all.values

    # Define keras model as a linear stack of layers. Add layers one at a time until we are happy with our network architecture.
    Verbose("Creating the sequential Keras model", True)
    model = Sequential()

    # The best network structure is found through a process of trial and error experimentation. Generally, you need a network large enough to capture the structure of the problem.
    # The Dense function defines each layer - how many neurons and mathematical function to use.
    for iLayer, n in enumerate(opts.neurons, 0):
        layer = "layer#%d" % (int(iLayer)+1)
        if iLayer == len(opts.neurons)-1:
            layer += " (output Layer)" # Layers of nodes between the input and output layers. There may be one or more of these layers.
        else:            
            layer += " (hidden layer)" # A layer of nodes that produce the output variables.
            
        Print("Adding %s, with %s%d neurons%s and activation function %s" % (ts + layer + ns, ls, n, ns, ls + opts.activation[iLayer] + ns), iLayer==0)
        if iLayer == 0:
            # Only first layer demands input_dim. For the rest it is implied.
            model.add( Dense(opts.neurons[iLayer], input_dim = nInputs) ) #, weights = [np.zeros([692, 50]), np.zeros(50)] OR bias_initializer='zeros',  or bias_initializer=initializers.Constant(0.1)
            model.add( Activation(opts.activation[iLayer]) )
            # model.add( Dense(opts.neurons[iLayer], input_dim = nInputs, activation = opts.activation[iLayer]) ) # her majerty Soti requested to break this into 2 lines
        else:
            if 0: #opts.neurons[iLayer] < nInputs:
                msg = "The number of  neurons (=%d) is less than the number of input variables (=%d). Please set the number of neurons to be at least the number of inputs!" % (opts.neurons[iLayer], nInputs)
                raise Exception(es + msg + ns)  
            model.add( Dense(opts.neurons[iLayer]))
            model.add( Activation(opts.activation[iLayer]) )
            # model.add( Dense(opts.neurons[iLayer], activation = opts.activation[iLayer]) ) 

    # Print a summary representation of your model. 
    Print("Printing model summary:", True)
    model.summary()
    
    # Get the number of parameters of the model
    nParams = model.count_params()
    msg = "Model has a total of %s%d%s parameters (neuron weights and biases) " % (hs, nParams, ns)
    Print(msg, True)
    total_count, trainable_count, non_trainable_count  = GetModelParams(model)
    nParams, nWeights, nBias = GetModelWeightsAndBiases(inputList, opts.neurons)
    opts.modelParams = total_count
    opts.modelParamsTrain = trainable_count
    opts.modelParamsNonTrainable = non_trainable_count
    opts.modelWeights = nWeights
    opts.modelBiases  = nBias

    # Get a dictionary containing the configuration of the model. The model can be reinstantiated from its config via: model = Model.from_config(config)
    if 0:
        config = model.get_config()
        for c in config:
            print "c = ", c
            print
    
    # Split data into input (X) and output (Y). (Note: dataset includes both signal and background sequentially)
    # Use only 2*nsignal rows => First nSignal rows is the signal. Another (equal) nSignal rows for the bkg. 
    # Therefore signal and bkg have exactly the same number
    X = dset_all[:2*nsignal,0:nInputs] # rows: 0 -> 2*signal, columns: 0 -> 19
    Y = dset_all[:2*nsignal,nInputs:]  # rows: 0 -> 2*signal, columns: 19 (isSignal = 0 or 1)
    X_signal     = dset_signal[:nsignal, 0:nInputs]
    X_background = dset_bkg[:nsignal, 0:nInputs]
    Print("Signal dataset has %s%d%s rows. Background dataset has %s%d%s rows" % (ss, len(X_signal), ns, es, len(X_background), ns), True)

    # Standardization of a dataset is a common requirement for many machine learning estimators: they might behave badly if the individual features do not
    # more or less look like standard normally distributed data (e.g. Gaussian with 0 mean and unit variance).
    # https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html
    if opts.standardise:
        msg = "Standardize features by removing the mean and scaling to unit variance"
        Print(hs + msg + ns, True)
        scaler = StandardScaler()

        # Compute the mean and std to be used for later scaling
        scaler.fit(X_signal)
        StandardScaler(copy=True, with_mean=True, with_std=True)
        Verbose("X_signal before standardisation: %s" % X_signal, True)

        # Fit to data, then transform it by performing standardization by centering and scaling (Mean and standard deviation were previously stored with StandardScaler())
        X_signal     = scaler.fit_transform(X_signal)
        X_background = scaler.fit_transform(X_background)
        Verbose("X_signal after standardisation: %s" % X_signal, True)

    # Plot the input variables
    if opts.plotInputs:
        Verbose("Plotting all %d input variables for signal and bacgkround" % (len(inputList)), True)
        for i, var in enumerate(inputList, 0):
            func.PlotInputs(dset_signal[:, i:i+1], dset_bkg[:, i:i+1], var, "%s/%s" % (opts.saveDir, "inputs"), opts.saveFormats)

    # Split the datasets (X= 19 inputs, Y=output variable). Test size 0.5 means half for training half for testing
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.5, random_state=opts.rndSeed, shuffle=True)
    opts.testSample  = len(X_test)
    opts.trainSample = len(X_train)

    # Early stop? Stop training when a monitored quantity has stopped improving.
    # Show patience of "50" epochs with a change in the loss function smaller than "min_delta" before stopping procedure
    # https://stackoverflow.com/questions/43906048/keras-early-stopping
    #earlystop = keras.callbacks.EarlyStopping(monitor='val_loss', min_delta=0.0, patience=10, verbose=1, mode='auto')
    earlystop = keras.callbacks.EarlyStopping(monitor='val_loss', min_delta=0.0001, patience=50)
    callbacks = [earlystop]

    # [Loss function is used to understand how well the network is working (compare predicted label with actual label via some function)]
    # Optimizer function is related to a function used to optimise the weights
    Print("Compiling the model with the loss function %s and optimizer %s " % (ls + opts.lossFunction + ns, ts + opts.optimizer + ns), True)
    model.compile(loss=opts.lossFunction, optimizer=opts.optimizer, metrics=['accuracy'])
    
    # Customise the optimiser settings?
    if 0: #opts.optimizer == "adam":  # does not work
        # Default parameters follow those provided in the original paper.
        # learning_rate: float >= 0. Learning rate.
        # beta_1: float, 0 < beta < 1. Generally close to 1.
        # beta_2: float, 0 < beta < 1. Generally close to 1.
        # amsgrad: boolean. Whether to apply the AMSGrad variant of this algorithm from the paper "On the Convergence of Adam and Beyond".
        keras.optimizers.Adam(learning_rate=0.001, beta_1=0.9, beta_2=0.999, amsgrad=False)

    # batch size equal to the training batch size (See https://machinelearningmastery.com/use-different-batch-sizes-training-predicting-python-keras/)
    if opts.batchSize == None:
        opts.batchSize = len(X_train)/2

    # Fit the model with our data
    # (An "epoch" is an arbitrary cutoff, generally defined as "one iteration of training on the whole dataset", 
    # used to separate training into distinct phases, which is useful for logging and periodic evaluation.)
    try:
        Verbose("Number of threads before fitting model is %s" % (ts + str(p.num_threads()) + ns), True)
        seqModel = model.fit(X_train,
                             Y_train,
                             validation_data=(X_test, Y_test),
                             epochs     = opts.epochs,    # a full pass over all of your training data
                             batch_size = opts.batchSize, # a set of N samples (https://stats.stackexchange.com/questions/153531/what-is-batch-size-in-neural-network)
                             shuffle    = False,
                             verbose    = 1, # 0=silent, 1=progress, 2=mention the number of epoch
                             callbacks  = callbacks,
                             sample_weight=sampleWeights
                             )
        Verbose("Number of threads after fitting model is %s" % (ts + str(p.num_threads()) + ns), True)

        # Retrieve  the training / validation loss / accuracy at each epoch
        trainLossList = seqModel.history['loss']
        trainAccList  = seqModel.history['acc']
        valLossList   = seqModel.history['val_loss']
        valAccList    = seqModel.history['val_acc']
        epochList     = range(opts.epochs)

    except KeyboardInterrupt: #(KeyboardInterrupt, SystemExit):
        msg = "Manually interrupted the training (keyboard interrupt)!"
        Print(es + msg + ns, True)

    # Write the model
    modelName = "model_trained.h5"
    model.save(os.path.join(opts.saveDir,  modelName) )
        
    # serialize model to JSON (contains arcitecture of model)
    model_json = model.to_json()
    with open(opts.saveDir + "/model_architecture.json", "w") as json_file:
        json_file.write(model_json)
        
    # serialize weights to HDF5
    model.save_weights(os.path.join(opts.saveDir, 'model_weights.h5'), overwrite=True)
    model.save(os.path.join(opts.saveDir, modelName))
        
    # write weights and architecture in txt file
    func.WriteModel(model, model_json, inputList, os.path.join(opts.saveDir, "model.txt") )

    # Produce method score (i.e. predict output value for given input dataset). Computation is done in batches.
    # https://stackoverflow.com/questions/49288199/batch-size-in-model-fit-and-model-predict
    Print("Generating output predictions for the input samples", True) # (e.g. Numpy array)
    pred_train  = model.predict(X_train     , batch_size=None, verbose=1, steps=None)
    pred_test   = model.predict(X_test      , batch_size=None, verbose=1, steps=None)
    pred_signal = model.predict(X_signal    , batch_size=None, verbose=1, steps=None)
    pred_bkg    = model.predict(X_background, batch_size=None, verbose=1, steps=None)    
    # pred_train = model.predict(x, batch_size=None, verbose=0, steps=None, callbacks=None, max_queue_size=10, workers=1, use_multiprocessing=False) # Keras version 2.2.5 or later (https://keras.io/models/model/)

    # Join a sequence of arrays (X and Y) along an existing axis (1). In other words, add the ouput variable (Y) to the input variables (X)
    XY_train = numpy.concatenate((X_train, Y_train), axis=1)
    XY_test  = numpy.concatenate((X_test , Y_test ), axis=1)

    # Pick events with output = 1
    Verbose("Select events/samples which have an output variable Y (last column) equal to 1 (i.e. prediction is combatible with signal)", True)
    x_train_S = XY_train[XY_train[:,nInputs] == 1]; x_train_S = x_train_S[:,0:nInputs]
    x_test_S  = XY_test[XY_test[:,nInputs] == 1];   x_test_S  = x_test_S[:,0:nInputs]

    Verbose("Select events/samples which have an output variable Y (last column) equal to 0 (i.e. prediction is NOT combatible with signal)", False)
    x_train_B = XY_train[XY_train[:,nInputs] == 0]; x_train_B = x_train_B[:,0:nInputs]
    x_test_B  = XY_test[XY_test[:,nInputs] == 0];   x_test_B  = x_test_B[:,0:nInputs]
    
    # Produce method score for signal (training and test) and background (training and test)
    pred_train_S =  model.predict(x_train_S, batch_size=None, verbose=1, steps=None)
    pred_train_B =  model.predict(x_train_B, batch_size=None, verbose=1, steps=None)
    pred_test_S  =  model.predict(x_test_S , batch_size=None, verbose=1, steps=None)
    pred_test_B  =  model.predict(x_test_B , batch_size=None, verbose=1, steps=None)

    # Inform user of early stop
    stopEpoch = earlystop.stopped_epoch
    if stopEpoch != 0 and stopEpoch < opts.epochs:
        msg = "Early stop occured after %d epochs!" % (stopEpoch)
        opts.epochs = stopEpoch
        Print(cs + msg + ns, True)

    # Create json files
    writeCfgFile(opts)
    writeGitFile(opts)
    jsonWr  = JsonWriter(saveDir=opts.saveDir, verbose=opts.verbose)

    # Plot selected output and save to JSON file for future use
    func.PlotAndWriteJSON(pred_signal , pred_bkg    , opts.saveDir, "Output"       , jsonWr, opts.saveFormats)
    func.PlotAndWriteJSON(pred_train  , pred_test   , opts.saveDir, "OutputPred"   , jsonWr, opts.saveFormats)
    func.PlotAndWriteJSON(pred_train_S, pred_train_B, opts.saveDir, "OutputTrain"  , jsonWr, opts.saveFormats)
    func.PlotAndWriteJSON(pred_test_S , pred_test_B , opts.saveDir, "OutputTest"   , jsonWr, opts.saveFormats)

    # For-loop: All branches to be plotted for various DNN score cuts (v. slow, especially for large number of "entrystop")
    for i, var in enumerate(inputList, 0):
        PrintFlushed("Variable %s (%d/%d)" % (var, i+1, len(inputList)), i==0)
        func.PlotAndWriteJSON(dset_signal[:, i:i+1], dset_bkg[:, i:i+1], opts.saveDir, var , jsonWr, opts.saveFormats, **GetKwargs(var))    

        #Verbose("Plotting variable %s for DNN score 0.1" % (var), True)
        #func.PlotAndWriteJSON_DNNscore(pred_signal, pred_bkg, 0.1, dset_signal[:, i:i+1], dset_bkg[:, i:i+1], opts.saveDir, var , jsonWr, opts.saveFormats, **GetKwargs(var))

        #Verbose("Plotting variable %s for DNN score 0.3" % (var), True)
        #func.PlotAndWriteJSON_DNNscore(pred_signal, pred_bkg, 0.3, dset_signal[:, i:i+1], dset_bkg[:, i:i+1], opts.saveDir, var , jsonWr, opts.saveFormats, **GetKwargs(var))

        Verbose("Plotting variable %s for DNN score 0.5" % (var), True)
        func.PlotAndWriteJSON_DNNscore(pred_signal, pred_bkg, 0.5, dset_signal[:, i:i+1], dset_bkg[:, i:i+1], opts.saveDir, var , jsonWr, opts.saveFormats, **GetKwargs(var))

        Verbose("Plotting variable %s for DNN score 0.7" % (var), True)
        func.PlotAndWriteJSON_DNNscore(pred_signal, pred_bkg, 0.7, dset_signal[:, i:i+1], dset_bkg[:, i:i+1], opts.saveDir, var , jsonWr, opts.saveFormats, **GetKwargs(var))

        Verbose("Plotting variable %s for DNN score 0.9" % (var), True)
        func.PlotAndWriteJSON_DNNscore(pred_signal, pred_bkg, 0.9, dset_signal[:, i:i+1], dset_bkg[:, i:i+1], opts.saveDir, var , jsonWr, opts.saveFormats, **GetKwargs(var))
    print

    # Plot overtraining test
    htrain_s, htest_s, htrain_b, htest_b = func.PlotOvertrainingTest(pred_train_S, pred_test_S, pred_train_B, pred_test_B, opts.saveDir, "OvertrainingTest", opts.saveFormats)

    # Plot summary plot (efficiency & singificance)
    func.PlotEfficiency(htest_s, htest_b, opts.saveDir, "Summary", opts.saveFormats)

    # Write efficiencies (signal and bkg)
    xVals_S, xErrs_S, effVals_S, effErrs_S  = func.GetEfficiency(htest_s)
    xVals_B, xErrs_B, effVals_B, effErrs_B  = func.GetEfficiency(htest_b)
    func.PlotTGraph(xVals_S, xErrs_S, effVals_S, effErrs_S, opts.saveDir, "EfficiencySig", jsonWr, opts.saveFormats)
    func.PlotTGraph(xVals_B, xErrs_B, effVals_B, effErrs_B, opts.saveDir, "EfficiencyBkg", jsonWr, opts.saveFormats)

    xVals, xErrs, sig_def, sig_alt = func.GetSignificance(htest_s, htest_b)
    func.PlotTGraph(xVals, xErrs, sig_def, effErrs_B, opts.saveDir, "SignificanceA", jsonWr, opts.saveFormats)
    func.PlotTGraph(xVals, xErrs, sig_alt, effErrs_B, opts.saveDir, "SignificanceB", jsonWr, opts.saveFormats)
    # Metrics
    xErr = [0.0 for i in range(0, opts.epochs)]
    yErr = [0.0 for i in range(0, opts.epochs)]
    func.PlotTGraph(epochList[0:opts.epochs], xErr, trainLossList[0:opts.epochs-1], yErr , opts.saveDir, "TrainLoss"    , jsonWr, opts.saveFormats, **GetKwargs("loss") )
    func.PlotTGraph(epochList[0:opts.epochs], xErr, trainAccList[0:opts.epochs-1] , yErr , opts.saveDir, "TrainAccuracy", jsonWr, opts.saveFormats, **GetKwargs("acc") )
    func.PlotTGraph(epochList[0:opts.epochs], xErr, valLossList[0:opts.epochs-1]  , yErr , opts.saveDir, "ValLoss"      , jsonWr, opts.saveFormats, **GetKwargs("loss") )
    func.PlotTGraph(epochList[0:opts.epochs], xErr, valAccList[0:opts.epochs-1]   , yErr , opts.saveDir, "ValAccuracy"  , jsonWr, opts.saveFormats, **GetKwargs("acc") )

    # Plot ROC curve
    gSig = func.GetROC(htest_s, htest_b)
    if 0:
        gBkg = func.GetROC(htest_b, htest_s)
        gDict = {"graph" : [gSig, gBkg], "name" : ["signal", "bkg"]}
    else:
        gDict = {"graph" : [gSig], "name" : [os.path.basename(opts.saveDir)]}
    style.setLogY(True)
    func.PlotROC(gDict, opts.saveDir, "ROC", opts.saveFormats)

    # Write the  resultsJSON file!
    jsonWr.write(opts.resultsJSON)
    
    # Print total time elapsed
    days, hours, mins, secs = GetTime(tStart)
    Print("Total elapsed time is %s days, %s hours, %s mins, %s secs" % (days[0], hours[0], mins[0], secs), True)
    
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
    STANDARDISE  = False
    DECORRELATE  = False
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
    PLOTINPUTS   = False

    # Define the available script options
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("--notBatchMode", dest="notBatchMode", action="store_true", default=NOTBATCHMODE, 
                      help="Disable batch mode (opening of X window) [default: %s]" % NOTBATCHMODE)

    parser.add_option("--decorrelate", dest="decorrelate", action="store_true", default=DECORRELATE,
                      help="Use the mass to calculate weights to decorrelated the mass from the training. This is done by reweighting the branches so that signal and background have similar mass distribution [default: %s]" % DECORRELATE)

    parser.add_option("--standardise", dest="standardise", action="store_true", default=STANDARDISE,
                      help="Standardizing a dataset involves rescaling the distribution of values so that the mean of observed values is 0 and the standard deviation is 1. This can be thought of as subtracting the mean value or centering the data. [default: %s]" % STANDARDISE)

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

    # Get the current date and time
    now    = datetime.now()
    nDay   = now.strftime("%d")
    nMonth = now.strftime("%h")
    nYear  = now.strftime("%Y")
    nTime  = now.strftime("%Hh%Mm%Ss") # w/ seconds
    nDate  = "%s%s%s_%s" % (nDay, nMonth, nYear, nTime)
    if opts.entrystop != None:
        sName  = "Keras_%s_%s" % (specs, str(opts.entrystop) + "Entrystop")
    if opts.decorrelate:
        sName += "_MassDecorrelated"
    if opts.standardise:
        sName += "_Standardised"
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

    # Inform user
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
    for r in table:
        Print(r, False)

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
    # Sanity checks
    if opts.optimizer not in optList:
        msg = "Unsupported loss function %s. Please select on of the following:%s\n\t%s" % (opts.optimizer, ss, "\n\t".join(optList))
        raise Exception(es + msg + ns)    

    # Determine number of layers
    opts.nLayers = len(opts.neurons)+1
    opts.nHiddenLayers = len(opts.neurons)-1
    # Sanity check
    if opts.nHiddenLayers < 1:
        msg = "The NN has %d hidden layers. It must have at least 1 hidden layer" % (opts.nHiddenLayers)
        raise Exception(es + msg + ns)    
    else:
        msg = "The NN has %s%d hidden layers%s (in addition to input and output layers)." % (hs, opts.nHiddenLayers, ns)
        Print(msg, True) 
        msg = "[NOTE: One hidden layer is sufficient for the large majority of problems. In general, the optimal size of the hidden layer is usually between the size of the input and size of the output layers.]"
        Verbose (msg, False)         

    # Get some basic information
    opts.keras      = keras.__version__
    opts.tensorflow = tf.__version__
    opts.hostname   = socket.gethostname()
    opts.python     = "%d.%d.%d" % (sys.version_info[0], sys.version_info[1], sys.version_info[2])
    opts.gitBranch  = subprocess.check_output(["git", "branch", "-a"])
    opts.gitStatus  = subprocess.check_output(["git", "status"])
    opts.gitDiff    = subprocess.check_output(["git", "diff"])

    # Call the main function
    Print("Hostname is %s" % (ls + opts.hostname + ns), True)
    Print("Using Keras %s" % (ss + opts.keras + ns), True)
    Print("Using Tensorflow %s" % (ss + tf.__version__ + ns), False)
    main(opts)

    Print("All output saved under directory %s" % (ls + opts.saveDir + ns), True)

    # Restore "stdout" to its original state and close the log file
    if opts.log:
        sys.stdout =  bak_stdout
        log_file.close()

    if opts.notBatchMode:
        raw_input("=== sequential.py: Press any key to quit ROOT ...")

#  LocalWords:  decorelate
