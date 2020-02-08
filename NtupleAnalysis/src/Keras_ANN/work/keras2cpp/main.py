#!/usr/bin/env python
''' 
DESCRIPTION:
This script loads a keras model and predicts the NN output value of the input samples


EXAMPLE:
./main.py --filename ../histograms-TT_19var.root --dir ../Keras_4Layers_32relu_32relu_32relu_1sigmoid_500Epochs_32BatchSize_500000Entrystop_19Inputs_07Feb2020_11h29m27s 

LAST USED:
./main.py --filename ../histograms-TT_19var.root --dir ../Keras_4Layers_32relu_32relu_32relu_1sigmoid_500Epochs_32BatchSize_500000Entrystop_19Inputs_07Feb2020_11h29m27s 
'''
#================================================================================================ 
# Import modules
#================================================================================================ 
import numpy
import pandas
import keras
import ROOT
import array
import math
import json
import sys
import os
import uproot
from keras.models import Sequential
from keras.layers import Dense, Dropout #https://machinelearningmastery.com/dropout-regularization-deep-learning-models-keras/
from keras.wrappers.scikit_learn import KerasRegressor
from keras.optimizers import Adam
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
# Regression predictions
from sklearn.datasets import make_regression
from sklearn.model_selection import train_test_split
from keras.models import load_model
from keras.models import model_from_json
from optparse import OptionParser

#import plot
#import tdrstyle
#import func

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' #Disable AVX/FMA Warning
# Do not display canvases
ROOT.gROOT.SetBatch(ROOT.kTRUE)
# Disable screen output info
ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")

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
# Function definition
#================================================================================================ 
def Print(msg, printHeader=False):
    fName = __file__.split("/")[-1]
    if printHeader==True:
        print "=== ", fName
        print "\t", msg
    else:
        print "\t", msg
    return

def main():
    ROOT.gStyle.SetOptStat(0)

    # Definitions
    filename = opts.filename
    Print("Opening ROOT file %s" %  (opts.filename), True)
    tfile    = ROOT.TFile.Open(filename)    
    sigTree  = "treeS"
    bkgTree  = "treeB"
    Nevts = 20

    #Signal and background branches
    signal     = uproot.open(filename)[sigTree]
    background = uproot.open(filename)[bkgTree]
    
    # Input list of the model to be used
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

    # Signal dataframe
    df_signal = signal.pandas.df(inputList, entrystop=opts.entries) # contains all TBranches (inputList) from the signal TTree 
    # Background dataframe
    df_background = background.pandas.df(inputList, entrystop=opts.entries) # contains all TBranches (inputList) from the background TTree 

    # Number of events to predict the output and plot the top-quark mass
    if opts.entries == None:
        if len(df_signal.index) > 0:
            opts.entries = min(len(df_signal.index), len(df_background.index))
        else:
            opts.entries = len(df_background.index)
        #Print("Number of events: %d" % (opts.entries), True)

    
    # Concatinate signal + background datasets
    df_all  = pandas.concat( [df_signal, df_background] )
    
    # Standardization of datasets? 
    if opts.standardise != None:
        msg  = "Standardising dataset features with the %sScaler%s" % (ls + opts.standardise, ns)
        #Print(msg, True)    
        '''
        if not df_signal.empty: 
            scaler_sig, df_signal = func.GetStandardisedDataFrame(df_signal, opts.inputList, scalerType=opts.standardise)
        if not df_background.empty: 
            scaler_bkg, df_background = func.GetStandardisedDataFrame(df_background, opts.inputList, scalerType=opts.standardise)
        scaler_all, df_all = func.GetStandardisedDataFrame(df_all, opts.inputList, scalerType=opts.standardise)
        '''
    # Get a Numpy representation of the DataFrames for signal and background datasets
    dset_signal     = df_signal.values
    dset_background = df_background.values
    dataset = df_all.values
    
    X_bkg      = dset_background[:opts.entries, 0:nInputs]
    X_sig      = dset_signal[:opts.entries, 0:nInputs]
    
    # Load & compile the model
    modelFile = os.path.join(opts.dir, 'model_trained.h5')
    Print("Loading model %s" % (ts + modelFile + ns), True)
    loaded_model = load_model(modelFile)
    loaded_model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['acc'])

    
    # Use the loaded model to generate output predictions for the input samples. Computation is done in batches.

    Print("Get the DNN score", True)
    Y_bkg = loaded_model.predict(X_bkg, verbose=0) 
    Y_sig = loaded_model.predict(X_sig, verbose=0) 
    
    Print("predicted signal: ", False)
    for i in range(Nevts):
        Print(Y_sig[i], False)
    Print("predicted background: ", False)
    for i in range(Nevts):
        Print(Y_bkg[i], False)

#================================================================================================
# Main
#================================================================================================
if __name__ == "__main__":
    '''
    https://docs.python.org/3/library/argparse. html

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
    FILENAME    = "histograms-TT_19var.root"
    WP          = 0.0
    DIR         = None
    SAVEDIR     = ""
    SAVEFORMATS = "pdf" #"png" does not work
    STANDARDISE = None
    SAVENAME    = None
    LOGY        = False
    ENTRIES     = None
    VERBOSE     = False
    SCALEBACK   = True #False

   # Define the available script options
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("--filename", dest="filename", type="string", default=FILENAME,
                      help="Input ROOT file containing the signal and backbground TTrees with the various TBranches *variables) [default: %s]" % FILENAME)

    parser.add_option("--wp", dest="wp", type=float, default=WP,
                      help="Neural Network output working point [default: %s]" % WP)

    parser.add_option("--standardise", dest="standardise", default=STANDARDISE,
                      help="Standardizing a dataset involves rescaling the distribution of INPUT values so that the mean of observed values is 0 and the standard deviation is 1 (e.g. StandardScaler) [default: %s]" % STANDARDISE)

    parser.add_option("--scaleBack", dest="scaleBack", action="store_true", default=SCALEBACK,
                      help="Scale back the data to the original representation (before the standardisation). i.e. Performing inverse-transform to all variables to get their original representation. [default: %s]" % SCALEBACK)

    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=VERBOSE,
                      help="Enable verbose mode (for debugging purposes mostly) [default: %s]" % VERBOSE)

    parser.add_option("--entries", dest="entries", type=int, default=ENTRIES,
                      help="Number of entries to be used in filling the mass histogram [default: %s]" % ENTRIES)

    parser.add_option("--dir", dest="dir", default=DIR,
                      help="Directory where the model training file is located (\"model_trained.h5\") [default: %s]" % DIR)

    parser.add_option("--saveDir", dest="saveDir", type="string", default=SAVEDIR,
                      help="Output name directory [default: %s]" % SAVEDIR)

    parser.add_option("-s", "--saveFormats", dest="saveFormats", default = SAVEFORMATS,
                      help="Save formats for all plots [default: %s]" % SAVEFORMATS)

    parser.add_option("--saveName", dest="saveName", type="string", default=SAVENAME,
                      help="Name of the output histogram")

    parser.add_option("--logY", dest="logY", action="store_true", default=LOGY,
                      help="Set logarithmic y-axis [default: %s]" % LOGY)

    (opts, parseArgs) = parser.parse_args()

    if opts.saveName == None:
        opts.saveName = "TopMass_DNN%s" % (str(opts.wp).replace(".", "p"))

    # Create save formats
    if "," in opts.saveFormats:
        opts.saveFormats = opts.saveFormats.split(",")
    else:
        opts.saveFormats = [opts.saveFormats]
    opts.saveFormats = [s for s in opts.saveFormats]

    if "TT" in os.path.basename(opts.filename):
        opts.dataset = "t#bar{t}"
    elif "QCD" in os.path.basename(opts.filename):
        opts.dataset = "QCD"
    else:
        opts.dataset = "dataset unknown"
    
    # Sanity check
    if opts.dir == None:
        raise Exception("No directory was defined where the model training file is located (\"model_trained.h5\")! Please define an input directory with the --dir option")

    # is the input standardised?
    cfgFile = os.path.join(opts.dir, "config.json")
    if os.path.exists(cfgFile):
        f = open(cfgFile, "r")
        config = json.load(f)
        f.close()
        if "standardised datasets" in config:
            if config["standardised datasets"] != "None":
                opts.standardise = (config["standardised datasets"])

    # Sanity checks 
    if opts.standardise == None:
        opts.scaleBack = False
    else:
        Print("Must standardise the input variables before applying the the deployed network", True)
        if opts.scaleBack:
            msg = "Performing inverse-transform (%s) to all variables to get their original representation" % (hs + "scaleBack" + ns)
            Print(msg, True)
            opts.scaleBack = True

    main()
