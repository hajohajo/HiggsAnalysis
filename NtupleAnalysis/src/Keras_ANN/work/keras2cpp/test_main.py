#!/usr/bin/env python
''' 
DESCRIPTION:
This script loads a keras model and predicts the NN output value of the input samples


EXAMPLE:
./test_main.py --filename ../histograms-TT_19var.root --dir ../Keras_4Layers_32relu_32relu_32relu_1sigmoid_500Epochs_32BatchSize_500000Entrystop_19Inputs_07Feb2020_11h29m27s 

LAST USED:
./test_main.py --filename ../histograms-TT_19var.root --dir ../Keras_4Layers_32relu_32relu_32relu_1sigmoid_500Epochs_32BatchSize_500000Entrystop_19Inputs_07Feb2020_11h29m27s 
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
from sklearn.externals import joblib
# Regression predictions
from keras.models import load_model
from optparse import OptionParser

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' #Disable AVX/FMA Warning
# Do not display canvases
ROOT.gROOT.SetBatch(ROOT.kTRUE)
# Disable screen output info
ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")

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
    #ROOT.gStyle.SetOptStat(0)

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

    # Signal and background dataframe (contain all TBranches (inputList) from the signal TTree)
    df_signal = signal.pandas.df(inputList, entrystop=opts.entries)
    df_background = background.pandas.df(inputList, entrystop=opts.entries)

    # Number of events to predict the output and plot the top-quark mass
    if opts.entries == None:
        if len(df_signal.index) > 0:
            opts.entries = min(len(df_signal.index), len(df_background.index))
        else:
            opts.entries = len(df_background.index)
        #Print("Number of events: %d" % (opts.entries), True)
    
    # Concatinate signal + background datasets
    df_all  = pandas.concat( [df_signal, df_background] )

    # Get a Numpy representation of the DataFrames for signal and background datasets
    dset_signal     = df_signal.values
    dset_background = df_background.values
    dataset = df_all.values
    
    X_bkg      = dset_background[:opts.entries, 0:nInputs]
    X_sig      = dset_signal[:opts.entries, 0:nInputs]
        
    # Standardization of datasets? 
    if opts.standardise != "None":
        msg  = "Standardising dataset features with the %sScaler" % (opts.standardise)
        Print(msg, True)
        scalerName = os.path.join(opts.dir, 'scaler.save')
        # Load the scaler files
        scaler = joblib.load(scalerName)
        
        #Transform inputs
        X_bkg = scaler.transform(X_bkg)
        X_sig = scaler.transform(X_sig)
        
    # Load & compile the model
    modelFile = os.path.join(opts.dir, 'model_trained.h5')
    Print("Loading model %s" % (modelFile), True)
    loaded_model = load_model(modelFile)
    loaded_model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['acc'])
    
    # Use the loaded model to generate output predictions for the input samples. Computation is done in batches.
    Print("Get the DNN score", True)
    Y_bkg = loaded_model.predict(X_bkg, verbose=0) 
    Y_sig = loaded_model.predict(X_sig, verbose=0) 
    
    # Create histograms
    h_sig = ROOT.TH1F("sig", '', 50, 0, 1)    
    h_bkg = ROOT.TH1F("bkg", '', 50, 0, 1)    

    #Print("predicted signal: ", False)
    for i in range(opts.entries):
        h_sig.Fill(Y_sig[i])
        #Print(Y_sig[i], False)
    #Print("predicted background: ", False)
    for i in range(opts.entries):
        h_bkg.Fill(Y_bkg[i])
        #Print(Y_bkg[i], False)

    outfile    = ROOT.TFile.Open("keras_py.root","RECREATE")
    outfile.cd()
    h_sig.Write()
    h_bkg.Write()
    outfile.Close()

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
    STANDARDISE = "None"
    SAVENAME    = None
    LOGY        = False
    ENTRIES     = None
    VERBOSE     = False

   # Define the available script options
    parser = OptionParser(usage="Usage: %prog [options]")

    parser.add_option("--filename", dest="filename", type="string", default=FILENAME,
                      help="Input ROOT file containing the signal and backbground TTrees with the various TBranches *variables) [default: %s]" % FILENAME)

    parser.add_option("--wp", dest="wp", type=float, default=WP,
                      help="Neural Network output working point [default: %s]" % WP)

    parser.add_option("--standardise", dest="standardise", default=STANDARDISE,
                      help="Standardizing a dataset involves rescaling the distribution of INPUT values so that the mean of observed values is 0 and the standard deviation is 1 (e.g. StandardScaler) [default: %s]" % STANDARDISE)

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

    main()