#!/usr/bin/env python
''' 
DESCRIPTION:
This script loads a set of mass-decorrelated top taggers with different lambda values and predicts their output. 
The output of the script, is the top-quark mass distribution of the selected top-candidates (truth-matched or 
unmatched) that pass a given working point. 


EXAMPLES:
./plotTopMass.py 
./plotMass.py --filename /uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-TT_19var.root --saveDir /publicweb/a/aattikis/Mass 
./plotMass.py --filename /uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-TT_19var.root --saveDir /publicweb/a/aattikis/TT --dir Keras_3Layers_50relu_50relu_1sigmoid_100Epochs_2000BatchSize_600000Entrystop_MassDecorrelated_27Nov2019_07h39m25s
./plotMass.py --filename /uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-TT_19var.root --dir Keras_3Layers_50relu_50relu_1sigmoid_100Epochs_2000BatchSize_600000Entrystop_MassDecorrelated_27Nov2019_07h39m25s --saveDir /publicweb/a/aattikis/ALEX --entries 100000
./plotMass.py --filename /uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-TT_19var_5Jets_1BJets.root --dir Keras_3Layers_50relu_50relu_1sigmoid_100Epochs_2000BatchSize_600000Entrystop_MassDecorrelated_27Nov2019_07h39m25s --saveDir /publicweb/a/aattikis/ALEX --entries 100000
./plotMass.py --filename /uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-TT_19var_6Jets_2BJets.root --dir Keras_3Layers_50relu_50relu_1sigmoid_100Epochs_2000BatchSize_600000Entrystop_MassDecorrelated_27Nov2019_07h39m25s --saveDir /publicweb/a/aattikis/ALEX --entries 100000
./plotMass.py --filename /uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-QCD_7Jets_3BJets.root --saveDir /publicweb/a/aattikis/Today --entries 500000 --dir Keras_3Layers_50relu_50relu_1sigmoid_10Epochs_2000BatchSize_600000Entrystop_MassDecorrelated_Standardised_27Nov2019_10h55m55s


LAST USED:
./plotMass.py --filename /uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-TT_19var_6Jets_2BJets.root --saveDir /publicweb/a/aattikis/Today --entries 500000 --dir Keras_3Layers_50relu_50relu_1sigmoid_10Epochs_2000BatchSize_600000Entrystop_MassDecorrelated_Standardised_27Nov2019_10h55m55s

./plotMass.py --filename /uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-QCD_7Jets_3BJets.root --saveDir /publicweb/a/aattikis/QCD --entries 200000 --dir Keras_3Layers_50relu_50relu_1sigmoid_10Epochs_2000BatchSize_600000Entrystop_MassDecorrelated_Standardised_27Nov2019_10h55m55s --wp 0.5


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

import plot
import tdrstyle
import func

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

def Verbose(msg, printHeader=True, verbose=False):
    if not opts.verbose:
        return
    Print(msg, printHeader)
    return

def main():
    ROOT.gStyle.SetOptStat(0)

    # Apply tdr style
    style = tdrstyle.TDRStyle() 
    style.setOptStat(False) 
    style.setGridX(False)
    style.setGridY(False)

    # Definitions
    filename  = opts.filename
    tfile     = ROOT.TFile.Open(filename)    
    
    #Signal and background branches
    signal     = uproot.open(filename)["treeS"]
    background = uproot.open(filename)["treeB"]
    
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
    nInputs = len(inputList)
    
    # Signal and background dataframes
    df_signal     = signal.pandas.df(inputList, entrystop=opts.entries)
    df_background = background.pandas.df(inputList, entrystop=opts.entries)

    # Number of events to predict the output and plot the top-quark mass
    if opts.entries == None:
        opts.entries = min(len(df_signal.index), len(df_background.index))
    Print("Number of events: %d" % (opts.entries), True)

    # Signal and background datasets
    dset_signal     = df_signal.values
    dset_background = df_background.values
    
    # Concatinate signal + background datasets
    df_list = [df_signal, df_background]
    df_all = pandas.concat(df_list)
    dataset = df_all.values
    
    # Target (top-quark mass) datasets (list of all mass values)
    if opts.dataset == "QCD":
        dset_target_sig = None
    else:
        dset_target_sig = signal.pandas.df(["TrijetMass"]).values
    dset_target_bkg = background.pandas.df(["TrijetMass"]).values
    dset_target_all = pandas.concat([signal.pandas.df(["TrijetMass"]), background.pandas.df(["TrijetMass"])]).values

    # Array of 19 columns (=# variables) and opts.entries)
    if opts.dataset == "QCD":
        X_sig = None
        target_sig = None
    else:
        X_sig = dset_signal[:opts.entries, 0:nInputs]
        target_sig = dset_target_sig[:opts.entries, :]

    X_bkg = dset_background[:opts.entries, 0:nInputs]
    target_bkg = dset_target_bkg[:opts.entries, :]
    
    # Create canvas
    colors = [ROOT.kAzure, ROOT.kOrange-2, ROOT.kMagenta, ROOT.kGreen+2, ROOT.kOrange+7, ROOT.kRed, ROOT.kRed, ROOT.kBlack]
    canvas = plot.CreateCanvas()
    canvas.cd()

    # Customise
    ymaxFactor = 1.1
    ymax = 0    
    histoList = []
    if (opts.logY):
        canvas.SetLogy()
        ymaxFactor = 2
        
    # load the models
    modelFile = os.path.join(opts.dir, 'model_trained.h5')
    Print("Loading model %s" % (hs + modelFile + ns), True)
    loaded_model = load_model(modelFile)

    # Compile the model
    loaded_model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['acc'])

    # Predict
    if opts.dataset == "QCD":
        Y_sig = None
    else:
        Y_sig = loaded_model.predict(X_sig, verbose=1)
    Y_bkg = loaded_model.predict(X_bkg, verbose=1)
    
    # Concatenate Y (predicted output) and target (top-quark mass)
    # Ymass 0 column:   DNN output
    # Ymass 1st column: top-quark mass
    if opts.dataset == "QCD":
        Ymass_sig = None
    else:
        Ymass_sig = numpy.concatenate((Y_sig, target_sig), axis=1)
    Ymass_bkg = numpy.concatenate((Y_bkg, target_bkg), axis=1)
    
    # Get selected top candidates (pass the output working point)
    if opts.dataset == "QCD":
        Ymass_sig_sel = None
        massSel_sig   = None
    else:
        Ymass_sig_sel = Ymass_sig[Ymass_sig[:,0] >= opts.wp] # Select samples with DNN output >= opts.wp (column 0)
        massSel_sig   = Ymass_sig_sel[:, 1]                  # Get the top-quark mass (col 1) for the selected samples.
    Ymass_bkg_sel = Ymass_bkg[Ymass_bkg[:,0] >= opts.wp] # Select samples with DNN output >= opts.wp (column 0)
    massSel_bkg   = Ymass_bkg_sel[:, 1]                  # Get the top-quark mass (col 1) for the selected samples.

    # Plot resutls
    nbins = 80 #100
    xmin  = 0
    xmax  = 800 #1000
    width = float(xmax)/nbins
    if opts.dataset == "QCD":
        hNameF = ''
    else:
        hNameF = 'non-matched'
    hNameT = 'truth-matched'
    histoT = ROOT.TH1F(hNameT, '', nbins, xmin, xmax)
    histoF = ROOT.TH1F(hNameF, '', nbins, xmin, xmax)
    Print("Selected entries: %d" % (len(massSel_bkg)), True)
    
    if opts.dataset != "QCD":
        for m in massSel_sig:
            histoT.Fill(m)
    for m in massSel_bkg:
        histoF.Fill(m)
    histoList.append(histoT)
    histoList.append(histoF)
    
    # Create legend
    leg = plot.CreateLegend(0.60, 0.70, 0.90, 0.90)
    leg.SetHeader(" DNN score #geq %s" % (opts.wp) )

    # For-loop: All histograms
    index = 0
    hList = []
    for i,h in enumerate(histoList, 0):

        if h.Integral() <= 0.0:
            continue
        else:
            index += 1

        Print("Plotting histogram #%d: %s" % (index, h.GetName()), index==1)

        # Normalize histograms to unity
        h.Scale(1.0/(h.Integral()))
        
        ymax = max(ymax, histoList[i].GetMaximum())
        h.GetXaxis().SetTitle("m_{top} (GeV)")
        h.GetYaxis().SetTitle("a.u. / %.0f GeV" % (width))
        plot.ApplyStyle(h, colors[i])
        if index==1 :
            h.SetFillColor(colors[i])
            h.SetFillStyle(1001)
            leg.AddEntry(h, "%s %s"  % (h.GetName(),  opts.dataset), "f")
            h.Draw("HIST")
        else:
            leg.AddEntry(h, "%s %s"  % (h.GetName(), opts.dataset), "l")
            h.Draw("HIST same")
        hList.append(h)
        
    # For-loop: All drawn Histograms
    for h in hList:
        h.SetMaximum(ymax*ymaxFactor)
    leg.Draw("same")
    
    # Additional text
    if 0:
        tex = plot.Text(" DNN score > %s" % (opts.wp), 0.82, 0.45)
        tex.Draw()
    
    # Draw line to indicate the real value of the top-quark mass
    graph = plot.CreateGraph([173., 173.], [0, ymax*ymaxFactor])
    graph.Draw("same")

    # Output directory
    dirName = plot.getDirName(opts.saveDir)

    # Save the histogram
    Print("Saving plot as %s" % (opts.saveName), True)
    plot.SavePlot(canvas, dirName, opts.saveName, saveFormats=opts.saveFormats, verbose=False)
    return

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
    SAVEDIR     = ""
    SAVEFORMATS = "pdf" #"png" does not work
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

    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=VERBOSE,
                      help="Enable verbose mode (for debugging purposes mostly) [default: %s]" % VERBOSE)

    parser.add_option("--entries", dest="entries", type=int, default=ENTRIES,
                      help="Number of entries to be used in filling the mass histogram [default: %s]" % ENTRIES)

    parser.add_option("--dir", dest="dir", default=None,
                      help="Directory where the model training file is located (\"model_trained.h5\") [default: %s]" % None)

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
    Verbose("Will save all output in %d format(s): %s" % (len(opts.saveFormats), ss + ", ".join(opts.saveFormats) + ns), True)

    if "TT" in os.path.basename(opts.filename):
        opts.dataset = "t#bar{t}"
    elif "QCD" in os.path.basename(opts.filename):
        opts.dataset = "QCD"
    else:
        opts.dataset = "dataset unknown"

    main()
