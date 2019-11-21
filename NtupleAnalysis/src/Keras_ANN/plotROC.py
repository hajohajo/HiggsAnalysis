#!/usr/bin/env python
''' 
DESCRIPTION:
This script loads a set of mass-decorrelated top taggers with different lambda values and predicts their output. 
The output of the script, is a comparison of the ROC curves of the given models (signal efficiency vs bkg efficiency)

EXAMPLES:
./plotROC.py
'''

import numpy
import pandas
import keras
import ROOT
import array
import math
import os
import uproot
from keras.models import Sequential
from keras.layers import Dense
from keras.wrappers.scikit_learn import KerasRegressor
from keras.optimizers import Adam
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
# Regression predictions
from sklearn.datasets import make_regression
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from keras.models import load_model
from sklearn.externals import joblib
from keras.models import model_from_json

import plot
import tdrstyle
import func

# Disable screen output info
ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' #Disable AVX/FMA Warning


def main():
    ROOT.gStyle.SetOptStat(0)

    style = tdrstyle.TDRStyle() 
    style.setOptStat(False) 
    style.setGridX(True)
    style.setGridY(True)

    # Definitions
    filename  = "histograms-TT_19var.root"
    debug   = 1
    nprint  = 100
    tfile    = ROOT.TFile.Open(filename)    
    
    #Signal and background branches
    signal     = uproot.open(filename)["treeS"]
    background = uproot.open(filename)["treeB"]
    
    # Input list
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
    
    #Signal and background dataframes
    df_signal     = signal.pandas.df(inputList)
    df_background = background.pandas.df(inputList)

    nEvts = len(df_signal.index)
    print "=== Number of signal events: ", nEvts
    
    #Signal and background datasets
    dset_signal     = df_signal.values
    dset_background = df_background.values
    
    # Concat signal, background datasets
    df_list = [df_signal, df_background]
    df_all = pandas.concat(df_list)
    dataset = df_all.values  
    dataset_target_all = pandas.concat([signal.pandas.df(["TrijetMass"]), background.pandas.df(["TrijetMass"])]).values
    dataset_target_bkg = background.pandas.df(["TrijetMass"]).values
    dataset_target_signal = signal.pandas.df(["TrijetMass"]).values

    X_signal     = dset_signal[:nEvts, 0:nInputs]
    X_background = dset_background[:nEvts, 0:nInputs]

    #nEvts = nEvts
    X = dset_background[:nEvts, 0:nInputs]
    
    
    colors = [ROOT.kBlue, ROOT.kMagenta, ROOT.kRed, ROOT.kOrange, ROOT.kYellow, ROOT.kGreen, ROOT.kCyan, ROOT.kViolet+5, ROOT.kPink+5, ROOT.kOrange+5, ROOT.kSpring+5, ROOT.kTeal+5]
    lines = [1, 2, 9, 3, 6, 7, 10]
    
    graphList = []
    nameList  = []

    lamValues = [0, 1, 5, 10, 20]#, 100, 500] # fixme! should be given as option
    for lam in lamValues:
        model_name = 'modelANN_32_Adam_0p0008_500_tanh_relu_msle_lam%s.h5' % (lam)
        if lam == 0:
            model_name = 'Model_relu_relu_relu_sigmoid.h5'
        model = load_model(model_name)
        print "Loading Model", model_name
        model.compile(loss='binary_crossentropy', optimizer='adam')
        Y_signal = model.predict(X_signal, verbose=0)
        Y_bkg    = model.predict(X_background, verbose=0)
        
        htrain_s, htest_s, htrain_b, htest_b = func.PlotOvertrainingTest(Y_signal, Y_signal, Y_bkg, Y_bkg, "plotROC", model_name, ["pdf"])
        graph = func.GetROC(htest_s, htest_b)
        graphList.append(graph.Clone(model_name))
        nameList.append('lam%s' % lam)
        print "graphList", len(graphList)
        graph.Delete()
        del htrain_s
        del htrain_b
        del htest_s
        del htest_b
    
    graph_roc = {"graph" : graphList, "name" : nameList}
    func.PlotROC(graph_roc, "TestROC", "plotROC", ["pdf"])

main()

# FIX:
# 1. Give the *h5 files as input
# 2. Give the lambda values as input
