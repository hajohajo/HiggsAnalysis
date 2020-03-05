#!/bin/bash

# EXAMPLES:
# source run.csh ../Keras_4Layers_32relu_32relu_32relu_1sigmoid_500Epochs_32BatchSize_500000Entrystop_19Inputs_07Feb2020_11h29m27s ../histograms-TT_19var.root None
# source run.csh ../Keras_4Layers_32relu_32relu_32relu_1sigmoid_500Epochs_32BatchSize_500000Entrystop_RobustScaler_ScaleBack_19Inputs_09Feb2020_08h41m23s_hadronicRobust ../histograms-TT_19var.root Robust   
# Parameters
set FOLDER     = $1
set ROOTFILE   = $2
set SCALER     = $3

set KERASMODEL = "../../../EventSelection/interface/KerasModel.h"

echo 'Copying KerasModel.h from ../../../EventSelection/interface/'
cp $KERASMODEL KerasModel.h

echo ''
echo 'NN output prediction using keras predict'
./test_main.py --filename $ROOTFILE --dir $FOLDER --standardise $SCALER

echo ''
echo 'NN output prediction using c++ KerasModel predict'

root -l -b -q test_main.cc\(\"$FOLDER\",\"$ROOTFILE\"\)

echo ''
echo 'Comparing the methods'
./plotComparisons.py --sName $FOLDER --standardise $SCALER

# Clean
echo ''
echo 'Removing the unneeded files'
rm -f KerasModel.h
rm -f keras_py.root
rm -f keras_cpp.root


