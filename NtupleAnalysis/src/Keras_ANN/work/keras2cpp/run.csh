#!/bin/bash

# Parameters
set FOLDER     = $1
set ROOTFILE   = $2
set KERASMODEL = "../../../EventSelection/interface/KerasModel.h"

echo $FOLDER
echo $ROOTFILE
echo $KERASMODEL

echo 'Copying KerasModel.h from ../../../EventSelection/interface/'
cp $KERASMODEL KerasModel.h

echo ''
echo 'NN output prediction using keras predict'
./test_main.py --filename $ROOTFILE --dir $FOLDER 

echo ''
echo 'NN output prediction using c++ KerasModel predict'

root -l -b -q test_main.cc\(\"$FOLDER\",\"$ROOTFILE\"\)

echo ''
# Clean
echo 'Remove unneeded files'
rm -f KerasModel.h
#rm -f $TEST_BIN

