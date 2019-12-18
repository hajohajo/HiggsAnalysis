#!/bin/csh   

#================================================================================================
# Ensure all script arguments are passed from command line
#================================================================================================
if ($#argv < 1) then
    echo "=== You must give at least 1 argument ($#argv provided):"
    echo "1=PSEUDO_MCRAB_DIR"
    echo "\n=== For example:"
    echo "./doPlots.csh <DIRS> <SAVE_DIR>"
    echo "./doPlots.csh results /publicweb/a/aattikis/results" # results contains multiple Keras_* directories
    echo
    exit 1
endif

#================================================================================================
# Define variables                                                                               
#================================================================================================
set INITIAL  = `echo $USER | cut -c1-1`
set DIRS     = ${1}
if ("${2}" == "") then
    set SAVE_DIR = /publicweb/a/aattikis/${DIRS}
else
    set SAVE_DIR = ${2}
endif

set FORMATS  = png #pdf,C,png 
set VERBOSE  = "" #--verbose

# echo "=== Variables (WPs)"
# ./plotResults.py -s ${FORMATS} --plotType var-wp --saveDir ${SAVE_DIR} --dirs ${DIRS} --yMin 0.5e-3 --yMaxFactor 2.0 ${VERBOSE} #--boldText 
# 
# echo "\n=== Variables"
# ./plotResults.py -s ${FORMATS} --plotType var --saveDir ${SAVE_DIR} --dirs ${DIRS} --yMin 0.5e-3 --yMaxFactor 2.0 ${VERBOSE} #--boldText
#  
# echo "\n=== Output"
# #./plotResults.py -s ${FORMATS} --plotType output --logY --saveDir ${SAVE_DIR} --dirs ${DIRS} --cutLineX 0.4 --xMin 0.0 --xMax 1.0 --yMin 0.5e3 ${VERBOSE}
# ./plotResults.py -s ${FORMATS} --plotType output --logY --saveDir ${SAVE_DIR} --dirs ${DIRS} --cutLineX 0.4 --xMin 0.0 --xMax 1.0 --yMin 0.5e3 ${VERBOSE}

echo "\n=== ROC" 
./plotResults.py -s ${FORMATS} --plotType roc --logY --saveDir ${SAVE_DIR} --dirs ${DIRS},Keras_BDTG2018_Aug2018  --cutLineX 0.65 --xMin 0.0 --refName BDTG ${VERBOSE}

# echo "\n=== Metrics (Loss Vs Epoch, Accuracy Vs Epoch, etc..)"
# ./plotResults.py -s ${FORMATS} --plotType metric --saveDir ${SAVE_DIR} --xMin 0.0 --dirs ${DIRS} --yMax 1.0 ${VERBOSE}

#echo "\n=== Significance"
#./plotResults.py -s ${FORMATS} --plotType significance --logY --saveDir ${SAVE_DIR} --dirs ${DIRS} --cutLineX 0.93 --xMin 0.0 ${VERBOSE} #--refName BDT 

#echo "\n=== Efficiency"
#./plotResults.py -s ${FORMATS} --plotType efficiency --logY --saveDir ${SAVE_DIR} --dirs ${DIRS} --cutLineX 0.93 --xMin 0.0 ${VERBOSE}
