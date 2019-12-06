#!/bin/csh   

#================================================================================================
# Ensure all script arguments are passed from command line
#================================================================================================
if ($#argv != 2) then
    echo "=== You must give exactly 1 argument:"
    echo "1=PSEUDO_MCRAB_DIR"
    echo "\n=== For example:"
    echo "./doPlots.csh <DIRS> <SAVE_DIR>"
    echo "./doPlots.csh results /publicweb/a/aattikis/results"
    echo
    exit 1
endif

#================================================================================================
# Define variables                                                                               
#================================================================================================
set INITIAL  = `echo $USER | cut -c1-1`
set DIRS     = ${1}
set SAVE_DIR = ${2}
set FORMATS  = png #pdf,C,png 
set VERBOSE  = "" #--verbose

echo "=== 1"
./plotResults.py -s ${FORMATS} --plotType var-wp --saveDir ${SAVE_DIR} --dirs ${DIRS} --yMin 0.5e-3 --yMaxFactor 2.0 ${VERBOSE} #--boldText 
#echo "\n=== 2"
#./plotResults.py -s ${FORMATS} --plotType var --saveDir ${SAVE_DIR} --dirs ${DIRS} --yMin 0.5e-3 --yMaxFactor 2.0 ${VERBOSE} #--boldText
# echo "\n=== 3"
# ./plotResults.py -s ${FORMATS} --plotType output --logY --saveDir ${SAVE_DIR} --dirs ${DIRS} --cutLineX 0.4 --xMin 0.0 --xMax 1.0 --yMin 0.5e3 ${VERBOSE}
# echo "\n=== 4" 
# ./plotResults.py -s ${FORMATS} --plotType roc --logY --saveDir ${SAVE_DIR} --dirs ${DIRS} --cutLineX 0.93 --xMin 0.0 --refName BDTG ${VERBOSE}
# echo "\n=== 5"
# #./plotResults.py -s ${FORMATS} --plotType metric --saveDir ${SAVE_DIR} --xMin 0.0 --xMax 30.0 --dirs ${DIRS} --yMax 1.0
# ./plotResults.py -s ${FORMATS} --plotType metric --saveDir ${SAVE_DIR} --xMin 0.0 --dirs ${DIRS} --yMax 1.0 ${VERBOSE}
# echo "\n=== 6"
# ./plotResults.py -s ${FORMATS} --plotType significance --logY --saveDir ${SAVE_DIR} --dirs ${DIRS} --cutLineX 0.93 --xMin 0.0 ${VERBOSE} #--refName BDT 
# echo "\n=== 7"
# ./plotResults.py -s ${FORMATS} --plotType efficiency --logY --saveDir ${SAVE_DIR} --dirs ${DIRS} --cutLineX 0.93 --xMin 0.0 ${VERBOSE}