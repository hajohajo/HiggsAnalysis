#!/bin/csh   

#================================================================================================
# Ensure all script arguments are passed from command line
#================================================================================================
if ($#argv < 1) then
    echo "=== You must give exactly 1 argument:"
    echo "1=PSEUDO_MCRAB_DIR"
    echo "\n=== For example:"
    echo "./doMass.csh <DIRS> <SAVE_DIR>"
    echo "./doMass.csh Keras_3Layers_50relu_50relu_1sigmoid_200Epochs_2000BatchSize_800000Entrystop_MassDecorrelated_Standardised_28Nov2019_07h40m57s /publicweb/a/aattikis/"
    echo
    exit 1
endif

#================================================================================================
# Define variables                                                                               
#================================================================================================
set INITIAL   = `echo $USER | cut -c1-1`
set DIR       = ${1}
if ("${2}" == "") then
    set SAVE_DIR = /publicweb/a/aattikis/${DIR}
else
    set SAVE_DIR = ${2}
endif
set SCALEBACK = "" #"--scaleBack"
set FORMATS   = pdf #pdf,C,png 
set VERBOSE   = "" #"--verbose"
#set ROOTFILE1 = /uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-TT_19var.root
set ROOTFILE1 = /uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-TT_19var_6Jets_2BJets.root
#set ROOTFILE1 = /uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-TT_19var_5Jets_1BJets.root
set ROOTFILE2 = /uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-QCD_7Jets_3BJets.root
set LOGY      = "--logY"
#set ENTRIES   = 50000
set ENTRIES   = 600000


#foreach WP (0.0 0.1)
foreach WP (0.0 0.1 0.3 0.5 0.7 0.9)
     ./plotMass.py -s ${FORMATS} --entries ${ENTRIES} --filename ${ROOTFILE1} --saveDir "${SAVE_DIR}" --dir ${DIR} --wp ${WP} --saveName "TopMass_TT_${WP}" ${VERBOSE} ${LOGY} ${SCALEBACK}
    ./plotMass.py -s ${FORMATS} --entries ${ENTRIES} --filename ${ROOTFILE2} --saveDir "${SAVE_DIR}" --dir ${DIR} --wp ${WP} --saveName "TopMass_QCD_${WP}" ${VERBOSE} ${SCALEBACK}
end    
