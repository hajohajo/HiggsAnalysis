#!/bin/csh   

#================================================================================================
# Ensure all script arguments are passed from command line
#================================================================================================
if ($#argv != 2) then
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
set INITIAL  = `echo $USER | cut -c1-1`
set DIR      = ${1}
set SAVE_DIR = ${2}
set FORMATS  = pdf #pdf,C,png 

foreach WP (0.0 0.2 0.4 0.6 0.8 0.9)
    echo "=== TT (WP = ${WP})"
    ./plotMass.py -s ${FORMATS} --entries 600000 --filename /uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-TT_19var_6Jets_2BJets.root --saveDir "${SAVE_DIR}/${DIR}" --dir ${DIR} --wp ${WP} --saveName "TopMass_TT_${WP}"

    echo "=== QCD (WP = ${WP})"
    ./plotMass.py -s ${FORMATS} --entries 600000 --filename /uscms_data/d3/aattikis/workspace/pseudomulticrab/Keras/TopTagger/histograms-QCD_7Jets_3BJets.root --saveDir "${SAVE_DIR}/${DIR}" --dir ${DIR} --wp ${WP} --saveName "TopMass_QCD_${WP}"
end    
