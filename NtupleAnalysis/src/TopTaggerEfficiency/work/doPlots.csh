#!/bin/csh
#================================
# LAST USED:
#================================
# Default Training:   DR < 0.3 && DPt/Pt < 0.32
# Other Training:     DR < 0.3
#
# BDT  0.40: ./doPlots.csh TopTaggerEfficiency_180717_081727_BDT_0p40_TopMassCut400_TrainingBJet40_DefaultTraining  TopTaggerEfficiency_180718_051601_BDT_0p40_TopMassCut400_TrainingBJet40_OtherTraining 0.40

#================================================================================================
# Ensure all script arguments are passed from command line
#================================================================================================
if ($#argv != 3) then
    echo "=== You must give exactly 3 arguments:"
    echo "1 = The pseudo-multicrab with the default training"
    echo "2 = The pseudo-multicrab with the other training"
    echo "3 = The MVA cut studied"
    echo "\n=== For example:"
    echo "./doPlots.csh TopTaggerEfficiency_180704_042603_BDT_0p40_DefaultTraining TopTaggerEfficiency_180704_072108_BDT_0p40_OtherTraining 0.40"
    exit 1
endif

#================================================================================================ 
# Define variables                                                                                                                   
#================================================================================================
set INITIAL = `echo $USER | cut -c1-1`
set PSEUDO_MCRAB_DIR_DEFAULT = ${1}
set PSEUDO_MCRAB_DIR_OTHER   = ${2}
set MVA_CUT = ${3}

echo "\n==== Running ttbar variations (shower scales, high pT radiation, mTop, parton shower, evtGen)"
./plot_EfficiencySystTop.py -m $PSEUDO_MCRAB_DIR_DEFAULT --type showerScales --mva $MVA_CUT
./plot_EfficiencySystTop.py -m $PSEUDO_MCRAB_DIR_DEFAULT --type highPtRadiation --mva $MVA_CUT
./plot_EfficiencySystTop.py -m $PSEUDO_MCRAB_DIR_DEFAULT --type mTop --mva $MVA_CUT
./plot_EfficiencySystTop.py -m $PSEUDO_MCRAB_DIR_DEFAULT --type partonShower --mva $MVA_CUT
./plot_EfficiencySystTop.py -m $PSEUDO_MCRAB_DIR_DEFAULT --type evtGen --mva $MVA_CUT

echo "\n==== Running default ttbar with different training (mathing definition)"
./plot_MatchingDefinition.py -m $PSEUDO_MCRAB_DIR_DEFAULT -s $PSEUDO_MCRAB_DIR_OTHER --mva $MVA_CUT

echo "\n==== Merging JSON files"
./mergeJSONs.py --mva $MVA_CUT

echo "\n==== Done"
