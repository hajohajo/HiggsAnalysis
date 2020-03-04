#!/bin/csh   

#================================================================================================
# Ensure all script arguments are passed from command line
#================================================================================================
if ($#argv != 1) then
    echo "=== You must give exactly 1 argument:"
    echo "1=PSEUDO_MCRAB_DIR"
    echo
    exit 1
endif

#================================================================================================
# Define variables                                                                               
#================================================================================================
set INITIAL = `echo $USER | cut -c1-1`
set PSEUDO_MCRAB_DIR = ${1}
#set FORMATS = "png,pdf,C"
set FORMATS = "png"
set DSETS = "QCD_bEnriched|Charged" #ttHJetToGG|ttHJetToTT"
#set DSETS = " QCD_bEnriched|ttHJetToGG_M125|ttHJetToNonbb_M125|ttHJetToTT_M125"
#set DSETS = " QCD_bEnriched|ttHJetToGG_M125|ttHJetToNonbb_M125|ttHJetToTT_M125|ttHJetTobb_M125"

#================================================================================================
# TH1
#================================================================================================
./plotDataMC.py --folder counters/weighted -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
./plotDataMC.py --folder "" -e $DSETS -s $FORMATS -m $PSEUDO_MCRAB_DIR --ratio

./plotDataMC.py --folder PUDependency -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
./plotDataMC.py --folder eSelection_Veto -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
./plotDataMC.py --folder muSelection_ -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
###./plotDataMC.py --folder metSelection_ -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
./plotDataMC.py --folder tauSelection_ -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
./plotDataMC.py --folder jetSelection_ -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
./plotDataMC.py --folder bjetSelection_ -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
./plotDataMC.py --folder topSelectionBDT_ -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
### ./plotDataMC.py --folder AngularCuts_Collinear -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
### ./plotDataMC.py --folder AngularCuts_BackToBack -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
 
./plotDataMC.py --folder ForDataDrivenCtrlPlots -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
./plotDataMC.py --folder ForDataDrivenCtrlPlotsEWKFakeTaus -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
./plotDataMC.py --folder ForDataDrivenCtrlPlotsEWKGenuineTaus -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR

### ./plotDataMC.py --folder ForDataDrivenCtrlPlotsEWKFakeB -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
### ./plotDataMC.py --folder ForDataDrivenCtrlPlotsEWKGenuineB -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR


#./plotMC_ValueVsMass.py --logY --yMin 1e-3 --yMax 1.2 --refCounter "passed trigger" -m $PSEUDO_MCRAB_DIR #-s $FORMATS

#================================================================================================
# TH2
#================================================================================================
./plotTH2.py --folder ForDataDrivenCtrlPlots --gridX --gridY --dataset TT --normalizeToOne --logZ -s png -m $PSEUDO_MCRAB_DIR
./plotTH2.py --folder ForDataDrivenCtrlPlots --gridX --gridY --dataset ChargedHiggs_HplusTB_HplusToHW_M300_mH200_2ta_NLO --normalizeToOne --logZ -s png -m $PSEUDO_MCRAB_DIR
./plotTH2.py --folder ForDataDrivenCtrlPlots --gridX --gridY --dataset ChargedHiggs_HplusTB_HplusToHW_M700_mH200_2ta_NLO --normalizeToOne --logZ -s png -m $PSEUDO_MCRAB_DIR

./plotTH2.py --folder AngularCuts_BackToBack --gridX --gridY --dataset TT --normalizeToOne --logZ -s png -m $PSEUDO_MCRAB_DIR
./plotTH2.py --folder AngularCuts_BackToBack --gridX --gridY --dataset ChargedHiggs_HplusTB_HplusToHW_M300_mH200_2ta_NLO --normalizeToOne --logZ -s png -m $PSEUDO_MCRAB_DIR
./plotTH2.py --folder AngularCuts_BackToBack --gridX --gridY --dataset ChargedHiggs_HplusTB_HplusToHW_M700_mH200_2ta_NLO --normalizeToOne --logZ -s png -m $PSEUDO_MCRAB_DIR
###./plotTH2.py --folder AngularCuts_BackToBack --gridX --gridY --dataset QCD --normalizeToOne --logZ -s png -m $PSEUDO_MCRAB_DIR
###./plotTH2.py --folder AngularCuts_BackToBack --gridX --gridY --dataset QCD-b --normalizeToOne --logZ -s png -m $PSEUDO_MCRAB_DIR

### Identical to AngularCuts_BackToBack folder
# ./plotTH2.py --folder AngularCuts_Collinear --gridX --gridY --dataset TT --normalizeToOne --logZ -s png -m $PSEUDO_MCRAB_DIR
# ./plotTH2.py --folder AngularCuts_Collinear --gridX --gridY --dataset ChargedHiggs_HplusTB_HplusToHW_M300_mH200_2ta_NLO --normalizeToOne --logZ -s png -m $PSEUDO_MCRAB_DIR
# ./plotTH2.py --folder AngularCuts_Collinear --gridX --gridY --dataset ChargedHiggs_HplusTB_HplusToHW_M700_mH200_2ta_NLO --normalizeToOne --logZ -s png -m $PSEUDO_MCRAB_DIR
