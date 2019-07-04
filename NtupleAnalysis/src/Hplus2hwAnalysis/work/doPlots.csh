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
set FORMATS = "png,pdf,C"
#set FORMATS = "png"
set DSETS = "#ZJetsToQQ_HT600toInf"

./plotDataMC_ControlPlots.py --folder "" -e $DSETS -s $FORMATS -m $PSEUDO_MCRAB_DIR --ratio
 ./plotDataMC_ControlPlots.py --folder ForDataDrivenCtrlPlotsEWKFakeB -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
 ./plotDataMC_ControlPlots.py --folder ForDataDrivenCtrlPlotsEWKGenuineB -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
 ./plotDataMC_ControlPlots.py --folder ForDataDrivenCtrlPlots -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
 ./plotDataMC_ControlPlots.py --folder AngularCuts_Collinear -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
 ./plotDataMC_ControlPlots.py --folder AngularCuts_BackToBack -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
 ./plotDataMC_ControlPlots.py --folder PUDependency -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
./plotDataMC_ControlPlots.py --folder counters/weighted -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
 ./plotDataMC_ControlPlots.py --folder eSelection_Veto -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
 ./plotDataMC_ControlPlots.py --folder muSelection_ -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
 ./plotDataMC_ControlPlots.py --folder tauSelection_ -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
 ./plotDataMC_ControlPlots.py --folder bjetSelection_ -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
 ./plotDataMC_ControlPlots.py --folder jetSelection_ -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
 ./plotDataMC_ControlPlots.py --folder topSelectionBDT_ -e $DSETS --ratio -s $FORMATS -m $PSEUDO_MCRAB_DIR
