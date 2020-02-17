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
set DSETS1 = "QCD|WJets"
set DSETS = ""

#================================================================================================
# TH1
#================================================================================================
./plotDataMC.py --folder counters/weighted -e $DSETS -s $FORMATS -m $PSEUDO_MCRAB_DIR --ratio
./plotDataMC.py --folder "" -e $DSETS1 -s $FORMATS -m $PSEUDO_MCRAB_DIR --ratio

./plotDataMC.py --folder PUDependency -e $DSETS -s $FORMATS -m $PSEUDO_MCRAB_DIR --ratio
./plotDataMC.py --folder eSelection_Veto -e $DSETS -s $FORMATS -m $PSEUDO_MCRAB_DIR --ratio
./plotDataMC.py --folder muSelection_ -e $DSETS -s $FORMATS -m $PSEUDO_MCRAB_DIR --ratio
./plotDataMC.py --folder tauSelection_ -e $DSETS -s $FORMATS -m $PSEUDO_MCRAB_DIR --ratio
./plotDataMC.py --folder jetSelection_ -e $DSETS -s $FORMATS -m $PSEUDO_MCRAB_DIR --ratio
./plotDataMC.py --folder metSelection_ -e $DSETS -s $FORMATS -m $PSEUDO_MCRAB_DIR --ratio
./plotDataMC.py --folder bjetSelection_ -e $DSETS -s $FORMATS -m $PSEUDO_MCRAB_DIR --ratio
./plotDataMC.py --folder ForDataDrivenCtrlPlots -e $DSETS -s $FORMATS -m $PSEUDO_MCRAB_DIR --ratio
./plotDataMC.py --folder ForDataDrivenCtrlPlotsEWKFakeTaus -e $DSETS -s $FORMATS -m $PSEUDO_MCRAB_DIR --ratio
./plotDataMC.py --folder ForDataDrivenCtrlPlotsEWKGenuineTaus -e $DSETS -s $FORMATS -m $PSEUDO_MCRAB_DIR --ratio
