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
set MYDIR   = ${1}
#set FORMATS = "png,pdf,C"
set FORMATS = "png"
set DSETS   = "QCD"
set SMODE   = "350to3000"
#set ANAME   = "TauFakeRate_ee"
set ANAME   = "TauFakeRate_mm"

#================================================================================================
# TH1
#================================================================================================
./plotDataMC.py --folder counters/weighted -e ${DSETS} -s ${FORMATS} -m $MYDIR --ratio -a ${ANAME} --searchMode ${SMODE}
# ./plotDataMC.py --folder "" -e ${DSETS} -s ${FORMATS} -m $MYDIR --ratio -a ${ANAME} --searchMode ${SMODE}
# #./plotDataMC.py --folder PUDependency -e ${DSETS} -s ${FORMATS} -m $MYDIR --ratio -a ${ANAME} --searchMode ${SMODE}
# ./plotDataMC.py --folder eSelection_Veto -e ${DSETS} -s ${FORMATS} -m $MYDIR --ratio -a ${ANAME} --searchMode ${SMODE}
# ./plotDataMC.py --folder muSelection_ -e ${DSETS} -s ${FORMATS} -m $MYDIR --ratio -a ${ANAME} --searchMode ${SMODE}
# ./plotDataMC.py --folder tauSelection_ -e ${DSETS} -s ${FORMATS} -m $MYDIR --ratio -a ${ANAME} --searchMode ${SMODE}
# ./plotDataMC.py --folder jetSelection_ -e ${DSETS} -s ${FORMATS} -m $MYDIR --ratio -a ${ANAME} --searchMode ${SMODE}
# ./plotDataMC.py --folder bjetSelection_ -e ${DSETS} -s ${FORMATS} -m $MYDIR --ratio -a ${ANAME} --searchMode ${SMODE}
# ./plotDataMC.py --folder metSelection_ -e ${DSETS} -s ${FORMATS} -m $MYDIR --ratio -a ${ANAME} --searchMode ${SMODE}
# ./plotDataMC.py --folder ForDataDrivenCtrlPlots -e ${DSETS} -s ${FORMATS} -m $MYDIR --ratio -a ${ANAME} --searchMode ${SMODE}
# ./plotDataMC.py --folder ForDataDrivenCtrlPlotsFakeTau -e ${DSETS} -s ${FORMATS} -m $MYDIR --ratio -a ${ANAME} --searchMode ${SMODE}
# ./plotDataMC.py --folder ForDataDrivenCtrlPlotsGenuineTau -e ${DSETS} -s ${FORMATS} -m $MYDIR --ratio -a ${ANAME} --searchMode ${SMODE}
