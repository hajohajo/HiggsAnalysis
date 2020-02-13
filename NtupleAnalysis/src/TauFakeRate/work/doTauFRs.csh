#!/bin/csh   

#================================================================================================
# Ensure all script arguments are passed from command line
#================================================================================================
if ($#argv != 1) then
    echo "=== You must give exactly 1 argument:"
    echo "1=PSEUDO_MCRAB_DIR"
    echo
    echo "=== For example:"
    echo "./doTauFRs.csh TauFakeRate_Attempt4_MuonPt40_AtLeast2Jets_08Feb2020"
    exit 1
endif

#================================================================================================
# Define variables                                                                               
#================================================================================================
set INITIAL = `echo $USER | cut -c1-1`
set MYDIR   = ${1} # pseudomulticrab directory name
#set FORMATS = "png"
set FORMATS = "png,pdf,C"
set DSETS   = "QCD|WJets"
set DATAERA = "Run2016"
set ANAME   = "TauFakeRate"
set SMODE   = "80to1000"
#set ANAME   = "Hplus2hwAnalysis_fake"
#set SMODE   = "350to3000"

#================================================================================================
# Data/MC plots
#================================================================================================
./plotDataMC.py --folder "" -e ${DSETS} -s {$FORMATS} -m ${MYDIR} --ratio --analysisName ${ANAME} --dataEra ${DATAERA} --searchMode ${SMODE}


#================================================================================================
# Fake Rate plots
#================================================================================================
./plotTauFakeRates.py -m ${MYDIR} --numHisto "tauPt_num_1pr" --denHisto "tauPt_den_1pr" -s ${FORMATS} --gridX --gridY --yMin 0.0 --yMax 0.6 -e ${DSETS} --analysisName ${ANAME} --dataEra ${DATAERA} --searchMode ${SMODE} --individualMC 
./plotTauFakeRates.py -m ${MYDIR} --numHisto "tauPt_num_2pr" --denHisto "tauPt_den_2pr" -s ${FORMATS} --gridX --gridY --yMin 0.0 --yMax 0.6 -e ${DSETS} --analysisName ${ANAME} --dataEra ${DATAERA} --searchMode ${SMODE} --individualMC
./plotTauFakeRates.py -m ${MYDIR} --numHisto "tauPt_num_3pr" --denHisto "tauPt_den_3pr" -s ${FORMATS} --gridX --gridY --yMin 0.0 --yMax 0.6 -e ${DSETS} --analysisName ${ANAME} --dataEra ${DATAERA} --searchMode ${SMODE} --individualMC
