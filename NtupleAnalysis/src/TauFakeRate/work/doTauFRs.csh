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
set FORMATS = "png,pdf,C"
set eDSETS   = "QCD|WJets"
set iDSETS   = "SingleE|SingleMuon|DYJetsToLL_M_50_ext|ST_|TT|WJetsToLNu_HT|WWTo2L2Nu|WZ|ZZ|ttH"
set DATAERA = "Run2016"
#set ANAME   = "TauFakeRate_ee"
set ANAME   = "TauFakeRate_mm"
#set ANAME   = "Hplus2hwAnalysis_fake"
set SMODE   = "350to3000"
#set GRID    = "--gridX --gridY"
set GRID    = ""


#================================================================================================
# Data/MC plots
#================================================================================================
#./plotDataMC.py --folder "" -e ${eDSETS} -s {$FORMATS} -m ${MYDIR} --ratio --analysisName ${ANAME} --dataEra ${DATAERA} --searchMode ${SMODE}
./plotDataMC.py --folder "" -i ${iDSETS} -s {$FORMATS} -m ${MYDIR} --ratio --analysisName ${ANAME} --dataEra ${DATAERA} --searchMode ${SMODE}


#================================================================================================
# Fake Rate plots
#================================================================================================
./plotTauFakeRates.py -a ${ANAME} -m ${MYDIR} --numHisto "tauPt_num_dm0_barrel" --denHisto "tauPt_den_dm0_barrel" -s ${FORMATS} ${GRID} --yMin 0.0 --yMax 0.8 -i ${iDSETS} --analysisName ${ANAME} --dataEra ${DATAERA} --searchMode ${SMODE} #--individualMC 
./plotTauFakeRates.py -a ${ANAME} -m ${MYDIR} --numHisto "tauPt_num_dm1_barrel" --denHisto "tauPt_den_dm1_barrel" -s ${FORMATS} ${GRID} --yMin 0.0 --yMax 0.8 -i ${iDSETS} --analysisName ${ANAME} --dataEra ${DATAERA} --searchMode ${SMODE} #--individualMC
./plotTauFakeRates.py -a ${ANAME} -m ${MYDIR} --numHisto "tauPt_num_dm10_barrel" --denHisto "tauPt_den_dm10_barrel" -s ${FORMATS} ${GRID} --yMin 0.0 --yMax 0.8 -i ${iDSETS} --analysisName ${ANAME} --dataEra ${DATAERA} --searchMode ${SMODE} #--individualMC

./plotTauFakeRates.py -a ${ANAME} -m ${MYDIR} --numHisto "tauPt_num_dm0_endcap" --denHisto "tauPt_den_dm0_endcap" -s ${FORMATS} ${GRID} --yMin 0.0 --yMax 0.8 -i ${iDSETS} --analysisName ${ANAME} --dataEra ${DATAERA} --searchMode ${SMODE} #--individualMC 
./plotTauFakeRates.py -a ${ANAME} -m ${MYDIR} --numHisto "tauPt_num_dm1_endcap" --denHisto "tauPt_den_dm1_endcap" -s ${FORMATS} ${GRID} --yMin 0.0 --yMax 0.8 -i ${iDSETS} --analysisName ${ANAME} --dataEra ${DATAERA} --searchMode ${SMODE} #--individualMC
./plotTauFakeRates.py -a ${ANAME} -m ${MYDIR} --numHisto "tauPt_num_dm10_endcap" --denHisto "tauPt_den_dm10_endcap" -s ${FORMATS} ${GRID} --yMin 0.0 --yMax 0.8 -i ${iDSETS} --analysisName ${ANAME} --dataEra ${DATAERA} --searchMode ${SMODE} #--individualMC

./plotTauFakeRates.py -a ${ANAME} -m ${MYDIR} --numHisto "tauPt_num_dm0" --denHisto "tauPt_den_dm0" -s ${FORMATS} ${GRID} --yMin 0.0 --yMax 0.8 -i ${iDSETS} --analysisName ${ANAME} --dataEra ${DATAERA} --searchMode ${SMODE} #--individualMC 
./plotTauFakeRates.py -a ${ANAME} -m ${MYDIR} --numHisto "tauPt_num_dm1" --denHisto "tauPt_den_dm1" -s ${FORMATS} ${GRID} --yMin 0.0 --yMax 0.8 -i ${iDSETS} --analysisName ${ANAME} --dataEra ${DATAERA} --searchMode ${SMODE} #--individualMC
./plotTauFakeRates.py -a ${ANAME} -m ${MYDIR} --numHisto "tauPt_num_dm10" --denHisto "tauPt_den_dm10" -s ${FORMATS} ${GRID} --yMin 0.0 --yMax 0.8 -i ${iDSETS} --analysisName ${ANAME} --dataEra ${DATAERA} --searchMode ${SMODE} #--individualMC
