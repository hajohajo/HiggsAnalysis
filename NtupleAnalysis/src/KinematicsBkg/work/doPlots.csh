#!/bin/csh   

#================================================================================================
# Ensure all script arguments are passed from command line
#================================================================================================
if ($#argv != 1) then
    echo "=== You must give exactly 1 argument:"
    echo "1=PSEUDO_MCRAB_DIR"
    echo "\n=== For example:"
    echo "./doPlots.csh <pseudo-mcrab>"
    echo
    exit 1
endif

#================================================================================================
# Define variables                                                                               
#================================================================================================
set INITIAL = `echo $USER | cut -c1-1`
set MYDIR   = ${1}

./plot_1d.py -n --url -e "ST_|WW|WZ|ZZ|mH150|mH200" -s png,pdf -m ${MYDIR}
./plot_2d.py --normalizeToOne --url  --gridX --gridY --logZ -s png,pdf -m ${MYDIR} -i "TT"
./plot_2d.py --normalizeToOne --url  --gridX --gridY --logZ -s png,pdf -m ${MYDIR} -i "M1500_mH150"
./plot_2d.py --normalizeToOne --url  --gridX --gridY --logZ -s png,pdf -m ${MYDIR} -i "M1500_mh125"
./plot_2d.py --normalizeToOne --url  --gridX --gridY --logZ -s png,pdf -m ${MYDIR} -i "M300_mH200"
./plot_2d.py --normalizeToOne --url  --gridX --gridY --logZ -s png,pdf -m ${MYDIR} -i "M350_mH150"
./plot_2d.py --normalizeToOne --url  --gridX --gridY --logZ -s png,pdf -m ${MYDIR} -i "M350_mh125"
