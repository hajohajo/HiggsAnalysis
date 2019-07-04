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
set DSETS   = "ST_"
#set FORMATS = "png"
set FORMATS = "png,pdf,C"

./plot_1d.py -n --url -e $DSETS -s $FORMATS -m ${MYDIR}
./plot_2d.py --normalizeToOne --url  --gridX --gridY --logZ -s $FORMATS -m ${MYDIR} -i "TT"
./plot_2d.py --normalizeToOne --url  --gridX --gridY --logZ -s $FORMATS -m ${MYDIR} -i "M300_mH200"
./plot_2d.py --normalizeToOne --url  --gridX --gridY --logZ -s $FORMATS -m ${MYDIR} -i "M700_mH200"
#./plot_2d.py --normalizeToOne --url  --gridX --gridY --logZ -s $FORMATS -m ${MYDIR} -i "M700_mH200"
#./plot_2d.py --normalizeToOne --url  --gridX --gridY --logZ -s $FORMATS -m ${MYDIR} -i "M300_mH200"
#./plot_2d.py --normalizeToOne --url  --gridX --gridY --logZ -s $FORMATS -m ${MYDIR} -i "M350_mH150"
#./plot_2d.py --normalizeToOne --url  --gridX --gridY --logZ -s $FORMATS -m ${MYDIR} -i "M350_mh125"
