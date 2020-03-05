# EXAMPLE:
# source run.csh -1
set dsets_exclude = "DY|ST|WW|WZ|ZZ|WJet|ZJet|JetHT|TTZToQQ|TTTT|TTZ|M_180|M_220|M_250|M_200|M_2000|M_1000|M_800|QCD_HT|QCD_b|M_10000|M_1500|M_2500|M_3000|M_5000|M_7000|TT|M_300|M_350|M_400|M_500|M_650"
set mcrab = "/uscms_data/d3/aattikis/workspace/multicrab/multicrab_Hplus2tbAnalysis_v8030_20180508T0644/"
#set dsets_include = "|"${2}"|"
set nEvts = ${1}
set replace = "|"
set output = `echo $dsets_exclude| sed "s/|TT|/$replace/g" `
set output = `echo $output| sed "s/|M_250|/$replace/g" `
set output = `echo $output| sed "s/|M_500|/$replace/g" `
set output = `echo $output| sed "s/|M_800|/$replace/g" `
set run="./run.py -m $mcrab -e $output --noPU"
echo $run
./run.py -m $mcrab -e $output -n $nEvts --noPU
