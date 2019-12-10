################################################
# NanoAOD skimming with trigger bit only
# 10.12.2019/SLehti
################################################

Copy NanoAOD_HplusTauNuAnalysisSkim.py for your analysis name
and change the triggers in the code.

Change the input file and test interacively:
python NanoAOD_HplusTauNuAnalysisSkim.py

Submit:
./multicrab.py --create -n NanoAOD_HplusTauNuAnalysisSkim.py

After DATA is 100% complete:
lumicalc.py <multicrabdir>
pileup.py <multicrabdir>
