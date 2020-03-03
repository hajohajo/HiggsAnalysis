nohup python run_muon.py /media/mlotti/Seagate\ Backup\ Plus\ Drive/data/multicrab_Hplus2hwAnalysis_v8030_20190611T1817/ &> muon_signal.out&
nohup python run_background.py /media/mlotti/Seagate\ Backup\ Plus\ Drive/data/multicrab_Hplus2hwAnalysis_v8030_20190611T1817/ &> muon_bkg.out&
nohup python run_ele.py /media/mlotti/Seagate\ Backup\ Plus\ Drive/data/multicrab_Hplus2hwAnalysis_v8030_20190611T1817/ &> ele_signal.out&
nohup python run_background_ele.py /media/mlotti/Seagate\ Backup\ Plus\ Drive/data/multicrab_Hplus2hwAnalysis_v8030_20190611T1817/ &> ele_bkg.out&
