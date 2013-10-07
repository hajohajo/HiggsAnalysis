#!/bin/sh

set -e

# Tag list modification history
# 10.11.2010/M.Kortelainen CMSSW_3_8_6
# 11.11.2010/M.Kortelainen the tag is needed to get the code to compile, moved one level up
# 12.11.2010/M.Kortelainen CMSSW_3_8_6 workaround added
# 17.11.2010/M.Kortelainen CMSSW_3_8_6 patch for OscarProducer
# 23.11.2010/M.Kortelainen CMSSW_3_8_6 updated MCEmbeddingTools tag
# 1.2.2011/M.Kortelainen CMSSW_3_9_7 Updated MCEmbeddingTools tag
# 22.3.2011/M.Kortelainen CMSSW_3_9_7 Added MuonAnalysis tag
# 14.4.2011/M.Kortelainen CMSSW_4_1_4 Updated MCEmbeddingTools tag
# 18.4.2011/M.Kortelainen CMSSW_4_1_4 Updated patMuonsWithTrigger_cff.py revision
# 11.5.2011/M.Kortelainen CMSSW_4_1_4 Added AnalysisDataFormats/EWK for edm::PtrVector<pat::Muon>
# 24.6.2011/M.Kortelainen CMSSW_4_2_4_patch1 Removed AnalysisDataFormats/EWK
# 28.6.2011/M.Kortelainen CMSSW_4_2_5 Updated MCEmbeddingTools tag
# 10.7.2011/M.Kortelainen CMSSW_4_2_5 Included AnalysisDataFormats/EWK back
# 13.8.2012/M.Kortelainen CMSSW_4_4_4 Fix compilation of SelectReplacementCandidates.cc
# 6.5.2013/M.Kortelainen CMSSW_4_4_5 Added MuScleFit correction class
# 7.10.2013/M.Kortelainen CMSSW_4_4_5 Add patch for the PAT helpers.py perf improvement

cvs co -r V00-00-13 TauAnalysis/MCEmbeddingTools
cvs up -r 1.7 TauAnalysis/MCEmbeddingTools/plugins/SelectReplacementCandidates.cc

cvs co -r V01-13-00 MuonAnalysis/MuonAssociators
cvs up -r 1.4 MuonAnalysis/MuonAssociators/python/patMuonsWithTrigger_cff.py

# We have to add protection for multiple runs in OscarProducer
addpkg SimG4Core/Application
patch -p0 < HiggsAnalysis/HeavyChHiggsToTauNu/test/tauEmbedding/OscarProducer.patch


addpkg AnalysisDataFormats/EWK
cvs up -r1.7 AnalysisDataFormats/EWK/src/classes.h
cvs up -r1.7 AnalysisDataFormats/EWK/src/classes_def.xml
cvs up -r1.3 AnalysisDataFormats/EWK/BuildFile.xml

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/MuScleFitCorrections2012
cvs co -r muscle_v4_2_0 -d MuScleFit/Calibration UserCode/scasasso/MuScleFit/Calibration 

# Improves the performance of massSearchReplaceAnyInputTag significantly
patch -p0 < HiggsAnalysis/HeavyChHiggsToTauNu/test/tauEmbedding/PAT_helpers.patch
