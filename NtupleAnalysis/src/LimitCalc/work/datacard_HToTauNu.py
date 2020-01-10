'''
DESCRIPTION:
-This is a datacard template for 2016 results. 
It can be used to generate datacards for H+ -> tb analysis, 
in the fully hadronic final state. 


USAGE:
./dcardGenerator.py -x dcardHplus2tb2017Datacard_v2.py -d [directory-containing-multicrab-named-SignalAnalysis_*]


EXAMPLES:
./dcardGenerator_v2.py -x dcardDefault_h2hw_mc.py -d limits2019/ --h2tb


LAST USED:
./dcardGenerator_v3.py --analysisType HToTauNu --datacard datacard_HToTauNu.py --dir 20180620_forUnblinding_RtauMore/ --barlowBeeston -v


HN Threads:
https://hypernews.cern.ch/HyperNews/CMS/get/HIG-18-015/11.html


REFERNCE CARD:
https://gitlab.cern.ch/HPlus/HiggsAnalysis/blob/HIG-18-015/NtupleAnalysis/src/LimitCalc/work/dcardDefault2016Datacard.py

'''
#================================================================================================  
# Imports
#================================================================================================  
import HiggsAnalysis.NtupleAnalysis.tools.systematics as systematics
import HiggsAnalysis.NtupleAnalysis.tools.aux as aux
import HiggsAnalysis.NtupleAnalysis.tools.ShellStyles as ShellStyles
from HiggsAnalysis.LimitCalc.InputClasses import ObservationInput
import sys
import re

#================================================================================================
# Shell Types
#================================================================================================
ss = ShellStyles.SuccessStyle()
ns = ShellStyles.NormalStyle()
ts = ShellStyles.NoteStyle()
hs = ShellStyles.HighlightAltStyle()
es = ShellStyles.ErrorStyle()

#================================================================================================
# Function Definition
#================================================================================================
def Print(msg, printHeader=False):
    fName = __file__.split("/")[-1]
    if printHeader==True:
        print "=== ", fName
        print "\t", msg
    else:
        print "\t", msg
    return

def PrintNuisancesTable(Nuisances, DataGroups):
    align   = "{:<3} {:<30} {:<10} {:<20} {:<15} {:<10} {:<40} {:<10}"
    hLine   = "="*150
    header  = align.format("#", "ID", "Distrib.", "Function", "Value (4f)", "Scaling", "Label", "# Datasets")
    table   = []
    table.append(hLine)
    table.append(header)
    table.append(hLine)

    # For-loop: All nuisances
    for i, n in enumerate(Nuisances, 1):
        
        datasetList = []
        for j, dg in enumerate(DataGroups, 1):
            if n.getId() in dg.getNuisances():
                datasetList.append(dg.getLabel())
        if isinstance(n.getArg("value"), float):
            value = "%.4f" % n.getArg("value")
        elif n.getId() == "lumi_13TeV":
            value = n.getArg("value").getUncertaintyMax()
        else:
            value = "N/A"

        # Create the row
        #row = align.format(i, n.getId(), n.getDistribution(), n.getFunction(), value, n.getArg("scaling"), n.getLabel(), ", ".join(datasetList) )
        row = align.format(i, n.getId(), n.getDistribution(), n.getFunction(), value, n.getArg("scaling"), n.getLabel(), len(datasetList) )
        table.append(row)
    table.append(hLine)
    table.append("")
    
    # For-loop: All table rows
    for i,row in enumerate(table, 1):
        Print(row, i==1)
    return


#================================================================================================  
# Options
#================================================================================================  
OptionTest                             = True # [default: False]
OptionPaper                            = True  # [default: True]   (Changes figure style to paper style)
OptionIncludeSystematics               = True  # [default: True]   (Shape systematics; Requires pseudo-multicrab produced with doSystematics=True) 
OptionShapeSystematics                 = False # [default: True]   (Shape systematics; Requires pseudo-multicrab produced with doSystematics=True) 
OptionDoControlPlots                   = True  # [default: True]   (Produce control plots defined at end of this file)
OptionGenuineTauBackgroundSource       = "MC"  # [Options: "DataDriven", "MC"]
OptionFakeTauMeasurementSource         = "MC"  # [default: "DataDriven"] (options: "DataDriven", "MC")
OptionBr                               = 1.0   # [default: 1.0]    (The Br(t->bH+) used in figures and tables)
OptionSqrtS                            = 13    # [default: 13]     (The sqrt(s) used in figures and tables)
OptionBlindThreshold                   = None  # [default: None]   (If signal exceeds this fraction of expected events, data is blinded in a given bin)
OptionCombineSingleColumnUncertainties = False # [default: False]  (Merge nuisances with quadratic sum using the TableProducer.py Only applied to nuisances with one column)
OptionDisplayEventYieldSummary         = False # [default: False]  (Print "Event yield summary", using the TableProducer.py. A bit messy at the moment)
OptionDoWithoutSignal                  = False # [default: False]  (Do the control plots without any signal present)
OptionLimitOnSigmaBr                   = True  # [default: True]   (Set to true for heavy H+)
OptionNumberOfDecimalsInSummaries      = 1     # [defaul: 1]       (Self explanatory)
OptionConvertFromShapeToConstantList   = []    # [default: []]     (Convert these nuisances from shape to constant; Makes limits run faster & converge more easily)
OptionSeparateShapeAndNormFromSystList = []    # [default: []]     (Separate in the following shape nuisances the shape and normalization components)
OptionMassShape                        = "shapeTransverseMass"
OptionPlotNamePrefix                   = None  #"Results" #[default: None] (Prefix for plots: "<Prefix>_M<mass>_<ControlPlotInput.Title>.<ext>". Default prefix is "DataDrivenCtrlPlot")

#================================================================================================  
# Definitions
#================================================================================================  
OptionPrintNuisances                   = False            # [default: True]
MassPoints                             = [80, 90, 100, 120, 140, 150, 155, 160, 180, 200, 220, 250, 300, 400, 500, 750, 800, 1000, 1500, 2000, 2500, 3000]
PlotLegendHeader                       = "H^{#pm}#rightarrow#tau^{#pm}_{h} #nu_{#tau}"
#SignalName                             = "ChargedHiggs_HplusTB_HplusToTauNu_M_%s"
SignalName                             = "HplusTB_M%s"
DataCardName                           = "HToTaunu"       # [default: Hplus2hw_13TeV]  (Used by TableProducer.py)
BlindAnalysis                          = False            # [default: True]   (Change only if "green light" for unblinding)
MinimumStatUncertainty                 = 0.5              # [default: 0.5]    (Minimum stat. uncert. to set to bins with zero events)
UseAutomaticMinimumStatUncertainty     = True             # [default: False]  (Do NOT use the MinimumStatUncertainty; determine value from lowest non-zero rate for each dataset   )
ToleranceForLuminosityDifference       = 0.05             # [default: 0.05]   (Tolerance for throwing error on luminosity difference; "0.01" means that a 1% is required) 
ToleranceForMinimumRate                = 0.0              # [default: 0.0]    (Tolerance for almost zero rate columns with smaller rate are suppressed) 
labelPrefix                            = ""               # [default: ""]     (Prefix for the labels of datacard columns; e.g. "CMS_Hptntj_", "CMS_H2tb_")
labelPostfix                           = ""               # [default: "_GenuineTau"] (Postfix for the labels of datacard columns; e.g. "TT" --> "TT_GenuineTau")
if OptionTest:
    MassPoints = [1000]


#================================================================================================  
# File-specific settings
#================================================================================================  
 # Get the binning for the shape histogram
ShapeHistogramsDimensions = systematics.getBinningForPlot("Mt", "HToHW")

# Counter and histogram path definitions
# fixme: https://gitlab.cern.ch/HPlus/HiggsAnalysis/blob/HIG-18-015/NtupleAnalysis/src/LimitCalc/work/dcardDefault_h2tb_2016_paper.py
# fixme: https://gitlab.cern.ch/HPlus/HiggsAnalysis/blob/master/NtupleAnalysis/src/LimitCalc/work/datacard_HToHW.py
histoPathInclusive        = "ForDataDrivenCtrlPlots"
histoPathGenuineTau       = histoPathInclusive #+ "EWKGenuineTaus" QCD inverted # genuine tau requirement already as part of event selection
histoPathFakeTau          = histoPathInclusive + "EWKFakeTaus"

# EWK Datasets should only be Genuibe-tau (FakeTau = QCD inclusive + EWK GenuineTaus)
if OptionFakeTauMeasurementSource == "DataDriven":
    histoPathEWK = histoPathGenuineTau
    dsetTypeEWK  = "GenuineTaus"
else:
    histoPathEWK = histoPathInclusive
    dsetTypeEWK  = "EWKMC"


# Observation definition (how to retrieve number of observed events)
Observation = ObservationInput(datasetDefinition="Data", shapeHistoName=OptionMassShape, histoPath=histoPathInclusive)

#================================================================================================  
# Nuisance Lists (Just the strings; The objects are defined later below)
#================================================================================================ 
myLumiSystematics       = ["lumi_13TeV"]
myPileupSystematics     = ["CMS_pileup"]
myTrgEffSystematics     = ["CMS_eff_trg_MC"]
myLeptonVetoSystematics = ["CMS_eff_e_veto", "CMS_eff_m", "CMS_eff_tau"]
myJetSystematics        = ["CMS_scale_j", "CMS_res_j"]
myBtagSystematics       = ["CMS_eff_b"]

# Define systematics dictionary (easier access)
mySystematics = {}
mySystematics["MC"]     = myLumiSystematics + myPileupSystematics + myTrgEffSystematics + myLeptonVetoSystematics + myJetSystematics + myBtagSystematics
mySystematics["Signal"] = mySystematics["MC"] #+ ["CMS_Hptn_mu_RF_Hptn","CMS_Hptn_pdf_Hptn"]
mySystematics["QCD"]    = mySystematics["MC"]
mySystematics["TT"]     = mySystematics["MC"] + ["QCDscale_ttbar", "pdf_ttbar", "mass_top"] + ["CMS_Hptn_mu_RF_top","CMS_Hptn_pdf_top"]
mySystematics["ttX"]    = mySystematics["MC"] + ["QCDscale_singleTop", "pdf_singleTop", "mass_top_forSingleTop"] + ["CMS_Hptn_mu_RF_top","CMS_Hptn_pdf_top"]
mySystematics["EWK"]    = mySystematics["MC"] + ["QCDscale_ewk", "pdf_ewk"] + ["CMS_Hptn_mu_RF_ewk","CMS_Hptn_pdf_ewk"]
mySystematics["WJets"]  = mySystematics["MC"] + ["QCDscale_ewk", "pdf_ewk"]

if not OptionIncludeSystematics:
    msg = "Disabled systematics for all datasets (Stat. only datacards)"
    # For-loop: All dataset-systematics pairs
    for i,dset in enumerate(mySystematics, 1):
        mySystematics[dset] = []
    Print(ShellStyles.ErrorStyle() + msg + ShellStyles.NormalStyle(), True)

#================================================================================================  
# DataGroups (= columns in datacard) 
#================================================================================================ 
from HiggsAnalysis.LimitCalc.InputClasses import DataGroup

# Create a signal tempate
signalTemplate   = DataGroup(datasetType="Signal", histoPath=histoPathInclusive, shapeHistoName=OptionMassShape)
signalDataGroups =  []
# For-loop: All mass points
for mass in MassPoints:
    hx=signalTemplate.clone()
    hx.setLabel("Hp" + str(mass) )
    hx.setLandSProcess(1)
    hx.setValidMassPoints([mass])
    hx.setNuisances(mySystematics["Signal"])
    #hx.setDatasetDefinition("ChargedHiggs_HplusTB_HplusToHW_M%s_mH200_2ta_NLO" % (mass))
    hx.setDatasetDefinition(SignalName % (mass))
    signalDataGroups.append(hx)

# Define all the background datasets
QCDandFakeTau = DataGroup(label  = labelPrefix + "QCDandFakeTau" + labelPostfix, #
               landsProcess      = 1,
               shapeHistoName    = OptionMassShape, 
               histoPath         = histoPathInclusive,
               datasetType       = "QCD inverted", 
               datasetDefinition = "QCDMeasurementMT",
               validMassPoints   = MassPoints,
               nuisances         = mySystematics["MC"]
               )

QCDMC = DataGroup(label            = labelPrefix + "QCD" + labelPostfix, #
               landsProcess      = 1,
               shapeHistoName    = OptionMassShape, 
               histoPath         = histoPathInclusive,
               datasetType       = "QCDMC",
               datasetDefinition = "QCD", # You must use the exact name given in plots.py (_datasetMerge dictionary)
               validMassPoints   = MassPoints,
               nuisances         = mySystematics["QCD"]
               )

TT = DataGroup(label             = labelPrefix + "TT" + labelPostfix, #
               landsProcess      = 2,
               shapeHistoName    = OptionMassShape, 
               histoPath         = histoPathInclusive,
               datasetType       = "EWKMC",
               datasetDefinition = "TT", # You must use the exact name given in plots.py (_datasetMerge dictionary)
               validMassPoints   = MassPoints,
               nuisances         = mySystematics["TT"]
               )

WJets = DataGroup(label          = labelPrefix + "WJets" + labelPostfix, #
               landsProcess      = 3,
               shapeHistoName    = OptionMassShape,
               histoPath         = histoPathInclusive,
               datasetType       = dsetTypeEWK,
               datasetDefinition = "WJets", # You must use the exact name given in plots.py (_datasetMerge dictionary)
               validMassPoints   = MassPoints,
               nuisances         = mySystematics["WJets"]
               )

TTX = DataGroup(label             = labelPrefix + "SingleTop" + labelPostfix, # labelPrefix + "ttX" + labelPostfix,
                landsProcess      = 4,
                shapeHistoName    = OptionMassShape,
                histoPath         = histoPathInclusive,
                datasetType       = dsetTypeEWK,
                datasetDefinition = "SingleTop", #"ttX", # You must use the exact name given in plots.py (_datasetMerge dictionary)
                validMassPoints   = MassPoints,
                nuisances         = mySystematics["ttX"]
                )

DYJets = DataGroup(label             = labelPrefix + "DYJets" + labelPostfix,
                   landsProcess      = 5,
                   shapeHistoName    = OptionMassShape,
                   histoPath         = histoPathInclusive,
                   datasetType       = dsetTypeEWK,
                   datasetDefinition = "DYJetsToLLHT", # You must use the exact name given in plots.py
                   validMassPoints   = MassPoints,
                   nuisances         = mySystematics["MC"]
                   )

Diboson = DataGroup(label             = labelPrefix + "Diboson" + labelPostfix,
                    landsProcess      = 6,
                    shapeHistoName    = OptionMassShape,
                    histoPath         = histoPathInclusive,
                    datasetType       = dsetTypeEWK,
                    datasetDefinition = "Diboson", # You must use the exact name given in plots.py
                    validMassPoints   = MassPoints,
                    nuisances         = mySystematics["MC"]
                   )

# Append datasets in order you want them to appear in the data-driven control plot stack
# The list \"DataGroups\" is required by the DataCardGenerator.py module
DataGroups = []
DataGroups.extend(signalDataGroups)
if 0:
    DataGroups.append(QCDandFakeTau)
else:
    DataGroups.append(QCDMC)
DataGroups.append(TT)
DataGroups.append(TTX)
DataGroups.append(WJets)
DataGroups.append(DYJets)
DataGroups.append(Diboson)

#================================================================================================  
# Shape Nuisance Parameters (aka Systematics)  (= rows in datacard) 
#================================================================================================ 
from HiggsAnalysis.LimitCalc.InputClasses import Nuisance

# Define all individual nuisances that can be potentially used (ShapeVariations require running with systematics flag! Defined in AnalysisBuilder.py)
JES_Shape        = Nuisance(id="CMS_scale_j"    , label="Jet Energy Scale (JES)", distr="shapeQ", function="ShapeVariation", systVariation="JES")
JER_Shape        = Nuisance(id="CMS_res_j"      , label="Jet Energy Resolution (JER)", distr="shapeQ", function="ShapeVariation", systVariation="JER")
bTagSF_Shape     = Nuisance(id="CMS_eff_b"      , label="b tagging", distr="shapeQ", function="ShapeVariation", systVariation="BTagSF")
topPt_Shape      = Nuisance(id="CMS_topreweight", label="Top pT reweighting", distr="shapeQ", function="ShapeVariation", systVariation="TopPt")
PU_Shape         = Nuisance(id="CMS_pileup"     , label="Pileup", distr="shapeQ", function="ShapeVariation", systVariation="PUWeight")
# NOTE: systVariation key is first declared in HiggsAnalysis/NtupleAnalysis/python/AnalysisBuilder.py

#================================================================================================  
# Constant Nuisance Parameters (aka Systematics)  (= rows in datacard) 
#================================================================================================ 
tt_scale_down        = systematics.getCrossSectionUncertainty("TTJets_scale").getUncertaintyDown()
tt_scale_up          = systematics.getCrossSectionUncertainty("TTJets_scale").getUncertaintyUp()
tt_pdf_down          = systematics.getCrossSectionUncertainty("TTJets_pdf").getUncertaintyDown()
tt_pdf_up            = systematics.getCrossSectionUncertainty("TTJets_pdf").getUncertaintyUp()
tt_mass_down         = systematics.getCrossSectionUncertainty("TTJets_mass").getUncertaintyDown()
tt_mass_up           = systematics.getCrossSectionUncertainty("TTJets_mass").getUncertaintyUp()
ttX_scale_down       = systematics.getCrossSectionUncertainty("SingleTop_scale").getUncertaintyDown()
ttX_pdf_down         = systematics.getCrossSectionUncertainty("SingleTop_pdf").getUncertaintyDown()
ewk_scale_down       = systematics.getCrossSectionUncertainty("Diboson_scale").getUncertaintyDown() #("EWK_scale").getUncertaintyDown()
ewk_scale_up         = systematics.getCrossSectionUncertainty("Diboson_scale").getUncertaintyUp()   #("EWK_scale").getUncertaintyUp()
ewk_pdf_down         = systematics.getCrossSectionUncertainty("Diboson_pdf").getUncertaintyDown()   #("EWK_pdf").getUncertaintyDown()
ewk_pdf_up           = systematics.getCrossSectionUncertainty("Diboson_pdf").getUncertaintyUp()     #("EWK_pdf").getUncertaintyUp()
lumi_2016            = systematics.getLuminosityUncertainty("2016")

# Default nuisances
lumi13TeV_Const = Nuisance(id="lumi_13TeV"         , label="Luminosity 13 TeV uncertainty", distr="lnN", function="Constant", value=lumi_2016)
trgMC_Const     = Nuisance(id="CMS_eff_trg_MC"     , label="Trigger MC efficiency (Approx.)", distr="lnN", function="Constant", value=0.05)
PU_Const        = Nuisance(id="CMS_pileup"         , label="Pileup (Approx.)", distr="lnN", function="Constant", value=0.05)
eVeto_Const     = Nuisance(id="CMS_eff_e_veto"     , label="e veto", distr="lnN", function="Ratio", numerator="passed e selection (Veto)", denominator="Primary vertex selection", scaling=0.02)
muVeto_Const    = Nuisance(id="CMS_eff_m"          , label="mu ID (Approx.)" , distr="lnN", function="Constant", value=0.10)
tauVeto_Const   = Nuisance(id="CMS_eff_tau"        , label="tau ID (Approx.)", distr="lnN", function="Constant", value=0.10)
bTagSF_Const    = Nuisance(id="CMS_eff_b"          , label="b tagging (Approx.)", distr="lnN", function="Constant", value=0.05)
JES_Const       = Nuisance(id="CMS_scale_j"        , label="Jet Energy Scale (JES) (Approx.)"     , distr="lnN", function="Constant", value=0.03)
JER_Const       = Nuisance(id="CMS_res_j"          , label="Jet Energy Resolution (JER) (Approx.)", distr="lnN", function="Constant", value=0.04)
topPt_Const     = Nuisance(id="CMS_topreweight"    , label="Top pT reweighting (Approx.)", distr="lnN", function="Constant", value=0.25)

# Cross section uncertainties
ttbar_scale_Const    = Nuisance(id="QCDscale_ttbar"       , label="QCD XSection uncertainties", distr="lnN", function="Constant", value=tt_scale_down, upperValue=tt_scale_up)
ttbar_pdf_Const      = Nuisance(id="pdf_ttbar"            , label="TTbar XSection pdf uncertainty", distr="lnN", function="Constant", value=tt_pdf_down, upperValue=tt_pdf_up)
ttbar_mass_Const     = Nuisance(id="mass_top"             , label="TTbar XSection top mass uncertainty", distr="lnN", function="Constant", value=tt_mass_down, upperValue=tt_mass_up) 
ttX_mass_Const       = Nuisance(id="mass_top_forSingleTop", label="Single top mass uncertainty", distr="lnN", function="Constant", value=0.022)
ttX_scale_Const      = Nuisance(id="QCDscale_singleTop"   , label="QCD XSection uncertainties", distr="lnN", function="Constant", value=ttX_scale_down)
ttX_pdf_Const        = Nuisance(id="pdf_singleTop"        , label="Single top XSection pdf ucnertainty", distr="lnN", function="Constant", value=ttX_pdf_down)
ewk_scale_Const      = Nuisance(id="QCDscale_ewk"         , label="EWK XSection uncertainties", distr="lnN", function="Constant", value=ewk_scale_down)
ewk_pdf_Const        = Nuisance(id="pdf_ewk"              , label="EWK XSection pdf uncertainty", distr="lnN", function="Constant", value=ewk_pdf_down)

#==== Acceptance uncertainties (QCDscale)
RF_QCDscale_top_const  = Nuisance(id="CMS_Hptn_mu_RF_top" , label="Scale acceptance uncertainty for top"   , distr="lnN", function="Constant",value=0.02)
RF_QCDscale_ewk_const  = Nuisance(id="CMS_Hptn_mu_RF_ewk" , label="Scale acceptance uncertainty for EWK"   , distr="lnN", function="Constant",value=0.05)
RF_QCDscale_Hptn_const = Nuisance(id="CMS_Hptn_mu_RF_Hptn", label="Scale acceptance uncertainty for signal", distr="lnN", function="Constant",value=0.048)
#RF_QCDscale_Hptn_const = Nuisance(id="CMS_Hptn_mu_RF_Hptn_heavy", label="QCDscale acceptance uncertainty for signal", distr="lnN", function="Constant",value=0.012)

#==== Acceptance uncertainties  (PDF)
RF_pdf_top_const  = Nuisance(id="CMS_Hptn_pdf_top", label="PDF acceptance uncertainty for top", distr="lnN", function="Constant",value=0.02,upperValue=0.0027)
RF_pdf_ewk_const  = Nuisance(id="CMS_Hptn_pdf_ewk", label="PDF acceptance uncertainty for EWK", distr="lnN", function="Constant",value=0.033,upperValue=0.046)
RF_pdf_Hptn_const = Nuisance(id="CMS_Hptn_pdf_Hptn", label="PDF acceptance uncertainty for signal", distr="lnN", function="Constant",value=0.004,upperValue=0.017)

# Fake-b nuisances
tf_FakeTau_Const          = Nuisance(id="CMS_Hptn_fake_t_transferfactor", label="Transfer Factor uncertainty", distr="lnN", function="Constant", value=0.10)
lumi13TeV_FakeTau_Const   = Nuisance(id="lumi_13TeV_forFakeTau"      , label="Luminosity 13 TeV uncertainty", distr="lnN", function="ConstantForFakeTau", value=lumi_2016)
trgMC_FakeTau_Const       = Nuisance(id="CMS_eff_trg_MC_forFakeTau"  , label="Trigger MC efficiency (Approx.)", distr="lnN", function="ConstantForFakeTau", value=0.05)
PU_FakeTau_Const          = Nuisance(id="CMS_pileup_forFakeTau"      , label="Pileup (Approx.)", distr="lnN", function="ConstantForFakeTau", value=0.05)
eVeto_FakeTau_Const       = Nuisance(id="CMS_eff_e_veto_forFakeTau"  , label="e veto", distr="lnN", function="Ratio", numerator="passed e selection (Veto)", denominator="passed PV", scaling=0.02) 
muVeto_FakeTau_Const      = Nuisance(id="CMS_eff_m_veto_forFakeTau"  , label="mu veto", distr="lnN", function="Ratio", numerator="passed mu selection (Veto)", denominator="passed e selection (Veto)", scaling=0.01)
tauVeto_FakeTau_Const     = Nuisance(id="CMS_eff_tau_veto_forFakeTau", label="tau veto", distr="lnN", function="Ratio", numerator="Passed tau selection (Veto)", denominator="passed mu selection (Veto)", scaling=0.01)
JES_FakeTau_Const         = Nuisance(id="CMS_scale_j_forFakeTau"     , label="Jet Energy Scale (JES) (Approx.)"     , distr="lnN", function="ConstantForFakeTau", value=0.03)
JER_FakeTau_Const         = Nuisance(id="CMS_res_j_forFakeTau"       , label="Jet Energy Resolution (JER) (Approx.)", distr="lnN", function="ConstantForFakeTau", value=0.04)
bTagSF_FakeTau_Const      = Nuisance(id="CMS_eff_b_forFakeTau"       , label="b tagging (Approx.)", distr="lnN", function="ConstantForFakeTau", value=0.05)
topPt_FakeTau_Const       = Nuisance(id="CMS_topreweight_forFakeTau" , label="Top pT reweighting (Approx.)", distr="lnN", function="ConstantForFakeTau", value=0.25)
ttbar_scale_FakeTau_Const = Nuisance(id="QCDscale_ttbar_forFakeTau"  , label="QCD XSection uncertainties", distr="lnN", function="ConstantForFakeTau", value=tt_scale_down, upperValue=tt_scale_up)
ttbar_pdf_FakeTau_Const   = Nuisance(id="pdf_ttbar_forFakeTau"       , label="TTbar XSection pdf uncertainty", distr="lnN", function="ConstantForFakeTau", value=tt_pdf_down, upperValue=tt_pdf_up)
ttbar_mass_FakeTau_Const  = Nuisance(id="mass_top_forFakeTau"        , label="TTbar XSection top mass uncertainty", distr="lnN", function="ConstantForFakeTau", value=tt_mass_down, upperValue=tt_mass_up) 
RF_QCDscale_FakeTau_const = Nuisance(id="CMS_Hptn_mu_RF_top_forFakeTau", label="Scale acceptance uncertainty for FakeTau" , distr="lnN", function="ConstantForFakeTau",value=0.02)
RF_pdf_FakeTau_const      = Nuisance(id="CMS_Hptn_pdf_top_forFakeTau"  , label="PDF acceptance uncertainty for FakeTau"      , distr="lnN", function="ConstantForFakeTau",value=0.02, upperValue=0.0027)


#================================================================================================ 
# Nuisance List (If a given nuisance "name" is used in any of the DataGroups it must be appended)
#================================================================================================ 
ReservedNuisances = []
Nuisances = []
Nuisances.append(lumi13TeV_Const)
if OptionShapeSystematics:
    Nuisances.append(PU_Shape)
    Nuisances.append(JES_Shape)
    Nuisances.append(JER_Shape)
    Nuisances.append(bTagSF_Shape) 
    Nuisances.append(tf_FakeTau_Shape)
else:
    Nuisances.append(PU_Const)
    Nuisances.append(PU_FakeTau_Const)
    Nuisances.append(JES_Const)
    Nuisances.append(JES_FakeTau_Const)
    Nuisances.append(JER_FakeTau_Const)
    Nuisances.append(JER_Const)
    Nuisances.append(bTagSF_Const) 
    Nuisances.append(bTagSF_FakeTau_Const)
    Nuisances.append(tf_FakeTau_Const)

# Common in Shapes/Constants
Nuisances.append(lumi13TeV_FakeTau_Const)
Nuisances.append(trgMC_Const)
Nuisances.append(trgMC_FakeTau_Const)

# Approximation 2: lepton and tau-veto neglected (negligible contribution)
Nuisances.append(eVeto_Const)
Nuisances.append(muVeto_Const)
Nuisances.append(tauVeto_Const)
# Nuisances.append(eVeto_FakeTau_Const)
# Nuisances.append(muVeto_FakeTau_Const)
# Nuisances.append(tauVeto_FakeTau_Const)

# Approximation 1: only ttbar xsect uncertainty applied to FakeTau, as ttbar dominates the EWK GenuineTau (but uncertainty is scaled according to 1-purity)
Nuisances.append(ttbar_scale_Const) 
Nuisances.append(ttbar_scale_FakeTau_Const) 
Nuisances.append(ttbar_pdf_Const)
Nuisances.append(ttbar_pdf_FakeTau_Const)
Nuisances.append(ttbar_mass_Const)
Nuisances.append(ttbar_mass_FakeTau_Const)
Nuisances.append(ttX_scale_Const)
Nuisances.append(ttX_pdf_Const)
Nuisances.append(ttX_mass_Const)
Nuisances.append(ewk_scale_Const)
Nuisances.append(ewk_pdf_Const)

# RF/QCD Scale
Nuisances.append(RF_QCDscale_top_const)
Nuisances.append(RF_QCDscale_ewk_const)
Nuisances.append(RF_QCDscale_Hptn_const)
#Nuisances.append(RF_QCDscale_Hptn_const)
Nuisances.append(RF_QCDscale_FakeTau_const)
Nuisances.append(RF_pdf_top_const)
Nuisances.append(RF_pdf_ewk_const)
Nuisances.append(RF_pdf_Hptn_const)
Nuisances.append(RF_pdf_FakeTau_const)

if OptionPrintNuisances:
    PrintNuisancesTable(Nuisances, DataGroups)

#================================================================================================ 
# Merge nuisances to same row (first item specifies the name for the row)
# This is for correlated uncertainties. It forces 2 nuisances to be on SAME datacard row
# For example, ttbar xs scale and singleTop pdf should be varied togethed (up or down) but always
# in phase.
# WARNING: This mostly (or solely?) applies for constants. not shape systematics!
# NOTE: When merging nuisances the resultant "merged" name will be the first element of the merging list
# MergeNuisances.append(["QCDscale_ttbar", "QCDscale_singleTop"]) #resultant name will be "QCDscale_ttbar"
#================================================================================================ 
MergeNuisances=[]

# Correlate ttbar and single top cross-section uncertainties
MergeNuisances.append(["QCDscale_ttbar", "QCDscale_ttbar_forFakeTau"])
MergeNuisances.append(["pdf_ttbar"  , "pdf_ttbar_forFakeTau"])

# Correlate FakeTau and GenuineTau uncerainties
MergeNuisances.append(["lumi_13TeV"       , "lumi_13TeV_forFakeTau"])
MergeNuisances.append(["CMS_eff_trg_MC"   , "CMS_eff_trg_MC_forFakeTau"])
MergeNuisances.append(["mass_top", "mass_top_forFakeTau", "mass_top_forSingleTop"])
MergeNuisances.append(["CMS_Hptn_mu_RF_top", "CMS_Hptn_mu_RF_top_forFakeTau"])
MergeNuisances.append(["CMS_Hptn_pdf_top", "CMS_Hptn_pdf_top_forFakeTau"])

if not OptionShapeSystematics:
    MergeNuisances.append(["CMS_pileup"   , "CMS_pileup_forFakeTau"])
    MergeNuisances.append(["CMS_eff_b"    , "CMS_eff_b_forFakeTau"])
    MergeNuisances.append(["CMS_scale_j"  , "CMS_scale_j_forFakeTau"])
    MergeNuisances.append(["CMS_res_j"    , "CMS_res_j_forFakeTau"])
    
#================================================================================================ 
# Convert shape systematics to constants if asked
#================================================================================================ 
from HiggsAnalysis.LimitCalc.InputClasses import convertFromSystVariationToConstant
nSysTotal     = len(Nuisances)
nSysToConvert = len(OptionConvertFromShapeToConstantList)
if nSysToConvert > 0:
    Print("Converting %s/%s shape systematics to constants if asked." % (nSysToConvert, nSysTotal), True)
    convertFromSystVariationToConstant(Nuisances, OptionConvertFromShapeToConstantList)

# Separate the shape nuisances and the shape and normalization components if asked
from HiggsAnalysis.LimitCalc.InputClasses import separateShapeAndNormalizationFromSystVariation
nSysShapeComponents = len(OptionSeparateShapeAndNormFromSystList)
if (nSysShapeComponents>0):
    Print("Separating %s/%s shape and normalization components" % (nSysShapeComponents, nSysTotal), True)
    separateShapeAndNormalizationFromSystVariation(Nuisances, OptionSeparateShapeAndNormFromSystList)

#================================================================================================ 
# Control plots
#================================================================================================ 
# https://twiki.cern.ch/twiki/bin/view/CMS/Internal/PubGuidelines#Figures_and_tables
uPt     = "GeV/c"
uMass   = "GeV/c^{2}"
uEnergy = "GeV"
if OptionPaper:
    uPt   = "GeV"
    uMass = "GeV"

from HiggsAnalysis.LimitCalc.InputClasses import ControlPlotInput
ControlPlots= []

hVertices_Std = ControlPlotInput(
    title     = "NVertices_AfterStandardSelections",
    histoName = "NVertices_AfterStandardSelections",
    details   = { "xlabel"             : "vertex multiplicity",
                  "ylabel"             : "< Events >",
                  "divideByBinWidth"   : True,
                  "unit"               : "",
                  "log"                : True,
                  "legendPosition"     : "NE",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": 0.5e-2, "ymaxfactor": 10, "xmax": 80.0} 
                  },
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot    
    )

hTauPt_Std = ControlPlotInput(
    title     = "SelectedTau_pT_AfterStandardSelections",
    histoName = "SelectedTau_pT_AfterStandardSelections",
    details   = { "xlabel": "Selected #tau ^{}p_{T}",
                  "ylabel": "Events/^{}#Deltap_{T}",
                  "divideByBinWidth": True,
                  "unit": "GeV/c",
                  "log": True,
                  "legendPosition": "NE",
                  "opts": {"ymin": 0.0009, "ymaxfactor": 25, "xmax": 500} 
                  },
    )

hTauLdgTrkPt_Std = ControlPlotInput(
    title     = "SelectedTau_ldgTrkPt_AfterStandardSelections",
    histoName = "SelectedTau_ldgTrkPt_AfterStandardSelections",
    details   = { "xlabel": "#tau leading track ^{}p_{T}",
                  "ylabel": "Events/^{}#Deltap_{T}",
                  "divideByBinWidth": True,
                  "unit": "GeV/c",
                  "log": True,
                  "ratioLegendPosition": "right",
                  "opts": {"ymin": 0.0009, "ymaxfactor": 10, "xmax": 500} 
                  },
    )

hTauEta_Std = ControlPlotInput(
    title      = "SelectedTau_eta_AfterStandardSelections",
    histoName  = "SelectedTau_eta_AfterStandardSelections",
    details    = { "xlabel": "Selected #tau #eta",
                   "ylabel": "Events",
                   "divideByBinWidth": False,
                   "unit": "",
                   "log": True,
                   "legendPosition": "SW",
                   "opts": {"ymin": 0.0009} 
                   },
    )

hTauPhi_Std = ControlPlotInput(
    title     = "SelectedTau_phi_AfterStandardSelections",
    histoName = "SelectedTau_phi_AfterStandardSelections",
    details   = { "xlabel": "Selected #tau #phi",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "{}^{o}",
                  "log": True,
                  "legendPosition": "SW",
                  "opts": {"ymin": 0.0009} 
                  },
    )

hTauRtauFullRange_Std = ControlPlotInput(
    title     = "SelectedTau_Rtau_FullRange_AfterStandardSelections",
    histoName = "SelectedTau_Rtau_AfterStandardSelections",
    details   = { "xlabel": "Selected #tau ^{}R_{#tau}",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "",
                  "log": True,
                  "legendPosition": "SE",
                  "ratioLegendPosition": "right",
                  "opts": {"ymin": 0.0009,  "xmax": 0.75} 
                  },
    )

hTauRtau_Std = ControlPlotInput(
    title     = "SelectedTau_Rtau_AfterStandardSelections",
    histoName = "SelectedTau_Rtau_AfterStandardSelections",
    details   = { "xlabel": "Selected #tau ^{}R_{#tau}",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "",
                  "log": True,
                  "legendPosition": "SE",
                  "ratioLegendPosition": "right",
                  "opts": {"ymin": 0.0009,  "xmin": 0.75} 
                  },
    )

hTauDM_Std = ControlPlotInput(
    title     = "SelectedTau_DecayMode_AfterStandardSelections",
    histoName = "SelectedTau_DecayMode_AfterStandardSelections",
    details   = { "xlabel": "Selected #tau Decay mode",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "",
                  "log": True,
                  "ratioLegendPosition": "right",
                  "opts": {"ymin": 0.9} 
                  },
    )

hTauNprongs_Std = ControlPlotInput(
    title     = "SelectedTau_Nprongs_AfterStandardSelections",
    histoName = "SelectedTau_Nprongs_AfterStandardSelections",
    details   = { "xlabel": "Selected #tau N_{prongs}",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "",
                  "log": True,
                  "ratioLegendPosition": "right",
                  "opts": {"ymin": 0.9} 
                  },
    )

hTauSource_Std = ControlPlotInput(
    title     = "SelectedTau_source_AfterStandardSelections",
    histoName = "SelectedTau_source_AfterStandardSelections",
    details   = { "xlabel": "",
                  "ylabel": "Events",
                  "xlabelsize": 10,
                  "divideByBinWidth": False,
                  "init": "",
                  "log": True,
                  "ratioLegendPosition": "right",
                  "opts": {"ymin": 0.9}
                  },
    )

hNJets_Std = ControlPlotInput(
    title     = "Njets_AfterStandardSelections",
    histoName = "Njets_AfterStandardSelections",
    details   = { "xlabel": "Number of selected jets",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "",
                  "log": True,
                  "opts": {"ymin": 0.9, "xmin": 3.0, "xmax": 14.0}
                  },
    flowPlotCaption  = "^{}#tau_{h}+#geq3j", # Leave blank if you don't want to include the item to the selection flow plot
    )

hJetPt_Std = ControlPlotInput(
    title     = "JetPt_AfterStandardSelections",
    histoName = "JetPt_AfterStandardSelections",
    details   = { "xlabel": "jet ^{}p_{T}",
                  "ylabel": "Events/^{}#Deltap_{T}",
                  "divideByBinWidth": True,
                  "unit": "GeV/c",
                  "log": True,
                  "opts": {"ymin": 0.0009, "ymaxfactor": 10, "xmax": 500} 
                  },
    )

hJetEta_Std = ControlPlotInput(
    title     = "JetEta_AfterStandardSelections",
    histoName = "JetEta_AfterStandardSelections",
    details   = { "xlabel": "jet #eta",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "",
                  "log": True,
                  "legendPosition": "SW",
                  "opts": {"ymin": 0.09} 
                  },
)

hCollAngularCutsMin = ControlPlotInput(
    title      = "CollinearAngularCutsMinimum",
    histoName  = "CollinearAngularCutsMinimum",
    details    = { "xlabel": "R_{coll}^{min}",
                   "ylabel": "Events",
                   "divideByBinWidth": False,
                   "unit": "{}^{o}",
                   "log": True,
                   "legendPosition": "SW",
                   "opts": {"ymin": 0.09} 
                   },
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    )

hNBjets = ControlPlotInput(
    title     = "BJetSelection",
    histoName = "NBjets",
    details   = { "xlabel": "Number of selected b jets",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "",
                  "log": True,
                  "opts": {"ymin": 0.09} 
                  },
    flowPlotCaption  = "#geq1 b tag", # Leave blank if you don't want to include the item to the selection flow plot
    )

hBJetPt = ControlPlotInput(
    title     = "BJetPt",
    histoName = "BJetPt",
    details   = { "xlabel": "b jet ^{}p_{T}",
                  "ylabel": "Events/^{}#Deltap_{T}",
                  "divideByBinWidth": True,
                  "unit": "GeV/c",
                  "log": True,
                  "opts": {"ymin": 0.0009, "ymaxfactor": 10, "xmax": 500} },
    )

hBJetEta = ControlPlotInput(
    title      = "BJetEta",
    histoName  = "BJetEta",
    details    = { "xlabel": "b jet #eta",
                   "ylabel": "Events",
                   "divideByBinWidth": False,
                   "unit": "",
                   "log": True,
                   "legendPosition": "SW",
                   "opts": {"ymin": 0.09} 
                   },
    )

hBtatDisc  = ControlPlotInput(
    title     = "BtagDiscriminator",
    histoName = "BtagDiscriminator",
    details   = { "xlabel": "b tag discriminator",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "",
                  "log": True,
                  "legendPosition": "SE",
                  "opts": {"ymin": 0.9} 
                  },
    )

hMET = ControlPlotInput(
    title     = "MET",
    histoName = "MET",
    details   = { "xlabel": "E_{T}^{miss}",
                  "ylabel": "Events/^{}#DeltaE_{T}^{miss}",
                  "divideByBinWidth": True,
                  "unit": "GeV",
                  "log": True,
                  "opts": {"ymin": 0.00009, "ymaxfactor": 10, "xmax": 400} 
                  },
    flowPlotCaption  = "^{}E_{T}^{miss}", # Leave blank if you don't want to include the item to the selection flow plot
    )

hMETPhi = ControlPlotInput(
    title     = "METPhi",
    histoName = "METPhi",
    details   = { "xlabel": "E_{T}^{miss} #phi",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "{}^{o}",
                  "log": True,
                  "legendPosition": "SW",
                  "opts": {"ymin": 0.09} 
                  }
    )

hBacktoBackAngularCutsMin = ControlPlotInput(
    title     = "BackToBackAngularCutsMinimum",
    histoName = "BackToBackAngularCutsMinimum",
    details   = { "xlabel": "^{}R_{bb}^{min}",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "{}^{o}",
                  "log": True,
                  "legendPosition": "SE",
                  "opts": {"ymin": 0.09} 
                  },
    flowPlotCaption  = "^{}R_{bb}^{min}", # Leave blank if you don't want to include the item to the selection flow plot
    )

hVertices_All = ControlPlotInput(
    title     = "NVertices_AfterAllSelections",
    histoName = "NVertices_AfterAllSelections",
    details   = { "xlabel"             : "vertex multiplicity",
                  "ylabel"             : "< Events >",
                  "divideByBinWidth"   : True,
                  "unit"               : "",
                  "log"                : True,
                  "legendPosition"     : "NE",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": 0.5e-2, "ymaxfactor": 10, "xmax": 80.0} 
                  },
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot    
    )

hTauPt_All = ControlPlotInput(
    title     = "SelectedTau_pT_AfterAllSelections",
    histoName = "SelectedTau_pT_AfterAllSelections",
    details   = { "xlabel": "Selected #tau ^{}p_{T}",
                  "ylabel": "Events/^{}#Deltap_{T}",
                  "divideByBinWidth": True,
                  "unit": "GeV/c",
                  "log": True,
                  "legendPosition": "NE",
                  "opts": {"ymin": 0.0009, "ymaxfactor": 25, "xmax": 500} 
                  },
    )

hTauLdgTrkPt_All = ControlPlotInput(
    title     = "SelectedTau_ldgTrkPt_AfterAllSelections",
    histoName = "SelectedTau_ldgTrkPt_AfterAllSelections",
    details   = { "xlabel": "#tau leading track ^{}p_{T}",
                  "ylabel": "Events/^{}#Deltap_{T}",
                  "divideByBinWidth": True,
                  "unit": "GeV/c",
                  "log": True,
                  "ratioLegendPosition": "right",
                  "opts": {"ymin": 0.0009, "ymaxfactor": 10, "xmax": 500} 
                  },
    )

hTauEta_All = ControlPlotInput(
    title      = "SelectedTau_eta_AfterAllSelections",
    histoName  = "SelectedTau_eta_AfterAllSelections",
    details    = { "xlabel": "Selected #tau #eta",
                   "ylabel": "Events",
                   "divideByBinWidth": False,
                   "unit": "",
                   "log": True,
                   "legendPosition": "SW",
                   "opts": {"ymin": 0.0009} 
                   },
    )

hTauPhi_All = ControlPlotInput(
    title     = "SelectedTau_phi_AfterAllSelections",
    histoName = "SelectedTau_phi_AfterAllSelections",
    details   = { "xlabel": "Selected #tau #phi",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "{}^{o}",
                  "log": True,
                  "legendPosition": "SW",
                  "opts": {"ymin": 0.0009} 
                  },
    )

hTauRtauFullRange_All = ControlPlotInput(
    title     = "SelectedTau_Rtau_FullRange_AfterAllSelections",
    histoName = "SelectedTau_Rtau_AfterAllSelections",
    details   = { "xlabel": "Selected #tau ^{}R_{#tau}",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "",
                  "log": True,
                  "legendPosition": "SE",
                  "ratioLegendPosition": "right",
                  "opts": {"ymin": 0.0009,  "xmax": 0.75} 
                  },
    )

hTauRtau_All = ControlPlotInput(
    title     = "SelectedTau_Rtau_AfterAllSelections",
    histoName = "SelectedTau_Rtau_AfterAllSelections",
    details   = { "xlabel": "Selected #tau ^{}R_{#tau}",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "",
                  "log": True,
                  "legendPosition": "SE",
                  "ratioLegendPosition": "right",
                  "opts": {"ymin": 0.0009,  "xmin": 0.75} 
                  },
    )

hTauDM_All = ControlPlotInput(
    title     = "SelectedTau_DecayMode_AfterAllSelections",
    histoName = "SelectedTau_DecayMode_AfterAllSelections",
    details   = { "xlabel": "Selected #tau Decay mode",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "",
                  "log": True,
                  "ratioLegendPosition": "right",
                  "opts": {"ymin": 0.9} 
                  },
    )

hTauNprongs_All = ControlPlotInput(
    title     = "SelectedTau_Nprongs_AfterAllSelections",
    histoName = "SelectedTau_Nprongs_AfterAllSelections",
    details   = { "xlabel": "Selected #tau N_{prongs}",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "",
                  "log": True,
                  "ratioLegendPosition": "right",
                  "opts": {"ymin": 0.9} 
                  },
    )

hTauSource_All = ControlPlotInput(
    title     = "SelectedTau_source_AfterAllSelections",
    histoName = "SelectedTau_source_AfterAllSelections",
    details   = { "xlabel": "",
                  "ylabel": "Events",
                  "xlabelsize": 10,
                  "divideByBinWidth": False,
                  "init": "",
                  "log": True,
                  "ratioLegendPosition": "right",
                  "opts": {"ymin": 0.9}
                  },
    )

hNJets_All = ControlPlotInput(
    title     = "Njets_AfterAllSelections",
    histoName = "Njets_AfterAllSelections",
    details   = { "xlabel": "Number of selected jets",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "",
                  "log": True,
                  "opts": {"ymin": 0.9} },
    flowPlotCaption  = "^{}#tau_{h}+#geq3j", # Leave blank if you don't want to include the item to the selection flow plot
    )

hJetPt_All = ControlPlotInput(
    title     = "JetPt_AfterAllSelections",
    histoName = "JetPt_AfterAllSelections",
    details   = { "xlabel": "jet ^{}p_{T}",
                  "ylabel": "Events/^{}#Deltap_{T}",
                  "divideByBinWidth": True,
                  "unit": "GeV/c",
                  "log": True,
                  "opts": {"ymin": 0.0009, "ymaxfactor": 10, "xmax": 500} 
                  },
    )

hJetEta_All = ControlPlotInput(
    title     = "JetEta_AfterAllSelections",
    histoName = "JetEta_AfterAllSelections",
    details   = { "xlabel": "jet #eta",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "",
                  "log": True,
                  "legendPosition": "SW",
                  "opts": {"ymin": 0.09} 
                  },
)

hCollAngularCutsMin_All = ControlPlotInput(
    title      = "CollinearAngularCutsMinimum_AfterAllSelections",
    histoName  = "CollinearAngularCutsMinimum_AfterAllSelections",
    details    = { "xlabel": "R_{coll}^{min}",
                   "ylabel": "Events",
                   "divideByBinWidth": False,
                   "unit": "{}^{o}",
                   "log": True,
                   "legendPosition": "SW",
                   "opts": {"ymin": 0.09} 
                   },
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot
    )

hNBjets_All = ControlPlotInput(
    title     = "BJetSelection_AfterAllSelections",
    histoName = "NBjets",
    details   = { "xlabel": "Number of selected b jets",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "",
                  "log": True,
                  "opts": {"ymin": 0.09} 
                  },
    flowPlotCaption  = "#geq1 b tag", # Leave blank if you don't want to include the item to the selection flow plot
    )

hBJetPt_All = ControlPlotInput(
    title     = "BJetPt_AfterAllSelections",
    histoName = "BJetPt_AfterAllSelections",
    details   = { "xlabel": "b jet ^{}p_{T}",
                  "ylabel": "Events/^{}#Deltap_{T}",
                  "divideByBinWidth": True,
                  "unit": "GeV/c",
                  "log": True,
                  "opts": {"ymin": 0.0009, "ymaxfactor": 10, "xmax": 500} },
    )

hBJetEta_All = ControlPlotInput(
    title      = "BJetEta_AfterAllSelections",
    histoName  = "BJetEta_AfterAllSelections",
    details    = { "xlabel": "b jet #eta",
                   "ylabel": "Events",
                   "divideByBinWidth": False,
                   "unit": "",
                   "log": True,
                   "legendPosition": "SW",
                   "opts": {"ymin": 0.09} 
                   },
    )

hBtatDisc_All  = ControlPlotInput(
    title     = "BtagDiscriminator_AfterAllSelections",
    histoName = "BtagDiscriminator_AfterAllSelections",
    details   = { "xlabel": "b tag discriminator",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "",
                  "log": True,
                  "legendPosition": "SE",
                  "opts": {"ymin": 0.9} 
                  },
    )

hMET_All = ControlPlotInput(
    title     = "MET_AfterAllSelections",
    histoName = "MET_AfterAllSelections",
    details   = { "xlabel": "E_{T}^{miss}",
                  "ylabel": "Events/^{}#DeltaE_{T}^{miss}",
                  "divideByBinWidth": True,
                  "unit": "GeV",
                  "log": True,
                  "opts": {"ymin": 0.00009, "ymaxfactor": 10, "xmax": 400} 
                  },
    flowPlotCaption  = "^{}E_{T}^{miss}", # Leave blank if you don't want to include the item to the selection flow plot
    )

hMETPhi_All = ControlPlotInput(
    title     = "METPhi_AfterAllSelections",
    histoName = "METPhi_AfterAllSelections",
    details   = { "xlabel": "E_{T}^{miss} #phi",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "{}^{o}",
                  "log": True,
                  "legendPosition": "SW",
                  "opts": {"ymin": 0.09} 
                  }
    )

hDeltaPhiTauMet_All = ControlPlotInput(
    title     = "DeltaPhiTauMet_AfterAllSelections",
    histoName = "DeltaPhiTauMet_AfterAllSelections",
    details   = {"xlabel": "#Delta#phi(#tau,E_{T}^{miss})",
                 "ylabel": "Events",
                 "divideByBinWidth": False,
                 "unit": "{}^{o}",
                 "log": True,
                 "legendPosition": "NE",
                 "opts": {"ymin": 0.9} 
                 }
    )

hBacktoBackAngularCutsMin_All = ControlPlotInput(
    title     = "BackToBackAngularCutsMinimum_AfterAllSelections",
    histoName = "BackToBackAngularCutsMinimum_AfterAllSelections",
    details   = { "xlabel": "^{}R_{bb}^{min}",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "{}^{o}",
                  "log": True,
                  "legendPosition": "SE",
                  "opts": {"ymin": 0.09} 
                  },
    flowPlotCaption  = "^{}R_{bb}^{min}", # Leave blank if you don't want to include the item to the selection flow plot
    )

hNJets_BtagSF = ControlPlotInput(
    title     = "Njets_AfterBtagSF",
    histoName = "Njets_AfterBtagSF",
    details   = { "xlabel": "Number of selected jets",
                  "ylabel": "Events",
                  "divideByBinWidth": False,
                  "unit": "",
                  "log": True,
                  "opts": {"ymin": 0.9, "xmin": 3.0, "xmax": 10.0} 
                  },
    )

JetPt_BtagSF = ControlPlotInput(
    title     = "JetPt_AfterBtagSF",
    histoName = "JetPt_AfterBtagSF",
    details   = { "xlabel": "jet ^{}p_{T}",
                  "ylabel": "Events/^{}#Deltap_{T}",
                  "divideByBinWidth": True,
                  "unit": "GeV/c",
                  "log": True,
                  "opts": {"ymin": 0.0009, "ymaxfactor": 10, "xmax": 500} 
                  },
    )

BJetPt_BtagSF = ControlPlotInput(
    title     = "BJetPt_AfterBtagSF",
    histoName = "BJetPt_AfterBtagSF",
    details   = { "xlabel": "b jet ^{}p_{T}",
                  "ylabel": "Events/^{}#Deltap_{T}",
                  "divideByBinWidth": True,
                  "unit": "GeV/c",
                  "log": True,
                  "opts": {"ymin": 0.0009, "ymaxfactor": 10, "xmax": 500} 
                  },
    )

hTransverseMass = ControlPlotInput(
    title       = "shapeTransverseMass",
    histoName   = "shapeTransverseMass", # OptionMassShape,
    details     = { "xlabel"             : "m_{T}",
                    "ylabel"             : "< Events / %s >" % (uMass),
                    "divideByBinWidth"   : True,
                    "unit"               : "%s" % (uMass),
                    "log"                : False,
                    "legendPosition"     : "NE",
                    "ratioLegendPosition": "right",
                    "opts"               : {"ymin": 0.0, "ymaxfactor": 1.1, "xmax": 500},
                    "opts2"              : {"ymin": 0.0, "ymax": 2.0}
                    },
    blindedRange=[81.0, 6000.0], # specify range min,max if blinding applies to this control plot
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot    
    # flowPlotCaption  = "m_{jjbb}", # Leave blank if you don't want to include the item to the selection flow plot    
    )

hTransverseMassLog = ControlPlotInput(
    title     = "TransverseMassLog",
    histoName = "shapeTransverseMass",
    details   = { "xlabel": "m_{T} (GeV)",
                  "cmsTextPosition": "right",
                  "ylabel": "< Events / bin >", "ylabelBinInfo": False,
                  "divideByBinWidth": True,
                  "log": True,
                  "opts": {"ymin": 1e-4, "ymaxfactor": 15, "xmax": 500},
                  "opts2": {"ymin": 0.0, "ymax": 2.0}
                  },
    blindedRange=[81, 6000], # specify range min,max if blinding applies to this control plot
    )

hTransverseMassLogx = ControlPlotInput(
    title     = "TransverseMassLogx",
    histoName = "shapeTransverseMass",
    details   = { "xlabel": "m_{T} (GeV)",
                  "cmsTextPosition": "right",
                  "ylabel": "< Events / bin >", "ylabelBinInfo": False,
                  "divideByBinWidth": True,
                  "logx": True,
                  "opts": {"ymin": 1e-4,"ymaxfactor": 1.3,"xmax": 5000},
                  "opts2": {"ymin": 0.0, "ymax": 2.0}
                  }, 
    blindedRange=[81, 6000], # specify range min,max if blinding applies to this control plot
    )

hTransverseMassLogxLog = ControlPlotInput(
    title     = "TransverseMassLogxLog",
    histoName = "shapeTransverseMass",
    details   = { "xlabel": "m_{T} (GeV)",
                  "cmsTextPosition": "right",
                  "ylabel": "< Events / bin >", "ylabelBinInfo": False,
                  "moveLegend": {"dx": -0.10, "dy": -0.12, "dh":0.1},
                  "ratioMoveLegend": {"dx": -0.06, "dy": -0.33},
                  "divideByBinWidth": True,
                  "log": True,
                  "logx": True,
                  "opts": {"xmin": 9, "ymin": 1e-4,"ymaxfactor": 10.0},
                  "opts2": {"ymin": 0.0, "ymax": 2.0}
                  }, 
    blindedRange=[81, 6000], # specify range min,max if blinding applies to this control plot
    )


#================================================================================================ 
# Create the list of control plots (NOTE: Set "OptionDoControlPlots" to True)
#================================================================================================ 
if 0:
    ControlPlots.append(hVertices_Std)
    ControlPlots.append(hTauPt_Std)
    ControlPlots.append(hTauLdgTrkPt_Std)
    ControlPlots.append(hTauEta_Std)
    ControlPlots.append(hTauPhi_Std)
    ControlPlots.append(hTauRtauFullRange_Std)
    ControlPlots.append(hTauRtau_Std)
    # ControlPlots.append(hTauDM_Std)
    ControlPlots.append(hTauNprongs_Std)
    # ControlPlots.append(hTauSource_Std)
    ControlPlots.append(hNJets_Std)
    ControlPlots.append(hJetPt_Std)
    ControlPlots.append(hJetEta_Std)
    ControlPlots.append(hCollAngularCutsMin)
    ControlPlots.append(hNBjets)
    ControlPlots.append(hBJetPt)
    ControlPlots.append(hBJetEta)
    ControlPlots.append(hBtatDisc) 
    ControlPlots.append(hMET)
    ControlPlots.append(hMETPhi)
    ControlPlots.append(hBacktoBackAngularCutsMin)
if 1:
    ControlPlots.append(hVertices_All)
    ControlPlots.append(hTauPt_All)
    ControlPlots.append(hTauLdgTrkPt_All)
    ControlPlots.append(hTauEta_All)
    ControlPlots.append(hTauPhi_All)
    ControlPlots.append(hTauRtauFullRange_All)
    ControlPlots.append(hTauRtau_All)
    # ControlPlots.append(hTauDM_All)
    ControlPlots.append(hTauNprongs_All)
    # ControlPlots.append(hTauSource_All)
    ControlPlots.append(hNJets_All)
    ControlPlots.append(hJetPt_All)
    ControlPlots.append(hJetEta_All)
    ControlPlots.append(hCollAngularCutsMin_All)
    ControlPlots.append(hNBjets_All)
    ControlPlots.append(hBJetPt_All)
    ControlPlots.append(hBJetEta_All)
    ControlPlots.append(hBtatDisc_All) 
    ControlPlots.append(hMET_All)
    ControlPlots.append(hMETPhi_All)
    ControlPlots.append(hDeltaPhiTauMet_All)
    ControlPlots.append(hBacktoBackAngularCutsMin_All)
if 0:
     ControlPlots.append(hNJets_BtagSF)
     ControlPlots.append(JetPt_BtagSF)
     ControlPlots.append(BJetPt_BtagSF)

ControlPlots.append(hTransverseMass)
ControlPlots.append(hTransverseMassLog)
ControlPlots.append(hTransverseMassLogx)
ControlPlots.append(hTransverseMassLogxLog)


if OptionTest:
    ControlPlots = []
    ControlPlots.append(hNJets_Std) # flow
    ControlPlots.append(hNBjets)    # flow
    ControlPlots.append(hMET)       # flow
    ControlPlots.append(hBacktoBackAngularCutsMin) #flow
    #ControlPlots.append(hDeltaPhiTauMet_All)
    ControlPlots.append(hTransverseMass)
