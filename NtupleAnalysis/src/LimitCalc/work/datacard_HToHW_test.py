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
./dcardGenerator_v2.py -x dcardDefault_h2hw_mc.py -d limits2019/ --h2hw


HN Threads:
https://hypernews.cern.ch/HyperNews/CMS/get/HIG-18-015/11.html

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
OptionTest                             = False  # [default: False]
OptionPaper                            = True  # [default: True]   (Changes figure style to paper style)
OptionIncludeSystematics               = True  # [default: True]   (Shape systematics; Requires pseudo-multicrab produced with doSystematics=True) 
OptionShapeSystematics                 = False # [default: True]   (Shape systematics; Requires pseudo-multicrab produced with doSystematics=True) 
OptionDoControlPlots                   = True  # [default: True]   (Produce control plots defined at end of this file)
OptionGenuineTauBackgroundSource       = "MC"  # [Options: "DataDriven", "MC"]
OptionFakeTauMeasurementSource         = "MC"  # [default: "DataDriven"] (options: "DataDriven", "MC")
OptionBr                               = 1.0   # [default: 1.0]    (The Br(t->bH+) used in figures and tables)
OptionSqrtS                            = 13    # [default: 13]     (The sqrt(s) used in figures and tables)
OptionBlindThreshold                   = 0.10  # [default: None]   (If signal exceeds this fraction of expected events, data is blinded in a given bin)
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
MassPoints                             = [300, 700]       # [default: [300, 700]] (Mass points to be considered)
PlotLegendHeader                       = "H^{#pm}#rightarrowW^{#pm}H^{0}, H^{0}#rightarrow #tau^{+} #tau^{-}"  #tau^{+}_{h} #tau^{-}_{h}"
SignalName                             = "ChargedHiggs_HplusTB_HplusToHW_M%s_mH200_2ta_NLO"
DataCardName                           = "HToHW"          # [default: Hplus2hw_13TeV]  (Used by TableProducer.py)
BlindAnalysis                          = False             # [default: True]   (Change only if "green light" for unblinding)
MinimumStatUncertainty                 = 0.5              # [default: 0.5]    (Minimum stat. uncert. to set to bins with zero events)
UseAutomaticMinimumStatUncertainty     = False            # [default: False]  (Do NOT use the MinimumStatUncertainty; determine value from lowest non-zero rate for each dataset   )
ToleranceForLuminosityDifference       = 0.05             # [default: 0.05]   (Tolerance for throwing error on luminosity difference; "0.01" means that a 1% is required) 
ToleranceForMinimumRate                = 0.0              # [default: 0.0]    (Tolerance for almost zero rate columns with smaller rate are suppressed) 
labelPrefix                            = ""               # [default: ""]     (Prefix for the labels of datacard columns; e.g. "CMS_Hptntj_", "CMS_H2tb_")
labelPostfix                           = ""               # [default: "_GenuineTau"] (Postfix for the labels of datacard columns; e.g. "TT" --> "TT_GenuineTau")
if OptionTest:
    MassPoints = [700]

#================================================================================================  
# File-specific settings
#================================================================================================  
 # Get the binning for the shape histogram
ShapeHistogramsDimensions = systematics.getBinningForPlot("Mt", "HToHW")

# Counter and histogram path definitions
histoPathInclusive   = "ForDataDrivenCtrlPlots"
histoPath_GenuineTau = "ForDataDrivenCtrlPlotsEWKGenuineTaus"
histoPath_FakeTau    = "ForDataDrivenCtrlPlotsEWKFakeTaus"

# Observation definition (how to retrieve number of observed events)
Observation = ObservationInput(datasetDefinition="Data", shapeHistoName=OptionMassShape, histoPath=histoPathInclusive)

#================================================================================================  
# Nuisance Lists (Just the strings; The objects are defined later below)
#================================================================================================ 
myLumiSystematics       = ["lumi_13TeV"]
myPileupSystematics     = ["CMS_pileup"]
myTopTagSystematics     = ["CMS_HPTB_toptagging"]
myTrgEffSystematics     = ["CMS_eff_trg_MC"]
myLeptonVetoSystematics = ["CMS_eff_e_veto", "CMS_eff_m", "CMS_eff_tau"]
myJetSystematics        = ["CMS_scale_j", "CMS_res_j"]
myBtagSystematics       = ["CMS_eff_b"]

# Define systematics dictionary (easier access)
mySystematics = {}
mySystematics["MC"]     = myLumiSystematics + myPileupSystematics + myTrgEffSystematics + myLeptonVetoSystematics + myJetSystematics + myBtagSystematics + myTopTagSystematics
mySystematics["Signal"] = mySystematics["MC"] + ["CMS_HPTB_mu_RF_HPTB","CMS_HPTB_pdf_HPTB"]
mySystematics["TT"]     = mySystematics["MC"] + ["QCDscale_ttbar", "pdf_ttbar", "mass_top"] + ["CMS_HPTB_mu_RF_top","CMS_HPTB_pdf_top"]
mySystematics["ttX"]    = mySystematics["MC"] + ["QCDscale_singleTop", "pdf_singleTop", "mass_top_forSingleTop"] + ["CMS_HPTB_mu_RF_top","CMS_HPTB_pdf_top"]
mySystematics["EWK"]    = mySystematics["MC"] + ["QCDscale_ewk", "pdf_ewk"] + ["CMS_HPTB_mu_RF_ewk","CMS_HPTB_pdf_ewk"]

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
TT = DataGroup(label             = labelPrefix + "TT" + labelPostfix, #
               landsProcess      = 2,
               shapeHistoName    = OptionMassShape, 
               histoPath         = histoPathInclusive,
               datasetType       = "EWKMC",
               datasetDefinition = "TT", # You must use the exact name given in plots.py (_datasetMerge dictionary)
               validMassPoints   = MassPoints,
               nuisances         = mySystematics["TT"]
               )

TTX = DataGroup(label                   = labelPrefix + "ttX" + labelPostfix,
                      landsProcess      = 3,
                      shapeHistoName    = OptionMassShape,
                      histoPath         = histoPathInclusive,
                      datasetType       = "EWKMC",
                      datasetDefinition = "ttX", # You must use the exact name given in plots.py (_datasetMerge dictionary)
                      validMassPoints   = MassPoints,
                      nuisances         = mySystematics["ttX"]
                      )

DYJets = DataGroup(label             = labelPrefix + "DYJets" + labelPostfix,
                   landsProcess      = 4,
                   shapeHistoName    = OptionMassShape,
                   histoPath         = histoPathInclusive,
                   datasetType       = "EWKMC",
                   datasetDefinition = "DYJetsToLLHT", # You must use the exact name given in plots.py
                   validMassPoints   = MassPoints,
                   nuisances         = mySystematics["MC"]
                   )

Diboson = DataGroup(label             = labelPrefix + "Diboson" + labelPostfix,
                    landsProcess      = 5,
                    shapeHistoName    = OptionMassShape,
                    histoPath         = histoPathInclusive,
                    datasetType       = "EWKMC",
                    datasetDefinition = "Diboson", # You must use the exact name given in plots.py
                    validMassPoints   = MassPoints,
                    nuisances         = mySystematics["MC"]
                   )

TT_GenuineTau = DataGroup(label             = labelPrefix + "TT" + labelPostfix, #
                          landsProcess      = 2,
                          shapeHistoName    = OptionMassShape, 
                          histoPath         = histoPath_GenuineTau,
                          datasetType       = "GenuineTau",
                          datasetDefinition = "TT", # You must use the exact name given in plots.py (_datasetMerge dictionary)
                          validMassPoints   = MassPoints,
                          nuisances         = mySystematics["TT"]
                          )

TTX_GenuineTau = DataGroup(label             = labelPrefix + "ttX" + labelPostfix,
                           landsProcess      = 3,
                           shapeHistoName    = OptionMassShape,
                           histoPath         = histoPath_GenuineTau, 
                           datasetType       = "GenuineTau",
                           datasetDefinition = "ttX", # You must use the exact name given in plots.py (_datasetMerge dictionary)
                           validMassPoints   = MassPoints,
                           nuisances         = mySystematics["ttX"]
                           )

DYJets_GenuineTau = DataGroup(label             = labelPrefix + "DYJets" + labelPostfix,
                              landsProcess      = 4,
                              shapeHistoName    = OptionMassShape,
                              histoPath         = histoPath_GenuineTau,
                              datasetType       = "GenuineTau",
                              datasetDefinition = "DYJetsToLLHT", # You must use the exact name given in plots.py
                              validMassPoints   = MassPoints,
                              nuisances         = mySystematics["MC"]
                              )

Diboson_GenuineTau = DataGroup(label             = labelPrefix + "Diboson" + labelPostfix,
                               landsProcess      = 5,
                               shapeHistoName    = OptionMassShape,
                               histoPath         = histoPath_GenuineTau,
                               datasetType       = "GenuineTau",
                               datasetDefinition = "Diboson", # You must use the exact name given in plots.py
                               validMassPoints   = MassPoints,
                               nuisances         = mySystematics["MC"]
                               )

TT_FakeTau = DataGroup(label             = labelPrefix + "TT" + labelPostfix, #
                       landsProcess      = 6,
                       shapeHistoName    = OptionMassShape, 
                       histoPath         = histoPath_FakeTau,
                       datasetType       = "FakeTau",
                       datasetDefinition = "TT", # You must use the exact name given in plots.py (_datasetMerge dictionary)
                       validMassPoints   = MassPoints,
                       nuisances         = mySystematics["TT"]
                       )

TTX_FakeTau = DataGroup(label             = labelPrefix + "ttX" + labelPostfix,
                        landsProcess      = 7,
                        shapeHistoName    = OptionMassShape,
                        histoPath         = histoPath_FakeTau, 
                        datasetType       = "FakeTau",
                        datasetDefinition = "ttX", # You must use the exact name given in plots.py (_datasetMerge dictionary)
                        validMassPoints   = MassPoints,
                        nuisances         = mySystematics["ttX"]
                        )

DYJets_FakeTau = DataGroup(label             = labelPrefix + "DYJets" + labelPostfix,
                           landsProcess      = 8,
                           shapeHistoName    = OptionMassShape,
                           histoPath         = histoPath_FakeTau,
                           datasetType       = "FakeTau",
                           datasetDefinition = "DYJetsToLLHT", # You must use the exact name given in plots.py
                           validMassPoints   = MassPoints,
                           nuisances         = mySystematics["MC"]
                           )

Diboson_FakeTau = DataGroup(label             = labelPrefix + "Diboson" + labelPostfix,
                            landsProcess      = 9,
                            shapeHistoName    = OptionMassShape,
                            histoPath         = histoPath_FakeTau,
                            datasetType       = "FakeTau",
                            datasetDefinition = "Diboson", # You must use the exact name given in plots.py
                            validMassPoints   = MassPoints,
                            nuisances         = mySystematics["MC"]
                            )

# Append datasets in order you want them to appear in the data-driven control plot stack
# The list \"DataGroups\" is required by the DataCardGenerator.py module
DataGroups = []
DataGroups.extend(signalDataGroups)
#DataGroups.append(TT)
#DataGroups.append(TTX)
#DataGroups.append(DYJets)
#DataGroups.append(Diboson)
#
DataGroups.append(TT_FakeTau)
DataGroups.append(TTX_FakeTau)
DataGroups.append(DYJets_FakeTau)
DataGroups.append(Diboson_FakeTau)
#
DataGroups.append(TT_GenuineTau)
DataGroups.append(TTX_GenuineTau)
DataGroups.append(DYJets_GenuineTau)
DataGroups.append(Diboson_GenuineTau)

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
topTag_Shape     = Nuisance(id="CMS_HPTB_toptagging", label="TopTag", distr="shapeQ", function="ShapeVariation", systVariation="TopTagSF")
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
eVeto_Const     = Nuisance(id="CMS_eff_e_veto"     , label="e veto", distr="lnN", function="Ratio", numerator="passed e selection (Veto)", denominator="passed PV", scaling=0.02)
muVeto_Const    = Nuisance(id="CMS_eff_m"          , label="mu ID (Approx.)" , distr="lnN", function="Constant", value=0.10)
tauVeto_Const   = Nuisance(id="CMS_eff_tau"        , label="tau ID (Approx.)", distr="lnN", function="Constant", value=0.10)
bTagSF_Const    = Nuisance(id="CMS_eff_b"          , label="b tagging (Approx.)", distr="lnN", function="Constant", value=0.05)
JES_Const       = Nuisance(id="CMS_scale_j"        , label="Jet Energy Scale (JES) (Approx.)"     , distr="lnN", function="Constant", value=0.03)
JER_Const       = Nuisance(id="CMS_res_j"          , label="Jet Energy Resolution (JER) (Approx.)", distr="lnN", function="Constant", value=0.04)
topPt_Const     = Nuisance(id="CMS_topreweight"    , label="Top pT reweighting (Approx.)", distr="lnN", function="Constant", value=0.25)
topTag_Const    = Nuisance(id="CMS_HPTB_toptagging", label="Top tagging (Approx.)", distr="lnN", function="Constant", value=0.10)

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
RF_QCDscale_top_const  = Nuisance(id="CMS_HPTB_mu_RF_top" , label="Scale acceptance uncertainty for top"   , distr="lnN", function="Constant",value=0.02)
RF_QCDscale_ewk_const  = Nuisance(id="CMS_HPTB_mu_RF_ewk" , label="Scale acceptance uncertainty for EWK"   , distr="lnN", function="Constant",value=0.05)
RF_QCDscale_HPTB_const = Nuisance(id="CMS_HPTB_mu_RF_HPTB", label="Scale acceptance uncertainty for signal", distr="lnN", function="Constant",value=0.048)
#RF_QCDscale_HPTB_const = Nuisance(id="CMS_HPTB_mu_RF_HPTB_heavy", label="QCDscale acceptance uncertainty for signal", distr="lnN", function="Constant",value=0.012)

#==== Acceptance uncertainties  (PDF)
RF_pdf_top_const  = Nuisance(id="CMS_HPTB_pdf_top", label="PDF acceptance uncertainty for top", distr="lnN", function="Constant",value=0.02,upperValue=0.0027)
RF_pdf_ewk_const  = Nuisance(id="CMS_HPTB_pdf_ewk", label="PDF acceptance uncertainty for EWK", distr="lnN", function="Constant",value=0.033,upperValue=0.046)
RF_pdf_HPTB_const = Nuisance(id="CMS_HPTB_pdf_HPTB", label="PDF acceptance uncertainty for signal", distr="lnN", function="Constant",value=0.004,upperValue=0.017)

# Fake-b nuisances
tf_FakeTau_Const          = Nuisance(id="CMS_HPTB_fakeB_transferfactor", label="Transfer Factor uncertainty", distr="lnN", function="Constant", value=0.10)
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
topTag_FakeTau_Const      = Nuisance(id="CMS_HPTB_toptagging_forFakeTau", label="Top tagging (Approx.)", distr="lnN", function="ConstantForFakeTau", value=0.20)
ttbar_scale_FakeTau_Const = Nuisance(id="QCDscale_ttbar_forFakeTau"  , label="QCD XSection uncertainties", distr="lnN", function="ConstantForFakeTau", value=tt_scale_down, upperValue=tt_scale_up)
ttbar_pdf_FakeTau_Const   = Nuisance(id="pdf_ttbar_forFakeTau"       , label="TTbar XSection pdf uncertainty", distr="lnN", function="ConstantForFakeTau", value=tt_pdf_down, upperValue=tt_pdf_up)
ttbar_mass_FakeTau_Const  = Nuisance(id="mass_top_forFakeTau"        , label="TTbar XSection top mass uncertainty", distr="lnN", function="ConstantForFakeTau", value=tt_mass_down, upperValue=tt_mass_up) 
RF_QCDscale_FakeTau_const = Nuisance(id="CMS_HPTB_mu_RF_top_forFakeTau", label="Scale acceptance uncertainty for FakeTau" , distr="lnN", function="ConstantForFakeTau",value=0.02)
RF_pdf_FakeTau_const      = Nuisance(id="CMS_HPTB_pdf_top_forFakeTau"  , label="PDF acceptance uncertainty for FakeTau"      , distr="lnN", function="ConstantForFakeTau",value=0.02, upperValue=0.0027)


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
    Nuisances.append(topTag_Shape)
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
    Nuisances.append(topTag_Const)
    Nuisances.append(topTag_FakeTau_Const)

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
Nuisances.append(RF_QCDscale_HPTB_const)
#Nuisances.append(RF_QCDscale_HPTB_const)
Nuisances.append(RF_QCDscale_FakeTau_const)
Nuisances.append(RF_pdf_top_const)
Nuisances.append(RF_pdf_ewk_const)
Nuisances.append(RF_pdf_HPTB_const)
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
MergeNuisances.append(["CMS_HPTB_mu_RF_top", "CMS_HPTB_mu_RF_top_forFakeTau"])
MergeNuisances.append(["CMS_HPTB_pdf_top", "CMS_HPTB_pdf_top_forFakeTau"])

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
yMin    = 0.5e-1 #0.5e-2
if OptionPaper:
    uPt   = "GeV"
    uMass = "GeV"

from HiggsAnalysis.LimitCalc.InputClasses import ControlPlotInput
ControlPlots= []

hMET = ControlPlotInput(
    title     = "Met", # MET_AfterAllSelections" (Used in determining binning from systematics.py)
    histoName = "MET_AfterAllSelections",
    details   = { "xlabel"             : "E_{T}^{miss}",
                  "ylabel"             : "< Events / %s >" % (uEnergy),
                  "divideByBinWidth"   : True,
                  "unit"               : "GeV",
                  "log"                : True,
                  "legendPosition"     : "NE",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": yMin, "ymaxfactor": 10}#, "xmax": 300.0}
                  },
    blindedRange=[100.0, 1000.0], # specify range min,max if blinding applies to this control plot      
    )

hHT = ControlPlotInput(
    title     = "HT", # (Used in determining binning from systematics.py)
    histoName = "HT_AfterAllSelections",
    details   = { "xlabel"             : "H_{T}",
                  "ylabel"             : "< Events / %s >" % (uEnergy),
                  "divideByBinWidth"   : True,
                  "unit"               : "GeV",
                  "log"                : True,
                  "legendPosition"     : "NE",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": yMin, "ymaxfactor": 10}#, "xmax": 3000.0} },
                  },
    blindedRange=[200.0, 1500.0], # specify range min,max if blinding applies to this control plot      
    )

hMHT = ControlPlotInput(
    title      = "MHT",
    histoName  = "MHT_AfterAllSelections",
    details    = { "xlabel"             : "MHT",                  
                   "ylabel"             : "< Events / %s >" % (uEnergy),
                   "divideByBinWidth"   : True,
                   "unit"               : "GeV",
                   "log"                : True,
                   "legendPosition"     : "NE",
                   "ratioLegendPosition": "right",
                   "opts"               : {"ymin": yMin, "ymaxfactor": 10}#, "xmax": 400.0} }
                   },
    )

hTopPt = ControlPlotInput(
    title     = "TopPt",
    histoName = "TopPt_AfterAllSelections",
    details   = { "xlabel"             : "p_{T}",
                  "ylabel"             : "< Events / %s >" % (uPt),
                  "divideByBinWidth"   : True,
                  "unit"               : "%s" % (uPt),
                  "log"                : True,
                  "legendPosition"     : "NE",
                  "ratioLegendPosition": "SE", #"right",
                  "opts"               : {"ymin": yMin, "ymaxfactor": 10, "xmax": 800.0}
                         },
    )

hTopMass = ControlPlotInput(
    title     = "TopMass",
    histoName = "TopMass_AfterAllSelections",
    details   = { "xlabel"             : "m_{jjb}",
                  "ylabel"             : "<Events / %s>" % (uMass),
                  "divideByBinWidth"   : True,
                  "unit"               : "%s" % (uMass),
                  "log"                : True,
                  "legendPosition"     : "NE",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": yMin, "ymaxfactor": 10}#, "xmax": 350.0} 
                  },
    )

hTopBjetPt = ControlPlotInput(
    title     = "TopBJetPt",
    histoName = "TopBjetPt_AfterAllSelections",
    details   = { "xlabel"             : "p_{T}",
                  "ylabel"             : "< Events / %s>" % (uPt),
                  "divideByBinWidth"   : True,
                  "unit"               : "%s" % (uPt),
                  "log"                : True,
                  "legendPosition"     : "NE",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": yMin, "ymaxfactor": 10}#, "xmax": 700.0} }
                  },
    )

hTopBjetEta = ControlPlotInput(
    title     = "TopBJetEta",
    histoName = "TopBjetEta_AfterAllSelections",
    details   = { "xlabel"             : "#eta",
                  "ylabel"             : "Events",
                  "divideByBinWidth"   : False,
                  "unit"               : "",
                  "log"                : True,
                  "legendPosition"     : "RM",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": yMin, "ymaxfactor": 10, "xmin": -2.5, "xmax": 2.5} 
                  },
    )

hTopBjetBdisc = ControlPlotInput(
    title     = "TopBJetBdisc",
    histoName = "TopBjetBdisc_AfterAllSelections",
    details   = { "xlabel"             : "CSVv2 discriminator",
                  "ylabel"             : "Events",
                  "divideByBinWidth"   : False,
                  "unit"               : "",
                  "log"                : True,
                  "legendPosition"     : "NW",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": yMin, "ymaxfactor": 10, "xmin": 0.7}
                  },
    )


hTopDijetPt = ControlPlotInput(
    title     = "DijetPt",
    histoName = "TopDijetPt_AfterAllSelections",
    details   = { "xlabel"             : "p_{T}",
                  "ylabel"             : "< Events / %s >" % (uPt),
                  "divideByBinWidth"   : True,
                  "unit"               : "%s" % (uPt),
                  "log"                : True,
                  "legendPosition"     : "NE",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": yMin, "ymaxfactor": 10}#, "xmax": 700.0} 
                  },
    )

hTopDijetMass = ControlPlotInput(
    title     = "DijetMass",
    histoName = "TopDijetMass_AfterAllSelections",
    details   = { "xlabel"             : "m_{jj}",
                  "ylabel"             : "< Events / %s>" % (uMass),
                  "divideByBinWidth"   : True,
                  "unit"               : "%s" % (uMass),
                  "log"                : True,
                  "legendPosition"     : "NE",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": yMin}
                         },
    )

hTransverseMass = ControlPlotInput(
    title       = "Mt",
    histoName   = OptionMassShape,
    details     = { "xlabel"             : "m_{T}",
                    "ylabel"             : "< Events / %s >" % (uMass),
                    "divideByBinWidth"   : True,
                    "unit"               : "%s" % (uMass),
                    "log"                : True,
                    "legendPosition"     : "NE",
                    "ratioLegendPosition": "right",
                    "opts"               : {"ymin": yMin, "ymaxfactor": 10}
                    },
    blindedRange=[0.0, 2500.0], # specify range min,max if blinding applies to this control plot
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot    
    # flowPlotCaption  = "m_{jjbb}", # Leave blank if you don't want to include the item to the selection flow plot    
    )

hVertices = ControlPlotInput(
    title     = "Vertices",
    histoName = "NVertices_AfterAllSelections",
    details   = { "xlabel"             : "vertex multiplicity",
                  "ylabel"             : "< Events >",
                  "divideByBinWidth"   : True,
                  "unit"               : "",
                  "log"                : True,
                  "legendPosition"     : "NE",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": yMin, "ymaxfactor": 10, "xmax": 105.0} 
                  },
    flowPlotCaption  = "", # Leave blank if you don't want to include the item to the selection flow plot    
    )

hNjets = ControlPlotInput(
    title     = "NJets",
    histoName = "Njets_AfterAllSelections",
    details   = { "xlabel"             : "jets multiplicity",
                  "ylabel"             : "Events",
                  "divideByBinWidth"   : False,
                  "unit"               : "",
                  "log"                : True,
                  "legendPosition"     : "NE",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": yMin, "ymaxfactor": 10}
                  },
    #blindedRange=[100.0, 400.0], # specify range min,max if blinding applies to this control plot      
    )

hNBjets = ControlPlotInput(
    title     = "NBjets",
    histoName = "NBjets_AfterAllSelections",
    details   = { "xlabel"             : "b-jets multiplicity",
                  "ylabel"             : "Events",
                  "divideByBinWidth"   : False,
                  "unit"               : "",
                  "log"                : True,
                  "legendPosition"     : "NE",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": yMin, "ymaxfactor": 10}
                  },
    #blindedRange=[100.0, 400.0], # specify range min,max if blinding applies to this control plot      
    )

hJetPt = ControlPlotInput(
    title     = "JetPt",
    histoName = "JetPt_AfterAllSelections",
    details   = { "xlabel"             : "p_{T}",
                  "ylabel"             : "< Events / %s >" % (uPt),
                  "divideByBinWidth"   : True,
                  "unit"               : "%s" % (uPt),
                  "log"                : True,
                  "legendPosition"     : "NE",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": yMin, "ymaxfactor": 10}
                  },
    #blindedRange=[100.0, 400.0], # specify range min,max if blinding applies to this control plot
    )

hJetEta = ControlPlotInput(
    title     = "JetEta",
    histoName = "JetEta_AfterAllSelections",
    details   = { "xlabel"             : "#eta",
                  "ylabel"             : "Events",
                  "divideByBinWidth"   : False,
                  "unit"               : "",
                  "log"                : True,
                  "legendPosition"     : "RM",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": yMin, "ymaxfactor": 10}}#,
    #blindedRange=[100.0, 400.0], # specify range min,max if blinding applies to this control plot      
    )

hBJetPt = ControlPlotInput(
    title     = "BJetPt",
    histoName = "BJetPt_AfterAllSelections",
    details   = { "xlabel"             : "p_{T}",
                  #"ylabel"             : "Events / #Deltap_{T}",
                  "ylabel"             : "< Events / %s >" % (uPt),
                  "divideByBinWidth"   : True,
                  "unit"               : "%s" % (uPt),
                  "log"                : True,
                  "legendPosition"     : "NE",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": yMin, "ymaxfactor": 10}}#,
    #blindedRange=[100.0, 400.0], # specify range min,max if blinding applies to this control plot      
    )

hBJetEta = ControlPlotInput(
    title     = "BJetEta",
    histoName = "BJetEta_AfterAllSelections",
    details   = { "xlabel"             : "#eta",
                  "ylabel"             : "Events",
                  "divideByBinWidth"   : False,
                  "unit"               : "",
                  "log"                : True,
                  "legendPosition"     : "RM",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": yMin, "ymaxfactor": 10}}#,
    #blindedRange=[100.0, 400.0], # specify range min,max if blinding applies to this control plot      
    )

hBtagDiscriminator = ControlPlotInput(
    title     = "BtagDisc",
    histoName = "BtagDiscriminator_AfterAllSelections",
    details   = { "xlabel"             : "CSVv2 discriminator",
                  "ylabel"             : "Events",
                  "divideByBinWidth"   : False,
                  "unit"               : "",
                  "log"                : True,
                  "legendPosition"     : "NW",
                  "ratioLegendPosition": "right",
                  "opts"               : {"ymin": yMin, "ymaxfactor": 10,  "xmin": 0.5}}
    #blindedRange=[100.0, 400.0], # specify range min,max if blinding applies to this control plot      
    )

#================================================================================================ 
# Create the list of control plots (NOTE: Set "OptionDoControlPlots" to True)
#================================================================================================ 
ControlPlots.append(hMET)
ControlPlots.append(hHT)
ControlPlots.append(hMHT)
#ControlPlots.append(hTopPt)
#ControlPlots.append(hTopMass)
#ControlPlots.append(hTopBjetPt)
#ControlPlots.append(hTopBjetEta)
#ControlPlots.append(hTopBjetBdisc)
#ControlPlots.append(hTopDijetPt)
#ControlPlots.append(hTopDijetMass)
ControlPlots.append(hVertices)
ControlPlots.append(hNjets)
ControlPlots.append(hJetPt)
ControlPlots.append(hJetEta)
ControlPlots.append(hNBjets)
ControlPlots.append(hBtagDiscriminator)
ControlPlots.append(hBJetPt)
ControlPlots.append(hBJetEta)
ControlPlots.append(hTransverseMass)

if OptionTest:
    ControlPlots = []
    #ControlPlots.append(hTopPt)
    #ControlPlots.append(hTopMass)
    ControlPlots.append(hNjets)
    ControlPlots.append(hMET)
    ControlPlots.append(hHT)
    #ControlPlots.append(hTransverseMass)
