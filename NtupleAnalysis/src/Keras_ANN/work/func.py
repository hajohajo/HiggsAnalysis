#================================================================================================
# Import modules
#================================================================================================   
import ROOT
import plot
import math
import array
import json
import pandas
import numpy 

from sklearn.utils.class_weight import compute_sample_weight
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.externals import joblib

#================================================================================================   
# Function definition
#================================================================================================   
def Print(msg, printHeader=False):
    fName = __file__.split("/")[-1]
    if printHeader==True:
        print "=== ", fName
        print "\t", msg
    else:
        print "\t", msg
    return

def Verbose(msg, printHeader=False, verbose=False):
    if not verbose:
        return
    fName = __file__.split("/")[-1]
    if printHeader==True:
        print "=== ", fName
        print "\t", msg
    else:
        print "\t", msg
    return


def split_list(a_list, firstHalf=True):
    half = len(a_list)//2
    if firstHalf:
        return a_list[:half]
    else:
        return a_list[half:]

def doSampleReweighing(df_sig, df_bkg, varName, _verbose=False, **kwargs):
    '''
     Optional numpy array of weights for the training & testing samples, used for weighting the loss function (during training only)
     Can be used to decorrelate a variable before training; alternative approach to adversarial neural network (ANN) for combatting mass sculpting

    df_sig is the signal dataFrame
    df_bkg is the bkg dataFrame
    
    '''
    msg = "De-correlating the variable \"%s\" by applying sample reweighting" % (varName)
    Verbose(msg, True, _verbose)

    # Definitions
    weights = {}
    events  = {}
    digi    = {}
    xMin    = kwargs["xMin"]
    xMax    = kwargs["xMax"]
    nBins   = kwargs["nBins"]
    myBins  = numpy.linspace(xMin, xMax, nBins)
 
    # Use Mass as the variable to be decorrelated
    events["all"] = pandas.concat( [df_sig, df_bkg] )[varName].values
    # events["all"] = pandas.concat( [split_list(df_sig), split_list(df_bkg)] )[varName].values
    events["sig"] = df_sig[varName].values
    events["bkg"] = df_bkg[varName].values
    #events["s/b"] = [float(s)/float(b) for s,b in zip(df_sig[varName].values, df_bkg[varName].values)]
    #events["s/b"] = [float(s)/float(b) for s,b in zip(events["sig"], events["bkg"])]
    
    # test-start
    hSig     = ROOT.TH1F("sig", "", nBins, xMin, xMax)
    hBkg     = ROOT.TH1F("bkg", "", nBins, xMin, xMax)
    hWeights = ROOT.TH1F("s/b", "", nBins, xMin, xMax)

    # For-loop: 
    for s,b in zip(events["sig"], events["bkg"]):
        # print "s = %s, b = %s" % (s, b)
        hSig.Fill(s)
        hBkg.Fill(b)
    hWeights = hSig.Clone("s/b")
    hWeights.Divide(hBkg)
    #for i in range(0, nBins+1):
    #    print "s = %s, b = %s, w = %s" % (hSig.GetBinContent(i), hBkg.GetBinContent(i), hWeights.GetBinContent(i)) 
    # test-end

    # Digitize will return numbers from 1 to len(bins) depending on which bin the event belongs to. However it won't handle
    # the situation if value is over the maximum bin edge, so you'll want to clip your values (or adjust the binning) accordingly
    # In other words; get the bin indices (convert values to bin indices according to the myBins variable)
    digi["all"] = numpy.digitize( numpy.clip(events["all"], xMin, xMax-1.0), bins=myBins, right=False )
    digi["sig"] = numpy.digitize( numpy.clip(events["sig"], xMin, xMax-1.0), bins=myBins, right=False )
    digi["bkg"] = numpy.digitize( numpy.clip(events["bkg"], xMin, xMax-1.0), bins=myBins, right=False )
    # test-start
    events["s/b"] = []
    for i, value in enumerate(events["bkg"], 0):
        bin    = hWeights.FindBin( value )
        weight = hWeights.GetBinContent(bin)
        events["s/b"].append(weight)
        #print "value = %s, weight = %s" % (value, weight)
    # test-end
    #digi["s/b"] = [float(s)/float(b) for s,b in zip(digi["sig"], digi["bkg"])]
    #digi["s/b"] = numpy.digitize( numpy.clip(events["s/b"], xMin, xMax-1.0), bins=myBins, right=False )

    # These are the weights you can give to keras.fit() function in parameter "sample_weight". The idea is to reweight the braches to get 
    # both signal and bkg to have similar mass (decouple mass from learning) 
    # The "balanced" mode uses the values of y to automatically adjust weights inversely proportional to variable frequencies in the input data as n_samples / (n_classes * np.bincount(y))
    # weights["sig"] = compute_sample_weight('balanced', y=digi["sig"]) # Flatten-out signal only
    weights["sig"] = numpy.ones(len(events["sig"]))  # an array of 1's of specific sig => leaves original distribution unchanged
    # weights["bkg"] = compute_sample_weight(class_weight='balanced', y=digi["bkg"])  # Flattens out bkg only (uniform distribution)
    #weights["bkg"] = numpy.asarray(events["bkg"], dtype=numpy.float32)  # smears a bit the bkg distribution
    weights["bkg"] = numpy.asarray(events["s/b"], dtype=numpy.float32)  
    #weights["bkg"] = numpy.asarray(digi["s/b"]) # brings core of bkg closer to signal
    if 1:
        weights["all"] = numpy.concatenate((weights["sig"], weights["bkg"]), axis=0) # Flatten-out both signal and bkg variable
    else:
        Print("WARNING! This option flattens out ALL variables resulting in no distinction between signal and background! And hence a useless DNN", True)
        weights["all"] = compute_sample_weight(class_weight='balanced', y=digi["all"]) 


    align = "{:>5} {:>15} {:>10} {:>15} {:>15} {:>15}"
    title = align.format("index", varName, "bin", "weight", "weight (s)", "weight (b)")
    hLine = "="*100
    Verbose(hLine, False, _verbose)
    Verbose("{:^100}".format("%s (%d bins, xMin = %.1f, xMax = %.1f)" % (varName, nBins, xMin, xMax) ))
    Verbose(title, False, _verbose)
    Verbose(hLine, False, _verbose)
    # For-loop: All weight information
    for i,w in enumerate(weights["all"], 0):
        s = "-"
        b = "-"
        v = events["all"][i]
        d = digi["all"][i]
        if i < len(weights["all"])/2:
            s = "%.2f" % weights["sig"][i]
            b = "%.2f" % weights["bkg"][i]
        msg = align.format("%d" % i, "%.2f" % v, "%d" % d, "%.2f" % w , "%s" % s, "%s" % b)
        Verbose(msg, False, _verbose)
    Verbose(hLine, False, _verbose)

    return weights

def GetOriginalDataFrame(scaler, df, inputList):
    '''
    https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html#sklearn.preprocessing.StandardScaler.inverse_transform
    https://machinelearningmastery.com/normalize-standardize-time-series-data-python/
    '''
    df_original = df.copy()
    features    = df_original[inputList]
    features    = scaler.inverse_transform(features.values)
    df_original[inputList] = features
    Verbose("Before:\n%s" % (df["TrijetMass"]), True)
    Verbose("After:\n%s" % (df_original["TrijetMass"]), True)
    return df_original


def GetStandardisedDataFrame(df, inputList, scalerType="robust"):
    '''
    Standardization of a dataset is a common requirement for many machine learning estimators: they might behave badly if the individual features do not
    more or less look like standard normally distributed data (e.g. Gaussian with 0 mean and unit variance).
    https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html
    N.B: When normalizing the input one must remember to additionally save the parameters for when you'll need them in the deployed network version to scale the inputs there as well!
         This is crucial since the network trained and learned in separating "standardised" input variables!
         Excellent example her ("Why to normalise?"): https://stackoverflow.com/questions/48284427/why-should-we-normalize-data-for-deep-learning-in-keras

    with_mean = True 
    If True, center the data before scaling. 

    with_std = True
    If True, scale the data to unit variance (or equivalently, unit standard deviation).


    The standard score of a sample x is calculated as:
    z = (x - u) / s
    where u is the mean of the training samples or zero if "with_mean=False", 
    and s is the standard deviation of the training samples or one if "with_std=False".

    Snippet of code taken from:
    https://stackoverflow.com/questions/38420847/apply-standardscaler-to-parts-of-a-data-set

    PLEASE READ: 
    https://jovianlin.io/feature-scaling/
    https://scikit-learn.org/stable/auto_examples/preprocessing/plot_all_scaling.html
    https://jovianlin.io/feature-scaling/
    '''    
    scalerTypres = ["standard", "robust", "minmax"]
    if scalerType.lower() == "standard":
        df_scaled = df.copy()
        features  = df_scaled[inputList]
        # Assumes that your data is normally distributed within each feature!     https://jovianlin.io/feature-scaling/
        scaler    = StandardScaler(copy=True, with_mean=True, with_std=True)
        features  = scaler.fit_transform(features.values)
        df_scaled[inputList] = features
    elif scalerType.lower() == "robust":
        df_scaled = df.copy()
        features  = df_scaled[inputList]
        # Scale features using statistics that are robust to outliers
        scaler    = RobustScaler(with_centering=True, with_scaling=True, quantile_range=(25.0, 75.0), copy=True)
        #scaler    = RobustScaler(with_centering=True, with_scaling=True, quantile_range=(30.0, 70.0), copy=True)
        #scaler    = RobustScaler(with_centering=True, with_scaling=True, quantile_range=(10.0, 90.0), copy=True)
        #scaler    = RobustScaler(with_centering=True, with_scaling=True, quantile_range=(5.0, 95.0), copy=True)
        features  = scaler.fit_transform(features.values)
        df_scaled[inputList] = features
    elif scalerType.lower() == "minmax": 
        df_scaled = df.copy()
        features  = df_scaled[inputList]
        # It shrinks the range such that it is now between 0 and 1 (or -1 to 1 if there exist negative values).
        scaler    = MinMaxScaler(feature_range=(0, 1), copy=True)
        features  = scaler.fit_transform(features.values)
        df_scaled[inputList] = features
    else:
        msg = "Unsupported scaler type \"%s\". Please select one of the following:%s" % (scalerType, ", ".join(scalerTypes) )
        raise Exception(msg)

    Verbose("Before:\n%s" % (df["TrijetMass"]), True)
    Verbose("After:\n%s" % (df_scaled["TrijetMass"]), True)
    return scaler, df_scaled


def convertHistoToGaph(histo, verbose=False):

    # Lists for values
    x     = []
    y     = []
    xerrl = []
    xerrh = []
    yerrl = []
    yerrh = []
    nBinsX = histo.GetNbinsX()
    
    for i in range (0, nBinsX+1):
        # Get values
        xVal  = histo.GetBinCenter(i)
        xLow  = histo.GetBinWidth(i)
        xHigh = xLow
        yVal  = histo.GetBinContent(i)
        yLow  = histo.GetBinError(i)
        yHigh = yLow
 
        # Store values
        x.append(xVal)
        xerrl.append(xLow)
        xerrh.append(xHigh)
        y.append(yVal)
        yerrl.append(yLow)
        yerrh.append(yHigh)
        
    # Create the TGraph with asymmetric errors
    tgraph = ROOT.TGraphAsymmErrors(len(x),
                                    array.array("d",x),
                                    array.array("d",y),
                                    array.array("d",xerrl),
                                    array.array("d",xerrh),
                                    array.array("d",yerrl),
                                    array.array("d",yerrh))
    if verbose:
        tgraph.GetXaxis().SetLimits(-0.05, 1.0)
    tgraph.SetName(histo.GetName())

    # Construct info table (debugging)
    table  = []
    align  = "{0:>6} {1:^10} {2:>10} {3:>10} {4:>10} {5:^3} {6:<10}"
    header = align.format("#", "x-", "x", "x+", "y", "+/-", "Error")
    hLine  = "="*70
    table.append("")
    table.append(hLine)
    table.append("{0:^70}".format(histo.GetName()))
    table.append(header)
    table.append(hLine)
 
    # For-loop: All values x-y and their errors
    for i, xV in enumerate(x, 0):
        row = align.format(i+1, "%.4f" % xerrl[i], "%.4f" %  x[i], "%.4f" %  xerrh[i], "%.5f" %  y[i], "+/-", "%.5f" %  yerrh[i])
        table.append(row)
    table.append(hLine)

    if 0:
        for i, line in enumerate(table, 1):
            print line
    return tgraph                              
        
def PlotOutput(Y_train, Y_test, saveDir, saveName, isSB, saveFormats):
    
    ROOT.gStyle.SetOptStat(0)

    # Create canvas
    canvas = plot.CreateCanvas()
    canvas.cd()
    canvas.SetLogy()

    # Create histograms
    h1 = ROOT.TH1F('train', '', 50, 0.0, 1.0)
    h2 = ROOT.TH1F('test' , '', 50, 0.0, 1.0)

    # Fill histograms
    for r in Y_train:
        h1.Fill(r)

    for r in Y_test:
        h2.Fill(r)

    if 0:
        h1.Scale(1./h1.Integral())
        h2.Scale(1./h2.Integral())

    ymax = max(h1.GetMaximum(), h2.GetMaximum())

    plot.ApplyStyle(h1, ROOT.kMagenta+1)
    plot.ApplyStyle(h2, ROOT.kGreen+2)

    for h in [h1,h2]:
        h.SetMinimum(100)
        h.SetMaximum(ymax*1.1)
        h.GetXaxis().SetTitle("Output")
        h.GetYaxis().SetTitle("Entries")
        h.Draw("HIST SAME")
    
    # What is this for? Ask soti
    graph = plot.CreateGraph([0.5, 0.5], [0, ymax*1.1])
    graph.Draw("same")
    
    # Create legend
    leg=plot.CreateLegend(0.6, 0.75, 0.9, 0.85)
    if isSB:
        leg.AddEntry(h1, "signal","l")
        leg.AddEntry(h2, "background","l")
    else:
        leg.AddEntry(h1, "train","l")
        leg.AddEntry(h2, "test","l")
    leg.Draw()

    plot.SavePlot(canvas, saveDir, saveName, saveFormats)
    canvas.Close()
    return

def PlotInputs(signal, bkg, var, saveDir, saveFormats, pType="sigVbkg", standardise=False, w1=numpy.array([]), w2=numpy.array([])):

    # Define plot type
    if pType == "sigVbkg":
        h1 = "signal"
        h2 = "background"
    elif pType == "testVtrain":
        h1 = "train data"
        h2 = "test data"        
    elif pType == "trainVtest":
        h1 = "test data"
        h2 = "train data"
    else:
        Print("Unknown plot type \"%s\". Using dummie names/legends for histograms" % (pType), True)
        h1 = "h1"
        h2 = "h2"
        
    ROOT.gStyle.SetOptStat(0)

    # Create canvas
    canvas = plot.CreateCanvas()
    canvas.cd()
    
    # Convert weights numpy arrays to list
    w1 = w1.tolist()
    w2 = w2.tolist()

    # Definitions 
    info  = plot.GetHistoInfo(var)
    nBins = info["nbins"]
    xMin  = info["xmin"]
    xMax  = info["xmax"]
    if standardise:
        xMin  = -5.0
        xMax  = +5.0
        nBins = 500

    # Create histograms
    hsignal = ROOT.TH1F(h1, '', nBins, xMin, xMax)
    hbkg    = ROOT.TH1F(h2, '', nBins, xMin, xMax)

    # Fill histogram (signal)
    for i, r in enumerate(signal, 0):
        if len(w1) == 0: 
            hsignal.Fill(r)
        else:
            hsignal.Fill(r, w1[i])
            
    # Fill histogram (bkg)    
    for j, r in enumerate(bkg, 0):
        if len(w2) == 0:
            hbkg.Fill(r)
        else:
            hbkg.Fill(r, w2[j])

    if hsignal.Integral() > 0:
        hsignal.Scale(1./hsignal.Integral())

    if hbkg.Integral() > 0:
        hbkg.Scale(1./hbkg.Integral())

    ymax = max(hsignal.GetMaximum(), hbkg.GetMaximum())

    plot.ApplyStyle(hsignal, ROOT.kAzure-3)
    hsignal.SetFillColor(ROOT.kAzure-3)

    plot.ApplyStyle(hbkg, ROOT.kRed +2 )

    for h in [hsignal, hbkg]:
        h.SetMaximum(ymax*1.2)
        h.GetXaxis().SetTitle(info["xlabel"])
        h.GetYaxis().SetTitle(info["ylabel"])

    hsignal.Draw("HIST")
    hbkg.Draw("HIST SAME")

    # Create legend
    leg=plot.CreateLegend(0.6, 0.75, 0.9, 0.85)
    leg.AddEntry(hsignal, h1, "f") 
    leg.AddEntry(hbkg   , h2, "l") 
    leg.Draw()

    plot.SavePlot(canvas, saveDir, var, saveFormats, False)
    canvas.Close()
    return


def PlotAndWriteJSON(signal, bkg, saveDir, saveName, jsonWr, saveFormats, **kwargs):

    resultsDict = {}
    resultsDict["signal"] = signal
    resultsDict["background"] = bkg

    normalizeToOne = False
    if "normalizeToOne" in kwargs:
        normalizeToOne = kwargs["normalizeToOne"]

   # Create canvas
    ROOT.gStyle.SetOptStat(0)
    canvas = plot.CreateCanvas()
    canvas.cd()

    hList  = []
    gList  = []
    yMin   = 100000
    yMax   = None
    xMin   = 0.0
    xMax   = 1.0
    nBins  = 50
    xTitle = "DNN output"
    yTitle = "Entries"
    log    = True

    if "log" in kwargs:
        log = kwargs["log"]
    
    if "xTitle" in kwargs:
        xTitle = kwargs["xTitle"]

    if "yTitle" in kwargs:
        yTitle = kwargs["yTitle"]
    if "xMin" in kwargs:
        xMin = kwargs['xMin']
    if "xMax" in kwargs:
        xMax = kwargs['xMax']
    if "yMin" in kwargs:
        yMin = kwargs['yMin']
    if "yMax" in kwargs:
        yMax = kwargs['yMax']
    if "nBins" in kwargs:
        xBins = kwargs['nBins']
    
    # For-loop: All results
    for i, key in enumerate(resultsDict.keys(), 0):
        # Print("Constructing histogram %s" % (key), i==0)
        h = ROOT.TH1F(key, '', nBins, xMin, xMax)
        
        # For-loop: All Entries
        for j, x in enumerate(resultsDict[key], 0):
            h.Fill(x)
            try:
                yMin_ = min(x[0], yMin)
            except:
                pass
                
        # Save maximum
        yMax_ = max(h.GetMaximum(), yMax)

        # Customise & append to list
        plot.ApplyStyle(h, i+1)

        if normalizeToOne:
            if h.Integral()>0.0:
                h.Scale(1./h.Integral())
        hList.append(h)

    if log:
        # print "yMin = %s, yMax = %s" % (yMin, yMax)
        canvas.SetLogy()

    if yMax == None:
        yMax = yMax_

    # For-loop: All histograms
    for i, h in enumerate(hList, 0):
        # h.SetMinimum(yMin_* 0.85)        
        # h.SetMaximum(yMax_* 1.15)
        h.SetMaximum(yMin)  # no guarantees when converted to TGraph!
        h.SetMaximum(yMax)  # no guarantees when converted to TGraph!
        h.GetXaxis().SetTitle(xTitle)
        h.GetYaxis().SetTitle(yTitle)
            
        if i==0:
            h.Draw("HIST")
        else:
            h.Draw("HIST SAME")
    
    # Create legend
    leg = plot.CreateLegend(0.70, 0.76, 0.90, 0.90)
    if "legHeader" in kwargs:
        leg.SetHeader(kwargs["legHeader"])

    # For-loop: All histograms
    for i, h in enumerate(hList, 0):
        if "legEntries" in kwargs:
            leg.AddEntry(h, kwargs["legEntries"][i], "l")
        else:
            leg.AddEntry(h, h.GetName(), "l")
    leg.Draw()

    plot.SavePlot(canvas, saveDir, saveName, saveFormats)
    canvas.Close()

    # Create TGraph
    for h in hList:
        gList.append(convertHistoToGaph(h))

    # Write the Tgraph into the JSON file
    for gr in gList:
        gName = "%s_%s" % (saveName, gr.GetName())
        jsonWr.addGraph(gName, gr)
    return

def PlotAndWriteJSON_DNNscore(sigOutput, bkgOutput, cutValue, signal, bkg, saveDir, saveName, jsonWr, saveFormats, **kwargs):

    resultsDict = {}
    resultsDict["signal"] = signal
    resultsDict["background"] = bkg

    normalizeToOne = False
    if "normalizeToOne" in kwargs:
        normalizeToOne = kwargs["normalizeToOne"]

   # Create canvas
    ROOT.gStyle.SetOptStat(0)
    canvas = plot.CreateCanvas()
    canvas.cd()

    hList  = []
    gList  = []
    yMin   = 100000
    yMax   = None
    xMin   = 0.0
    xMax   = 1.0
    nBins  = 50
    xTitle = "DNN output"
    yTitle = "Entries"
    log    = True

    if "log" in kwargs:
        log = kwargs["log"]
    if "xTitle" in kwargs:
        xTitle = kwargs["xTitle"]
    if "yTitle" in kwargs:
        yTitle = kwargs["yTitle"]
    if "xMin" in kwargs:
        xMin = kwargs['xMin']
    if "xMax" in kwargs:
        xMax = kwargs['xMax']
    if "nBins" in kwargs:
        nBins = kwargs['nBins']
    if "yMax" in kwargs:
        yMax = kwargs["yMax"]
    if "yMin" in kwargs:
        yMin = kwargs["yMin"]
    
    # For-loop: 
    for i, key in enumerate(resultsDict.keys(), 0):

        # print "nBins = %d, xMin = %.2f, xMax = %.2f" % (nBins, xMin, xMax)
        h = ROOT.TH1F(key, '', nBins, xMin, xMax)
        for j, x in enumerate(resultsDict[key], 0):
            
            score = None
            if key == "signal":
                score = sigOutput[j]
            elif key == "background":
                score = bkgOutput[j]
            else:
                raise Exception("This should not be reached")
    
            # Check if DNN score satisfies cut
            if score < cutValue:
                Verbose("%d) DNN score for %s is %.3f" % (j, key, score), False)
                continue
            else:
                Verbose("%d) DNN score for %s is %.3f" % (j, key, score), False)

            # Fill the histogram
            h.Fill(x)
            try:
                yMin_ = min(x[0], yMin)
            except:
                pass
                
        # Save maximum
        yMax_ = max(h.GetMaximum(), yMax)

        # Customise & append to list
        plot.ApplyStyle(h, i+1)

        if normalizeToOne:
            if h.Integral()>0.0:
                h.Scale(1./h.Integral())
        hList.append(h)

    if yMax == None:
        yMax = yMax_

    if log:
        canvas.SetLogy()

    # For-loop: All histograms
    for i, h in enumerate(hList, 0):
        #h.SetMinimum(yMin*0.85)        
        #h.SetMaximum(yMax*1.15)
        h.SetMinimum(yMin)
        h.SetMaximum(yMax)

        h.GetXaxis().SetTitle(xTitle)
        h.GetYaxis().SetTitle(yTitle)
            
        if i==0:
            h.Draw("HIST")
        else:
            h.Draw("HIST SAME")
    
    # Create legend
    leg = plot.CreateLegend(0.6, 0.75, 0.9, 0.85)
    for h in hList:
        leg.AddEntry(h, h.GetName(),"l")
    leg.Draw()
    
    postfix = "_%s%s" % ("WP", str(cutValue).replace(".", "p"))
    plot.SavePlot(canvas, saveDir, saveName + postfix, saveFormats)
    canvas.Close()

    # Create TGraph
    for h in hList:
        gList.append(convertHistoToGaph(h))

    # Write the Tgraph into the JSON file
    for gr in gList:
        gName = "%s_%s" % (saveName, gr.GetName()) + postfix
        jsonWr.addGraph(gName, gr)
    return

def PlotTGraph(xVals, xErrs, yVals, yErrs, saveDir, saveName, jsonWr, saveFormats, **kwargs):

    # Create a TGraph object
    graph = plot.GetGraph(xVals, yVals, xErrs, xErrs, yErrs, yErrs)
 
    # Create a TCanvas object
    ROOT.gStyle.SetOptStat(0)
    canvas = plot.CreateCanvas()
    canvas.cd()
    #canvas.SetLogy()
    
    # Create the TGraph with asymmetric errors                                                                                                                                                              
    tgraph = ROOT.TGraphAsymmErrors(len(xVals),
                                    array.array("d", xVals),
                                    array.array("d", yVals),
                                    array.array("d", xErrs),
                                    array.array("d", xErrs),
                                    array.array("d", yErrs),
                                    array.array("d", yErrs))
    tgraph.SetName(saveName)
    legName = None
    if "efficiency" in saveName.lower():
        legName = "efficiency"
        if saveName.lower().endswith("sig"):
            plot.ApplyStyle(tgraph, ROOT.kBlue)
        else:
            plot.ApplyStyle(tgraph, ROOT.kRed)
    elif "significance" in saveName.lower():
        legName = "significance"
        if saveName.lower().endswith("sig"):
            plot.ApplyStyle(tgraph, ROOT.kGreen)
        else:
            plot.ApplyStyle(tgraph, ROOT.kGreen+3)
    else:
        plot.ApplyStyle(tgraph, ROOT.kOrange)
        
    # Draw the TGraph
    if "xTitle" in kwargs:
        tgraph.GetXaxis().SetTitle(kwargs["xTitle"])
    if "yTitle" in kwargs:
        tgraph.GetYaxis().SetTitle(kwargs["yTitle"])
    if "xMin" in kwargs and "xMax" in kwargs:
        tgraph.GetXaxis().SetLimits(kwargs["xMin"], kwargs["xMax"])
    else:
        tgraph.GetXaxis().SetLimits(-0.05, 1.0)
    tgraph.Draw("AC")
        
    # Create legend
    leg = plot.CreateLegend(0.60, 0.70, 0.85, 0.80)
    if legName != None:
        leg.AddEntry(tgraph, legName, "l")
        leg.Draw()

    plot.SavePlot(canvas, saveDir, saveName, saveFormats)
    canvas.Close()

    # Write the Tgraph into the JSON file
    jsonWr.addGraph(saveName, tgraph)
    return

def GetEfficiency(histo):

    # Initialize sigma variables
    nbins    = histo.GetNbinsX()
    intErrs  = ROOT.Double(0.0)
    intVals  = histo.IntegralAndError(0, nbins+1, intErrs, "")
    xVals    = []
    xErrs    = []
    yVals    = []
    yErrs    = []
    yTmp     = ROOT.Double(0.0)
    if intVals == 0.0:
        Print("WARNING! The integral of histogram \"%s\" is zero!" % (histo.GetName()), True)
    
    # For-loop: All bins
    for i in range(0, nbins+1):
        xVal = histo.GetBinCenter(i)
        if xVal < 0.0:
            continue
        xErr = histo.GetBinWidth(i)*0.5
        intBin = histo.IntegralAndError(i, nbins+1, yTmp, "")
        if intVals > 0:
            yVals.append(intBin/intVals)
        else:
            yVals.append(0.0)
        xVals.append(xVal)
        xErrs.append(xErr)
        if intVals > 0:
            yErrs.append(yTmp/intVals)
        else:
            yVals.append(0.0)

    return xVals, xErrs, yVals, yErrs


def CalcEfficiency(htest_s, htest_b):

    # Initialize sigma variables
    nbins    = htest_s.GetNbinsX()
    sigmaAll = ROOT.Double(0.0)
    sigmaSel = ROOT.Double(0.0)
    All_s    = htest_s.IntegralAndError(0, nbins+1, sigmaAll, "")
    All_b    = htest_b.IntegralAndError(0, nbins+1, sigmaAll, "")
    eff_s    = []
    eff_b    = [] 
    xvalue   = []
    error    = []
    
    # For-loop: All bins
    for i in range(0, nbins+1):
        Sel_s = htest_s.IntegralAndError(i, nbins+1, sigmaSel, "")
        Sel_b = htest_b.IntegralAndError(i, nbins+1, sigmaSel, "")

        if (All_s <= 0):
            All_s = 1
            Sel_s = 0
        if (All_b <= 0):
            All_b = 1
            Sel_b = 0

        eff_s.append(Sel_s/All_s)
        eff_b.append(Sel_b/All_b)
        error.append(0)
        xvalue.append(htest_s.GetBinCenter(i))
    
    #print "%d: %s" % (len(xvalue), xvalue)
    #print "%d: %s" % (len(eff_s), eff_s)
    return xvalue, eff_s, eff_b, error

def CalcSignificance(htest_s, htest_b):
    nbins = htest_s.GetNbinsX()
    h_signif0=ROOT.TH1F('signif0', '', nbins, 0.0, 1.)
    h_signif1=ROOT.TH1F('signif1', '', nbins, 0.0, 1.)
    
    sigmaSel_s = ROOT.Double(0.0)
    sigmaSel_b = ROOT.Double(0.0)
    
    for i in range(0, nbins+1):
        # Get selected events                                                                                                                                                                       
        sSel = htest_s.IntegralAndError(i, nbins+1, sigmaSel_s, "")
        bSel = htest_b.IntegralAndError(i, nbins+1, sigmaSel_b, "")
        # Calculate Significance
        _sign0 = 0
        if (sSel+bSel > 0):
            _sign0 = sSel/math.sqrt(sSel+bSel)

        _sign1 = 2*(math.sqrt(sSel+bSel) - math.sqrt(bSel))
        h_signif0.Fill(htest_s.GetBinCenter(i), _sign0)        
        h_signif1.Fill(htest_s.GetBinCenter(i), _sign1)        
    return h_signif0, h_signif1

def GetSignificance(htest_s, htest_b):
    nbinsX     = htest_s.GetNbinsX()
    sigmaSel_s = ROOT.Double(0.0)
    sigmaSel_b = ROOT.Double(0.0)
    xVals      = []
    xErrs      = []
    signif_def = [] # S/sqrt(S+B) - same definition as TMVA
    signif_alt = [] # 2[sqrt(S+B) -sqrt(B)]

    # For-loop: All histogram bins
    for i in range(0, nbinsX+1):
        xVal = htest_s.GetBinCenter(i)
        if xVal < 0.0:
            continue
        xErr = htest_s.GetBinWidth(i)*0.5

        # Get selected events                                                                                                                                                                       
        sSel = htest_s.IntegralAndError(i, nbinsX+1, sigmaSel_s, "")
        bSel = htest_b.IntegralAndError(i, nbinsX+1, sigmaSel_b, "")

        # Calculate Significance
        _sign0 = 0
        if (sSel+bSel > 0):
            _sign0 = sSel/math.sqrt(sSel+bSel)
        _sign1 = 2*(math.sqrt(sSel+bSel) - math.sqrt(bSel))

        # Append values
        xVals.append(xVal)
        xErrs.append(xErr)
        signif_def.append(_sign0)
        signif_alt.append(_sign1)
    return xVals, xErrs, signif_def, signif_alt

def PlotEfficiency(htest_s, htest_b, saveDir, saveName, saveFormats):
    ROOT.gStyle.SetOptStat(0)
    canvas = plot.CreateCanvas()
    canvas.cd()
    canvas.SetLeftMargin(0.145)
    canvas.SetRightMargin(0.11)

    # Calculate signal and background efficiency vs output
    xvalue, eff_s, eff_b, error = CalcEfficiency(htest_s, htest_b)
    graph_s = plot.GetGraph(xvalue, eff_s, error, error, error, error)
    graph_b = plot.GetGraph(xvalue, eff_b, error, error, error, error)
    
    plot.ApplyStyle(graph_s, ROOT.kBlue)
    plot.ApplyStyle(graph_b, ROOT.kRed)
    
    # Calculate significance vs output
    h_signif0, h_signif1 = CalcSignificance(htest_s, htest_b)
    
    plot.ApplyStyle(h_signif0, ROOT.kGreen)
    plot.ApplyStyle(h_signif1, ROOT.kGreen+3)
    
    #=== Get maximum of significance
    maxSignif0 = h_signif0.GetMaximum()
    maxSignif1 = h_signif1.GetMaximum()
    maxSignif = max(maxSignif0, maxSignif1)
    
    # Normalize significance
    h_signifScaled0 = h_signif0.Clone("signif0")
    if maxSignif > 0:
        h_signifScaled0.Scale(1./float(maxSignif))

    h_signifScaled1 = h_signif1.Clone("signif1")
    if maxSignif > 0:
        h_signifScaled1.Scale(1./float(maxSignif))

    #Significance: Get new maximum
    ymax = max(h_signifScaled0.GetMaximum(), h_signifScaled1.GetMaximum())

    for obj in [graph_s, graph_b, h_signifScaled0, h_signifScaled1]:
        obj.GetXaxis().SetTitle("Output")
        obj.GetYaxis().SetTitle("Efficiency")
        obj.SetMaximum(ymax*1.1)
        obj.SetMinimum(0)
    #Draw    
    h_signifScaled0.Draw("HIST")
    h_signifScaled1.Draw("HIST SAME")
    graph_s.Draw("PL SAME")
    graph_b.Draw("PL SAME")

    graph = plot.CreateGraph([0.5, 0.5], [0, ymax*1.1])
    graph.Draw("same")

    #Legend
    leg=plot.CreateLegend(0.50, 0.25, 0.85, 0.45)    
    leg.AddEntry(graph_s, "Signal Efficiency", "l")
    leg.AddEntry(graph_b, "Bkg Efficiency", "l")
    leg.AddEntry(h_signifScaled0, "S/#sqrt{S+B}", "l")
    leg.AddEntry(h_signifScaled1, "2#times(#sqrt{S+B} - #sqrt{B})", "l")
    leg.Draw()

    # Define Right Axis (Significance)
    signifColor = ROOT.kGreen+2
    rightAxis = ROOT.TGaxis(1, 0, 1, 1.1, 0, 1.1*maxSignif, 510, "+L")
    rightAxis.SetLineColor ( signifColor )
    rightAxis.SetLabelColor( signifColor )
    rightAxis.SetTitleColor( signifColor )
    rightAxis.SetTitleOffset(1.25)
    rightAxis.SetLabelOffset(0.005)
    rightAxis.SetLabelSize(0.04)
    rightAxis.SetTitleSize(0.045)
    rightAxis.SetTitle( "Significance" )
    rightAxis.Draw()
        
    plot.SavePlot(canvas, saveDir, saveName, saveFormats)
    canvas.Close()
    return

def GetROC(htest_s, htest_b):
    '''
    Get ROC curve (signal efficiency vs bkg efficiency)
    '''
    # Calculate signal and background efficiency vs output
    xvalue, eff_s, eff_b, error = CalcEfficiency(htest_s, htest_b)
    graph_roc = plot.GetGraph(eff_s, eff_b, error, error, error, error)
    return graph_roc

def PlotROC(graphMap, saveDir, saveName, saveFormats):
    '''
    Plot ROC curves (signal efficiency vs bkg efficiency)
    '''
    ROOT.gStyle.SetOptStat(0)
    canvas = plot.CreateCanvas()
    canvas.cd()

    leg = plot.CreateLegend(0.25, 0.25, 0.60, 0.45) #0.50, 0.25, 0.85, 0.45

    # For-loop: All graphs
    for i, k in enumerate(graphMap["graph"], 0):
        gr = graphMap["graph"][i]
        gr_name = graphMap["name"][i]
        plot.ApplyStyle(gr, i+2)
        gr.GetXaxis().SetTitle("Signal Efficiency")
        gr.GetYaxis().SetTitle("Misidentification rate")
        if i == 0:
            gr.Draw("apl")
        else:
            gr.Draw("pl same")
        leg.AddEntry(gr, gr_name, "l")

    leg.Draw("same")
    
    plot.SavePlot(canvas, saveDir, saveName, saveFormats)
    canvas.Close()
    return

def PlotOvertrainingTest(Y_train_S, Y_test_S, Y_train_B, Y_test_B, saveDir, saveName, saveFormats):
    def ApplyStyle(histo):
        if "_s" in histo.GetName():
            color = ROOT.kBlue
        else:
            color = ROOT.kRed
            
        plot.ApplyStyle(h, color)
        histo.SetMarkerSize(0.5)
        histo.GetXaxis().SetTitle("Output")
        histo.GetYaxis().SetTitle("Entries")
        histo.SetMinimum(10)
        return
        
    def GetLegendStyle(histoName):
        if "test" in histoName:
            legStyle = "f"
            if "_s" in histoName:
                legText = "signal (test)"
            else:
                legText = "background (test)"
        elif "train" in histoName:
            legStyle = "p"
            if "_s" in histoName:
                legText = "signal (train)"
            else:
                legText = "background (train)"
        return legText, legStyle

    def DrawStyle(histoName):
        if "train" in histoName:
            _style = "P"
        else:
            _style = "HIST"
        return _style

    ROOT.gStyle.SetOptStat(0)
    canvas = plot.CreateCanvas()
    canvas.cd()
    canvas.SetLogy()

    hList     = []
    DataList = [Y_train_S, Y_test_S, Y_train_B, Y_test_B]
    ymax     = 0
    nbins    = 500
    
    # Create the histograms
    htrain_s = ROOT.TH1F('train_s', '', nbins, 0.0, 1.0)
    htest_s  = ROOT.TH1F('test_s' , '', nbins, 0.0, 1.0)
    htrain_b = ROOT.TH1F('train_b', '', nbins, 0.0, 1.0)
    htest_b  = ROOT.TH1F('test_b' , '', nbins, 0.0, 1.0)

    # Append to list
    hList.append(htrain_s)
    hList.append(htest_s)
    hList.append(htrain_b)
    hList.append(htest_b)
    
    for i in range(len(DataList)):
        for r in DataList[i]:
            hList[i].Fill(r)            

    # Clone the histograms
    htrain_s1 = htrain_s.Clone("train_s")
    htrain_b1 = htrain_b.Clone("train_b")
    htest_s1  = htest_s.Clone("test_s")
    htest_b1  = htest_b.Clone("test_b")
    drawStyle = "HIST SAME"
    leg=plot.CreateLegend(0.55, 0.68, 0.85, 0.88)

    for h in hList:
        h.Rebin(10)
        # Legend
        legText, legStyle = GetLegendStyle(h.GetName())
        leg.AddEntry(h, legText, legStyle)
        ApplyStyle(h)

    ymax = max(htrain_s.GetMaximum(), htest_s.GetMaximum(), htrain_b.GetMaximum(), htest_b.GetMaximum())
    for h in hList:
        h.SetMaximum(ymax*2)        
        h.Draw(DrawStyle(h.GetName())+" SAME")

    #graph = plot.CreateGraph([0.5, 0.5], [0, ymax*2])
    #graph.Draw("same")
    #leg.Draw()
    
    # Save & close canvas
    plot.SavePlot(canvas, saveDir, saveName, saveFormats)
    canvas.Close()
    return htrain_s1, htest_s1, htrain_b1, htest_b1


def WriteModel(model, model_json, inputList, output, verbose=False):
    '''
    Write model weights and architecture in txt file
    '''
    arch = json.loads(model_json)
    with open(output, 'w') as fout:
        #Store input variable names
        fout.write('inputs ' + str(len(inputList)) + '\n')
        for var in inputList:
            fout.write(var + '\n')
                   
        # Store number of layers
        fout.write( 'layers ' + str(len(model.layers)) + '\n')
        layers = []
        
        # Iterate over each layer
        for index, l in enumerate(arch["config"]):

            # Store type of layer
            fout.write(l['class_name'] + '\n')
            #layers += [l['class_name']]

            # Convolution2D layer
            if l['class_name'] == 'Convolution2D':
                # Get weights of layer
                W = model.layers[index].get_weights()[0]
                fout.write(str(W.shape[0]) + ' ' + str(W.shape[1]) + ' ' + str(W.shape[2]) + ' ' + str(W.shape[3]) + ' ' + l['config']['border_mode'] + '\n')

                for i in range(W.shape[0]):
                    for j in range(W.shape[1]):
                        for k in range(W.shape[2]):
                            fout.write(str(W[i,j,k]) + '\n')
                fout.write(str(model.layers[index].get_weights()[1]) + '\n')

            # Activation layer
            if l['class_name'] == 'Activation':
                # Store activation function
                fout.write(l['config']['activation'] + '\n')

            # MaxPooling2D layer
            if l['class_name'] == 'MaxPooling2D':
                fout.write(str(l['config']['pool_size'][0]) + ' ' + str(l['config']['pool_size'][1]) + '\n')

            # Dense layer
            if l['class_name'] == 'Dense':
                # Store number of inputs, outputs for each layer
                W = model.layers[index].get_weights()[0]                
                fout.write(str(W.shape[0]) + ' ' + str(W.shape[1]) + '\n')
                
                for w in W:
                    # Store weights
                    fout.write(str(w) + '\n')
                # Store bias values (shifts the activation function : output[i] = (Sum(weights[i,j]*inputs[j]) + bias[i]))
                biases = model.layers[index].get_weights()[1]
                fout.write(str(biases) + '\n')

        if verbose:
            Print('Writing model in file %s' % (output), True)
        return
