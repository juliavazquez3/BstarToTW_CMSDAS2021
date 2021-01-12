from TIMBER.Tools.Common import OpenJSON
import math, ROOT, collections
from collections import OrderedDict

def getNormFactor(setname,year,configPath,genEventCount):
    # Config loading - will have cuts, xsec, and lumi
    if isinstance(configPath,str): config = OpenJSON(configPath)
    else: config = configPath
    cuts = config['CUTS'][year]
    lumi = config['lumi'+str(year)]

    # Deal with unique ttbar cases
    if setname == 'ttbar' and year == '16':
        setname = 'ttbar'
    elif setname == 'ttbar' and (year == '17' or year == '18'):
        setname = 'ttbar-allhad'

    if setname in config['XSECS'].keys(): 
        xsec = config['XSECS'][setname]
    else:
        raise KeyError('Key "%s" does not exist in config["XSECS"]'%setname)

    norm = (xsec*lumi)/genEventCount

    return norm

def CompareShapesWithSoverB(outfilename,year,prettyvarname,bkgs={},signals={},names={},colors={},scale=True,stackBkg=True):
    '''Create a plot that compares the shapes of backgrounds versus signal.
       Backgrounds will be stacked together and signals will be plot separately.
       Total background and signals are scaled to 1 if scale = True. Inputs organized 
       as dicts so that keys can match across dicts (ex. bkgs and bkgNames).

       Added sub pad with signal/sqrt(background) calculation

    Args:
        outfilename (string): Path where plot will be saved.
        prettyvarname (string): What will be assigned to as the axis title.
        bkgs ({string:TH1}, optional): . Defaults to {}.
        signals ({string:TH1], optional): [description]. Defaults to {}.
        names ({string:string}, optional): Formatted version of names for backgrounds and signals to appear in legend. Keys must match those in bkgs and signal. Defaults to {}. 
        colors ({string:int}, optional): TColor code for backgrounds and signals to appear in plot. Keys must match those in bkgs and signal. Defaults to {}.
        scale (bool, optional): Scales everything to unity if true. Defaults to True.
        stackBkg (bool, optional): Stack backgrounds together. Defaults to True
    '''
    # Prep everything as with CompareShapes() but first check we're stacking backgrounds or this won't work
    if not stackBkg:
        raise ValueError('Cannot run without backgrounds stacked or s/sqrt(b) will not be valid.')

    # Initialize
    c = ROOT.TCanvas('c','c',800,700)
    legend = ROOT.TLegend(0.6,0.72,0.87,0.88)
    legend.SetBorderSize(0)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetOptStat(0)
    tot_bkg_int = 0
    if stackBkg:
        bkgStack = ROOT.THStack('Totbkg','Total Bkg - '+prettyvarname)
        bkgStack.SetTitle(';%s;%s'%(prettyvarname,'A.U.'))
         # Add bkgs to integral
        for bkey in bkgs.keys():
            tot_bkg_int += bkgs[bkey].Integral()
    else:
        for bkey in bkgs.keys():
            bkgs[bkey].SetTitle(';%s;%s'%(prettyvarname,'A.U.'))

    if colors == None:
        colors = {'signal':ROOT.kBlue,'qcd':ROOT.kYellow,'ttbar':ROOT.kRed,'multijet':ROOT.kYellow}
        
    if scale:
        # Scale bkgs to total integral
        for bkey in bkgs.keys():
            if stackBkg: bkgs[bkey].Scale(1.0/tot_bkg_int)
            else: bkgs[bkey].Scale(1.0/bkgs[bkey].Integral())
        # Scale signals
        for skey in signals.keys():
            signals[skey].Scale(1.0/signals[skey].Integral())

    # Now add bkgs to stack, setup legend, and draw!
    colors_in_legend = []
    procs = OrderedDict() 
    procs.update(bkgs)
    procs.update(signals)
    for pname in procs.keys():
        h = procs[pname]
        # Legend names
        if pname in names.keys(): leg_name = names[pname]
        else: leg_name = pname
        # If bkg, set fill color and add to stack
        if pname in bkgs.keys():
            h.SetFillColorAlpha(colors[pname],0.2 if not stackBkg else 1)
            h.SetLineWidth(0) 
            if stackBkg: bkgStack.Add(h)
            if colors[pname] not in colors_in_legend:
                legend.AddEntry(h,leg_name,'f')
                colors_in_legend.append(colors[pname])
                
        # If signal, set line color
        else:
            h.SetLineColor(colors[pname])
            h.SetLineWidth(2)
            if colors[pname] not in colors_in_legend:
                legend.AddEntry(h,leg_name,'l')
                colors_in_legend.append(colors[pname])

    if stackBkg:
        maximum =  max(bkgStack.GetMaximum(),signals.values()[0].GetMaximum())*1.4
        bkgStack.SetMaximum(maximum)
    else:
        maximum = max(bkgs.values()[0].GetMaximum(),signals.values()[0].GetMaximum())*1.4
        for p in procs.values():
            p.SetMaximum(maximum)
    
    #
    # Build sub pads
    #
    c.cd()
    main = ROOT.TPad('c_main','c_main',0, 0.3, 1, 1)
    SoverB = ROOT.TPad('c_sub','c_sub',0, 0, 1, 0.3)

    main.SetBottomMargin(0.0)
    main.SetLeftMargin(0.1)
    main.SetRightMargin(0.05)
    main.SetTopMargin(0.1)

    SoverB.SetLeftMargin(0.1)
    SoverB.SetRightMargin(0.05)
    SoverB.SetTopMargin(0)
    SoverB.SetBottomMargin(0.3)

    main.Draw()
    SoverB.Draw()

    ################
    # Draw on main #
    ################
    main.cd()
    if len(bkgs.keys()) > 0:
        if stackBkg:
            bkgStack.Draw('hist')
            bkgStack.GetXaxis().SetTitleOffset(1.1)
            bkgStack.Draw('hist')
        else:
            for bkg in bkgs.values():
                bkg.GetXaxis().SetTitleOffset(1.1)
                bkg.Draw('same hist')
    for h in signals.values():
        h.Draw('same hist')
    legend.Draw()

    #
    # Make the S/sqrt(B) sub pad with dedicated function
    #
    s_over_b,line_pos = MakeSoverB(bkgStack,signals.values()[0])
    SoverB.cd()
    s_over_b.GetYaxis().SetTitle('S/#sqrt{B}')
    s_over_b.GetXaxis().SetTitle(prettyvarname)
    s_over_b.SetTitle('')
    s_over_b.SetLineColorAlpha(ROOT.kBlack,1)
    s_over_b.SetLineWidth(2)
    s_over_b.SetFillColorAlpha(ROOT.kWhite,0)
    s_over_b.GetYaxis().SetLabelSize(0.08)
    s_over_b.GetYaxis().SetTitleSize(0.08)
    s_over_b.GetYaxis().SetNdivisions(306)
    s_over_b.GetXaxis().SetLabelSize(0.09)
    s_over_b.GetXaxis().SetTitleSize(0.09)
    s_over_b.GetYaxis().SetTitleOffset(0.4)
    s_over_b.Draw('hist')
    if line_pos:
        line = ROOT.TLine(line_pos,s_over_b.GetMinimum(),line_pos,s_over_b.GetMaximum())
        line.SetLineColor(ROOT.kRed)
        line.SetLineStyle(10)
        line.SetLineWidth(2)
        line.Draw('same')

    # Canvas wide settings
    c.cd()
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.08)
    c.SetRightMargin(0.11)
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Preliminary simulation"
    CMS_lumi.lumi_sqrtS = "13 TeV"
    CMS_lumi.cmsTextSize = 0.6
    CMS_lumi.CMS_lumi(c, year, 11)

    c.Print(outfilename,'png')

def MakeSoverB(stack_of_bkgs,signal):
    '''Makes the SoverB distribution and returns it.
    Assumes that signal and stack_of_bkgs have same binning.

    S/sqrt(B) is defined from the cumulative distribution
    of the histogram. In other words, S = total amount of 
    signal kept by a cut and B = total amount of backgroud
    kept by a cut. So the cumulative distribution for each
    must be calculated and then the ratio of signal to square
    root of background is taken of those. 

    There is a question then of which direction to integrate
    for a distribution. For something like tau32, you want to 
    integrate "forward" from 0 up since a signal-like tau32 cut
    is defined as keeping less than the cut value. For a 
    machine learning algorithm score (like DeepCSV), one needs
    to integrate "backward" since the cut is defined as selecting
    signal-like events when keeping values greater than the cut.

    The script will automatically find which of these to do and if 
    the signal peak is not at edge of the space. If it is not at
    the edge, it will find the 
    signal peak and build the cumulative distributions backwards
    to the left of the peak and forwards to the right of the peak.

    Args:
        pad (TPad): TPad that's already built
        stack_of_bkgs (THStack): Stack of backgrounds, already normalized
            together, and as a sum normalized to 1.
        signal (TH1): One histogram for signal. Can only calculate
            s/sqrt(b) one signal at a time.

    Returns:
        None
    '''
    # Check where signal peak is relative to distribution
    total_bkgs = stack_of_bkgs.GetStack().Last()
    nbins = total_bkgs.GetNbinsX()
    peak_bin = signal.GetMaximumBin()

    if total_bkgs.GetXaxis().GetXmin() == 0:
        if peak_bin == nbins:
            forward = False
        elif peak_bin == 1:
            forward = True
        else:
            forward = True
        peak_bin = False
        print 'Not a mass distribution. Forward = %s'%forward
    # If peak is non-zero, do background cumulative scan to left of peak
    # and forward scan to right  
    else:
        forward = None
        print 'Mass-like distribution.'
        # Clone original distirbution, set new range around peak, get cumulative
        bkg_int_low  = MakeCumulative(total_bkgs,1,       peak_bin,forward=False)
        bkg_int_high = MakeCumulative(total_bkgs,peak_bin,nbins+1, forward=True)

        sig_int_low  = MakeCumulative(signal,1,       peak_bin,forward=False)
        sig_int_high = MakeCumulative(signal,peak_bin,nbins+1, forward=True)

        # Make empty versions of original histograms
        bkg_int = total_bkgs.Clone()
        bkg_int.Reset()
        sig_int = signal.Clone()
        sig_int.Reset()     

        bkg_int.Add(bkg_int_low)
        bkg_int.Add(bkg_int_high)
        sig_int.Add(sig_int_low)
        sig_int.Add(sig_int_high)

    if forward != None:
        # if forward == False:
        #     total_bkgs.GetXaxis().SetRange(0,total_bkgs.GetNbinsX())
        #     signal.GetXaxis().SetRange(0,signal.GetNbinsX())
        bkg_int = MakeCumulative(total_bkgs,1, total_bkgs.GetNbinsX()+1,forward)
        sig_int = MakeCumulative(signal,    1, signal.GetNbinsX()+1,    forward)
        
    # Clone and empty one for binning structure
    s_over_b = bkg_int.Clone()
    s_over_b.Reset()

    # Build s/sqrt(b) per-bin
    for ix in range(1,nbins+1):
        if bkg_int.GetBinContent(ix) != 0:
            val = sig_int.GetBinContent(ix)/math.sqrt(bkg_int.GetBinContent(ix))
            s_over_b.SetBinContent(ix,val)
        else:
            s_over_b.SetBinContent(ix,0)
            print ('WARNING: Background is empty for bin %s'%ix)
        
    peak_bin_edge = False
    if peak_bin != False:
        peak_bin_edge = bkg_int.GetBinLowEdge(peak_bin)

    return s_over_b, peak_bin_edge

def MakeCumulative(hist,low,high,forward=True):
    out = hist.Clone(hist.GetName()+'_cumul')
    out.Reset()
    prev = 0
    if forward: to_scan = range(low,high)
    else: to_scan = range(high-1,low-1,-1)
    for ix in to_scan:
        val = prev + hist.GetBinContent(ix)
        out.SetBinContent(ix,val)
        prev = val
    return out
