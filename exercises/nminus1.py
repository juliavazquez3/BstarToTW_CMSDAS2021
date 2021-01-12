'''
   Apply simple kinematic selection and plot variables for N-1 selections
   for signal and background MC and compare.
'''
import ROOT, collections,sys,os, time
sys.path.append('./')
from optparse import OptionParser
from collections import OrderedDict

from TIMBER.Analyzer import analyzer, HistGroup, VarGroup, CutGroup
from TIMBER.Tools.Common import CompileCpp, openJSON
from TIMBER.Tools.Plot import *
import helpers

ROOT.gROOT.SetBatch(True) 

# CL options
parser = OptionParser()
parser.add_option('-y', '--year', metavar='YEAR', type='string', action='store',
                default   =   '',
                dest      =   'year',
                help      =   'Year (16,17,18)')
parser.add_option('--select', metavar='BOOL', action='store_true',
                default   =   False,
                dest      =   'select',
                help      =   'Whether to run the selection. If False, will attempt to recycle previous run histograms.')
(options, args) = parser.parse_args()

###########################################
# Establish some global variables for use #
###########################################
plotdir = 'plots/' # this is where we'll save your plots                                                                                                                                                 
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

rootfile_path = 'root://cmsxrootd.fnal.gov///store/user/cmsdas/2021/long_exercises/BstarTW/rootfiles'
config = openJSON('bstar_config.json')
cuts = config['CUTS'][options.year]

CompileCpp("TIMBER/Framework/include/common.h")
CompileCpp('bstar.cc') # has the c++ functions we need when looping of the RDataFrame

# Sets we want to process and some nice naming for our plots
signal_names = ['signalLH%s'%(mass) for mass in [2000]]#range(1400,4200,600)]
bkg_names = ['singletop_tW','singletop_tWB','ttbar','QCDHT700','QCDHT1000','QCDHT1500','QCDHT2000']
names = {
    "singletop_tW":"single top (tW)",
    "singletop_tWB":"single top (tW)",
    "ttbar":"t#bar{t}",
    "QCDHT700":"QCD",
    "QCDHT1000":"QCD",
    "QCDHT1500":"QCD",
    "QCDHT2000":"QCD",
    "QCD":"QCD",
    "singletop":"single top (tW)"
}
for sig in signal_names:
    names[sig] = "b*_{LH} %s (GeV)"%(sig[-4:])
# ... and some nice colors
colors = {}
for p in signal_names+bkg_names:
    if 'signal' in p:
        colors[p] = ROOT.kCyan-int((int(p[-4:])-1400)/600)
    elif 'ttbar' in p:
        colors[p] = ROOT.kRed
    elif 'singletop' in p:
        colors[p] = ROOT.kBlue
    elif 'QCD' in p:
        colors[p] = ROOT.kYellow
colors["QCD"] = ROOT.kYellow
colors["singletop"] = ROOT.kBlue
# ... and names
prettynames = {
    'deltaY':'|#Delta y|',
    'tau21':'#tau_{2} / #tau_{1}',
    'mW':'m_{W} [GeV]',
    'tau32':'#tau_{3} / #tau_{2}',
    'subjet_btag':'Sub-jet DeepCSV',
    'mtop':'m_{top} [GeV]'
}

# Flags - https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
flags = ["Flag_goodVertices",
        "Flag_globalSuperTightHalo2016Filter", 
        "Flag_HBHENoiseFilter", 
        "Flag_HBHENoiseIsoFilter",
        "Flag_EcalDeadCellTriggerPrimitiveFilter",
        "Flag_BadPFMuonFilter",
        "Flag_ecalBadCalibReducedMINIAODFilter"
    ]
# Triggers
if options.year == '16': 
    triggers = ["HLT_PFHT800","HLT_PFHT900","HLT_PFJet450"]
else: 
    triggers = ["HLT_PFHT1050","HLT_PFJet500","HLT_AK8PFJet380_TrimMass30","HLT_AK8PFJet400_TrimMass30"]

#########################################
# Define function for actual processing #
#########################################
def select(setname,year):
    ROOT.ROOT.EnableImplicitMT(2) # Just use two threads - no need to kill the interactive nodes

    # Initialize TIMBER analyzer
    file_path = '%s/%s_bstar%s.root' %(rootfile_path,setname, year)
    a = analyzer(file_path)

    # Determine normalization weight
    if not a.isData: 
        norm = helpers.getNormFactor(setname,year,config,a.genEventCount)
    else: 
        norm = 1.

    # Book actions on the RDataFrame
    a.Cut('filters',a.GetFlagString(flags))
    a.Cut('trigger',a.GetTriggerString(triggers))
    a.Define('jetIdx','hemispherize(FatJet_phi, FatJet_jetId)') # need to calculate if we have two jets (with Id) that are back-to-back
    a.Cut('nFatJets_cut','nFatJet > max(jetIdx[0],jetIdx[1])') # If we don't do this, we may try to access variables of jets that don't exist! (leads to seg fault)
    a.Cut("hemis","(jetIdx[0] != -1)&&(jetIdx[1] != -1)") # cut on that calculation
    a.Cut('pt_cut','FatJet_pt[jetIdx[0]] > 400 && FatJet_pt[jetIdx[1]] > 400')
    a.Cut('eta_cut','abs(FatJet_eta[jetIdx[0]]) < 2.4 && abs(FatJet_eta[jetIdx[1]]) < 2.4')
    a.Define('norm',str(norm))

    #################################
    # Build some variables for jets #
    #################################
    # Wtagging decision logic
    # Returns 0 for no tag, 1 for lead tag, 2 for sublead tag, and 3 for both tag (which is physics-wise equivalent to 2)
    wtag_str = "1*Wtag(FatJet_tau2[jetIdx[0]]/FatJet_tau1[jetIdx[0]],0,{0}, FatJet_msoftdrop[jetIdx[0]],65,105) + 2*Wtag(FatJet_tau2[jetIdx[1]]/FatJet_tau1[jetIdx[1]],0,{0}, FatJet_msoftdrop[jetIdx[1]],65,105)".format(cuts['tau21'])

    jets = VarGroup('jets')
    jets.Add('wtag_bit',    wtag_str)
    jets.Add('top_bit',     '(wtag_bit & 2)? 0: (wtag_bit & 1)? 1: -1') # (if wtag==3 or 2 (subleading w), top_index=0) else (if wtag==1, top_index=1) else (-1)
    jets.Add('top_index',   'top_bit >= 0 ? jetIdx[top_bit] : -1')
    jets.Add('w_index',     'top_index == 0 ? jetIdx[1] : top_index == 1 ? jetIdx[0] : -1')
    # Calculate some new comlumns that we'd like to cut on (that were costly to do before the other filtering)
    jets.Add("lead_vect",   "analyzer::TLvector(FatJet_pt[jetIdx[0]],FatJet_eta[jetIdx[0]],FatJet_phi[jetIdx[0]],FatJet_msoftdrop[jetIdx[0]])")
    jets.Add("sublead_vect","analyzer::TLvector(FatJet_pt[jetIdx[1]],FatJet_eta[jetIdx[1]],FatJet_phi[jetIdx[1]],FatJet_msoftdrop[jetIdx[1]])")
    jets.Add("deltaY",      "abs(lead_vect.Rapidity()-sublead_vect.Rapidity())")
    jets.Add("mtw",         "analyzer::invariantMass(lead_vect,sublead_vect)")
    
    #########
    # N - 1 #
    #########
    plotting_vars = VarGroup('plotting_vars') # assume leading is top and subleading is W
    plotting_vars.Add("mtop",        "FatJet_msoftdrop[jetIdx[0]]")
    plotting_vars.Add("mW",          "FatJet_msoftdrop[jetIdx[1]]")
    plotting_vars.Add("tau32",       "FatJet_tau3[jetIdx[0]]/FatJet_tau2[jetIdx[0]]")
    plotting_vars.Add("subjet_btag", "max(SubJet_btagDeepB[FatJet_subJetIdx1[jetIdx[0]]],SubJet_btagDeepB[FatJet_subJetIdx2[jetIdx[0]]])")
    plotting_vars.Add("tau21",       "FatJet_tau2[jetIdx[1]]/FatJet_tau1[jetIdx[1]]")

    N_cuts = CutGroup('Ncuts') # cuts
    N_cuts.Add("deltaY_cut",      "deltaY<1.6")
    N_cuts.Add("mtop_cut",        "(mtop > 105.)&&(mtop < 220.)")
    N_cuts.Add("mW_cut",          "(mW > 65.)&&(mW < 105.)")
    N_cuts.Add("tau32_cut",       "(tau32 > 0.0)&&(tau32 < %s)"%(cuts['tau32']))
    N_cuts.Add("subjet_btag_cut", "(subjet_btag > %s)&&(subjet_btag < 1.)"%(cuts['sjbtag']))
    N_cuts.Add("tau21_cut",       "(tau21 > 0.0)&&(tau21 < %s)"%(cuts['tau21']))

    # Organize N-1 of tagging variables when assuming top is always leading
    nodeToPlot = a.Apply([jets,plotting_vars])
    nminus1Nodes = a.Nminus1(nodeToPlot,N_cuts) # constructs N nodes with a different N-1 selection for each
    nminus1Hists = HistGroup('nminus1Hists')
    binning = {
        'mtop': [25,50,300],
        'mW': [25,30,270],
        'tau32': [20,0,1],
        'tau21': [20,0,1],
        'subjet_btag': [20,0,1],
        'deltaY': [20,0,2.0]
    }
    # Add hists to group and write out
    for nkey in nminus1Nodes.keys():
        if nkey == 'full': continue
        var = nkey.replace('_cut','').replace('minus_','')
        hist_tuple = (var,var,binning[var][0],binning[var][1],binning[var][2])
        hist = nminus1Nodes[nkey].DataFrame.Histo1D(hist_tuple,var,'norm')
        hist.GetValue()
        nminus1Hists.Add(var,hist)

    # Return the group
    return nminus1Hists

# Runs when calling `python ex4.py`
if __name__ == "__main__":
    start_time = time.time()
    histgroups = {} # all of the HistGroups we want to track
    varnames = []
    for setname in signal_names+bkg_names:
        print ('Selecting for %s...'%setname)
        rootfile_name = "rootfiles/%s_%s_Nminus1.root"%(setname,options.year)
    
        # Perform selection and write out histograms if using --select
        if options.select:
            histgroup = select(setname,options.year)
            outfile = ROOT.TFile.Open(rootfile_name,'RECREATE')
            outfile.cd()
            histgroup.Do('Write') # This will call TH1.Write() for all of the histograms in the group
            outfile.Close()
            del histgroup # Now that they are saved out, drop from memory
            
        # Open histogram files that we saved
        print ('Opening '+rootfile_name)
        infile = ROOT.TFile.Open(rootfile_name)
        # ... raise exception if we forgot to run with --select!
        if infile == None:
            raise TypeError(rootfile_name)
        # Put histograms back into HistGroups
        histgroups[setname] = HistGroup(setname)
        for key in infile.GetListOfKeys(): # loop over histograms in the file
            keyname = key.GetName()
            inhist = infile.Get(key.GetName()) # get it from the file
            inhist.SetDirectory(0) # set the directory so hist is stored in memory and not as reference to TFile (this way it doesn't get tossed by python garbage collection when infile changes)
            histgroups[setname].Add(keyname,inhist) # add to our group
            if keyname not in varnames:
                varnames.append(keyname)
            
    # For each variable to plot...
    for varname in varnames:
        if varname == 'deltaY': continue # deltaY optimization requires cuts on mtw to make sense so skipping
        plot_filename = plotdir+'/%s_%s_Nminus1.png'%(varname,options.year)

        # Setup ordered dictionaries so processes plot in the order we specify
        bkg_hists,signal_hists = OrderedDict(),OrderedDict()
        # Get the backgrounds
        for bkg in bkg_names: 
            histgroups[bkg][varname].SetTitle('%s N-1 20%s'%(varname, options.year)) # empty title
            if 'QCD' in bkg: # Add the QCD HT bins together for one QCD sample
                if 'QCD' not in bkg_hists.keys():
                    bkg_hists['QCD'] = histgroups[bkg][varname].Clone('QCD_'+varname)
                else:
                    bkg_hists['QCD'].Add(histgroups[bkg][varname])
            elif 'singletop' in bkg: # Add the single top channels together for one single top sample
                if 'singletop' not in bkg_hists.keys():
                    bkg_hists['singletop'] = histgroups[bkg][varname].Clone('singletop_'+varname)
                else:
                    bkg_hists['singletop'].Add(histgroups[bkg][varname])
            else: # add everything else normally (ttbar)
                bkg_hists[bkg] = histgroups[bkg][varname]#all_hists[bkg] = histgroups[bkg][varname]#
                
        # Add the signals
        for sig in signal_names: signal_hists[sig] = histgroups[sig][varname]#all_hists[sig] = histgroups[sig][varname]#

        # Plot everything together!
        print ('Plotting %s'%(plot_filename))
        helpers.CompareShapesWithSoverB(plot_filename,options.year,prettynames[varname],bkgs=bkg_hists,signals=signal_hists,colors=colors,names=names,stackBkg=True)

    print ("Total time: "+str((time.time()-start_time)/60.) + ' min')
