'''
   Apply simple kinematic selection and plot substructure variables
   for signal and background MC and compare.
'''
import ROOT, collections,sys,os
sys.path.append('./')
from optparse import OptionParser
from collections import OrderedDict

from TIMBER.Analyzer import analyzer, HistGroup
from TIMBER.Tools.Common import CompileCpp
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
config = 'bstar_config.json' # holds luminosities and cross sections

# common c++ functions that we will need when looping of the RDataFrame
CompileCpp("TIMBER/Framework/include/common.h") 
CompileCpp('bstar.cc') 

# define sample sets that we want to process, label them and define colors
# here we are only going to work with a single signal dataset
signal_names = ['signalLH2000']
names = {'signalLH2000': "b*_{LH} 2000 (GeV)"}
colors = {'signalLH2000': ROOT.kBlue}

# define some filters that we will use later: here are MET filter names and Trigger path names
# MET Flags - https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
flags = ["Flag_goodVertices",
        "Flag_globalSuperTightHalo2016Filter", 
        "Flag_HBHENoiseFilter", 
        "Flag_HBHENoiseIsoFilter",
        "Flag_EcalDeadCellTriggerPrimitiveFilter",
        "Flag_BadPFMuonFilter"
        #"Flag_ecalBadCalibReducedMINIAODFilter"  # Still work in progress flag, may not be used
    ]
# Triggers
if options.year == '16': 
    triggers = ["HLT_PFHT800","HLT_PFHT900","HLT_PFJet450"]
else: 
    triggers = ["HLT_PFHT1050","HLT_PFJet500","HLT_AK8PFJet380_TrimMass30","HLT_AK8PFJet400_TrimMass30"]

# Variables we want to plot 
# These need to be constructed as variables in the RDataFrame
# We will start with a simple one: the pT of the leading jet, and a more convoluted one: the number of loose b-jets

varnames = {
        'nbjet_loose':'loosebjets',
        'lead_jetPt':'p_{T}^{jet0}',
        'lead_softdrop_mass':'softdrop'
    }


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
        # For MC we need to apply the xsec * lumi / NumberOfGeneratedEvents weight
        # This function is a helper defined here: https://github.com/cmantill/BstarToTW_CMSDAS2021/blob/master/helpers.py#L5-L18
        norm = helpers.getNormFactor(setname,year,config,a.genEventCount)
    else: 
        norm = 1.
        
    # Book actions on the RDataFrame

    # First - we will cut on the filters we specified above
    a.Cut('filters',a.GetFlagString(flags))
    a.Cut('trigger',a.GetTriggerString(triggers))

    # Second - we need to calculate if we have two jets (with Id) that are back-to-back                                                                   
    # The following function will check for jets in opposite hemispheres (of phi) that also pass a jetId
    # it is defined here: https://github.com/cmantill/BstarToTW_CMSDAS2021/blob/master/bstar.cc#L17-L66
    # so first we *define* jetIdx as the index of these two jets back-to-back - ordered by pT
    a.Define('jetIdx','hemispherize(FatJet_phi, FatJet_jetId)') 
    
    # Third - we will perform a selection:
    # by requiring at least two fat-jets (step 1) that are back to back (step 2) and that have a minimum pT of 400 (step 3)
    # some of these functions used below such as max() and Sum() are defined in RDataFrame - see the cheatsheet: https://root.cern/doc/master/classROOT_1_1RDataFrame.html#cheatsheet
    a.Cut('nFatJets_cut','nFatJet > max(jetIdx[0],jetIdx[1])') # (step 1) If we don't do this, we may try to access variables of jets that don't exist! (leads to seg fault)
    a.Cut("hemis","(jetIdx[0] != -1)&&(jetIdx[1] != -1)") # (step 2) we cut on the variable we just defined - so that both jet indices exist and are different that the default value -1 
    a.Cut('pt_cut','FatJet_pt[jetIdx[0]] > 400 && FatJet_pt[jetIdx[1]] > 400') # (step 3) 

    # Now we are ready to define our first variable to plot: lead_jetPt
    a.Define('lead_jetPt','FatJet_pt[jetIdx[0]]')





    #ADD SOFT DROP MASS
    a.Define('lead_softdrop_mass','FatJet_msoftdrop[0]')
    a.Cut('softdrop_cut','lead_softdrop_mass > 50')



    # To define our second variable, the number of loose b-jets, let's define the b-tagging working points
    # These [loose, medium, tight] working points are for the DeepCSV variable (ranging between 0 and 1) - saved in NanoAOD as Jet_btagDeepB:
    bcut = []
    if year == '16' :
        bcut = [0.2217,0.6321,0.8953]
    elif year == '17' :
        bcut = [0.1522,0.4941,0.8001]
    elif year == '18' :
        bcut = [0.1241,0.4184,0.7571]
    # Then, we use the Sum function of RDataFrame to count the number of AK4Jets with DeepCSV score larger than the loose WP
    a.Define('nbjet_loose','Sum(Jet_btagDeepB > '+str(bcut[0])+')') # DeepCSV loose WP

    # Finally let's define the normalization weight of the sample as one variable as well
    a.Define('norm',str(norm))

    # A nice functionality of TIMBER is to print all the selections that we have done:
    a.PrintNodeTree(plotdir+'/signal_tree.dot',verbose=True)

    # Now we are ready to save histograms (in a HistGroup)
    out = HistGroup("%s_%s"%(setname,year))
    for varname in varnames.keys():
        histname = '%s_%s_%s'%(setname,year,varname)
        # Arguments for binning that you would normally pass to a TH1 (histname, histname, number of bins, min bin, max bin)
        if "nbjet" in varname :
            hist_tuple = (histname,histname, 10,0,10)
        elif "Pt" in varname :
            hist_tuple = (histname,histname,30,400,2000)
        elif "softdrop" in varname :
            hist_tuple = (histname,histname,30,0,300)
        hist = a.GetActiveNode().DataFrame.Histo1D(hist_tuple,varname,'norm') # Project dataframe into a histogram (hist name/binning tuple, variable to plot from dataframe, weight)
        hist.GetValue() # This gets the actual TH1 instead of a pointer to the TH1
        out.Add(varname,hist) # Add it to our group

    # Return the group
    return out

# Runs when calling `python example.py`
if __name__ == "__main__":
    histgroups = {} # all of the HistGroups we want to track
    for setname in signal_names:
        print ('Selecting for %s...'%setname)
    
        # Perform selection and write out histograms if using --select
        if options.select:
            histgroup = select(setname,options.year)
            outfile = ROOT.TFile.Open("rootfiles/%s_%s.root"%(setname,options.year),'RECREATE')
            outfile.cd()
            histgroup.Do('Write') # This will call TH1.Write() for all of the histograms in the group
            outfile.Close()
            del histgroup # Now that they are saved out, drop from memory
            
        # Open histogram files that we saved
        infile = ROOT.TFile.Open("rootfiles/%s_%s.root"%(setname,options.year))
        # ... raise exception if we forgot to run with --select!
        if infile == None:
            raise TypeError("rootfiles/%s_%s.root does not exist"%(setname,options.year))
        # Put histograms back into HistGroups
        histgroups[setname] = HistGroup(setname)
        for key in infile.GetListOfKeys(): # loop over histograms in the file
            keyname = key.GetName()
            if setname not in keyname: continue # skip if it's not for the current set we are on
            varname = keyname[len(setname+'_'+options.year)+1:] # get the variable name (ex. lead_tau32)
            inhist = infile.Get(key.GetName()) # get it from the file
            inhist.SetDirectory(0) # set the directory so hist is stored in memory and not as reference to TFile (this way it doesn't get tossed by python garbage collection when infile changes)
            histgroups[setname].Add(varname,inhist) # add to our group
            print('add var  name for ',setname)

    # For each variable to plot...
    for varname in varnames.keys():
        plot_filename = plotdir+'/%s_%s.png'%(varname,options.year)

        # Setup ordered dictionaries so processes plot in the order we specify
        #signal_hists = []
        #for sig in signal_names: signal_hists.append(histgroups[sig][varname])
        signal_hists = OrderedDict()
        for sig in signal_names: signal_hists[sig] = histgroups[sig][varname]

        # Plot everything together!
        # Using TIMBER's functionaly of Plotting https://github.com/lcorcodilos/TIMBER/blob/master/TIMBER/Tools/Plot.py#L293
        #EasyPlots(plot_filename, signal_hists, colors=[2], ytitle='Events', datastyle='hist')
        # This one is more pretty: https://github.com/lcorcodilos/TIMBER/blob/master/TIMBER/Tools/Plot.py#L10
        CompareShapes(plot_filename, options.year, varnames[varname], bkgs={}, signals=signal_hists,scale=True,stackBkg=False,doSoverB=False)
