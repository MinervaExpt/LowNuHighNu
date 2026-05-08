import ROOT
from functions import *
from plottingClasses import *
from errorMaps import *

##set ROOT to batch mode
ROOT.gROOT.SetBatch()

#Load and implement Phil's plot style header file
ROOT.gROOT.ProcessLine(".L myPlotStyle.h")
ROOT.myPlotStyle()

# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
ROOT.TH1.AddDirectory(False)

plotter = PlotUtils.MnvPlotter()
#plotter.ApplyStyle(8)
# Manually override default error summary groups
plotter.error_summary_group_map.clear()
for group in error_bands:
  for error in error_bands[group]:
    plotter.error_summary_group_map[group].push_back(error)

# Set the number of columns to use in legends
plotter.SetLegendNColumns(2)

#ROOT.gStyle.SetPalette(54) # I think this needs to be done again after instantiating MnvPlotter, because I guess MnvPlotter sets the palette as well
plotter.SetROOT6Palette(87)
ROOT.gStyle.SetNumberContours(200)

############################################################################################## Preamble above
#############################################################################################################
#############################################################################################################

## Written by ChatGPT 
def getPOT(dataSwitch,filePath):
  meta_tree = WrapChain("Meta")
  meta_tree.Add(filePath)
  meta_tree.bind("POT_Used",'d')
  meta_tree.bind("POT_Total",'d')
  
  pot_used = 0.0
  pot_total = 0.0
 
  print("DEBUG: meta_tree.GetEntries() = {0}".format(meta_tree.GetEntries())) 
  for entry in range(meta_tree.GetEntries()):
    meta_tree.GetEntry(entry)
    pot_used += meta_tree["POT_Used"]
    pot_total += meta_tree["POT_Total"]

  if dataSwitch == 'data':
    return pot_used
  else:
    return [pot_used,pot_total] 


################################################################ Define functions above

SIGNAL_DEFINITIONS = [
  'inclusive',
  'lowNu',
  'highNu'
]

FLUX_COMPONENTS = [
  'dataRate',
  'effNumerator',
  'effDenominator'
]

thesis_file_path_analysisHists = "/exp/minerva/data/users/finer/highNu/analysisHists/thesisProcessing_noFluxConstraint_yesMuonCVShift_2020-06-13.root"
thesis_file_analysisHists = ROOT.TFile(thesis_file_path_analysisHists)

## Digging further upstream in the thesis processing to gain access to a variant of the
## effDenom which adds kinematic cuts, and likely explains for the present tension b/w
## modern efficiency and thesis-era efficiency
thesis_file_path_mergedHists = "/exp/minerva/data/users/finer/highNu/mergedHists/thesisProcessing_noFluxConstraint_yesMuonCVShift_2020-06-13/highNuHists_MC_MuonKludged_minervame1D.root"
thesis_file_mergedHists = ROOT.TFile(thesis_file_path_mergedHists)
thesis_file_path_mergedHists_data = "/exp/minerva/data/users/finer/highNu/mergedHists/thesisProcessing_noFluxConstraint_yesMuonCVShift_2020-06-13/highNuHists_Data_MuonKludged_AnaTuple_minervame1D.root "
thesis_file_mergedHists_data = ROOT.TFile(thesis_file_path_mergedHists_data)
mcPOT_used,mcPOT_total = getPOT('mc',thesis_file_path_mergedHists)
thesis_pot_me1D_mc = mcPOT_used
thesis_pot_me1D_data = getPOT('data',thesis_file_path_mergedHists_data)
scaleFactor = thesis_pot_me1D_data/thesis_pot_me1D_mc
print("scaleFactor: {0}".format(scaleFactor))
##

#HISTFILE_NAME = "processedHists_2025-03-11"
#OUTFILE_NAME = "study_thesisComparison_2025-04-29"
#HISTFILE_NAME = "processedHists_2025-12-10"
#OUTFILE_NAME = "study_thesisComparison_2025-12-16"
HISTFILE_NAME = "processedHists_2026-02-19"
OUTFILE_NAME = "study_thesisComparison_2026-02-19"

#HISTDIR_NAME = "/exp/minerva/data/users/finer/MATAna/2025-05_development"
#HISTDIR_NAME = "/exp/minerva/data/users/finer/MATAna/2025-12_development"
HISTDIR_NAME = "/exp/minerva/data/users/finer/MATAna/2026-02_development"
histFileLocation = "{0}/{1}.root".format(HISTDIR_NAME,HISTFILE_NAME)

histFile = ROOT.TFile(histFileLocation)

#PLOTDIR_ROOT = "/exp/minerva/data/users/finer/MATAna/2025-05_development/plots"
#PLOTDIR_ROOT = "/exp/minerva/data/users/finer/MATAna/2025-12_development/plots"
PLOTDIR_ROOT = "/exp/minerva/data/users/finer/MATAna/2026-02_development/plots"
plotDir = "{0}/{1}".format(PLOTDIR_ROOT,OUTFILE_NAME)

if not os.path.isdir(plotDir):
  print("Making plot directory {0}".format(plotDir))
  os.system( "mkdir %s" % plotDir )

horizontalAxis_lowerBound = 2
horizontalAxis_upperBound = 22 
#horizontalAxis_upperBound = 60 

##################################################################### Common definitions, paths defined above
#############################################################################################################
#############################################################################################################

for playlist,playlist_thesis in [('ME1D','minervame1D')]:

  # Extract relevant hists from histFile
  for sigDef in SIGNAL_DEFINITIONS:

    ## eff
    exec("eff_{0} = histFile.Get('{1}/eff_{0}_{1}')".format(sigDef,playlist))
    exec("thesis_eff_{0} = thesis_file_analysisHists.Get('{1}/eff_{0}_{1}')".format(sigDef,playlist_thesis))
 
    ## eff num
    exec("effNumerator_{0} = histFile.Get('{1}/effnum_enu_ehad_{0}').ProjectionX()".format(sigDef,playlist))
    exec("thesis_effNumerator_{0} = thesis_file_analysisHists.Get('{1}/h_ENu_VS_nu_{0}_reco').ProjectionX()".format(sigDef,playlist_thesis))
  
    ## eff denom
    exec("effDenominator_{0} = histFile.Get('{1}/effdenom_enu_ehad_{0}').ProjectionX()".format(sigDef,playlist))
    exec("thesis_effDenominator_{0} = thesis_file_analysisHists.Get('{1}/h_ENu_VS_nu_{0}_truth').ProjectionX()".format(sigDef,playlist_thesis))
    exec("thesis_effDenominator_{0}_addKinematicCuts = thesis_file_mergedHists.Get('h_ENu_VS_nu_{0}_truth_addKinematicCuts').ProjectionX()".format(sigDef))
    exec("thesis_effDenominator_{0}_addKinematicCuts.Scale(scaleFactor*10**-4,\"width\")".format(sigDef)) ## This one isn't POT-normalized yet because it's being pulled from an upstage part of the analysis chain
  
  # Plot flux components, efficiency, and error summaries
  for sigDef in SIGNAL_DEFINITIONS:
    for component,yesAdd in zip(['eff','effNumerator','effDenominator'],[False,False,True]):

      # Create local, scaled copies of the flux components for each signal definition
      exec('local_{0} = {0}_{1}.GetCVHistoWithError()'.format(component,sigDef))
      exec('local_{0}.Scale(10**-4,"width")'.format(component))
      exec('local_thesis_{0} = thesis_{0}_{1}.GetCVHistoWithError()'.format(component,sigDef))
      exec('local_thesis_{0}.Scale(10**-4,"width")'.format(component))
      if(yesAdd):
        exec('local_thesis_{0}_addKinematicCuts = thesis_{0}_{1}_addKinematicCuts.GetCVHistoWithError()'.format(component,sigDef))
        exec('local_thesis_{0}_addKinematicCuts.SetMarkerColor(ROOT.kGreen+2)'.format(component,sigDef))
        exec('local_thesis_{0}_addKinematicCuts.SetLineColor(ROOT.kGreen+2)'.format(component,sigDef))
      with makeEnv_TCanvas('{0}/{1}_{2}.png'.format(plotDir,component,sigDef)):
        exec("setPlotSpecs_{0}(local_{0})".format(component))
        exec("local_{0}.GetXaxis().SetRangeUser(2,20)".format(component))
        exec("local_{0}.Draw('PE')".format(component))
      with makeEnv_TCanvas('{0}/{1}_{2}_thesis.png'.format(plotDir,component,sigDef)):
        exec("setPlotSpecs_{0}(local_thesis_{0})".format(component))
        exec("local_thesis_{0}.GetXaxis().SetRangeUser(2,20)".format(component))
        exec("local_thesis_{0}.Draw('PE')".format(component))
        if(yesAdd):
          exec("local_thesis_{0}_addKinematicCuts.Draw('PE,SAME')".format(component))
      with makeEnv_TCanvas('{0}/{1}_{2}_comparison.png'.format(plotDir,component,sigDef)):
        #exec("setPlotSpecs_{0}(local_thesis_{0})".format(component))
        #exec("local_thesis_{0}.GetXaxis().SetRangeUser(2,20)".format(component))
        exec("local_thesis_{0}.Draw('PE')".format(component))
        exec("local_{0}.SetMarkerColor(ROOT.kRed)".format(component))
        exec("local_{0}.SetLineColor(ROOT.kRed)".format(component))
        exec("local_{0}.Draw('same')".format(component))
        if(yesAdd):
          exec("local_thesis_{0}_addKinematicCuts.Draw('PE,SAME')".format(component))
      with makeEnv_TCanvas('{0}/{1}_{2}_ratio.png'.format(plotDir,component,sigDef)):
        exec("local_ratio_{0} = local_thesis_{0}.Clone(\"Thesis/Current {0}\")".format(component))
        exec("setPlotSpecs_thesisComparisonRatio(local_ratio_{0})".format(component))
        for i in range(2,16):
          exec("thesis_bin = local_thesis_{0}.GetBinContent(i)".format(component))
          exec("thesis_bin_error = local_thesis_{0}.GetBinError(i)".format(component))
          exec("current_bin = local_{0}.GetBinContent(i)".format(component))
          exec("current_bin_error = local_{0}.GetBinError(i)".format(component))
          exec("local_ratio_{0}.SetBinContent(i,thesis_bin/current_bin)".format(component))
          if not (sigDef=="highNu" and i==2): # highNu distribution has zero content/error below 3 GeV, which causes a division by zero error
            exec("local_ratio_{0}.SetBinError(i,math.sqrt(thesis_bin/current_bin*((thesis_bin_error/thesis_bin)**2+(current_bin_error/current_bin)**2)))".format(component))
        exec("local_ratio_{0}.Draw('PE')".format(component))
        exec("local_ratio_{0}.SetMarkerColor(ROOT.kViolet+2)".format(component))
        exec("local_ratio_{0}.SetLineColor(ROOT.kViolet+2)".format(component))
        if(yesAdd):
          exec("local_ratio_{0}_addKinematicCuts = local_thesis_{0}_addKinematicCuts.Clone(\"Thesis addKinematicCuts/Current {0}\")".format(component))
          for i in range(2,16):
            exec("thesis_bin_addKC = local_thesis_{0}_addKinematicCuts.GetBinContent(i)".format(component))
            exec("thesis_bin_addKC_error = local_thesis_{0}_addKinematicCuts.GetBinError(i)".format(component))
            exec("current_bin = local_{0}.GetBinContent(i)".format(component))
            exec("current_bin_error = local_{0}.GetBinError(i)".format(component)) 
            exec("local_ratio_{0}_addKinematicCuts.SetBinContent(i,thesis_bin_addKC/current_bin)".format(component))
            if not (sigDef=="highNu" and i==2): # highNu distribution has zero content/error below 3 GeV, which causes a division by zero error
              exec("local_ratio_{0}_addKinematicCuts.SetBinError(i,math.sqrt(thesis_bin_addKC/current_bin*((thesis_bin_addKC_error/thesis_bin_addKC)**2+(current_bin_error/current_bin)**2)))".format(component))
          exec("local_ratio_{0}_addKinematicCuts.SetMarkerColor(ROOT.kGreen+2)".format(component))
          exec("local_ratio_{0}_addKinematicCuts.SetLineColor(ROOT.kGreen+2)".format(component))
          exec("local_ratio_{0}_addKinematicCuts.Draw('PE,SAME')".format(component))

    #  with makeEnv_TCanvas('{0}/fluxComponents/errorSummary_{1}_{2}_{3}.png'.format(plotDir,component,sigDef,LEMEString)):
    #    exec("setPlotSpecs_{0}({0}Hist_{1}_{2})".format(component,sigDef,LEMEString))
    #    exec("localDrawErrorSummary(plotter,{0}Hist_{1}_{2})".format(component,sigDef,LEMEString))
    #  with makeEnv_TCanvas('{0}/fluxComponents/correlationMatrix_{1}_{2}_{3}.png'.format(plotDir,component,sigDef,LEMEString)):
    #    exec("tmp = localDrawCorrelationMatrix({0}Hist_{1}_{2})".format(component,sigDef,LEMEString))
    #    tmp.Draw("colz")
    #with makeEnv_TCanvas('{0}/eff/eff_{1}_{2}.png'.format(plotDir,sigDef,LEMEString)):
    #  exec("local_eff_{0}_{1} = eff_{0}_{1}.GetCVHistoWithError()".format(sigDef,LEMEString))
    #  exec("setPlotSpecs_eff(local_eff_{0}_{1})".format(sigDef,LEMEString))
    #  exec("local_eff_{0}_{1}.Draw('PE')".format(sigDef,LEMEString))
    #with makeEnv_TCanvas('{0}/eff/errorSummary_eff_{1}_{2}.png'.format(plotDir,sigDef,LEMEString)):
    #  exec("setPlotSpecs_eff(eff_{0}_{1})".format(sigDef,LEMEString))
    #  exec("localDrawErrorSummary(plotter,eff_{0}_{1})".format(sigDef,LEMEString))
    #with makeEnv_TCanvas('{0}/eff/correlationMatrix_eff_{1}_{2}.png'.format(plotDir,sigDef,LEMEString)):
    #  exec("tmp = localDrawCorrelationMatrix(eff_{0}_{1})".format(sigDef,LEMEString))
    #  tmp.Draw("colz")

 
