import ROOT
#from functions import *
from plottingClasses import *
from errorMaps import *
#from binning import *
#from errorMaps_fracDiff import *
import os,sys

##set ROOT to batch mode
ROOT.gROOT.SetBatch()

##Load and implement Phil's plot style header file
ROOT.gROOT.ProcessLine(".L $PLOTROOT/myPlotStyle.h")
ROOT.myPlotStyle()

# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
ROOT.TH1.AddDirectory(False)

plotter = ROOT.PlotUtils.MnvPlotter()
# Manually override default error summary groups
plotter.error_summary_group_map.clear()
for group in error_bands:
  for error in error_bands[group]:
    plotter.error_summary_group_map[group].push_back(error)

# I think this needs to be done again after instantiating MnvPlotter, because I guess MnvPlotter sets the palette as well
ROOT.gStyle.SetPalette(54)
#plotter.SetROOT6Palette(87)
ROOT.gStyle.SetNumberContours(200)

############################################################################################## Preamble above
#############################################################################################################
#############################################################################################################

#############################################################################################################
### Input file location #####################################################################################
#############################################################################################################

#histFileLocation = "../runStudyTemplate_0000_ME1B_2023-11-17.root"
#histFileLocation = "/minerva/app/users/finer/MATAna/LowNuHighNu/runStudyTemplate_0000_ME1B_2023-11-30_static.root"
histFileLocation = "/minerva/app/users/finer/MATAna/LowNuHighNu/runStudyTemplate_0000_ME1A_2023-12-06.root"
histFile = ROOT.TFile(histFileLocation)

#############################################################################################################
### Make output directory ###################################################################################
#############################################################################################################

#plotSubdir = "2023-11-17_highNu_RecoilEnergyStudy"
#plotSubdir = "2023-11-30_highNu_RecoilEnergyStudy"
plotSubdir = "2023-12-07_highNu_RecoilEnergyStudy"

PLOTDIR_ROOT = "/minerva/data/users/finer/highNu/analysisPlots/highNu_studies"
plotDir = "{0}/{1}".format(PLOTDIR_ROOT,plotSubdir)
if not os.path.isdir(plotDir):
  print "Making plot directory {0}".format(plotDir)
  os.system( "mkdir %s" % plotDir )

################
########## Plots
################

## Get histograms from input file
mHist_ehad_migration = histFile.Get("migration_vmat_ehad")
tHist_ehad_migration = mHist_ehad_migration.GetCVHistoWithError()
mHist_ehad_migration_lessOD = histFile.Get("migration_vmat_ehad_lessOD")
tHist_ehad_migration_lessOD = mHist_ehad_migration_lessOD.GetCVHistoWithError()

tHist_ehad_migration_diff = tHist_ehad_migration.Clone("migration_vmat_ehad_diff")
tHist_ehad_migration_diff.Add(tHist_ehad_migration_lessOD,-1.0)

for histCat in ["reco","reco_lessOD","residual","residual_zoom","residual_lessOD","residual_lessOD_zoom"]:
  for slice_label in ["03","05","10","20"]:

    ## ## Kludge because of typo in event loop script
    ## if histCat == "residual" and slice_label == "03": 
    ##   exec("mHist_ehad_{0}_{1} = histFile.Get(\"ehad_reco_{0}_{1}\")".format(histCat,slice_label))
    ## else:
    exec("mHist_ehad_{0}_{1} = histFile.Get(\"ehad_{0}_{1}\")".format(histCat,slice_label))
    exec("tHist_ehad_{0}_{1} = mHist_ehad_{0}_{1}.GetCVHistoWithError()".format(histCat,slice_label))

## Create hists for slices of truth and reco projections
for str_label,low_bound,up_bound in [["03",2,3],["04",3,4],["05",4,5],["06",5,6],["10",9,10],
                                     ["11",10,11],["20",19,20],["21",20,21]]:

  exec("tHist_ehad_migration_reco_{0} = tHist_ehad_migration.ProjectionX(\"migration_ehad_reco_{0}\",{1},{2})".format(str_label,low_bound,up_bound))
  exec("tHist_ehad_migration_true_{0} = tHist_ehad_migration.ProjectionY(\"migration_ehad_true_{0}\",{1},{2})".format(str_label,low_bound,up_bound))

## Plot migration matrices
with makeEnv_TCanvas('{0}/ehad_migration.png'.format(plotDir)):

  tHist_ehad_migration.GetXaxis().SetTitle("Reconstructed E_{had} (MeV)")
  tHist_ehad_migration.GetYaxis().SetTitle("True E_{had} (MeV)")
  tHist_ehad_migration.SetTitle("Migration using ID+OD calorimetry")
  tHist_ehad_migration.Draw("colz")

with makeEnv_TCanvas('{0}/ehad_migration_lessOD.png'.format(plotDir)):

  tHist_ehad_migration_lessOD.GetXaxis().SetTitle("Reconstructed E_{had} (MeV)")
  tHist_ehad_migration_lessOD.GetYaxis().SetTitle("True E_{had} (MeV)")
  tHist_ehad_migration_lessOD.SetTitle("Migration using ID-only calorimetry")
  tHist_ehad_migration_lessOD.Draw("colz")

ROOT.gStyle.SetPalette(57)
with makeEnv_TCanvas('{0}/ehad_migration_diff.png'.format(plotDir)):

  tHist_ehad_migration_diff.GetXaxis().SetTitle("Reconstructed E_{had} (MeV)")
  tHist_ehad_migration_diff.GetYaxis().SetTitle("True E_{had} (MeV)")
  tHist_ehad_migration_diff.SetTitle("Migration using ID+OD calorimetry minus migration using ID-only calorimetry")
  tHist_ehad_migration_diff.GetZaxis().SetRangeUser(-2500,2500)
  tHist_ehad_migration_diff.Draw("colz")

## Plot slices of truth and reco projections
for label_first_hist,label_second_hist,low_bin,center_bin,up_bin in [["03","04","0.2","0.3","0.4"],
                                                                    ["05","06","0.3","0.4","0.5"],
                                                                    ["10","11","0.9","1.0","1.1"],
                                                                    ["20","21","1.9","2.0","2.1"]]:

  for true_reco_switch in ["reco","true"]:

    x_axis_label = "Reconstructed" if true_reco_switch == "reco" else "True"
    anti_x_axis_label = "True"     if true_reco_switch == "reco" else "Reconstructed"

    with makeEnv_TCanvas('{0}/ehad_migration_{1}_{2}.png'.format(plotDir,true_reco_switch,label_first_hist)):
    
      exec("tHist_ehad_migration_{0}_{1}.GetXaxis().SetTitle(\"{2} E_{{had}} (MeV)\")".format(true_reco_switch,label_first_hist,x_axis_label))
      exec("tHist_ehad_migration_{0}_{1}.Draw()".format(true_reco_switch,label_first_hist))
    
      exec("tHist_ehad_migration_{0}_{1}.SetLineStyle(10)".format(true_reco_switch,label_second_hist))
      exec("tHist_ehad_migration_{0}_{1}.SetLineColor(ROOT.kCyan+2)".format(true_reco_switch,label_second_hist))
      exec("tHist_ehad_migration_{0}_{1}.SetMarkerColor(ROOT.kCyan+2)".format(true_reco_switch,label_second_hist))
      exec("tHist_ehad_migration_{0}_{1}.Draw(\"same\")".format(true_reco_switch,label_second_hist))
    
      leg = ROOT.TLegend(0.55,0.7,0.84,0.9)
      exec("leg.AddEntry(tHist_ehad_migration_{0}_{1},\"{2} < {3} E_{{had}} (GeV) < {4}\",\"lep\")".format(true_reco_switch,label_first_hist,low_bin,anti_x_axis_label,center_bin))
      exec("leg.AddEntry(tHist_ehad_migration_{0}_{1},\"{2} < {3} E_{{had}} (GeV) < {4}\",\"lep\")".format(true_reco_switch,label_second_hist,center_bin,anti_x_axis_label,up_bin))
      leg.SetBorderSize(0)
      leg.Draw()

for slice_label,low_bin,up_bin in [["03","0.2","0.4"],["05","0.4","0.6"],["10","0.9","1.1"],["20","1.9","2.1"]]:

  with makeEnv_TCanvas("{0}/ehad_reco_alt_{1}.png".format(plotDir,slice_label)):

    exec("tHist_ehad_reco_{0}.GetXaxis().SetTitle(\"Reconstructed E_{{ehad}} (MeV)\")".format(slice_label))
    exec("tHist_ehad_reco_{0}.Draw()".format(slice_label))

    exec("tHist_ehad_reco_lessOD_{0}.SetLineColor(ROOT.kRed)".format(slice_label))
    exec("tHist_ehad_reco_lessOD_{0}.SetMarkerColor(ROOT.kRed)".format(slice_label))
    exec("tHist_ehad_reco_lessOD_{0}.Draw(\"same\")".format(slice_label))

    leg = ROOT.TLegend(0.55,0.7,0.84,0.9)
    exec("leg.AddEntry(tHist_ehad_reco_{0},\"{1} < True E_{{had}} (GeV) < {2}\",\"lep\")".format(slice_label,low_bin,up_bin))
    exec("leg.AddEntry(tHist_ehad_reco_lessOD_{0},\"{1} < True E_{{had}} (GeV) < {2}; excluding OD activity\",\"lep\")".format(slice_label,low_bin,up_bin))
    leg.SetBorderSize(0)
    leg.Draw()

  ## Failed attempt at making log-y version of above plots
  ## with makeEnv_TCanvas("{0}/ehad_reco_alt_{1}_logy.png".format(plotDir,slice_label,True)):

  ##   #exec("tHist_ehad_reco_{0}.GetXaxis().SetTitle(\"Reconstructed E_{{ehad}} (MeV)\")".format(slice_label))
  ##   exec("tHist_ehad_reco_{0}.Draw()".format(slice_label))

  ##   #exec("tHist_ehad_reco_lessOD_{0}.SetLineColor(ROOT.kRed)".format(slice_label))
  ##   #exec("tHist_ehad_reco_lessOD_{0}.SetMarkerColor(ROOT.kRed)".format(slice_label))
  ##   exec("tHist_ehad_reco_lessOD_{0}.Draw(\"same\")".format(slice_label))

  ##   leg = ROOT.TLegend(0.55,0.7,0.84,0.9)
  ##   exec("leg.AddEntry(tHist_ehad_reco_{0},\"{1} < True E_{{had}} (GeV) < {2}\",\"lep\")".format(slice_label,low_bin,up_bin))
  ##   exec("leg.AddEntry(tHist_ehad_reco_lessOD_{0},\"{1} < True E_{{had}} (GeV) < {2}; excluding OD activity\",\"lep\")".format(slice_label,low_bin,up_bin))
  ##   leg.SetBorderSize(0)
  ##   leg.Draw()

  with makeEnv_TCanvas("{0}/ehad_residual_{1}.png".format(plotDir,slice_label)):

    exec("tHist_ehad_residual_lessOD_{0}.GetXaxis().SetRangeUser(-20,80)".format(slice_label))
    exec("tHist_ehad_residual_lessOD_{0}.GetXaxis().SetTitle(\"E_{{ehad}} #frac{{Reco - True}}{{True}}\")".format(slice_label))
    exec("tHist_ehad_residual_lessOD_{0}.SetLineColor(ROOT.kRed)".format(slice_label))
    exec("tHist_ehad_residual_lessOD_{0}.SetMarkerColor(ROOT.kRed)".format(slice_label))
    exec("tHist_ehad_residual_lessOD_{0}.Draw()".format(slice_label))

    exec("tHist_ehad_residual_{0}.Draw(\"same\")".format(slice_label))

    leg = ROOT.TLegend(0.55,0.7,0.84,0.9)
    exec("leg.AddEntry(tHist_ehad_residual_{0},\"{1} < True E_{{had}} (GeV) < {2}\",\"lep\")".format(slice_label,low_bin,up_bin))
    exec("leg.AddEntry(tHist_ehad_residual_lessOD_{0},\"{1} < True E_{{had}} (GeV) < {2}; excluding OD activity\",\"lep\")".format(slice_label,low_bin,up_bin))
    leg.SetBorderSize(0)
    leg.Draw()

  with makeEnv_TCanvas("{0}/ehad_residual_zoom_{1}.png".format(plotDir,slice_label)):

    #exec("tHist_ehad_residual_lessOD_zoom_{0}.GetXaxis().SetRangeUser(-20,80)".format(slice_label))
    exec("tHist_ehad_residual_lessOD_zoom_{0}.GetXaxis().SetTitle(\"E_{{ehad}} #frac{{Reco - True}}{{True}}\")".format(slice_label))
    exec("tHist_ehad_residual_lessOD_zoom_{0}.SetLineColor(ROOT.kRed)".format(slice_label))
    exec("tHist_ehad_residual_lessOD_zoom_{0}.SetMarkerColor(ROOT.kRed)".format(slice_label))
    exec("tHist_ehad_residual_lessOD_zoom_{0}.Draw()".format(slice_label))

    exec("tHist_ehad_residual_zoom_{0}.Draw(\"same\")".format(slice_label))

    leg = ROOT.TLegend(0.55,0.7,0.84,0.9)
    exec("leg.AddEntry(tHist_ehad_residual_zoom_{0},\"{1} < True E_{{had}} (GeV) < {2}\",\"lep\")".format(slice_label,low_bin,up_bin))
    exec("leg.AddEntry(tHist_ehad_residual_lessOD_zoom_{0},\"{1} < True E_{{had}} (GeV) < {2}; excluding OD activity\",\"lep\")".format(slice_label,low_bin,up_bin))
    leg.SetBorderSize(0)
    leg.Draw()

