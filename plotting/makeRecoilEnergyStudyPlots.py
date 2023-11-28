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
#ROOT.gROOT.ProcessLine(".L myPlotStyle.h")
#ROOT.myPlotStyle()

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

histFileLocation = "../runStudyTemplate_0000_ME1B_2023-11-17.root"
histFile = ROOT.TFile(histFileLocation)

#############################################################################################################
### Make output directory ###################################################################################
#############################################################################################################

plotSubdir = "2023-11-17_highNu_RecoilEnergyStudy"

PLOTDIR_ROOT = "/minerva/data/users/finer/highNu/analysisPlots/highNu_studies"
plotDir = "{0}/{1}".format(PLOTDIR_ROOT,plotSubdir)
if not os.path.isdir(plotDir):
  print "Making plot directory {0}".format(plotDir)
  os.system( "mkdir %s" % plotDir )

################
########## Plots
################

ehad_migration = histFile.Get("migration_vmat_ehad")
p_ehad_migration = ehad_migration.GetCVHistoWithError()

p_ehad_migration_reco_03 = p_ehad_migration.ProjectionX("migration_ehad_reco_03",2,3)
p_ehad_migration_reco_04 = p_ehad_migration.ProjectionX("migration_ehad_reco_04",3,4)
p_ehad_migration_reco_05 = p_ehad_migration.ProjectionX("migration_ehad_reco_05",4,5)
p_ehad_migration_reco_06 = p_ehad_migration.ProjectionX("migration_ehad_reco_06",5,6)
p_ehad_migration_reco_10 = p_ehad_migration.ProjectionX("migration_ehad_reco_10",9,10)
p_ehad_migration_reco_11 = p_ehad_migration.ProjectionX("migration_ehad_reco_11",10,11)
p_ehad_migration_reco_20 = p_ehad_migration.ProjectionX("migration_ehad_reco_20",19,20)
p_ehad_migration_reco_21 = p_ehad_migration.ProjectionX("migration_ehad_reco_21",20,21)

p_ehad_migration_true_03 = p_ehad_migration.ProjectionY("migration_ehad_true_03",2,3)
p_ehad_migration_true_04 = p_ehad_migration.ProjectionY("migration_ehad_true_04",3,4)
p_ehad_migration_true_05 = p_ehad_migration.ProjectionY("migration_ehad_true_05",4,5)
p_ehad_migration_true_06 = p_ehad_migration.ProjectionY("migration_ehad_true_06",5,6)
p_ehad_migration_true_10 = p_ehad_migration.ProjectionY("migration_ehad_true_10",9,10)
p_ehad_migration_true_11 = p_ehad_migration.ProjectionY("migration_ehad_true_11",10,11)
p_ehad_migration_true_20 = p_ehad_migration.ProjectionY("migration_ehad_true_20",19,20)
p_ehad_migration_true_21 = p_ehad_migration.ProjectionY("migration_ehad_true_21",20,21)

with makeEnv_TCanvas('{0}/ehad_migration.png'.format(plotDir)):

  p_ehad_migration.Draw("colz")

with makeEnv_TCanvas('{0}/ehad_migration_reco_03.png'.format(plotDir)):

  p_ehad_migration_reco_03.GetXaxis().SetTitle("Reconstructed E_{had} (MeV)")
  p_ehad_migration_reco_03.Draw()

  p_ehad_migration_reco_04.SetLineColor(ROOT.kRed)
  p_ehad_migration_reco_04.SetMarkerColor(ROOT.kRed)
  p_ehad_migration_reco_04.Draw("same")

  leg = ROOT.TLegend(0.55,0.7,0.84,0.9)
  leg.AddEntry(p_ehad_migration_reco_03,"0.2 < True E_{had} (GeV) < 0.3","lep")
  leg.AddEntry(p_ehad_migration_reco_04,"0.3 < True E_{had} (GeV) < 0.4","lep")
  leg.SetBorderSize(0)
  leg.Draw()

with makeEnv_TCanvas('{0}/ehad_migration_reco_03_frac.png'.format(plotDir)):

  p_ehad_migration_reco_03_frac = p_ehad_migration_reco_03.Clone("p_ehad_migration_reco_03_frac")
  p_ehad_migration_reco_04_frac = p_ehad_migration_reco_04.Clone("p_ehad_migration_reco_04_frac")
  p_ehad_migration_reco_03_frac.Add(p_ehad_migration_reco_03_frac,-300)
  p_ehad_migration_reco_04_frac.Add(p_ehad_migration_reco_04_frac,-300)

  p_ehad_migration_reco_03_frac.GetXaxis().SetTitle("Reconstructed E_{had} (MeV) [Frac diff]")
  p_ehad_migration_reco_03_frac.Draw()

  p_ehad_migration_reco_04_frac.SetLineColor(ROOT.kRed)
  p_ehad_migration_reco_04_frac.SetMarkerColor(ROOT.kRed)
  p_ehad_migration_reco_04_frac.Draw("same")

  leg = ROOT.TLegend(0.55,0.7,0.84,0.9)
  leg.AddEntry(p_ehad_migration_reco_03_frac,"0.2 < True E_{had} (GeV) < 0.3","lep")
  leg.AddEntry(p_ehad_migration_reco_04_frac,"0.3 < True E_{had} (GeV) < 0.4","lep")
  leg.SetBorderSize(0)
  leg.Draw()

with makeEnv_TCanvas('{0}/ehad_migration_reco_05.png'.format(plotDir)):

  p_ehad_migration_reco_05.GetXaxis().SetTitle("Reconstructed E_{had} (MeV)")
  p_ehad_migration_reco_05.Draw()

  p_ehad_migration_reco_06.SetLineColor(ROOT.kRed)
  p_ehad_migration_reco_06.SetMarkerColor(ROOT.kRed)
  p_ehad_migration_reco_06.Draw("same")

  leg = ROOT.TLegend(0.55,0.7,0.84,0.9)
  leg.AddEntry(p_ehad_migration_reco_04,"0.3 < True E_{had} (GeV) < 0.4","lep")
  leg.AddEntry(p_ehad_migration_reco_05,"0.4 < True E_{had} (GeV) < 0.5","lep")
  leg.SetBorderSize(0)
  leg.Draw()

with makeEnv_TCanvas('{0}/ehad_migration_reco_10.png'.format(plotDir)):

  p_ehad_migration_reco_10.GetXaxis().SetTitle("Reconstructed E_{had} (MeV)")
  p_ehad_migration_reco_10.Draw()

  p_ehad_migration_reco_11.SetLineColor(ROOT.kRed)
  p_ehad_migration_reco_11.SetMarkerColor(ROOT.kRed)
  p_ehad_migration_reco_11.Draw("same")

  leg = ROOT.TLegend(0.55,0.7,0.84,0.9)
  leg.AddEntry(p_ehad_migration_reco_10,"0.9 < True E_{had} (GeV) < 1.0","lep")
  leg.AddEntry(p_ehad_migration_reco_11,"1.0 < True E_{had} (GeV) < 1.1","lep")
  leg.SetBorderSize(0)
  leg.Draw()

with makeEnv_TCanvas('{0}/ehad_migration_reco_20.png'.format(plotDir)):

  p_ehad_migration_reco_20.GetXaxis().SetTitle("Reconstructed E_{had} (MeV)")
  p_ehad_migration_reco_20.Draw()

  p_ehad_migration_reco_21.SetLineColor(ROOT.kRed)
  p_ehad_migration_reco_21.SetMarkerColor(ROOT.kRed)
  p_ehad_migration_reco_21.Draw("same")

  leg = ROOT.TLegend(0.55,0.7,0.84,0.9)
  leg.AddEntry(p_ehad_migration_reco_20,"1.9 < True E_{had} (GeV) < 2.0","lep")
  leg.AddEntry(p_ehad_migration_reco_21,"2.0 < True E_{had} (GeV) < 2.1","lep")
  leg.SetBorderSize(0)
  leg.Draw()

with makeEnv_TCanvas('{0}/ehad_migration_true_03.png'.format(plotDir)):

  p_ehad_migration_true_03.GetXaxis().SetTitle("True E_{had} (MeV)")
  p_ehad_migration_true_03.Draw()

  p_ehad_migration_true_04.SetLineColor(ROOT.kRed)
  p_ehad_migration_true_04.SetMarkerColor(ROOT.kRed)
  p_ehad_migration_true_04.Draw("same")

  leg = ROOT.TLegend(0.55,0.7,0.84,0.9)
  leg.AddEntry(p_ehad_migration_true_03,"0.2 < Reco E_{had} (GeV) < 0.3","lep")
  leg.AddEntry(p_ehad_migration_true_04,"0.3 < Reco E_{had} (GeV) < 0.4","lep")
  leg.SetBorderSize(0)
  leg.Draw()

with makeEnv_TCanvas('{0}/ehad_migration_true_05.png'.format(plotDir)):

  p_ehad_migration_true_05.GetXaxis().SetTitle("True E_{had} (MeV)")
  p_ehad_migration_true_05.Draw()

  p_ehad_migration_true_06.SetLineColor(ROOT.kRed)
  p_ehad_migration_true_06.SetMarkerColor(ROOT.kRed)
  p_ehad_migration_true_06.Draw("same")

  leg = ROOT.TLegend(0.55,0.7,0.84,0.9)
  leg.AddEntry(p_ehad_migration_true_04,"0.3 < Reco E_{had} (GeV) < 0.4","lep")
  leg.AddEntry(p_ehad_migration_true_05,"0.4 < Reco E_{had} (GeV) < 0.5","lep")
  leg.SetBorderSize(0)
  leg.Draw()

with makeEnv_TCanvas('{0}/ehad_migration_true_10.png'.format(plotDir)):

  p_ehad_migration_true_10.GetXaxis().SetTitle("True E_{had} (MeV)")
  p_ehad_migration_true_10.Draw()

  p_ehad_migration_true_11.SetLineColor(ROOT.kRed)
  p_ehad_migration_true_11.SetMarkerColor(ROOT.kRed)
  p_ehad_migration_true_11.Draw("same")

  leg = ROOT.TLegend(0.55,0.7,0.84,0.9)
  leg.AddEntry(p_ehad_migration_true_10,"0.9 < Reco E_{had} (GeV) < 1.0","lep")
  leg.AddEntry(p_ehad_migration_true_11,"1.0 < Reco E_{had} (GeV) < 1.1","lep")
  leg.SetBorderSize(0)
  leg.Draw()

with makeEnv_TCanvas('{0}/ehad_migration_true_20.png'.format(plotDir)):

  p_ehad_migration_true_20.GetXaxis().SetTitle("True E_{had} (MeV)")
  p_ehad_migration_true_20.Draw()

  p_ehad_migration_true_21.SetLineColor(ROOT.kRed)
  p_ehad_migration_true_21.SetMarkerColor(ROOT.kRed)
  p_ehad_migration_true_21.Draw("same")

  leg = ROOT.TLegend(0.55,0.7,0.84,0.9)
  leg.AddEntry(p_ehad_migration_true_20,"1.9 < Reco E_{had} (GeV) < 2.0","lep")
  leg.AddEntry(p_ehad_migration_true_21,"2.0 < Reco E_{had} (GeV) < 2.1","lep")
  leg.SetBorderSize(0)
  leg.Draw()

