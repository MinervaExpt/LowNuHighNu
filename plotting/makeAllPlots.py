import ROOT
from functions import *
from plottingClasses import *
from errorMaps import *

##set ROOT to batch mode
ROOT.gROOT.SetBatch()

#Load and implement Phil's plot style header file
ROOT.gROOT.ProcessLine(".L ../style/myPlotStyle.h")
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

#HISTFILE_NAME = "Checkpoint_ZExpansionFixedSyncHistos_noFluxConstraint_noMuonCVShift_2020-05-05"
#HISTFILE_NAME = "ResponseSystematicsTest_noFluxConstraint_noMuonCVShift_2020-05-07"
#HISTFILE_NAME = "ResponseSystematicsOnlyInME_noFluxConstraint_noMuonCVShift_2020-05-11"
#HISTFILE_NAME = "CCQE3DFits_FirstAttempt_noFluxConstraint_noMuonCVShift_2020-05-14"
#HISTFILE_NAME = 'CCQE3DFits_SecondAttempt_noFluxConstraint_noMuonCVShift_2020-05-15'
#HISTFILE_NAME = 'CCQE3DFits_FourthAttempt_noFluxConstraint_noMuonCVShift_2020-05-20'
#HISTFILE_NAME = 'alternativeHighNuBinning1_noFluxConstraint_noMuonCVShift_2020-05-22'
#HISTFILE_NAME = 'alternativeHighNuBinningVeryFine_noFluxConstraint_noMuonCVShift_2020-05-22'
#HISTFILE_NAME = 'nowWithMigrationMatrices_noFluxConstraint_noMuonCVShift_2020-05-25'
#HISTFILE_NAME = 'nowWithMigrationMatrices_noFluxConstraint_yesMuonCVShift_2020-05-25'
#HISTFILE_NAME = 'newBinning_noFluxConstraint_noMuonCVShift_2020-05-27'
#HISTFILE_NAME = 'newBinning_noFluxConstraint_yesMuonCVShift_2020-05-27'
HISTFILE_NAME = 'thesisProcessing_noFluxConstraint_noMuonCVShift_2020-06-18'
#HISTFILE_NAME = 'thesisProcessing_noFluxConstraint_yesMuonCVShift_2020-06-13'

HISTDIR_NAME = "/minerva/data/users/finer/highNu/analysisHists"
histFileLocation = "{0}/{1}.root".format(HISTDIR_NAME,HISTFILE_NAME)

histFile = ROOT.TFile(histFileLocation)

PLOTDIR_ROOT = "/minerva/data/users/finer/highNu/analysisPlots"
plotDir = "{0}/{1}".format(PLOTDIR_ROOT,HISTFILE_NAME)

if not os.path.isdir(plotDir):
  print "Making plot directory {0}".format(plotDir)
  os.system( "mkdir %s" % plotDir )

horizontalAxis_lowerBound = 2
horizontalAxis_upperBound = 22 
#horizontalAxis_upperBound = 60 

##################################################################### Common definitions, paths defined above
#############################################################################################################
#############################################################################################################

#############################################################################################################
### PPFX Plots ##############################################################################################
#############################################################################################################

for fakeFold in [1]:

  ppfxFlux_LE_Nom = histFile.Get('flux/flux_PPFX_LE_originalBinning')
  ppfxFlux_LE_Nom_rebinned = histFile.Get('flux/flux_PPFX_LE_originalBinning_rebinned')
  ppfxFlux_LE_Nom_fineBins = histFile.Get('flux/flux_PPFX_LE_newOriginalBinning')
  ppfxFlux_LE_rebinned = histFile.Get('flux/flux_PPFX_LE')
  
  ppfxFlux_ME_originalBinning = histFile.Get('flux/flux_PPFX_ME_originalBinning')
  ppfxFlux_ME_lowNu_rebinned = histFile.Get('flux/flux_PPFX_ME_lowNuBinning')
  ppfxFlux_ME_highNu_rebinned = histFile.Get('flux/flux_PPFX_ME')
  
  # Create local copy of PPFX fluxes that we'll scale
  ppfxFlux_LE_Nom_scaled = ppfxFlux_LE_Nom.Clone()
  ppfxFlux_LE_Nom_scaled.Scale(10**6)
  ppfxFlux_LE_Nom_rebinned_scaled = ppfxFlux_LE_Nom_rebinned.Clone()
  ppfxFlux_LE_Nom_rebinned_scaled.Scale(10**6)
  ppfxFlux_LE_Nom_fineBins_scaled = ppfxFlux_LE_Nom_fineBins.Clone()
  ppfxFlux_LE_Nom_fineBins_scaled.Scale(10**6)
  ppfxFlux_LE_rebinned_scaled = ppfxFlux_LE_rebinned.Clone()
  ppfxFlux_LE_rebinned_scaled.Scale(10**6)
  
  ppfxFlux_ME_originalBinning_scaled = ppfxFlux_ME_originalBinning.Clone()
  ppfxFlux_ME_originalBinning_scaled.Scale(10**6)
  ppfxFlux_ME_lowNu_rebinned_scaled = ppfxFlux_ME_lowNu_rebinned.Clone("ppfxFlux_ME_lowNu_rebinned_scaled")
  ppfxFlux_ME_lowNu_rebinned_scaled.Scale(10**6)
  ppfxFlux_ME_highNu_rebinned_scaled = ppfxFlux_ME_highNu_rebinned.Clone("ppfxFlux_ME_highNu_rebinned_scaled")
  ppfxFlux_ME_highNu_rebinned_scaled.Scale(10**6)

  with makeEnv_TCanvas('{0}/flux/ppfxFlux_LE.png'.format(plotDir)):
    ppfxFlux_LE_Nom_fineBins_scaled.SetLineColor(ROOT.kGray)
    ppfxFlux_LE_Nom_fineBins_scaled.SetMarkerColor(ROOT.kGray)
    ppfxFlux_LE_Nom_fineBins_scaled.SetMarkerSize(0.5)
    ppfxFlux_LE_Nom_fineBins_scaled.Draw()
    ppfxFlux_LE_Nom_fineBins_scaled.GetXaxis().SetRangeUser(0,40)
    ppfxFlux_LE_Nom_scaled.SetMarkerSize(0.5)
    ppfxFlux_LE_Nom_scaled.Draw("same")
    #ppfxFlux_LE_Nom_scaled.GetYaxis().SetRangeUser(0,10**-5) # zoom in on tail
    ppfxFlux_LE_rebinned_scaled.SetLineColor(ROOT.kRed)
    ppfxFlux_LE_rebinned_scaled.SetMarkerColor(ROOT.kRed)
    ppfxFlux_LE_rebinned_scaled.Draw("same")
    leg = declareLegend(3,1,"UR-Flux")
    leg.AddEntry(ppfxFlux_LE_Nom_scaled,"original LE PPFX flux","p")
    leg.AddEntry(ppfxFlux_LE_Nom_fineBins_scaled,"LE PPFX flux,","p")
    leg.AddEntry(None,"regenerated","")
    leg.AddEntry(ppfxFlux_LE_rebinned_scaled,"LE PPFX flux, ","lep")
    leg.AddEntry(None,"regenerated, rebinned","")
    leg.Draw()
  
  with makeEnv_TCanvas('{0}/flux/ppfxFlux_LE_ratios.png'.format(plotDir)):
    # Divide each of the fluxes by the nominal PPFX flux
    setPlotSpecs_fluxRatio(ppfxFlux_LE_rebinned)
    ppfxFlux_LE_rebinned.Divide(ppfxFlux_LE_rebinned,ppfxFlux_LE_Nom_rebinned)
    ppfxFlux_LE_Nom_rebinned.Divide(ppfxFlux_LE_Nom_rebinned,ppfxFlux_LE_Nom_rebinned)
  
    ppfxFlux_LE_rebinned.GetYaxis().SetRangeUser(0,2)
    ppfxFlux_LE_rebinned.SetMarkerColor(ROOT.kRed)
    ppfxFlux_LE_rebinned.SetLineColor(ROOT.kRed)
    
    ppfxFlux_LE_rebinned.Draw()
    ppfxFlux_LE_Nom_rebinned.Draw("same")
    
    box = ROOT.TBox(18,0,22,2)
    box.SetFillColor(ROOT.kGray)
    box.SetFillStyle(3001)
    box.Draw()
  
    leg = declareLegend(2,1,"UL")
    leg.AddEntry(ppfxFlux_LE_Nom_rebinned,"original LE PPFX flux,","l")
    leg.AddEntry(None,"rebinned","")
    leg.AddEntry(ppfxFlux_LE_rebinned,"LE PPFX flux, ","l")
    leg.AddEntry(None,"regenerated, rebinned","")
    leg.Draw()
 
  for sigDef in ["lowNu","highNu"]: 
    with makeEnv_TCanvas('{0}/flux/ppfxFlux_ME.png'.format(plotDir)):
      ppfxFlux_ME_originalBinning_scaled.SetLineColor(ROOT.kGray)
      ppfxFlux_ME_originalBinning_scaled.SetMarkerColor(ROOT.kGray)
      ppfxFlux_ME_originalBinning_scaled.SetMarkerSize(0.5)
      ppfxFlux_ME_originalBinning_scaled.Draw()
      ppfxFlux_ME_originalBinning_scaled.GetXaxis().SetRangeUser(0,40)
      exec("ppfxFlux_ME_{0}_rebinned_scaled.SetMarkerSize(0.5)".format(sigDef))
      exec("ppfxFlux_ME_{0}_rebinned_scaled.Draw(\"same\")".format(sigDef))
      leg = declareLegend(2,1,"UR-Flux")
      leg.AddEntry(ppfxFlux_ME_originalBinning_scaled,"original ME PPFX flux","p")
      exec("leg.AddEntry(ppfxFlux_ME_{0}_rebinned_scaled,\"ME PPFX flux,\",\"p\")".format(sigDef))
      leg.AddEntry(None,"rebinned","")
      leg.Draw()

  with makeEnv_TCanvas('{0}/flux/errorSummary_ppfxFlux_ME.png'.format(plotDir)):
    setPlotSpecs_flux(ppfxFlux_ME_highNu_rebinned)
    localDrawErrorSummary(plotter,ppfxFlux_ME_highNu_rebinned)
  with makeEnv_TCanvas('{0}/flux/correlationMatrix_ppfxFlux_ME.png'.format(plotDir)):
    tmp = localDrawCorrelationMatrix(ppfxFlux_ME_highNu_rebinned)
    tmp.Draw("colz")

#############################################################################################################
### Flux Component Plots ####################################################################################
#############################################################################################################

for LEMEString in ['LE','ME']:

  # Extract all hists from histFile
  for sigDef in SIGNAL_DEFINITIONS:
    exec("eff_{0}_{1} = histFile.Get('all-{1}/eff_{0}_{1}')".format(sigDef,LEMEString))
    exec("eff_geometric_{0}_{1} = histFile.Get('all-{1}/eff_geometric_{0}_{1}')".format(sigDef,LEMEString))
    for component in FLUX_COMPONENTS:
      exec("{0}Hist_{1}_{2} = histFile.Get('all-{2}/{0}Hist2D_{1}_{2}').ProjectionX()".format(component,sigDef,LEMEString))
      exec("{0}Hist2D_{1}_{2} = histFile.Get('all-{2}/{0}Hist2D_{1}_{2}')".format(component,sigDef,LEMEString))
  
  # Plot flux components, efficiency, and error summaries
  for sigDef in SIGNAL_DEFINITIONS:
    for component in FLUX_COMPONENTS:
      # Create local, scaled copies of the flux components for each signal definition
      exec('{0}Hist_{1}_{2}_scaled = {0}Hist_{1}_{2}.GetCVHistoWithError()'.format(component,sigDef,LEMEString))
      exec('{0}Hist_{1}_{2}_scaled.Scale(10**-4,"width")'.format(component,sigDef,LEMEString))
      with makeEnv_TCanvas('{0}/fluxComponents/{1}_{2}_{3}.png'.format(plotDir,component,sigDef,LEMEString)):
        exec("setPlotSpecs_{0}({0}Hist_{1}_{2}_scaled)".format(component,sigDef,LEMEString))
        exec("{0}Hist_{1}_{2}_scaled.Draw('PE')".format(component,sigDef,LEMEString))
      with makeEnv_TCanvas('{0}/fluxComponents/errorSummary_{1}_{2}_{3}.png'.format(plotDir,component,sigDef,LEMEString)):
        exec("setPlotSpecs_{0}({0}Hist_{1}_{2})".format(component,sigDef,LEMEString))
        exec("localDrawErrorSummary(plotter,{0}Hist_{1}_{2})".format(component,sigDef,LEMEString))
      with makeEnv_TCanvas('{0}/fluxComponents/correlationMatrix_{1}_{2}_{3}.png'.format(plotDir,component,sigDef,LEMEString)):
        exec("tmp = localDrawCorrelationMatrix({0}Hist_{1}_{2})".format(component,sigDef,LEMEString))
        tmp.Draw("colz")
    with makeEnv_TCanvas('{0}/eff/eff_{1}_{2}.png'.format(plotDir,sigDef,LEMEString)):
      exec("local_eff_{0}_{1} = eff_{0}_{1}.GetCVHistoWithError()".format(sigDef,LEMEString))
      exec("setPlotSpecs_eff(local_eff_{0}_{1})".format(sigDef,LEMEString))
      exec("local_eff_{0}_{1}.Draw('PE')".format(sigDef,LEMEString))
    with makeEnv_TCanvas('{0}/eff/errorSummary_eff_{1}_{2}.png'.format(plotDir,sigDef,LEMEString)):
      exec("setPlotSpecs_eff(eff_{0}_{1})".format(sigDef,LEMEString))
      exec("localDrawErrorSummary(plotter,eff_{0}_{1})".format(sigDef,LEMEString))
    with makeEnv_TCanvas('{0}/eff/correlationMatrix_eff_{1}_{2}.png'.format(plotDir,sigDef,LEMEString)):
      exec("tmp = localDrawCorrelationMatrix(eff_{0}_{1})".format(sigDef,LEMEString))
      tmp.Draw("colz")
  
#############################################################################################################
### Plot Total and Kinematic Efficiencies ###################################################################
#############################################################################################################

for LEMEString in ['LE','ME']:
  for sigDef in SIGNAL_DEFINITIONS:

    with makeEnv_TCanvas('{0}/eff/eff_geometric_{1}_{2}.png'.format(plotDir,sigDef,LEMEString)):
      exec("local_eff_geometric_{0}_{1} = eff_geometric_{0}_{1}.GetCVHistoWithError()".format(sigDef,LEMEString))
      exec("local_eff_{0}_{1}.Draw(\"PE\")".format(sigDef,LEMEString))
      exec("local_eff_geometric_{0}_{1}.SetMarkerColor(ROOT.kRed)".format(sigDef,LEMEString))
      exec("local_eff_geometric_{0}_{1}.SetLineColor(ROOT.kRed)".format(sigDef,LEMEString))
      exec("local_eff_geometric_{0}_{1}.Draw(\"PE,same\")".format(sigDef,LEMEString))
      leg = ROOT.TLegend(0.45,0.2,0.85,0.3)
      setPlotSpecs_legend(leg)
      exec("leg.AddEntry(local_eff_geometric_{0}_{1},\"Restricted phase space \",\"l\")".format(sigDef,LEMEString))
      exec("leg.AddEntry(local_eff_{0}_{1},\"Full phase space\",\"l\")".format(sigDef,LEMEString))
      leg.Draw()

    with makeEnv_TCanvas('{0}/eff/errorSummary_eff_geometric_{1}_{2}.png'.format(plotDir,sigDef,LEMEString)):
      exec("setPlotSpecs_eff(eff_geometric_{0}_{1})".format(sigDef,LEMEString))
      exec("localDrawErrorSummary(plotter,eff_geometric_{0}_{1})".format(sigDef,LEMEString))
    with makeEnv_TCanvas('{0}/eff/correlationMatrix_eff_geometric_{1}_{2}.png'.format(plotDir,sigDef,LEMEString)):
      exec("tmp = localDrawCorrelationMatrix(eff_geometric_{0}_{1})".format(sigDef,LEMEString))
      tmp.Draw("colz")

#############################################################################################################
### Plot Data/MC Event Rate Ratios ##########################################################################
#############################################################################################################

lineAt1 = ROOT.TLine(1.5,1.,25,1.)
lineAt1.SetLineColor(ROOT.kGray+2)
lineAt1.SetLineWidth(3)
lineAt1.SetLineStyle(9)

for LEMEString in ['LE','ME']:
  for sigDef in SIGNAL_DEFINITIONS:

    # Extract hists
    exec("dataMCEventRateRatio_{0}_{1} = histFile.Get('dataMCRatios/dataMCEventRateRatio_{0}_{1}')".format(sigDef,LEMEString))
    exec("dataMCEventRateRatio_{0}_{1}_areaNorm = histFile.Get('dataMCRatios/dataMCEventRateRatio_{0}_{1}_areaNorm')".format(sigDef,LEMEString))
    
    # Create local copy of hists
    exec('dataMCEventRateRatio_{0}_{1}_scaled = dataMCEventRateRatio_{0}_{1}.Clone()'.format(sigDef,LEMEString))
    exec('dataMCEventRateRatio_{0}_{1}_areaNorm_scaled = dataMCEventRateRatio_{0}_{1}_areaNorm.Clone()'.format(sigDef,LEMEString))
    
    # Plot hists
    with makeEnv_TCanvas('{0}/dataMCRatios/dataMCEventRateRatio_{1}_{2}.png'.format(plotDir,sigDef,LEMEString)):
      exec("setPlotSpecs_dataMCRatio(dataMCEventRateRatio_{0}_{1}_scaled)".format(sigDef,LEMEString))
      exec("dataMCEventRateRatio_{0}_{1}_scaled.GetYaxis().SetRangeUser(0.7,1.3)".format(sigDef,LEMEString))
      exec("dataMCEventRateRatio_{0}_{1}_scaled.Draw('PE')".format(sigDef,LEMEString))
      lineAt1.Draw() 
 
    with makeEnv_TCanvas('{0}/dataMCRatios/dataMCEventRateRatio_{1}_{2}_areaNorm.png'.format(plotDir,sigDef,LEMEString)):
      exec("setPlotSpecs_dataMCRatio(dataMCEventRateRatio_{0}_{1}_areaNorm_scaled)".format(sigDef,LEMEString))
      exec("dataMCEventRateRatio_{0}_{1}_areaNorm_scaled.GetYaxis().SetRangeUser(0.7,1.3)".format(sigDef,LEMEString))
      exec("dataMCEventRateRatio_{0}_{1}_areaNorm_scaled.Draw('PE')".format(sigDef,LEMEString))
      lineAt1.Draw() 
  
#############################################################################################################
### Lu-style Low Nu XS,Flux, Inclusive XS Plots #############################################################
#############################################################################################################

for fakeFold in [1]:

  # Fetch Lu's cross section histogram
  xSection_lowNu_LE_Lu = histFile.Get('xSections/xSection_lowNu_LE_Lu')
  
  # Create local copy of Lu's flux that we'll scale
  xSection_lowNu_LE_Lu_scaled = xSection_lowNu_LE_Lu.Clone() 
  xSection_lowNu_LE_Lu_scaled.Scale(10**38) # Remove units for nicer plotting
  # Plot Lu's flux
  with makeEnv_TCanvas('{0}/xSections/xSection_lowNu_LE_Lu.png'.format(plotDir)):
    setPlotSpecs_xSection(xSection_lowNu_LE_Lu_scaled)
    xSection_lowNu_LE_Lu_scaled.Draw()
  with makeEnv_TCanvas('{0}/xSections/errorSummary_xSection_lowNu_LE_Lu.png'.format(plotDir)):
    localDrawErrorSummary(plotter,xSection_lowNu_LE_Lu_scaled)
  with makeEnv_TCanvas('{0}/xSections/correlationMatrix_xSection_lowNu_LE_Lu.png'.format(plotDir)):
    tmp = localDrawCorrelationMatrix(xSection_lowNu_LE_Lu_scaled)
    tmp.Draw("colz")
  
  # Fetch Lu's flux histogram
  flux_lowNu_LE_Lu = histFile.Get('flux/flux_lowNu_LE_Lu')
  
  # Create local copy of fluxes that we'll scale
  flux_lowNu_LE_Lu_scaled = flux_lowNu_LE_Lu.Clone() 
  flux_lowNu_LE_Lu_scaled.Scale(10**6,"width") # Remove units for nicer plotting
  
  # Plot lowNu Fluxes
  with makeEnv_TCanvas('{0}/flux/flux_lowNu_LE_Lu.png'.format(plotDir)):
    setPlotSpecs_flux(flux_lowNu_LE_Lu_scaled)
    flux_lowNu_LE_Lu_scaled.GetYaxis().SetRangeUser(0,100)
    flux_lowNu_LE_Lu_scaled.Draw()
  with makeEnv_TCanvas('{0}/flux/errorSummary_flux_lowNu_LE_Lu.png'.format(plotDir)):
    localDrawErrorSummary(plotter,flux_lowNu_LE_Lu_scaled)
  with makeEnv_TCanvas('{0}/flux/correlationMatrix_flux_lowNu_LE_Lu.png'.format(plotDir)):
    tmp = localDrawCorrelationMatrix(flux_lowNu_LE_Lu_scaled)
    tmp.Draw("colz")
  
  xSection_inclusive_LE_Lu = histFile.Get('xSections/xSection_inclusive_LE_Lu')
  
  # Create local copy of inclusive xSections that we'll scale
  xSection_inclusive_LE_Lu_scaled = xSection_inclusive_LE_Lu.Clone() 
  xSection_inclusive_LE_Lu_scaled.Scale(10**38) # Remove units for nicer plotting
  
  # Plot inclusive xSections
  with makeEnv_TCanvas('{0}/xSections/xSection_inclusive_LE_Lu.png'.format(plotDir)):
    setPlotSpecs_xSection(xSection_inclusive_LE_Lu_scaled)
    xSection_inclusive_LE_Lu_scaled.Draw('PE')
  with makeEnv_TCanvas('{0}/xSections/errorSummary_xSection_inclusive_LE_Lu.png'.format(plotDir)):
    localDrawErrorSummary(plotter,xSection_inclusive_LE_Lu_scaled)
  with makeEnv_TCanvas('{0}/xSections/correlationMatrix_xSection_inclusive_LE_Lu.png'.format(plotDir)):
    tmp = localDrawCorrelationMatrix(xSection_inclusive_LE_Lu_scaled)
    tmp.Draw("colz")

#############################################################################################################
### Nu cut Breakdown Plots ##################################################################################
#############################################################################################################

def plotNuCutBreakout( plotPath , plotType , plotsByNuCut , yMax , lowNuOrInclusive , legPos='UR',addLu=False , addNOMAD=False ):
  for nuCut in range(1,5):
    exec('localPlot_nuCut_{0} = plotsByNuCut[{1}].GetCVHistoWithError()'.format(nuCut,nuCut-1))
  with makeEnv_TCanvas(plotPath):
    exec('setPlotSpecs_{0}(localPlot_nuCut_1)'.format(plotType))
    localPlot_nuCut_1.Draw()
    exec('localPlot_nuCut_1.GetYaxis().SetRangeUser(0,{0})'.format(yMax))
    for nuCut,color in zip(range(2,5),[ROOT.kRed,ROOT.kBlue,ROOT.kGreen+2]):
      exec('localPlot_nuCut_{0}.SetMarkerColor({1})'.format(nuCut,color))
      exec('localPlot_nuCut_{0}.SetLineColor({1})'.format(nuCut,color))
      exec('localPlot_nuCut_{0}.Draw("same")'.format(nuCut))
    nEntriesInLegend = 4
    addLu = False ## For thesis style
    if addLu:
      # add Lu's lowNu cross section
      xSection_lowNu_LE_Lu_scaled.SetMarkerColor(ROOT.kMagenta+3)
      xSection_lowNu_LE_Lu_scaled.SetLineColor(ROOT.kMagenta+3)
      xSection_lowNu_LE_Lu_scaled.SetLineStyle(2)
      xSection_lowNu_LE_Lu_scaled.Draw("same")
      nEntriesInLegend += 1
    if addNOMAD:  
      # add NOMAD data point
      if addNOMAD == 'scaled':
        dataPoint_NOMAD = ROOT.TGraphErrors(1,array('d',[17]),array('d',[0.699*17.]),array('d',[5]),array('d',[0.025*17.]))
      else:
        dataPoint_NOMAD = ROOT.TGraphErrors(1,array('d',[17]),array('d',[0.699]),array('d',[5]),array('d',[0.025]))
      dataPoint_NOMAD.SetMarkerColor(ROOT.kWhite)
      dataPoint_NOMAD.SetLineColor(ROOT.kMagenta+2)
      dataPoint_NOMAD.Draw("same")
      nEntriesInLegend += 1
    # add legend
    leg = declareLegend(nEntriesInLegend,1,legPos)
    if addLu: leg.AddEntry(xSection_lowNu_LE_Lu_scaled,"#sigma^{Lu}_{low-#nu} ","l")
    if addNOMAD: leg.AddEntry(dataPoint_NOMAD,"#sigma^{NOMAD}_{inclusive} ","lep")
    sigDef = 'low-#nu' if lowNuOrInclusive == 'lowNu' else 'inclusive'
    for nuCut,nuVal in zip(range(1,5),[0.3,0.5,1.0,2.0]):
      exec('leg.AddEntry(localPlot_nuCut_{0},"#sigma^{{#nu<{1}}}_{{{2}}} ","lep")'.format(nuCut,nuVal,sigDef))
    leg.Draw()

for LEMEString in ['LE','ME']:

  # Fetch flux components, efficiency, xSection histograms
  for nuCut in range(1,5):
    for sigDef in SIGNAL_DEFINITIONS:
      exec("eff_{2}_{0}_nuCut_{1} = histFile.Get('all-{0}/eff_{2}_{0}_nuCut_{1}')".format(LEMEString,nuCut,sigDef))
      for component in FLUX_COMPONENTS:
        exec("{0}Hist_{1}_{2}_nuCut_{3} = histFile.Get('all-{2}/{0}Hist_{1}_{2}_nuCut_{3}')".format(component,sigDef,LEMEString,nuCut))
    for sigDef in ['xSection_lowNu_{0}_GENIE'.format(LEMEString),'xSection_lowNu_{0}_GENIE_RAW'.format(LEMEString),'xSection_lowNu_{0}'.format(LEMEString),'xSection_inclusive_{0}'.format(LEMEString),'xSectionPerE_inclusive_{0}'.format(LEMEString)]:
      exec("{0}_nuCut_{1} = histFile.Get('xSections/{0}_nuCut_{1}')".format(sigDef,nuCut))
      exec("{0}_nuCut_{1}_normalized = histFile.Get('xSections/{0}_nuCut_{1}_normalized')".format(sigDef,nuCut))
    exec("flux_lowNu_{0}_nuCut_{1} = histFile.Get('flux/flux_lowNu_{0}_nuCut_{1}')".format(LEMEString,nuCut))
    exec("flux_lowNu_normalized_{0}_nuCut_{1} = histFile.Get('flux/flux_lowNu_normalized_{0}_nuCut_{1}')".format(LEMEString,nuCut))
  
  # For debugging
  exec("ppfxFlux_{0}_tmp = histFile.Get('flux/flux_PPFX_{0}_notBinWidthNormalized')".format(LEMEString))
  with makeEnv_TCanvas('{0}/xSections/fluxUniverses_PPFX_{1}.png'.format(plotDir,LEMEString)):
    exec("tmp = ppfxFlux_{0}_tmp.GetVertErrorBand(\"Flux\")".format(LEMEString))
    tmp.GetXaxis().SetRangeUser(0,10)
    tmp.DrawAll("",True)
  with makeEnv_TCanvas('{0}/xSections/fluxUniverses_PPFX_{1}_CV.png'.format(plotDir,LEMEString)):
    exec("tmp = ppfxFlux_{0}_tmp.GetCVHistoWithError()".format(LEMEString))
    tmp.GetXaxis().SetRangeUser(0,10)
    tmp.Draw()
  
  for nuCut in range(1,5):
    # Create local, scaled copies of the flux components for various nu cuts
    for sigDef in SIGNAL_DEFINITIONS:
      for component in FLUX_COMPONENTS:
        print 'Attempting: {0}Hist_{1}_{2}_nuCut_{3}_scaled = {0}Hist_{1}_{2}_nuCut_{3}.Clone()'.format(component,sigDef,LEMEString,nuCut)
        exec('{0}Hist_{1}_{2}_nuCut_{3}_scaled = {0}Hist_{1}_{2}_nuCut_{3}.Clone()'.format(component,sigDef,LEMEString,nuCut))
        exec('{0}Hist_{1}_{2}_nuCut_{3}_scaled.Scale(10**-4,"width")'.format(component,sigDef,LEMEString,nuCut))
    # Create local, scaled copies of the xSections for various nu cuts
    #for sigDef in ['xSection_lowNu_LE_GENIE','xSection_lowNu_LE_GENIE_RAW','xSection_inclusive_LE','xSectionPerE_inclusive_LE']:
    for sigDef in ['xSection_lowNu_{0}_GENIE'.format(LEMEString),'xSection_inclusive_{0}'.format(LEMEString),'xSectionPerE_inclusive_{0}'.format(LEMEString)]:
      # Remove units for nicer plotting
      exec('{0}_nuCut_{1}_scaled = {0}_nuCut_{1}.Clone()'.format(sigDef,nuCut))
      exec('{0}_nuCut_{1}_normalized_scaled = {0}_nuCut_{1}_normalized.Clone()'.format(sigDef,nuCut))
      exec('{0}_nuCut_{1}_scaled.Scale(10**38)'.format(sigDef,nuCut)) 
      exec('{0}_nuCut_{1}_normalized_scaled.Scale(10**38)'.format(sigDef,nuCut)) 
      # Plot each xSection separately
      with makeEnv_TCanvas('{0}/xSections/{1}_nuCut_{2}.png'.format(plotDir,sigDef,nuCut)):
        exec('setPlotSpecs_xSection({0}_nuCut_{1}_scaled)'.format(sigDef,nuCut))
        exec('{0}_nuCut_{1}_scaled.Draw()'.format(sigDef,nuCut))
      with makeEnv_TCanvas('{0}/xSections/errorSummary_{1}_nuCut_{2}.png'.format(plotDir,sigDef,nuCut)):
        #exec('localDrawErrorSummary(plotter,{0}_nuCut_{1}_scaled,"TR",True,True,0.00001,False,"Flux")'.format(sigDef,nuCut))
        #exec('localDrawErrorSummary(plotter,{0}_nuCut_{1}_scaled,"TR",True,True,0.00001,False,"MuonReconstruction")'.format(sigDef,nuCut))
        #exec('localDrawErrorSummary(plotter,{0}_nuCut_{1}_scaled,"TR",False)'.format(sigDef,nuCut))
        #exec('localDrawErrorSummary(plotter,{0}_nuCut_{1}_scaled,"TR",True)'.format(sigDef,nuCut))
        exec('localDrawErrorSummary(plotter,{0}_nuCut_{1}_scaled)'.format(sigDef,nuCut))
      with makeEnv_TCanvas('{0}/xSections/correlationMatrix_{1}_nuCut_{2}.png'.format(plotDir,sigDef,nuCut)):
        exec('tmp = localDrawCorrelationMatrix({0}_nuCut_{1}_scaled)'.format(sigDef,nuCut))
        tmp.Draw("colz")
      with makeEnv_TCanvas('{0}/xSections/{1}_nuCut_{2}_normalized.png'.format(plotDir,sigDef,nuCut)):
        exec('setPlotSpecs_xSection({0}_nuCut_{1}_normalized_scaled)'.format(sigDef,nuCut))
        exec('{0}_nuCut_{1}_normalized_scaled.Draw()'.format(sigDef,nuCut))
      with makeEnv_TCanvas('{0}/xSections/errorSummary_{1}_normalized_nuCut_{2}.png'.format(plotDir,sigDef,nuCut)):
        exec('localDrawErrorSummary(plotter,{0}_nuCut_{1}_normalized_scaled)'.format(sigDef,nuCut))
      with makeEnv_TCanvas('{0}/xSections/correlationMatrix_{1}_normalized_nuCut_{2}.png'.format(plotDir,sigDef,nuCut)):
        exec('tmp = localDrawCorrelationMatrix({0}_nuCut_{1}_normalized_scaled)'.format(sigDef,nuCut))
        tmp.Draw("colz")
      #Commented out for a run that excluded flux universes. Is this necessary, anyways?
      #with makeEnv_TCanvas('{0}/xSections/fluxUniverses_{1}_nuCut_{2}.png'.format(plotDir,sigDef,nuCut)):
      #  exec('tmp = {0}_nuCut_{1}_scaled.GetVertErrorBand("Flux")'.format(sigDef,nuCut))
      #  exec('tmp.GetXaxis().SetRangeUser(0,10)'.format(sigDef,nuCut))
      #  tmp.DrawAll("",True)
  
    exec('flux_lowNu_{0}_nuCut_{1}_scaled = flux_lowNu_{0}_nuCut_{1}.Clone()'.format(LEMEString,nuCut))
    exec('flux_lowNu_{0}_nuCut_{1}_scaled.Scale(10**6,"width")'.format(LEMEString,nuCut))
    exec('flux_lowNu_normalized_{0}_nuCut_{1}_scaled = flux_lowNu_normalized_{0}_nuCut_{1}.Clone()'.format(LEMEString,nuCut))
    exec('flux_lowNu_normalized_{0}_nuCut_{1}_scaled.Scale(10**6,"width")'.format(LEMEString,nuCut))
  
  # Plot efficiency for each nu cut separately
  for nuCut in range(1,5):
    # Plot efficiency
    with makeEnv_TCanvas('{0}/eff/eff_lowNu_{1}_nuCut_{2}.png'.format(plotDir,LEMEString,nuCut)):
      exec("setPlotSpecs_eff(eff_lowNu_{0}_nuCut_{1})".format(LEMEString,nuCut))
      #exec("eff_lowNu_{0}_nuCut_{1}.GetXaxis().SetRangeUser(0,30)".format(LEMEString,nuCut))
      exec("eff_lowNu_{0}_nuCut_{1}.Draw('PE')".format(LEMEString,nuCut))
    with makeEnv_TCanvas('{0}/eff/errorSummary_eff_lowNu_{1}_nuCut_{2}.png'.format(plotDir,LEMEString,nuCut)):
      exec("localDrawErrorSummary(plotter,eff_lowNu_{0}_nuCut_{1})".format(LEMEString,nuCut))
    #with makeEnv_TCanvas('{0}/eff/correlationMatrix_eff_lowNu_{1}_nuCut_{2}.png'.format(plotDir,LEMEString,nuCut)):
    #  exec("tmp = localDrawCorrelationMatrix(eff_lowNu_{0}_nuCut_{1})".format(LEMEString,nuCut))
    #  tmp.Draw("colz")
  
  # Plot un-normalized lowNu and inclusive xSections by nu cut
  plotNuCutBreakout('{0}/xSections/xSection_lowNu_{1}_GENIE_all.png'.format(plotDir,LEMEString),'xSection',[eval('xSection_lowNu_{0}_GENIE_nuCut_{1}_scaled'.format(LEMEString,nuCut)) for nuCut in range(1,5)],2,'lowNu',legPos='UL',addLu=True)
  plotNuCutBreakout('{0}/xSections/xSection_inclusive_{1}_all.png'.format(plotDir,LEMEString),'xSection',[eval('xSection_inclusive_{0}_nuCut_{1}_scaled'.format(LEMEString,nuCut)) for nuCut in range(1,5)],18,'inclusive',legPos='UL',addNOMAD='scaled')
  plotNuCutBreakout('{0}/xSections/xSectionPerE_inclusive_{1}_all.png'.format(plotDir,LEMEString),'xSectionPerE',[eval('xSectionPerE_inclusive_{0}_nuCut_{1}_scaled'.format(LEMEString,nuCut)) for nuCut in range(1,5)],1.6,'inclusive',legPos='UL',addNOMAD=True)
  
  # Plot normalized lowNu and inclusive xSections by nu cut
  plotNuCutBreakout('{0}/xSections/xSection_lowNu_{1}_GENIE_all-normalized.png'.format(plotDir,LEMEString),'xSection',[eval('xSection_lowNu_{0}_GENIE_nuCut_{1}_normalized_scaled'.format(LEMEString,nuCut)) for nuCut in range(1,5)],2,'lowNu',legPos='UL',addLu=True)
  plotNuCutBreakout('{0}/xSections/xSection_inclusive_{1}_all-normalized.png'.format(plotDir,LEMEString),'xSection',[eval('xSection_inclusive_{0}_nuCut_{1}_normalized_scaled'.format(LEMEString,nuCut)) for nuCut in range(1,5)],18,'inclusive',legPos='UL',addNOMAD='scaled')
  plotNuCutBreakout('{0}/xSections/xSectionPerE_inclusive_{1}_all-normalized.png'.format(plotDir,LEMEString),'xSectionPerE',[eval('xSectionPerE_inclusive_{0}_nuCut_{1}_normalized_scaled'.format(LEMEString,nuCut)) for nuCut in range(1,5)],1.6,'inclusive',legPos='UL',addNOMAD=True)
  
  # for debuggging
  for nuCut in range(1,5):
    with makeEnv_TCanvas('{0}/fluxComponents/errorSummary_effNumerator_lowNu_{1}_nuCut_{2}.png'.format(plotDir,LEMEString,nuCut)):
      exec("localDrawErrorSummary(plotter,effNumeratorHist_lowNu_{0}_nuCut_{1}_scaled)".format(LEMEString,nuCut))
    with makeEnv_TCanvas('{0}/fluxComponents/correlationMatrix_effNumerator_lowNu_{1}_nuCut_{2}.png'.format(plotDir,LEMEString,nuCut)):
      exec("tmp = localDrawCorrelationMatrix(effNumeratorHist_lowNu_{0}_nuCut_{1}_scaled)".format(LEMEString,nuCut))
      tmp.Draw("colz")
    with makeEnv_TCanvas('{0}/fluxComponents/errorSummary_effDenominator_lowNu_{1}_nuCut_{2}.png'.format(plotDir,LEMEString,nuCut)):
      exec("localDrawErrorSummary(plotter,effDenominatorHist_lowNu_{0}_nuCut_{1}_scaled)".format(LEMEString,nuCut))
    with makeEnv_TCanvas('{0}/fluxComponents/correlationMatrix_effDenominator_lowNu_{1}_nuCut_{2}.png'.format(plotDir,LEMEString,nuCut)):
      exec("tmp = localDrawCorrelationMatrix(effDenominatorHist_lowNu_{0}_nuCut_{1}_scaled)".format(LEMEString,nuCut))
      tmp.Draw("colz")
  
  # Plot flux components by nu cut for all signal definitions
  for sigDef in SIGNAL_DEFINITIONS:
    dataRateBreakoutMax = 50 if LEMEString == "ME" else 5
    plotNuCutBreakout('{0}/fluxComponents/dataRate_{1}_{2}_allNuCuts.png'.format(plotDir,sigDef,LEMEString),'dataRate',[eval('dataRateHist_{1}_{0}_nuCut_{2}_scaled'.format(LEMEString,sigDef,nuCut)) for nuCut in range(1,5)],dataRateBreakoutMax,sigDef)
    plotNuCutBreakout('{0}/fluxComponents/effNumerator_{1}_{2}_allNuCuts.png'.format(plotDir,sigDef,LEMEString),'effNumerator',[eval('effNumeratorHist_{1}_{0}_nuCut_{2}_scaled'.format(LEMEString,sigDef,nuCut)) for nuCut in range(1,5)],5,sigDef)
    plotNuCutBreakout('{0}/fluxComponents/effDenominator_{1}_{2}_allNuCuts.png'.format(plotDir,sigDef,LEMEString),'effDenominator',[eval('effDenominatorHist_{1}_{0}_nuCut_{2}_scaled'.format(LEMEString,sigDef,nuCut)) for nuCut in range(1,5)],15,sigDef)
    plotNuCutBreakout('{0}/eff/eff_{1}_{2}_allNuCuts.png'.format(plotDir,sigDef,LEMEString),'eff',[eval('eff_{1}_{0}_nuCut_{2}'.format(LEMEString,sigDef,nuCut)) for nuCut in range(1,5)],1,sigDef,legPos='LR')
  
  # Plot flux by nu cut
  fluxUnnormalizedBreakoutMax = 180 if LEMEString == "ME" else 100
  plotNuCutBreakout('{0}/flux/flux_lowNu_{1}_allNuCuts.png'.format(plotDir,LEMEString),'flux',[eval('flux_lowNu_{0}_nuCut_{1}_scaled'.format(LEMEString,nuCut)) for nuCut in range(1,5)],fluxUnnormalizedBreakoutMax,'lowNu')
  plotNuCutBreakout('{0}/flux/flux_lowNu_{1}_allNuCuts-normalized.png'.format(plotDir,LEMEString),'flux',[eval('flux_lowNu_normalized_{0}_nuCut_{1}_scaled'.format(LEMEString,nuCut)) for nuCut in range(1,5)],fluxUnnormalizedBreakoutMax,'lowNu')

  # Plot Additional Error Breakdowns for each nuCut, for completeness
  for nuCut in range(1,5):
    for sigDef in SIGNAL_DEFINITIONS:
      for component in FLUX_COMPONENTS:

        with makeEnv_TCanvas('{0}/fluxComponents/errorSummary_{1}_{2}_{3}_nuCut_{4}.png'.format(plotDir,component,sigDef,LEMEString,nuCut)):
          exec("setPlotSpecs_ENu({0}Hist_{1}_{2}_nuCut_{3})".format(component,sigDef,LEMEString,nuCut))
          exec("localDrawErrorSummary(plotter,{0}Hist_{1}_{2}_nuCut_{3})".format(component,sigDef,LEMEString,nuCut))
        with makeEnv_TCanvas('{0}/fluxComponents/correlationMatrix_{1}_{2}_{3}_nuCut_{4}.png'.format(plotDir,component,sigDef,LEMEString,nuCut)):
          exec("tmp = localDrawCorrelationMatrix({0}Hist_{1}_{2}_nuCut_{3})".format(component,sigDef,LEMEString,nuCut))
          tmp.Draw("colz")

    with makeEnv_TCanvas('{0}/flux/errorSummary_flux_lowNu_{1}_nuCut_{2}.png'.format(plotDir,LEMEString,nuCut)):
      exec("setPlotSpecs_flux(flux_lowNu_{0}_nuCut_{1}_scaled)".format(LEMEString,nuCut))
      exec("localDrawErrorSummary(plotter,flux_lowNu_{0}_nuCut_{1}_scaled)".format(LEMEString,nuCut))
    with makeEnv_TCanvas('{0}/flux/correlationMatrix_flux_lowNu_{1}_nuCut_{2}.png'.format(plotDir,LEMEString,nuCut)):
      exec("tmp = localDrawCorrelationMatrix(flux_lowNu_{0}_nuCut_{1}_scaled)".format(LEMEString,nuCut))
      tmp.Draw("colz")
    with makeEnv_TCanvas('{0}/flux/errorSummary_flux_lowNu_normalized_{1}_nuCut_{2}.png'.format(plotDir,LEMEString,nuCut)):
      exec("setPlotSpecs_flux(flux_lowNu_normalized_{0}_nuCut_{1}_scaled)".format(LEMEString,nuCut))
      exec("localDrawErrorSummary(plotter,flux_lowNu_normalized_{0}_nuCut_{1}_scaled)".format(LEMEString,nuCut))
    with makeEnv_TCanvas('{0}/flux/correlationMatrix_flux_lowNu_normalized_{1}_nuCut_{2}.png'.format(plotDir,LEMEString,nuCut)):
      exec("tmp = localDrawCorrelationMatrix(flux_lowNu_normalized_{0}_nuCut_{1}_scaled)".format(LEMEString,nuCut))
      tmp.Draw("colz")

#############################################################################################################
### LE,ME Low Nu XSection Plots #############################################################################
#############################################################################################################

for LEMEString in ['LE','ME']:

  exec("xSection_lowNu_{0} = histFile.Get('xSections/xSection_lowNu_{0}')".format(LEMEString))
  
  # Create local copy of lowNu xSection that we'll scale
  exec("xSection_lowNu_{0}_scaled = xSection_lowNu_{0}.Clone()".format(LEMEString)) 
  exec("xSection_lowNu_{0}_scaled.Scale(10**38)".format(LEMEString)) # Remove units for nicer plotting
  
  # Plot lowNu xSection
  with makeEnv_TCanvas('{0}/xSections/xSection_lowNu_{1}.png'.format(plotDir,LEMEString)):
    exec("setPlotSpecs_xSection(xSection_lowNu_{0}_scaled)".format(LEMEString))
    exec("xSection_lowNu_{0}_scaled.GetYaxis().SetRangeUser(0,2)".format(LEMEString))
    exec("xSection_lowNu_{0}_scaled.Draw('PE')".format(LEMEString))
  with makeEnv_TCanvas('{0}/xSections/errorSummary_xSection_lowNu_{1}.png'.format(plotDir,LEMEString)):
    exec("localDrawErrorSummary(plotter,xSection_lowNu_{0}_scaled)".format(LEMEString))
  with makeEnv_TCanvas('{0}/xSections/correlationMatrix_xSection_lowNu_{1}.png'.format(plotDir,LEMEString)):
    exec("tmp = localDrawCorrelationMatrix(xSection_lowNu_{0}_scaled)".format(LEMEString))
    tmp.Draw("colz")

#############################################################################################################
### LE,ME Low Nu Flux Plots #################################################################################
#############################################################################################################

for LEMEString in ['LE','ME']:

  exec("flux_lowNu_{0} = histFile.Get('flux/flux_lowNu_{0}')".format(LEMEString))
  
  # Create local copy of flux that we'll scale
  exec("setPlotSpecs_flux(flux_lowNu_{0})".format(LEMEString))
  exec('flux_lowNu_{0}_scaled = flux_lowNu_{0}.Clone()'.format(LEMEString))
  exec('flux_lowNu_{0}_scaled.Scale(10**6,"width")'.format(LEMEString)) # Remove units for nicer plotting
  exec('local_flux_lowNu_{0} = flux_lowNu_{0}_scaled.GetCVHistoWithError()'.format(LEMEString)) # Remove units for nicer plotting
  
  # Plot lowNu Flux
  with makeEnv_TCanvas('{0}/flux/flux_lowNu_{1}.png'.format(plotDir,LEMEString)):
    #exec("setPlotSpecs_flux(flux_lowNu_{0}_scaled)".format(LEMEString))
    fluxMax = 180 if LEMEString == "ME" else 100
    exec("local_flux_lowNu_{0}.GetYaxis().SetRangeUser(0,fluxMax)".format(LEMEString))
    exec("local_flux_lowNu_{0}.Draw()".format(LEMEString))
  with makeEnv_TCanvas('{0}/flux/errorSummary_flux_lowNu_{1}.png'.format(plotDir,LEMEString)):
    exec("localDrawErrorSummary(plotter,flux_lowNu_{0}_scaled)".format(LEMEString))
  with makeEnv_TCanvas('{0}/flux/correlationMatrix_flux_lowNu_{1}.png'.format(plotDir,LEMEString)):
    exec("tmp = localDrawCorrelationMatrix(flux_lowNu_{0}_scaled)".format(LEMEString))
    tmp.Draw("colz")

#############################################################################################################
### LE,ME Inclusive, high Nu XS Plots #######################################################################
#############################################################################################################

for LEMEString in ['LE','ME']:
  for sigDef in ['highNu','inclusive']:

    exec("xSection_{0}_{1} = histFile.Get('xSections/xSection_{0}_{1}')".format(sigDef,LEMEString))
    
    # Create local copy of inclusive xSections that we'll scale
    exec("setPlotSpecs_xSection(xSection_{0}_{1})".format(sigDef,LEMEString))
    exec("xSection_{0}_{1}_scaled = xSection_{0}_{1}.GetCVHistoWithError()".format(sigDef,LEMEString)) 
    exec("xSection_{0}_{1}_scaled.Scale(10**38)".format(sigDef,LEMEString)) # Remove units for nicer plotting
    
    # Plot inclusive xSections
    with makeEnv_TCanvas('{0}/xSections/xSection_{1}_{2}.png'.format(plotDir,sigDef,LEMEString)):
      exec("xSection_{0}_{1}_scaled.Draw('PE')".format(sigDef,LEMEString))
    with makeEnv_TCanvas('{0}/xSections/errorSummary_xSection_{1}_{2}.png'.format(plotDir,sigDef,LEMEString)):
      exec("localDrawErrorSummary(plotter,xSection_{0}_{1})".format(sigDef,LEMEString))
    with makeEnv_TCanvas('{0}/xSections/correlationMatrix_xSection_{1}_{2}.png'.format(plotDir,sigDef,LEMEString)):
      exec("tmp = localDrawCorrelationMatrix(xSection_{0}_{1})".format(sigDef,LEMEString))
      tmp.Draw("colz")

#############################################################################################################
### ME high Nu Flux Plot ####################################################################################
#############################################################################################################

for fakeFold in [1]:

  flux_highNu_ME = histFile.Get('flux/flux_highNu_ME')
  
  # Create local copy of flux that we'll scale
  flux_highNu_ME_scaled = flux_highNu_ME.Clone() 
  flux_highNu_ME_scaled.Scale(10**6,"width") # Remove units for nicer plotting
  flux_highNu_ME_scaled.Scale(10**4) # Change units from cm^{-2} to m^{-2}
  
  with makeEnv_TCanvas('{0}/flux/flux_highNu_ME.png'.format(plotDir)):
    setPlotSpecs_flux(flux_highNu_ME_scaled)
    flux_highNu_ME_scaled.Draw()
  
  with makeEnv_TCanvas('{0}/flux/errorSummary_flux_highNu_ME.png'.format(plotDir)):
    #flux_highNu_ME_scaled.GetYaxis().SetMaxDigits(2)
    #ROOT.TGaxis.SetMaxDigits(3)
    localDrawErrorSummary(plotter,flux_highNu_ME_scaled)
    #localDrawErrorSummary(plotter,flux_highNu_ME_scaled,"TR",True,True,0.00001,False,"Genie_FSI")
  with makeEnv_TCanvas('{0}/flux/correlationMatrix_flux_highNu_ME.png'.format(plotDir)):
    tmp = localDrawCorrelationMatrix(flux_highNu_ME_scaled)
    tmp.Draw("colz")

#############################################################################################################
### Chi^2 Extraction ########################################################################################
#############################################################################################################

## LE low-nu
chi2_LE_lowNu_matrix = ROOT.TMatrixD(16,16)
## Exclude Muon_Energy_MINOS from calculation of chi2, because the main point of the chi2 is to see if the pmu-shifted or unshifted result fits better
flux_lowNu_LE_chi2 = flux_lowNu_LE_scaled.Clone("flux_lowNu_LE_chi2")
flux_lowNu_LE_chi2.PopVertErrorBand("Muon_Energy_MINOS")
ppfxFlux_LE_rebinned_chi2 = ppfxFlux_LE_rebinned_scaled.Clone("ppfxFlux_LE_rebinned_chi2")
ppfxFlux_LE_rebinned_chi2.PopVertErrorBand("Muon_Energy_MINOS")

## options:    <        ...      >(       hist1                    hist2            mcScale,useDataErrorMatrix,useOnlyShapeErrors,useModelStat,Chi2ByBin
chi2_LE_lowNu = plotter.Chi2DataMC(flux_lowNu_LE_chi2,ppfxFlux_LE_rebinned_chi2,1.0,True,False,False,chi2_LE_lowNu_matrix)
print 'chi2_LE_lowNu: ', chi2_LE_lowNu

with makeEnv_TCanvas('{0}/covariance/Chi2ContributionMatrix_lowNuFlux_LE.png'.format(plotDir)):
  chi2_LE_lowNu_matrix.Draw("colz")

## ME low-nu,high-nu
for sigDef in ['lowNu','highNu']:

  exec("chi2_ME_{0}_matrix = ROOT.TMatrixD(16,16)".format(sigDef))
  ## Exclude Muon_Energy_MINOS from calculation of chi2, because the main point of the chi2 is to see if the pmu-shifted or unshifted result fits better
  exec("flux_{0}_ME_chi2 = flux_{0}_ME_scaled.Clone(\"flux_{0}_ME_chi2\")".format(sigDef))
  exec("flux_{0}_ME_chi2.PopVertErrorBand(\"Muon_Energy_MINOS\")".format(sigDef))
  exec("ppfxFlux_ME_{0}_rebinned_chi2 = ppfxFlux_ME_{0}_rebinned_scaled.Clone(\"ppfxFlux_ME_{0}_rebinned_chi2\")".format(sigDef))
  exec("ppfxFlux_ME_{0}_rebinned_chi2.PopVertErrorBand(\"Muon_Energy_MINOS\")".format(sigDef))

  exec("chi2_ME_{0} = plotter.Chi2DataMC(flux_{0}_ME_chi2,ppfxFlux_ME_{0}_rebinned_chi2,1.0,True,False,False,chi2_ME_{0}_matrix)".format(sigDef))
  exec("print 'chi2_ME_{0}: ',chi2_ME_{0}".format(sigDef))

  with makeEnv_TCanvas("{0}/covariance/Chi2ContributionMatrix_{1}Flux_ME.png".format(plotDir,sigDef)):
    exec("chi2_ME_{0}_matrix.Draw(\"colz\")".format(sigDef))

## ME low-nu vs high-nu
#chi2_ME_lowNu_highNu_matrix = ROOT.TMatrixD(15,15)
#chi2_ME_lowNU_highNu = plotter.Chi2DataMC(flux_lowNu_ME_scaled,flux_highNu_ME_scaled,1.0,True,False,False,chi2_ME_lowNu_highNu_matrix)
#print 'chi2_ME_lowNu_highNu: ' , chi2_ME_lowNu_highNu

#############################################################################################################
### Extracted Fluxes vs PPFX Comparison Plots ###############################################################
#############################################################################################################

## I fucked this up while thesis editing. Come back to it and make treatment of LE and ME symmetric? Why wouldn't I?
### with makeEnv_TCanvas('{0}/LE_lowNu-VS-PPFX.png'.format(plotDir)):
### 
###   ppfxFlux_LE_originalBinning_scaled.Draw()
###   ppfxFlux_LE_originalBinning_scaled.GetXaxis().SetRangeUser(horizontalAxis_lowerBound,horizontalAxis_upperBound)
###   ppfxFlux_LE_originalBinning_scaled.GetYaxis().SetRangeUser(0,150)
###   tmp = ppfxFlux_LE_lowNu_rebinned_scaled.GetCVHistoWithError()
###   tmp.SetMarkerColor(ROOT.kRed)
###   tmp.SetLineColor(ROOT.kRed)
###   tmp.Draw("same")
###   local_flux_lowNu_LE = flux_lowNu_LE_scaled.GetCVHistoWithError()
###   local_flux_lowNu_LE.Draw("same")
### 
###   leg = declareLegend(3,1,"UR-Flux2")
###   leg.AddEntry(local_flux_lowNu_LE,"LE low-#nu flux","l")
###   leg.AddEntry(ppfxFlux_LE_originalBinning_scaled,"LE PPFX flux","p")
###   leg.AddEntry(tmp,"LE PPFX flux, ","lep")
###   leg.AddEntry(None,"rebinned","")
###   leg.Draw()
### 
###   texbox = declareChi2TextLine(chi2_LE_lowNu)
###   texbox.Draw() 
### 
###   #local_flux_lowNu_LE = flux_lowNu_LE_scaled.GetCVHistoWithError()
###   #local_flux_lowNu_LE.Draw()
###   #ppfxFlux_LE_Nom_fineBins_scaled.Draw("same")
###   #ppfxFlux_LE_rebinned_scaled.Draw("same")
###   #leg.AddEntry(ppfxFlux_LE_Nom_fineBins_scaled,"LE PPFX flux","p")
###   #leg.AddEntry(ppfxFlux_LE_rebinned_scaled,"LE PPFX flux, ","lep")
### 
with makeEnv_TCanvas('{0}/LE_lowNu-VS-PPFX_ratio.png'.format(plotDir)):

  # Make local copies of MnvH1Ds
  local_ppfxFlux_LE_MnvH1D = ppfxFlux_LE_rebinned_scaled.Clone("local_ppfxFlux_LE_MnvH1D") 
  local_ppfxFlux_LE_MnvH1D_lessError = ppfxFlux_LE_rebinned_scaled.Clone("local_ppfxFlux_LE_MnvH1D_lessError") 
  local_ppfxFlux_LE_MnvH1D_lessError.PopVertErrorBand("Flux") # Pop the actual flux errors, then create a dummy replacement of the same name
  local_ppfxFlux_LE_MnvH1D_lessError.AddVertErrorBandAndFillWithCV("Flux",100)
  local_lowNuFlux_LE_MnvH1D = flux_lowNu_LE_scaled.Clone("local_lowNuFlux_LE_MnvH1D")

  #Divide
  local_lowNuFlux_LE_MnvH1D.Divide(local_lowNuFlux_LE_MnvH1D,local_ppfxFlux_LE_MnvH1D_lessError)
  local_ppfxFlux_LE_MnvH1D.Divide(local_ppfxFlux_LE_MnvH1D,local_ppfxFlux_LE_MnvH1D_lessError)
 
  # Extract CV TH1Ds from MnvH1Ds 
  local_lowNuFlux_LE = local_lowNuFlux_LE_MnvH1D.GetCVHistoWithError()
  local_ppfxFlux_LE = local_ppfxFlux_LE_MnvH1D.GetCVHistoWithError()

  setPlotSpecs_fluxRatio(local_ppfxFlux_LE)

  local_ppfxFlux_LE.GetYaxis().SetRangeUser(0.8,1.2)
  local_ppfxFlux_LE.SetFillColor(plotter.mc_error_color)
  local_ppfxFlux_LE.SetFillStyle(plotter.mc_error_style)
  local_ppfxFlux_LE.SetMarkerColor(ROOT.kRed)
  local_ppfxFlux_LE.SetMarkerStyle(1)
  local_ppfxFlux_LE.SetLineColor(ROOT.kRed)
    
  local_ppfxFlux_LE.Draw("E2")
  lineAtOne = ROOT.TLine(horizontalAxis_lowerBound,1,horizontalAxis_upperBound,1)
  lineAtOne.SetLineColor(ROOT.kRed)
  lineAtOne.Draw("same")
  local_lowNuFlux_LE.Draw("same")
  
  leg = declareLegend(2,1,"UR-Flux2")
  leg.AddEntry(local_lowNuFlux_LE,"LE low-#nu flux","lep")
  leg.AddEntry(local_ppfxFlux_LE,"LE PPFX flux, ","lf")
  leg.AddEntry(None,"rebinned","")
  leg.Draw()

  texbox = declareChi2TextLine(chi2_LE_lowNu)
  texbox.Draw() 

for sigDef in ['low','high']:
  with makeEnv_TCanvas('{0}/ME_{1}Nu-VS-PPFX.png'.format(plotDir,sigDef)):

    ppfxFlux_ME_originalBinning_scaled.Draw()
    ppfxFlux_ME_originalBinning_scaled.GetXaxis().SetRangeUser(horizontalAxis_lowerBound,horizontalAxis_upperBound)
    ppfxFlux_ME_originalBinning_scaled.GetYaxis().SetRangeUser(0,150)
    exec("tmp = ppfxFlux_ME_{0}Nu_rebinned_scaled.GetCVHistoWithError()".format(sigDef)) 
    tmp.SetMarkerColor(ROOT.kRed)
    tmp.SetLineColor(ROOT.kRed)
    tmp.Draw("same")
    exec('local_flux_{0}Nu_ME = flux_{0}Nu_ME_scaled.GetCVHistoWithError()'.format(sigDef))
    exec('local_flux_{0}Nu_ME.Draw("same")'.format(sigDef))

    leg = declareLegend(3,1,"UR-Flux2")
    exec('leg.AddEntry(local_flux_{0}Nu_ME,"ME {0}-#nu flux","l")'.format(sigDef))
    leg.AddEntry(ppfxFlux_ME_originalBinning_scaled,"ME PPFX flux","p")
    leg.AddEntry(tmp,"ME PPFX flux, ","lep")
    leg.AddEntry(None,"rebinned","")
    leg.Draw()

    exec('texbox = declareChi2TextLine(chi2_ME_{0}Nu)'.format(sigDef))
    texbox.Draw() 

  with makeEnv_TCanvas('{0}/ME_{1}Nu-VS-PPFX_ratio.png'.format(plotDir,sigDef)):

    # Make local copies of MnvH1Ds
    exec("local_ppfxFlux_ME_MnvH1D = ppfxFlux_ME_{0}Nu_rebinned_scaled.Clone(\"local_ppfxFlux_ME_MnvH1D\")".format(sigDef))
    exec("local_ppfxFlux_ME_MnvH1D_lessError = ppfxFlux_ME_{0}Nu_rebinned_scaled.Clone(\"local_ppfxFlux_ME_MnvH1D_lessError\")".format(sigDef)) # Need a second copy that I can pop the Flux error out of
    local_ppfxFlux_ME_MnvH1D_lessError.PopVertErrorBand("Flux") # Pop the actual flux errors, then create a dummy replacement of the same name
    local_ppfxFlux_ME_MnvH1D_lessError.AddVertErrorBandAndFillWithCV("Flux",100)
    exec('local_{0}NuFlux_ME_MnvH1D = flux_{0}Nu_ME_scaled.Clone("local_{0}NuFlux_ME_MnvH1D")'.format(sigDef))

    # Divide
    exec('local_{0}NuFlux_ME_MnvH1D.Divide(local_{0}NuFlux_ME_MnvH1D,local_ppfxFlux_ME_MnvH1D_lessError)'.format(sigDef))
    local_ppfxFlux_ME_MnvH1D.Divide(local_ppfxFlux_ME_MnvH1D,local_ppfxFlux_ME_MnvH1D_lessError)

    # Extract CV TH1Ds from MnvH1Ds
    exec('local_{0}NuFlux_ME = local_{0}NuFlux_ME_MnvH1D.GetCVHistoWithError()'.format(sigDef))
    local_ppfxFlux_ME = local_ppfxFlux_ME_MnvH1D.GetCVHistoWithError() 

    setPlotSpecs_fluxRatio(local_ppfxFlux_ME)
 
    local_ppfxFlux_ME.SetFillColor(plotter.mc_error_color)
    local_ppfxFlux_ME.SetFillStyle(plotter.mc_error_style)
    local_ppfxFlux_ME.SetMarkerColor(ROOT.kRed)
    local_ppfxFlux_ME.SetMarkerStyle(1)
    local_ppfxFlux_ME.SetLineColor(ROOT.kRed)
    #local_ppfxFlux_ME.SetLineStyle(1)

    local_ppfxFlux_ME.Draw("E2")  
    #local_ppfxFlux_ME.Draw("same")
    lineAtOne = ROOT.TLine(horizontalAxis_lowerBound,1,horizontalAxis_upperBound,1)
    lineAtOne.SetLineColor(ROOT.kRed)
    lineAtOne.Draw("same")
    exec('local_{0}NuFlux_ME.Draw("same")'.format(sigDef))
    
    leg = declareLegend(2,1,"UR-Flux2")
    exec('leg.AddEntry(local_{0}NuFlux_ME,"ME {0}-#nu flux","lep")'.format(sigDef))
    leg.AddEntry(local_ppfxFlux_ME,"ME PPFX flux, ","lf")
    leg.AddEntry(None,"rebinned","")
    leg.Draw()
    
    exec('texbox = declareChi2TextLine(chi2_ME_{0}Nu)'.format(sigDef))
    texbox.Draw() 

## This has become difficult because of binning differences
     
with makeEnv_TCanvas('{0}/ME_lowNu-VS-ME_highNu_ratio.png'.format(plotDir)):

  #local_lowNuFlux_ME = flux_lowNu_ME_scaled.GetCVHistoWithError()
  #local_highNuFlux_ME = flux_highNu_ME_scaled.GetCVHistoWithError()
  #setPlotSpecs_fluxRatio(local_lowNuFlux_ME)
  #local_highNuFlux_ME.Divide(local_highNuFlux_ME,local_lowNuFlux_ME)
  #local_lowNuFlux_ME.Divide(local_lowNuFlux_ME,local_lowNuFlux_ME)
  #
  #local_lowNuFlux_ME.GetYaxis().SetRangeUser(0.8,1.2)
  #local_lowNuFlux_ME.SetMarkerColor(ROOT.kRed)
  #local_lowNuFlux_ME.SetLineColor(ROOT.kRed)
  #  
  #local_lowNuFlux_ME.Draw()
  #local_highNuFlux_ME.Draw("same")
  #
  #leg = declareLegend(2,1,"UR-Flux2")
  #leg.AddEntry(local_lowNuFlux_ME,"ME low-#nu flux","l")
  #leg.AddEntry(local_highNuFlux_ME,"ME high-#nu flux","l")
  #leg.Draw()

  # Make local copies of MnvH1Ds
  local_ppfxFlux_ME_lowNu_MnvH1D = ppfxFlux_ME_lowNu_rebinned_scaled.Clone("local_ppfxFlux_ME_lowNu_MnvH1D") 
  local_ppfxFlux_ME_highNu_MnvH1D = ppfxFlux_ME_highNu_rebinned_scaled.Clone("local_ppfxFlux_ME_highNu_MnvH1D") 
  local_ppfxFlux_ME_lowNu_MnvH1D_lessError = ppfxFlux_ME_lowNu_rebinned_scaled.Clone("local_ppfxFlux_ME_lowNu_MnvH1D_lessError") # Need a second copy that I can pop the Flux error out of
  local_ppfxFlux_ME_lowNu_MnvH1D_lessError.PopVertErrorBand("Flux") # Pop the actual flux errors, then create a dummy replacement of the same name
  local_ppfxFlux_ME_lowNu_MnvH1D_lessError.AddVertErrorBandAndFillWithCV("Flux",100)
  local_ppfxFlux_ME_highNu_MnvH1D_lessError = ppfxFlux_ME_highNu_rebinned_scaled.Clone("local_ppfxFlux_ME_highNu_MnvH1D_lessError") # Need a second copy that I can pop the Flux error out of
  local_ppfxFlux_ME_highNu_MnvH1D_lessError.PopVertErrorBand("Flux") # Pop the actual flux errors, then create a dummy replacement of the same name
  local_ppfxFlux_ME_highNu_MnvH1D_lessError.AddVertErrorBandAndFillWithCV("Flux",100)
  local_lowNuFlux_ME_MnvH1D = flux_lowNu_ME_scaled.Clone("local_lowNuFlux_ME_MnvH1D")
  local_highNuFlux_ME_MnvH1D = flux_highNu_ME_scaled.Clone("local_highNuFlux_ME_MnvH1D")

  # Divide
  local_lowNuFlux_ME_MnvH1D.Divide(local_lowNuFlux_ME_MnvH1D,local_ppfxFlux_ME_lowNu_MnvH1D_lessError)
  local_highNuFlux_ME_MnvH1D.Divide(local_highNuFlux_ME_MnvH1D,local_ppfxFlux_ME_highNu_MnvH1D_lessError)
  local_ppfxFlux_ME_highNu_MnvH1D.Divide(local_ppfxFlux_ME_highNu_MnvH1D,local_ppfxFlux_ME_highNu_MnvH1D_lessError)

  # Extract CV TH1Ds from MnvH1Ds
  local_lowNuFlux_ME = local_lowNuFlux_ME_MnvH1D.GetCVHistoWithError()
  local_highNuFlux_ME = local_highNuFlux_ME_MnvH1D.GetCVHistoWithError()
  local_ppfxFlux_ME = local_ppfxFlux_ME_highNu_MnvH1D.GetCVHistoWithError() 

  setPlotSpecs_fluxRatio(local_ppfxFlux_ME)
 
  local_ppfxFlux_ME.SetFillColor(plotter.mc_error_color)
  local_ppfxFlux_ME.SetFillStyle(plotter.mc_error_style)
  local_ppfxFlux_ME.SetMarkerColor(ROOT.kRed)
  local_ppfxFlux_ME.SetMarkerStyle(1)
  local_ppfxFlux_ME.SetLineColor(ROOT.kRed)
  local_ppfxFlux_ME.Draw("E2")  

  lineAtOne = ROOT.TLine(horizontalAxis_lowerBound,1,horizontalAxis_upperBound,1)
  lineAtOne.SetLineColor(ROOT.kRed)
  lineAtOne.Draw("same")

  local_lowNuFlux_ME.Draw("same")

  local_highNuFlux_ME.SetMarkerColor(ROOT.kBlue)
  local_highNuFlux_ME.SetLineColor(ROOT.kBlue)
  local_highNuFlux_ME.SetLineStyle(9)
  local_highNuFlux_ME.Draw("same")
  
  leg = declareLegend(2,1,"UR-Flux2")
  leg.AddEntry(local_lowNuFlux_ME,"ME low-#nu flux","lep")
  leg.AddEntry(local_highNuFlux_ME,"ME high-#nu flux","lep")
  leg.AddEntry(local_ppfxFlux_ME,"ME PPFX flux, ","lf")
  leg.AddEntry(None,"rebinned","")
  leg.Draw()
  
  #texbox = declareChi2TextLine(chi2_ME_lowNu_highNu)
  #texbox.Draw() 

# Make local copies of MnvH1Ds
local_lowNuHighNuFluxRatio = flux_lowNu_ME_scaled.Clone("local_lowNuHighNuFluxRatio")
local_lowNuFlux_ME_MnvH1D = flux_lowNu_ME_scaled.Clone("local_lowNuFlux_ME_MnvH1D")
local_highNuFlux_ME_MnvH1D = flux_highNu_ME_scaled.Clone("local_highNuFlux_ME_MnvH1D")

## This is tough to do right now because of the different binning between ME low- and high-nu fluxes
if False:
  # Divide
  local_lowNuHighNuFluxRatio.Divide(local_lowNuFlux_ME_MnvH1D,local_highNuFlux_ME_MnvH1D)
  # Test
  local_lowNuHighNuFluxRatio.PopVertErrorBand("Flux")
  
  with makeEnv_TCanvas('{0}/ME_lowNu-VS-ME_highNu.png'.format(plotDir)):
  
    setPlotSpecs_fluxRatio(local_lowNuHighNuFluxRatio)
    local_lowNuHighNuFluxRatio.Draw()
   
    lineAtOne = ROOT.TLine(horizontalAxis_lowerBound,1,horizontalAxis_upperBound,1)
    lineAtOne.SetLineColor(ROOT.kRed)
    lineAtOne.Draw("same")
  
    leg = declareLegend(2,1,"UR-Flux2")
    leg.AddEntry(local_lowNuHighNuFluxRatio,"#frac{ME low-#nu flux}{ME high-#nu flux}","lep")
    leg.Draw()
    
  with makeEnv_TCanvas('{0}/errorSummary_ME_lowNu-VS-ME_highNu.png'.format(plotDir)):
    localDrawErrorSummary(plotter,local_lowNuHighNuFluxRatio)
  with makeEnv_TCanvas('{0}/correlationMatrix_ME_lowNu-VS-ME_highNu.png'.format(plotDir)):
    tmp = localDrawCorrelationMatrix(local_lowNuHighNuFluxRatio)
    tmp.Draw("colz")

  # Make local copies of MnvH1Ds
  local_lowNuHighNuXSectionRatio = xSection_inclusive_LE_scaled.Clone("local_lowNuHighNuXSectionRatio")
  local_lowNuXSection_ME_MnvH1D = xSection_inclusive_LE_scaled.Clone("local_lowNuXSection_ME_MnvH1D")
  local_highNuXSection_ME_MnvH1D = xSection_inclusive_ME_scaled.Clone("local_highNuXSection_ME_MnvH1D")
  
  ## Temporary - this will be fixed when I respin LE hists
  # Add missing error bands of each other
  local_lowNuXSection_ME_MnvH1D.AddMissingErrorBandsAndFillWithCV(local_highNuXSection_ME_MnvH1D)
  local_highNuXSection_ME_MnvH1D.AddMissingErrorBandsAndFillWithCV(local_lowNuXSection_ME_MnvH1D)
  
  # Divide
  local_lowNuHighNuXSectionRatio.Divide(local_lowNuXSection_ME_MnvH1D,local_highNuXSection_ME_MnvH1D)
  
  with makeEnv_TCanvas('{0}/LE_lowNu-VS-ME_lowNu.png'.format(plotDir)):
  
    setPlotSpecs_fluxRatio(local_lowNuHighNuXSectionRatio)
    local_lowNuHighNuXSectionRatio.Draw()
   
    lineAtOne = ROOT.TLine(horizontalAxis_lowerBound,1,horizontalAxis_upperBound,1)
    lineAtOne.SetLineColor(ROOT.kRed)
    lineAtOne.Draw("same")
  
    leg = declareLegend(2,1,"UR-Flux2")
    leg.AddEntry(local_lowNuHighNuXSectionRatio,"#frac{ME low-#nu flux}{ME high-#nu flux}","lep")
    leg.Draw()
    
  with makeEnv_TCanvas('{0}/errorSummary_LE_lowNu-VS-ME_lowNu.png'.format(plotDir)):
    localDrawErrorSummary(plotter,local_lowNuHighNuXSectionRatio)
  with makeEnv_TCanvas('{0}/correlationMatrix_LE_lowNu-VS-ME_lowNu.png'.format(plotDir)):
    tmp = localDrawCorrelationMatrix(local_lowNuHighNuXSectionRatio)
    tmp.Draw("colz")

#############################################################################################################
### Covariance Matrix Tests #################################################################################
#############################################################################################################

with makeEnv_TCanvas('{0}/covariance/covarianceMatrix_lowNuFlux_ME.png'.format(plotDir)) as canvas:
  tmp = flux_lowNu_ME.GetTotalErrorMatrix(True,True)
  tmp.Draw("colz")
  canvas.canvas.SetLogz()

with makeEnv_TCanvas('{0}/covariance/covarianceMatrix_highNuFlux_ME.png'.format(plotDir)):
  tmp = flux_highNu_ME.GetTotalErrorMatrix(True,True)
  tmp.Draw("colz")

#############################################################################################################
### Correlation Matrix Tests #################################################################################
#############################################################################################################

# Correlation matrices should have a different color palette
#plotter.SetROOT6Palette(87)

for LEMEString in ["LE","ME"]:

  with makeEnv_TCanvas('{0}/covariance/correlationMatrix_lowNuFlux_{1}.png'.format(plotDir,LEMEString)) as canvas:
    exec("tmp = flux_lowNu_{0}.GetTotalCorrelationMatrix(True,True)".format(LEMEString))
    tmp.Draw("colz")

  if not LEMEString =="LE": 
    with makeEnv_TCanvas('{0}/covariance/correlationMatrix_highNuFlux_{1}.png'.format(plotDir,LEMEString)):
      exec("tmp = flux_highNu_{0}.GetTotalCorrelationMatrix(True,True)".format(LEMEString))
      tmp.Draw("colz")

    with makeEnv_TCanvas('{0}/covariance/correlationMatrix_ppfxFlux_{1}.png'.format(plotDir,LEMEString)):
      exec("tmp = ppfxFlux_{0}_highNu_rebinned.GetTotalCorrelationMatrix(True,True)".format(LEMEString))
      tmp.Draw("colz")
 
  with makeEnv_TCanvas('{0}/covariance/correlationMatrix_inclusiveXSection_{1}.png'.format(plotDir,LEMEString)):
    exec("tmp = xSection_inclusive_{0}.GetTotalCorrelationMatrix(True,True)".format(LEMEString))
    tmp.Draw("colz")

# Revert to standard color palette
#ROOT.gStyle.SetPalette(54)

#############################################################################################################
### Migration Matrices ######################################################################################
#############################################################################################################

for LEMEString in ['LE','ME']:
  for sigDef in SIGNAL_DEFINITIONS:

    exec("migrationMatrix_Enu_{0}_{1} = histFile.Get('migrationMatrices/migrationMatrix_Enu_{0}_{1}')".format(sigDef,LEMEString))

    with makeEnv_TCanvas('{0}/migrationMatrices/migrationMatrix_Enu_{1}_{2}.png'.format(plotDir,sigDef,LEMEString)) as canvas:
      exec("setPlotSpecs_Migration(migrationMatrix_Enu_{0}_{1})".format(sigDef,LEMEString))
      exec("rowNormalize(migrationMatrix_Enu_{0}_{1})".format(sigDef,LEMEString))
      exec("migrationMatrix_Enu_{0}_{1}.Draw(\"colz\")".format(sigDef,LEMEString))
      #canvas.canvas.SetLogz()

#############################################################################################################
### Emu by Range and Curve ##################################################################################
#############################################################################################################

for LEMEString in ['LE','ME']:
  for sigDef in SIGNAL_DEFINITIONS:

    exec("EMu_rangeOnly_{0}_{1} = histFile.Get('all-{1}/EMu_rangeOnly_{0}_{1}')".format(sigDef,LEMEString))
    exec("EMu_rangeOnly_{0}_{1}.Scale(1.0,\"width\")".format(sigDef,LEMEString))
    exec("EMu_rangeAndCurve_{0}_{1} = histFile.Get('all-{1}/EMu_rangeAndCurve_{0}_{1}')".format(sigDef,LEMEString))
    exec("EMu_rangeAndCurve_{0}_{1}.Scale(1.0,\"width\")".format(sigDef,LEMEString))

    with makeEnv_TCanvas('{0}/reco/EMu_reco_{1}_{2}.png'.format(plotDir,sigDef,LEMEString)) as canvas:
      exec("setPlotSpecs_EMu(EMu_rangeOnly_{0}_{1})".format(sigDef,LEMEString))
      exec("EMu_rangeOnly_{0}_{1}_local = EMu_rangeOnly_{0}_{1}.GetCVHistoWithError()".format(sigDef,LEMEString))
      exec("EMu_rangeAndCurve_{0}_{1}_local = EMu_rangeAndCurve_{0}_{1}.GetCVHistoWithError()".format(sigDef,LEMEString))
      plotMax = 1.2*10**6 if LEMEString == "ME" else 1*10**5
      exec("EMu_rangeOnly_{0}_{1}_local.SetMaximum(plotMax)".format(sigDef,LEMEString))
      exec("EMu_rangeOnly_{0}_{1}_local.Draw()".format(sigDef,LEMEString))
      exec("EMu_rangeAndCurve_{0}_{1}_local.SetMarkerColor(ROOT.kRed)".format(sigDef,LEMEString))
      exec("EMu_rangeAndCurve_{0}_{1}_local.SetLineColor(ROOT.kRed)".format(sigDef,LEMEString))
      exec("EMu_rangeAndCurve_{0}_{1}_local.Draw(\"same\")".format(sigDef,LEMEString))

      leg = ROOT.TLegend(0.325,0.7,0.95,0.9)
      setPlotSpecs_legend(leg)
      exec("leg.AddEntry(EMu_rangeOnly_{0}_{1}_local,\"Reco by range only\",\"lep\")".format(sigDef,LEMEString))
      exec("leg.AddEntry(EMu_rangeAndCurve_{0}_{1}_local,\"Reco by range and curvature\",\"lep\")".format(sigDef,LEMEString))
      leg.Draw()

#############################################################################################################
### LE,ME 2D Inclusive Plots, for illustration purposes #####################################################
#############################################################################################################

#nuLines= [ROOT.TLine(i,j,i,k) for i,j,k in zip([3,7,12],[0.3,0.5,1.0],[0.5,1.0,2.0])]+[ROOT.TLine(i,j,k,j) for i,j,k in zip([1,3,7,12],[0.3,0.5,1.0,2.0],[3,7,12,horizontalAxis_upperBound])]
nuLines= [ROOT.TLine(i,j,i,k) for i,j,k in zip([3,7,12],[0.3,0.5,1.0],[0.5,1.0,2.0])]+[ROOT.TLine(i,j,k,j) for i,j,k in zip([horizontalAxis_lowerBound,3,7,12],[0.3,0.5,1.0,2.0],[3,7,12,horizontalAxis_upperBound])]
for line in nuLines:
  line.SetLineColor(ROOT.kRed+1)
  line.SetLineWidth(3)
amitNuLine = ROOT.TLine(horizontalAxis_lowerBound,.8,horizontalAxis_upperBound,.8)
amitNuLine.SetLineColor(ROOT.kGreen+1)
amitNuLine.SetLineWidth(3)

for LEMEString in ['LE','ME']:
  exec("dataRateHist2D_inclusive_{0} = histFile.Get('all-{0}/dataRateHist2D_inclusive_{0}')".format(LEMEString))
  exec("dataRateHist2D_inclusive_{0}.Scale(1,\"width\") # bin-width normalize".format(LEMEString))
  exec("nEntries_{0} = dataRateHist2D_inclusive_{0}.GetEntries()".format(LEMEString))

#scaleFactor = nEntries_LE/nEntries_ME
#dataRateHist2D_inclusive_ME.Scale(scaleFactor)
#print 'scaleFactor: ' , 1/scaleFactor

# Write out cumulative 2D inclusive hist
for LEMEString in ['LE','ME']:

  with makeEnv_TCanvas('{0}/2D/dataRate_2DInclusive_{1}.png'.format(plotDir,LEMEString)) as canvas:
    canvas.canvas.SetLeftMargin(0.12);
    canvas.canvas.SetRightMargin(0.20);
    exec("setPlotSpecs_2D(dataRateHist2D_inclusive_{0})".format(LEMEString))
    exec("dataRateHist2D_inclusive_{0}.GetZaxis().SetTitle('# {0} Data Events/GeV^{{2}}')".format(LEMEString))
    exec("dataRateHist2D_inclusive_{0}.GetZaxis().SetTitleOffset(1.2)".format(LEMEString))
    exec("dataRateHist2D_inclusive_{0}.SetMinimum(10**0)".format(LEMEString))
    #exec("dataRateHist2D_inclusive_{0}.SetMaximum(10**5)".format(LEMEString))
    exec("dataRateHist2D_inclusive_{0}.SetMaximum(10**6)".format(LEMEString))
    exec("dataRateHist2D_inclusive_{0}.GetYaxis().SetTitleOffset(0.8)".format(LEMEString))
    exec("dataRateHist2D_inclusive_{0}.Draw('colz')".format(LEMEString))
    ROOT.setCustomRochesterPalette()
    canvas.canvas.SetLogz()
    #plotter.WritePreliminary("TC",0.035,0,0.115,True) 
    for line in nuLines:
      line.Draw()
    t1 = ROOT.TLatex(18,2.1,"#it{high-#nu}")
    t1.SetTextSize(0.05)
    t1.SetTextColor(ROOT.kRed+1)
    t1.Draw()
    t2 = ROOT.TLatex(18,1.7,"#it{low-#nu}")
    t2.SetTextSize(0.05)
    t2.SetTextColor(ROOT.kRed+1)
    t2.Draw()
    #if LEMEString == 'ME':
    #  t3 = ROOT.TLatex(4,2.5,"scaled down by factor of 13")
    #  t3.SetTextSize(0.03)
    #  t3.SetTextAngle(90)
    #  t3.SetTextColor(ROOT.kWhite)
    #  t3.Draw()
    #  t4 = ROOT.TLatex(15,0.9,"#it{Amit's low-#nu}")
    #  t4.SetTextSize(0.05)
    #  t4.SetTextColor(ROOT.kGreen+1)
    #  #t4.Draw()
    #  #amitNuLine.Draw()

#############################################################################################################
### Inclusive XSection vs World Data ########################################################################
#############################################################################################################

## Pull MC prediction from separate file
fakeDataHistFile = ROOT.TFile("/minerva/data/users/finer/highNu/analysisHists/test_thesisProcessing_noFluxConstraint_noMuonCVShift_2020-06-18.root")
fakeDataHist = fakeDataHistFile.Get("xSections/xSection_inclusive_ME")

with makeEnv_TCanvas("{0}/xSections/xSectionPerE_inclusive_world.png".format(plotDir)):  

  binCenterVals = histFile.Get('xSections/h_binCenterVals_ME')

  ## LE data
  xSectionPerE_inclusive_LE = xSection_inclusive_LE.Clone('xSectionPerE_inclusive_LE')
  xSectionPerE_inclusive_LE.Divide(xSectionPerE_inclusive_LE,binCenterVals)
  xSectionPerE_inclusive_LE.Scale(10**38)
 
  xSectionPerE_inclusive_LE.SetMarkerColor(ROOT.kViolet+2)
  xSectionPerE_inclusive_LE.SetLineColor(ROOT.kViolet+2)

  ## ME data
  xSectionPerE_inclusive_ME = xSection_inclusive_ME.Clone('xSectionPerE_inclusive_ME')
  xSectionPerE_inclusive_ME.Divide(xSectionPerE_inclusive_ME,binCenterVals)
  xSectionPerE_inclusive_ME.Scale(10**38)

  ## GENIE
  #xSectionPerE_inclusive_GENIE = effDenominatorHist_inclusive_ME.Clone('xSectionPerE_inclusive_GENIE')
  xSectionPerE_inclusive_GENIE = fakeDataHist.Clone('xSectionPerE_inclusive_GENIE')
  xSectionPerE_inclusive_GENIE.Divide(xSectionPerE_inclusive_GENIE,binCenterVals)
  xSectionPerE_inclusive_GENIE.Scale(10**38)

  xSectionPerE_inclusive_GENIE.SetMarkerColor(plotter.good_colors[6])
  xSectionPerE_inclusive_GENIE.SetLineColor(plotter.good_colors[6])

  local_xSectionPerE_inclusive_LE = xSectionPerE_inclusive_LE.GetCVHistoWithError()
  local_xSectionPerE_inclusive_ME = xSectionPerE_inclusive_ME.GetCVHistoWithError()
  local_xSectionPerE_inclusive_GENIE = xSectionPerE_inclusive_GENIE.GetCVHistoWithError()

  local_xSectionPerE_inclusive_LE.GetYaxis().SetTitle('#sigma/E (10^{-38}cm^{2}/nucleon/GeV)')
  local_xSectionPerE_inclusive_LE.GetYaxis().SetRangeUser(0.5,1.25)
  #local_xSectionPerE_inclusive_LE.GetYaxis().SetRangeUser(0.6,0.8)
  local_xSectionPerE_inclusive_LE.GetXaxis().SetRangeUser(0,50)
  local_xSectionPerE_inclusive_ME.GetXaxis().SetRangeUser(0,50)
  local_xSectionPerE_inclusive_LE.Draw()
  local_xSectionPerE_inclusive_ME.Draw("same")
  local_xSectionPerE_inclusive_GENIE.Draw("same hist")
  
  ###################################################### MINOS Inclusive
  ######################################################################

  MINOS_bins =      [3,4,5,7,9,12,15,18,22,26,30,36,42,50]
  MINOS_xSections = [0.748,0.711,0.708,0.722,0.699,0.691,0.708,0.689,0.683,0.686,0.675,0.675,0.676]
  MINOS_errors =    [0.061,0.033,0.032,0.045,0.043,0.028,0.020,0.016,0.015,0.016,0.016,0.018,0.019]

  MINOS_InclusiveXSection = ROOT.TH1D( 'MINOS_InclusiveXSection' , 'MINOS_InclusiveXSection' , len(MINOS_bins)-1 , array('d',MINOS_bins))
  for i in range(len(MINOS_bins)-1):
    MINOS_InclusiveXSection.SetBinContent(i+1,MINOS_xSections[i])
    MINOS_InclusiveXSection.SetBinError(i+1,MINOS_errors[i])

  MINOS_InclusiveXSection.SetMarkerColor(plotter.good_colors[3]) 
  MINOS_InclusiveXSection.SetLineColor(plotter.good_colors[3]) 
  MINOS_InclusiveXSection.Draw("same")

  ######################################################## NOMAD Inclusive 
  ######################################################################

  NOMAD_bins =      [2.5,6,7,8,9,10,11,12,13,14,15,17.5,20,22.5,25,27.5,30,35,40,45,50]
  NOMAD_xSections = [0.786,0.763,0.722,0.701,0.716,0.706,0.705,0.697,0.700,0.698,0.698,0.700,0.699,0.694,0.694,0.694,0.677,0.681,0.675,0.682]
  NOMAD_errors =    [0.037,0.038,0.036,0.034,0.034,0.026,0.025,0.025,0.025,0.025,0.025,0.025,0.024,0.024,0.025,0.025,0.026,0.026,0.028,0.027]

  NOMAD_InclusiveXSection = ROOT.TH1D( 'NOMAD_InclusiveXSection' , 'NOMAD_InclusiveXSection' , len(NOMAD_bins)-1 , array('d',NOMAD_bins))
  for i in range(len(NOMAD_bins)-1):
    NOMAD_InclusiveXSection.SetBinContent(i+1,NOMAD_xSections[i])
    NOMAD_InclusiveXSection.SetBinError(i+1,NOMAD_errors[i])

  NOMAD_InclusiveXSection.SetMarkerColor(plotter.good_colors[5]) 
  NOMAD_InclusiveXSection.SetLineColor(plotter.good_colors[5]) 
  NOMAD_InclusiveXSection.Draw("same")

  ######################################################## T2K Inclusive 
  ######################################################################

  T2K_bins =      [0.7,1.5,2.5,4.1]
  T2K_xSections = [1.10/1.1,2.07/2.0,2.29/3.3]
  T2K_errors =    [0.15/1.1,0.27/2.0,0.45/3.3]

  T2K_InclusiveXSection = ROOT.TH1D( 'T2K_InclusiveXSection' , 'T2K_InclusiveXSection' , len(T2K_bins)-1 , array('d',T2K_bins))
  for i in range(len(T2K_bins)-1):
    T2K_InclusiveXSection.SetBinContent(i+1,T2K_xSections[i])
    T2K_InclusiveXSection.SetBinError(i+1,T2K_errors[i])

  T2K_InclusiveXSection.SetMarkerColor(plotter.good_colors[2]) 
  T2K_InclusiveXSection.SetLineColor(plotter.good_colors[2]) 
  T2K_InclusiveXSection.Draw("same,E0X0")

  ############################################################### Redraw
  ######################################################################

  ## Redraw things I want on top
  local_xSectionPerE_inclusive_LE.Draw("same")
  local_xSectionPerE_inclusive_ME.Draw("same")

  ############################################################### Legend
  ######################################################################

  leg = ROOT.TLegend(0.55,0.7,0.95,0.9)
  setPlotSpecs_legend(leg)
  leg.AddEntry(local_xSectionPerE_inclusive_LE,"MINERvA LE","lep")
  leg.AddEntry(local_xSectionPerE_inclusive_ME,"MINERvA ME","lep")
  leg.AddEntry(local_xSectionPerE_inclusive_GENIE,"MnvGENIE-v1","l")
  leg.AddEntry(MINOS_InclusiveXSection,"MINOS","lep")
  leg.AddEntry(NOMAD_InclusiveXSection,"NOMAD","lep")
  leg.AddEntry(T2K_InclusiveXSection,"T2K","ep")
  leg.Draw()

with makeEnv_TCanvas("{0}/xSections/xSectionPerE_inclusive_ratio.png".format(plotDir)):  

  binCenterVals = histFile.Get('xSections/h_binCenterVals_ME')

  ## GENIE
  xSectionPerE_inclusive_ratio_GENIE = fakeDataHist.Clone('xSectionPerE_inclusive_ratio_GENIE')
  xSectionPerE_inclusive_ratio_GENIE_lessErrors = fakeDataHist.Clone('xSectionPerE_inclusive_ratio_GENIE_lessErrors')
  xSectionPerE_inclusive_ratio_GENIE.Divide(xSectionPerE_inclusive_ratio_GENIE,xSectionPerE_inclusive_ratio_GENIE_lessErrors)

  local_xSectionPerE_inclusive_ratio_GENIE = xSectionPerE_inclusive_ratio_GENIE.GetCVHistoWithError()
  local_xSectionPerE_inclusive_ratio_GENIE.GetYaxis().SetTitle('#sigma/MnvGENIE-v1')
  local_xSectionPerE_inclusive_ratio_GENIE.GetYaxis().SetRangeUser(0.75,1.25)
  local_xSectionPerE_inclusive_ratio_GENIE.GetXaxis().SetRangeUser(0,50)

  #local_xSectionPerE_inclusive_ratio_GENIE.SetFillColor(plotter.mc_error_color)
  #local_xSectionPerE_inclusive_ratio_GENIE.SetFillStyle(plotter.mc_error_style)
  local_xSectionPerE_inclusive_ratio_GENIE.SetMarkerColor(ROOT.kRed)
  local_xSectionPerE_inclusive_ratio_GENIE.SetMarkerStyle(1)
  local_xSectionPerE_inclusive_ratio_GENIE.SetLineColor(ROOT.kRed)
  #local_xSectionPerE_inclusive_ratio_GENIE.Draw()  

  ## ME data
  xSectionPerE_inclusive_ratio_ME = xSection_inclusive_ME.Clone('xSectionPerE_inclusive_ratio_ME')
  xSectionPerE_inclusive_ratio_ME.Divide(xSectionPerE_inclusive_ratio_ME,xSectionPerE_inclusive_ratio_GENIE_lessErrors)
  
  local_xSectionPerE_inclusive_ratio_ME = xSectionPerE_inclusive_ratio_ME.GetCVHistoWithError()
  local_xSectionPerE_inclusive_ratio_ME.GetYaxis().SetTitle('#sigma Data/MnvGENIE-v1')
  local_xSectionPerE_inclusive_ratio_ME.GetYaxis().SetRangeUser(0.75,1.25)
  local_xSectionPerE_inclusive_ratio_ME.GetXaxis().SetRangeUser(0,50)

  local_xSectionPerE_inclusive_ratio_ME.Draw()
  
  lineAtOne = ROOT.TLine(0,1,50,1)
  lineAtOne.SetLineColor(ROOT.kRed)
  lineAtOne.Draw("same")

  # Just to get this stupid thing to draw on top of the lineAtOne...
  local_xSectionPerE_inclusive_ratio_ME.Draw("same")

  ## ME data
  xSectionPerE_inclusive_ratio_LE = xSection_inclusive_LE.Clone('xSectionPerE_inclusive_ratio_LE')
  xSectionPerE_inclusive_ratio_LE.Divide(xSectionPerE_inclusive_ratio_LE,xSectionPerE_inclusive_ratio_GENIE_lessErrors)

  local_xSectionPerE_inclusive_ratio_LE = xSectionPerE_inclusive_ratio_LE.GetCVHistoWithError()
  local_xSectionPerE_inclusive_ratio_LE.GetXaxis().SetRangeUser(0,50)

  #local_xSectionPerE_inclusive_ratio_LE.Draw("same")
  
  ############################################################### Legend
  ######################################################################

  leg = ROOT.TLegend(0.55,0.7,0.95,0.9)
  setPlotSpecs_legend(leg)
  #leg.AddEntry(local_xSectionPerE_inclusive_ratio_LE,"MINERvA LE","lep")
  leg.AddEntry(local_xSectionPerE_inclusive_ratio_ME,"MINERvA ME","lep")
  leg.AddEntry(local_xSectionPerE_inclusive_ratio_GENIE,"MnvGENIE-v1","lf")
  leg.Draw()
