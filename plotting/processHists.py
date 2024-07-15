import ROOT
from functions import *
from plottingClasses import *
from binning import *

##set ROOT to batch mode
ROOT.gROOT.SetBatch()

# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)

targetUtils = PlotUtils.TargetUtils()
ROOT.gStyle.SetPalette(54)

################################################################ Preamble above

def undoBinWidthNormalization(hist):
  tmpHist = hist.GetCVHistoWithStatError()
  binWidths = tmpHist.Clone("notImportant")
  for i in range(1,hist.GetNbinsX()):
    binWidths.SetBinContent(i,hist.GetBinWidth(i))
  hist.MultiplySingle(hist,binWidths)
  return hist

def AddNOMADErrorBand( inHist , applyWgts = False ):
  # This is designed for a call with applyWgts to be made only after a call w/o applyWgts has already been made on the MnvH1D 
  NOMAD_uncertainty = 0.036
  uncertaintyWgts = [1.-NOMAD_uncertainty,1.+NOMAD_uncertainty] if applyWgts else [1.,1.]

  #print("Adding NOMAD error band to " , inHist.GetName())
  #print("Applying these weights: " , uncertaintyWgts)

  if not applyWgts: inHist.AddVertErrorBandAndFillWithCV("NOMAD_normalization",2)
  inHist.GetVertErrorBand("NOMAD_normalization").GetHist(0).Scale(uncertaintyWgts[0])
  inHist.GetVertErrorBand("NOMAD_normalization").GetHist(1).Scale(uncertaintyWgts[1])

def getPOT( dataSwitch , filePath ):
  meta_tree = WrapChain( "Meta" )
  meta_tree.Add(filePath)
  pot_used = 0.0
  pot_total = 0.0
  for entry in range(meta_tree.GetEntries()):
    meta_tree.LoadTree(entry)
    pot_used += meta_tree.POT_Used
    pot_total += meta_tree.POT_Total

  if dataSwitch == 'data':
    return pot_used
  else:
    return [pot_used,pot_total] 

################################################################ Define functions above

HISTDIR = "/exp/minerva/data/users/finer/MATAna/2024-06_development/root"

HISTDIR_ROOT_OUTPUT = "/exp/minerva/data/users/finer/MATAna/2024-06_development"
histOutputFilePath = "{0}/processedHists_2024-06-24.root".format(HISTDIR_ROOT_OUTPUT)

################################################################ Specify indir, outpath above

# Make reference hist that has the binning of this analysis
referenceHist = PlotUtils.MnvH1D( 'h_flux_reference' , 'h_flux_reference' , nBins_nuE_LE , array('d',bins_nuE_GeV))

## original PPFX LE
###################
# Declare instance of FluxReweighter
fluxReweighter_LE_originalBinning = PlotUtils.FluxReweighter(14,False,PlotUtils.FluxReweighter.minerva1,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv5,100)
# Get flux with unaltered binning
flux_PPFX_LE_originalBinning = fluxReweighter_LE_originalBinning.GetFluxReweighted(14)
flux_PPFX_LE_originalBinning.SetName('flux_PPFX_LE_originalBinning')
# Get flux with analysis binning 
flux_PPFX_LE_analysisBinning = fluxReweighter_LE_originalBinning.GetRebinnedFluxReweighted(14,referenceHist)
flux_PPFX_LE_analysisBinning.SetName('flux_PPFX_LE_analysisBinning')

## REVISIT THIS?? RDF 2024-03-26
##
## # regenerated PPFX LE
## #fluxReweighter_LE = PlotUtils.FluxReweighter(14,False,PlotUtils.FluxReweighter.minerva1_Rob,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv5)
## fluxReweighter_LE = PlotUtils.FluxReweighter(14,False,PlotUtils.FluxReweighter.minerva1_Dan,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv5,100)
## flux_PPFX_LE_newOriginalBinning = fluxReweighter_LE.GetFluxReweighted(14)
## flux_PPFX_LE_newOriginalBinning.SetName('flux_PPFX_LE_newOriginalBinning')
## flux_PPFX_LE = fluxReweighter_LE.GetRebinnedFluxReweighted(14,referenceHist)
## flux_PPFX_LE.SetName('flux_PPFX_LE')
## # This is the one that will actually be used; undo bin-width normalization
## flux_PPFX_LE_notBinWidthNormalized = undoBinWidthNormalization(fluxReweighter_LE.GetRebinnedFluxReweighted(14,referenceHist))
## flux_PPFX_LE_notBinWidthNormalized.SetName("flux_PPFX_LE_notBinWidthNormalized")

# original PPFX ME w/ nue constraint
###################
#Declare instance of FluxReweighter
fluxReweighter_ME_withNuEConstraint = PlotUtils.FluxReweighter(14,True,PlotUtils.FluxReweighter.minervame1D1M1NWeightedAve,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6,100)
# Get flux with unaltered binning
flux_PPFX_ME_yesNuEConstraint_originalBinning = fluxReweighter_ME_withNuEConstraint.GetFluxReweighted(14)
flux_PPFX_ME_yesNuEConstraint_originalBinning.SetName('flux_PPFX_ME_yesNuEConstraint_originalBinning')
# Get flux with analysis binning 
flux_PPFX_ME_yesNuEConstraint_analysisBinning = fluxReweighter_ME_withNuEConstraint.GetRebinnedFluxReweighted(14,referenceHist)
flux_PPFX_ME_yesNuEConstraint_analysisBinning.SetName('flux_PPFX_ME_yesNuEConstraint_analysisBinning')

# original PPFX ME, no nue constraint
###################
#Declare instance of FluxReweighter
fluxReweighter_ME = PlotUtils.FluxReweighter(14,False,PlotUtils.FluxReweighter.minervame1D1M1NWeightedAve,PlotUtils.FluxReweighter.gen2thin,PlotUtils.FluxReweighter.g4numiv6,100)
# Get flux with unaltered binning; undo bin-width normalization
flux_PPFX_ME_originalBinning_notBinWidthNormalized = undoBinWidthNormalization(fluxReweighter_ME.GetRebinnedFluxReweighted(14,referenceHist))
flux_PPFX_ME_originalBinning_notBinWidthNormalized.SetName("flux_PPFX_ME_originalBinning_notBinWidthNormalized")
# Get flux with analysis binning; undo bin-width normalization
# THIS IS THE FLUX THAT GETS USED FOR THE ANALYSIS
flux_PPFX_ME_analysisBinning_notBinWidthNormalized = undoBinWidthNormalization(fluxReweighter_ME.GetRebinnedFluxReweighted(14,referenceHist))
flux_PPFX_ME_analysisBinning_notBinWidthNormalized.SetName("flux_PPFX_ME_analysisBinning_notBinWidthNormalized")

################################################################ Get fluxes above

PLAYLISTS_LE = [
  #'minerva1',
  #'minerva7',
  #'minerva9',
  #'minerva13',
  #'2p2h'
]

PLAYLISTS_ME = [
  #'minervame1A',
  #'minervame1B',
  #'minervame1C',
  #'minervame1D',
  #'minervame1E',
  #'minervame1F',
  #'minervame1G',
  #'minervame1L',
  #'minervame1M',
  #'minervame1N',
  #'minervame1O',
  #'minervame1P'
  #'ME1A',
  'ME1B',
  'ME1C',
  #'ME1D',
  #'ME1E',
  #'ME1F',
  #'ME1G',
  'ME1L',
  #'ME1M',
]

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

################################################################ Define lists above

histOutputFile = ROOT.TFile(histOutputFilePath,"recreate")
histOutputFile.cd()

histOutputDir_xSections = histOutputFile.mkdir('xSections')
histOutputDir_flux = histOutputFile.mkdir('flux')
histOutputDir_allLE = histOutputFile.mkdir('all-LE')
histOutputDir_allME = histOutputFile.mkdir('all-ME')
histOutputDir_dataMCRatios = histOutputFile.mkdir('dataMCRatios')
histOutputDir_migrationMatrices = histOutputFile.mkdir('migrationMatrices')

## See above RDF 2024-03-26
## # I don't quite remember why this is, but something to do with getting rid of the fluxReweighter object once we don't explicitly need it...
## del fluxReweighter_LE

nTargets = targetUtils.GetTrackerNNucleons(106,False) # 106 planes, not sure what the False is about
print('nTargets: ' , nTargets)
totalDataPOT_LE = 0.0
totalDataPOT_ME = 0.0

#############################################################################################################
### Extract all hists, make local copies, take care of scaling (for LE and ME) ##############################
#############################################################################################################
for PLAYLISTS,isME in zip([PLAYLISTS_LE,PLAYLISTS_ME],[False,True]):
  for playlist in PLAYLISTS: 
    # Make new directory in output file for this playlist
    histOutputDir_playlist = histOutputFile.mkdir(playlist)
   
    print('playlist: ' , playlist)
    print('HISTDIR: ' , HISTDIR)
    print('isME: ' , isME)

    # Define file locations
    #fileString = "XSecInputs_1110_{0}_2024-05-07.root".format(playlist)
    dateString = "2024-06-21" if playlist == "ME1C" else "2024-06-20"
    fileString = "XSecInputs_1110_{0}_{1}.root".format(playlist,dateString)
    histsFileLocation = "{0}/{1}".format(HISTDIR,fileString)
    histsFile = ROOT.TFile(histsFileLocation)

    print('Opening hists file: ' , histsFileLocation)

    #############################################################################
    ## Old POT methodology, some aspects may yet be needed
    ## Specifically: 2p2h data exception, used/total corrective factor in LE
    ##     
    ## # Get POT from Meta for data, MC
    ## mcPOT_used,mcPOT_total = getPOT('mc',mcHistsFileLocation)
    ## mcPOT_ratio = mcPOT_used/mcPOT_total # POT-counting-bug correction factor 
    ## if not playlist == '2p2h': # There is no 2p2h data
    ##   dataPOT = getPOT('data',dataHistsFileLocation)
    ##   print('POT of playlist {0}: {1}'.format(playlist,dataPOT))
    ##   if isME:  totalDataPOT_ME += dataPOT
    ##   else:     totalDataPOT_LE += dataPOT
    ##   scaleFactor = dataPOT/mcPOT_used
    ## else:
    ##   mcPOT_used_2p2h = mcPOT_used # Save these for later
    ##   mcPOT_ratio_2p2h = mcPOT_ratio
    #############################################################################
   
    mcPOT_hist = histsFile.Get("mc_pot")
    mcPOT = mcPOT_hist.GetBinContent(1)
    mcPOT_ratio = 1.0 # revisit this point in context of correcting the LE POT-counting-bug

    dataPOT_hist = histsFile.Get("data_pot")
    dataPOT = dataPOT_hist.GetBinContent(1)

    print('mc POT of playlist {0}: {1}'.format(playlist,mcPOT))
    print('data POT of playlist {0}: {1}'.format(playlist,dataPOT))

    if isME:  totalDataPOT_ME += dataPOT
    else:     totalDataPOT_LE += dataPOT
    scaleFactor = dataPOT/mcPOT
 
    # Extract flux components for each signal definition
    for sigDef in SIGNAL_DEFINITIONS:
      key = '{0}_{1}'.format(sigDef,playlist)
      if not playlist == '2p2h': # There is no 2p2h data 
        exec("dataRateHist2D_{0} = histsFile.Get('selection_data_enu_ehad_{1}')".format(key,sigDef))
        exec("dataRateHist_{0} = histsFile.Get('selection_data_enu_ehad_{1}').ProjectionX()".format(key,sigDef))
      exec("effNumeratorHist2D_{0} = histsFile.Get('effnum_enu_ehad_{1}')".format(key,sigDef))
      exec("effNumeratorHist_{0} = histsFile.Get('effnum_enu_ehad_{1}').ProjectionX()".format(key,sigDef))
      exec("effDenominatorHist2D_{0} = histsFile.Get('effdenom_enu_ehad_{1}')".format(key,sigDef))
      exec("effDenominatorHist_{0} = histsFile.Get('effdenom_enu_ehad_{1}').ProjectionX()".format(key,sigDef))
      ##exec("effDenominatorHist2D_addKinematicCuts_{0} = histsFile.Get('h_ENu_VS_nu_{1}_truth_addKinematicCuts')".format(key,sigDef))
      ##exec("effDenominatorHist_addKinematicCuts_{0} = histsFile.Get('h_ENu_VS_nu_{1}_truth_addKinematicCuts').ProjectionX()".format(key,sigDef))
      exec("migrationMatrix_Enu_{0} = histsFile.Get('migration_enu')".format(key,sigDef))
      ##exec("EMu_rangeOnly_{0} = histsFile.Get('h_hist_EMu_rangeOnly_{1}')".format(key,sigDef))
      ##exec("EMu_rangeAndCurve_{0} = histsFile.Get('h_hist_EMu_rangeAndCurve_{1}')".format(key,sigDef))

      # Scale MC to the POT of data, skipping 2p2h for now
      if not playlist == '2p2h':
        exec("effNumeratorHist_{0}.Scale(scaleFactor)".format(key))
        exec("effNumeratorHist2D_{0}.Scale(scaleFactor)".format(key))
        exec("effDenominatorHist_{0}.Scale(scaleFactor*mcPOT_ratio)".format(key)) # Truth distribution doesn't have data-quality pre-selection that reco MC does. This is the Eroica 'POT-Counting' bug, which is corrected for by scaling the truth distribution by the ratio of POT_Used to POT_Total
        exec("effDenominatorHist2D_{0}.Scale(scaleFactor*mcPOT_ratio)".format(key))
        exec("migrationMatrix_Enu_{0}.Scale(scaleFactor*mcPOT_ratio)".format(key))
      ##   ##exec("effDenominatorHist_addKinematicCuts_{0}.Scale(scaleFactor*mcPOT_ratio)".format(key))
      ##   ##exec("effDenominatorHist2D_addKinematicCuts_{0}.Scale(scaleFactor*mcPOT_ratio)".format(key))
      ##   exec("EMu_rangeOnly_{0}.Scale(scaleFactor*mcPOT_ratio)".format(key))
      ##   exec("EMu_rangeAndCurve_{0}.Scale(scaleFactor*mcPOT_ratio)".format(key))

      # Write raw flux components to output file
      for component in FLUX_COMPONENTS:
        if component == 'dataRate' and playlist == '2p2h': continue # There is no 2p2h data 
        exec("writeHist({0}Hist2D_{1},histOutputDir_playlist)".format(component,key))
      exec("writeHist(migrationMatrix_Enu_{0},histOutputDir_playlist)".format(key))

      # Construct efficiency
      exec("eff_{0} = effNumeratorHist_{0}.Clone('eff_{0}')".format(key))
      exec("eff_{0}.Divide(effNumeratorHist_{0},effDenominatorHist_{0})".format(key))

      # Write efficiency to output file
      exec('writeHist(eff_{0},histOutputDir_playlist)'.format(key))

#############################################################################################################
### Combine LE,ME playlists into single objects for all LE,ME ###############################################
#############################################################################################################
#for PLAYLISTS,LEMEString,isME in zip([PLAYLISTS_LE,PLAYLISTS_ME],['LE','ME'],[False,True]):
if True:
  PLAYLISTS = PLAYLISTS_ME
  LEMEString = "ME"
  isME = True
### Swap out above when running with LE

  firstPlaylist = PLAYLISTS[0]

  # Create all-LE(ME) container hists
  for sigDef in SIGNAL_DEFINITIONS:
    for component in FLUX_COMPONENTS: 
      exec("{0}Hist_{1}_{2} = {0}Hist_{1}_{3}.Clone('{0}Hist_{1}_{2}')".format(component,sigDef,LEMEString,firstPlaylist))
      exec("{0}Hist2D_{1}_{2} = {0}Hist2D_{1}_{3}.Clone('{0}Hist2D_{1}_{2}')".format(component,sigDef,LEMEString,firstPlaylist))
    exec("migrationMatrix_Enu_{0}_{1} = migrationMatrix_Enu_{0}_{2}.Clone('migrationMatrix_Enu_{0}_{1}')".format(sigDef,LEMEString,firstPlaylist))
    #exec("effDenominatorHist2D_addKinematicCuts_{0}_{1} = effDenominatorHist2D_addKinematicCuts_{0}_{2}.Clone('effDenominatorHist2D_addKinematicCuts_{0}_{1}')".format(sigDef,LEMEString,firstPlaylist))
    #exec("effDenominatorHist_addKinematicCuts_{0}_{1} = effDenominatorHist_addKinematicCuts_{0}_{2}.Clone('effDenominatorHist_addKinematicCuts_{0}_{1}')".format(sigDef,LEMEString,firstPlaylist))
    #exec("EMu_rangeOnly_{0}_{1} = EMu_rangeOnly_{0}_{2}.Clone('EMu_rangeOnly_{0}_{1}')".format(sigDef,LEMEString,firstPlaylist))
    #exec("EMu_rangeAndCurve_{0}_{1} = EMu_rangeAndCurve_{0}_{2}.Clone('EMu_rangeAndCurve_{0}_{1}')".format(sigDef,LEMEString,firstPlaylist))

  # Add all playlists to container hists
  for playlist in PLAYLISTS:
    if playlist == firstPlaylist or playlist == '2p2h': continue # We'll come back to 2p2h (LE only), which has to be handled separately, and we've already added the first playlist to our cumulative hists
    print('adding playlist {0} to cumulative {1} hists.'.format(playlist,LEMEString))
    for sigDef in SIGNAL_DEFINITIONS:
      for component in FLUX_COMPONENTS: 
        exec("{0}Hist_{1}_{2}.Add({0}Hist_{1}_{3})".format(component,sigDef,LEMEString,playlist))
        exec("{0}Hist2D_{1}_{2}.Add({0}Hist2D_{1}_{3})".format(component,sigDef,LEMEString,playlist))
      exec("migrationMatrix_Enu_{0}_{1}.Add(migrationMatrix_Enu_{0}_{2})".format(sigDef,LEMEString,playlist))
      #exec("effDenominatorHist2D_addKinematicCuts_{0}_{1}.Add(effDenominatorHist2D_addKinematicCuts_{0}_{2})".format(sigDef,LEMEString,playlist))
      #exec("effDenominatorHist_addKinematicCuts_{0}_{1}.Add(effDenominatorHist_addKinematicCuts_{0}_{2})".format(sigDef,LEMEString,playlist))
      #exec("EMu_rangeOnly_{0}_{1}.Add(EMu_rangeOnly_{0}_{2})".format(sigDef,LEMEString,playlist))
      #exec("EMu_rangeAndCurve_{0}_{1}.Add(EMu_rangeAndCurve_{0}_{2})".format(sigDef,LEMEString,playlist))

  exec("totalDataPOT = totalDataPOT_{0}".format(LEMEString))
  print('totalDataPOT_{0}: {1}'.format(LEMEString,totalDataPOT))

  if not isME: # The below is LE only: Add 2p2h into cumulative LE objects
    scaleFactor_2p2h = totalDataPOT/mcPOT_used_2p2h
    print('scaleFactor_2p2h: ' , scaleFactor_2p2h)
    print('mcPOT_ratio_2p2h: ' , mcPOT_ratio_2p2h)
    
    print('adding 2p2h to cumulative LE hists.')
    for sigDef in SIGNAL_DEFINITIONS:
      # POT scale 2p2h samples
      exec("effNumeratorHist_{0}_2p2h.Scale(scaleFactor_2p2h)".format(sigDef))
      exec("effNumeratorHist2D_{0}_2p2h.Scale(scaleFactor_2p2h)".format(sigDef))
      exec("effDenominatorHist_{0}_2p2h.Scale(scaleFactor_2p2h*mcPOT_ratio_2p2h)".format(sigDef))
      exec("effDenominatorHist2D_{0}_2p2h.Scale(scaleFactor_2p2h*mcPOT_ratio_2p2h)".format(sigDef))
      exec("migrationMatrix_Enu_{0}_2p2h.Scale(scaleFactor_2p2h*mcPOT_ratio_2p2h)".format(sigDef))
      exec("effDenominatorHist2D_addKinematicCuts_{0}_2p2h.Scale(scaleFactor_2p2h*mcPOT_ratio_2p2h)".format(sigDef))
      exec("effDenominatorHist_addKinematicCuts_{0}_2p2h.Scale(scaleFactor_2p2h*mcPOT_ratio_2p2h)".format(sigDef))
      exec("EMu_rangeOnly_{0}_2p2h.Scale(scaleFactor_2p2h*mcPOT_ratio_2p2h)".format(sigDef))
      exec("EMu_rangeAndCurve_{0}_2p2h.Scale(scaleFactor_2p2h*mcPOT_ratio_2p2h)".format(sigDef))
      # Add 2p2h sample to the rest of the LE sample
      for component in ['effNumerator','effDenominator']:
        exec("{0}Hist_{1}_LE.Add({0}Hist_{1}_2p2h)".format(component,sigDef))
        exec("{0}Hist2D_{1}_LE.Add({0}Hist2D_{1}_2p2h)".format(component,sigDef))
      exec("migrationMatrix_Enu_{0}_LE.Add(migrationMatrix_Enu_{0}_2p2h)".format(sigDef))
      exec("effDenominatorHist2D_addKinematicCuts_{0}_LE.Add(effDenominatorHist2D_addKinematicCuts_{0}_2p2h)".format(sigDef))
      exec("effDenominatorHist_addKinematicCuts_{0}_LE.Add(effDenominatorHist_addKinematicCuts_{0}_2p2h)".format(sigDef))
      exec("EMu_rangeOnly_{0}_LE.Add(EMu_rangeOnly_{0}_2p2h)".format(sigDef))
      exec("EMu_rangeAndCurve_{0}_LE.Add(EMu_rangeAndCurve_{0}_2p2h)".format(sigDef))

  for sigDef in SIGNAL_DEFINITIONS:

    # Construct efficiency
    exec("eff_{0}_{1} = effNumeratorHist_{0}_{1}.Clone('eff_{0}_{1}')".format(sigDef,LEMEString))
    exec("eff_{0}_{1}.Divide(effNumeratorHist_{0}_{1},effDenominatorHist_{0}_{1})".format(sigDef,LEMEString))
    #exec("eff_geometric_{0}_{1} = effNumeratorHist_{0}_{1}.Clone('eff_geometric_{0}_{1}')".format(sigDef,LEMEString))
    #exec("eff_geometric_{0}_{1}.Divide(effNumeratorHist_{0}_{1},effDenominatorHist_addKinematicCuts_{0}_{1})".format(sigDef,LEMEString))

    # Add normalization error band to all MnvH1Ds
    # In principle this could be done at the initial histogram-generation stage, but because the band is blank at this stage, it seems tricky to deal with there
    for component in FLUX_COMPONENTS: 
      exec("AddNOMADErrorBand({0}Hist_{1}_{2})".format(component,sigDef,LEMEString))
      exec("AddNOMADErrorBand({0}Hist2D_{1}_{2})".format(component,sigDef,LEMEString))
      print("Adding NOMAD error band to {0}Hist_{1}_{2}".format(component,sigDef,LEMEString))
    exec("AddNOMADErrorBand(eff_{0}_{1})".format(sigDef,LEMEString))
    #exec("AddNOMADErrorBand(eff_geometric_{0}_{1})".format(sigDef,LEMEString))
    exec("AddNOMADErrorBand(migrationMatrix_Enu_{0}_{1})".format(sigDef,LEMEString))
    #exec("AddNOMADErrorBand(effDenominatorHist2D_addKinematicCuts_{0}_{1})".format(sigDef,LEMEString))
    #exec("AddNOMADErrorBand(effDenominatorHist_addKinematicCuts_{0}_{1})".format(sigDef,LEMEString))
    #exec("AddNOMADErrorBand(EMu_rangeOnly_{0}_{1})".format(sigDef,LEMEString))
    #exec("AddNOMADErrorBand(EMu_rangeAndCurve_{0}_{1})".format(sigDef,LEMEString))

    # Write components of flux to output file
    for component in FLUX_COMPONENTS:
      exec('writeHist({0}Hist2D_{1}_{2},histOutputDir_all{2})'.format(component,sigDef,LEMEString))
    # Write efficiency to output file
    exec('writeHist(eff_{0}_{1},histOutputDir_all{1})'.format(sigDef,LEMEString))
    #exec('writeHist(eff_geometric_{0}_{1},histOutputDir_all{1})'.format(sigDef,LEMEString))
    exec("writeHist(migrationMatrix_Enu_{0}_{1},histOutputDir_migrationMatrices)".format(sigDef,LEMEString))
    #exec("writeHist(effDenominatorHist2D_addKinematicCuts_{0}_{1},histOutputDir_all{1})".format(sigDef,LEMEString))
    #exec("writeHist(effDenominatorHist_addKinematicCuts_{0}_{1},histOutputDir_all{1})".format(sigDef,LEMEString))
    #exec("writeHist(EMu_rangeOnly_{0}_{1},histOutputDir_all{1})".format(sigDef,LEMEString))
    #exec("writeHist(EMu_rangeAndCurve_{0}_{1},histOutputDir_all{1})".format(sigDef,LEMEString))
 
#############################################################################################################
### Extract Data/MC Event Rate Ratios #######################################################################
#############################################################################################################

#for LEMEString in ['LE','ME']:
if True: 
  LEMEString = "ME"
### Swap out above when running with LE

  exec("temp1 = dataRateHist_inclusive_{0}.GetVertErrorBand(\"Flux\").GetNHists()".format(LEMEString))
  print("universe in dataRateHist_inclusive_{0}: {1}".format(LEMEString,temp1))
  exec("temp2 = effNumeratorHist_inclusive_{0}.GetVertErrorBand(\"Flux\").GetNHists()".format(LEMEString))
  print("universe in effNumeratorHist_inclusive_{0}: {1}".format(LEMEString,temp2))

  for sigDef in SIGNAL_DEFINITIONS:
  
    # Construct efficiency
    exec("dataMCEventRateRatio_{0}_{1} = dataRateHist_{0}_{1}.Clone('dataMCEventRateRatio_{0}_{1}')".format(sigDef,LEMEString))
    exec("dataMCEventRateRatio_{0}_{1}_areaNorm = dataRateHist_{0}_{1}.Clone('dataMCEventRateRatio_{0}_{1}_areaNorm')".format(sigDef,LEMEString))
    exec("MCEventRate_{0}_{1} = effNumeratorHist_{0}_{1}.Clone('MCEventRate_{0}_{1}')".format(sigDef,LEMEString))
    exec("areaNormScaleFactor = dataRateHist_{0}_{1}.Integral()/MCEventRate_{0}_{1}.Integral()".format(sigDef,LEMEString))
    exec("MCEventRate_{0}_{1}.Scale(areaNormScaleFactor)".format(sigDef,LEMEString)) 
    exec("dataMCEventRateRatio_{0}_{1}.Divide(dataRateHist_{0}_{1},effNumeratorHist_{0}_{1})".format(sigDef,LEMEString))
    exec("dataMCEventRateRatio_{0}_{1}_areaNorm.Divide(dataRateHist_{0}_{1},MCEventRate_{0}_{1})".format(sigDef,LEMEString))

    # Write data/mc event rate ratio hist to output file
    exec('writeHist(dataMCEventRateRatio_{0}_{1},histOutputDir_dataMCRatios)'.format(sigDef,LEMEString))
    exec('writeHist(dataMCEventRateRatio_{0}_{1}_areaNorm,histOutputDir_dataMCRatios)'.format(sigDef,LEMEString))

#############################################################################################################
### Lu-specific Low-nu extraction  ##########################################################################
#############################################################################################################

# Put Lu's lowNu xSection into an MnvH1D for future comparisons
lowNuXSection_Lu_lessUnits = [0.286,0.531,0.527,0.522,0.987,0.985,1.75,1.73,1.73]
lowNuXSection_Lu_list = [i*10**-38 for i in lowNuXSection_Lu_lessUnits] # convert to cm^-2

xSection_lowNu_LE_Lu = PlotUtils.MnvH1D( 'xSection_lowNu_LE_Lu' , 'xSection_lowNu_LE_Lu' , nBins_nuE_LE , array('d',bins_nuE_LE))

for i in range(3,12):
  local_xSectionNum = lowNuXSection_Lu_list[i-3]
  xSection_lowNu_LE_Lu.SetBinContent(i,local_xSectionNum)

xSection_lowNu_LE_Lu.AddMissingErrorBandsAndFillWithCV(effNumeratorHist_highNu_ME)
writeHist(xSection_lowNu_LE_Lu,histOutputDir_xSections)

#############################################################################################################
### Give PPFX Flux Hists Correct Error Bands and Write to Output File #######################################
#############################################################################################################

# PPFX fluxes needd to be given the correct set of error bands to know about
# Using effNumeratorHist_highNu_ME as the reference hist for no particular reason

#flux_PPFX_LE_notBinWidthNormalized.AddMissingErrorBandsAndFillWithCV(effNumeratorHist_highNu_LE) # Using effNumeratorHist_highNu_LE as the reference hist for no particular reason
flux_PPFX_ME_yesNuEConstraint_originalBinning.AddMissingErrorBandsAndFillWithCV(effNumeratorHist_highNu_ME) 
flux_PPFX_ME_yesNuEConstraint_analysisBinning.AddMissingErrorBandsAndFillWithCV(effNumeratorHist_highNu_ME) 
flux_PPFX_ME_originalBinning_notBinWidthNormalized.AddMissingErrorBandsAndFillWithCV(effNumeratorHist_highNu_ME) 
flux_PPFX_ME_analysisBinning_notBinWidthNormalized.AddMissingErrorBandsAndFillWithCV(effNumeratorHist_highNu_ME) # Using effNumeratorHist_highNu_ME as the reference hist for no particular reason

# There are some legacy systematics universes propagated into the fluxes that aren't used any more which can be removed
#flux_PPFX_LE_notBinWidthNormalized.PopVertErrorBand("Flux_BeamFocus")
#flux_PPFX_LE_notBinWidthNormalized.PopVertErrorBand("ppfx1_Total")
flux_PPFX_ME_yesNuEConstraint_originalBinning.PopVertErrorBand("Flux_BeamFocus")
flux_PPFX_ME_yesNuEConstraint_originalBinning.PopVertErrorBand("ppfx1_Total")
flux_PPFX_ME_yesNuEConstraint_analysisBinning.PopVertErrorBand("Flux_BeamFocus")
flux_PPFX_ME_yesNuEConstraint_analysisBinning.PopVertErrorBand("ppfx1_Total")
flux_PPFX_ME_originalBinning_notBinWidthNormalized.PopVertErrorBand("Flux_BeamFocus")
flux_PPFX_ME_originalBinning_notBinWidthNormalized.PopVertErrorBand("ppfx1_Total")
flux_PPFX_ME_analysisBinning_notBinWidthNormalized.PopVertErrorBand("Flux_BeamFocus")
flux_PPFX_ME_analysisBinning_notBinWidthNormalized.PopVertErrorBand("ppfx1_Total")

# Write fluxes to output file
#writeHist(flux_PPFX_LE_notBinWidthNormalized,histOutputDir_flux)
writeHist(flux_PPFX_LE_originalBinning,histOutputDir_flux)
writeHist(flux_PPFX_LE_analysisBinning,histOutputDir_flux)
writeHist(flux_PPFX_ME_yesNuEConstraint_originalBinning,histOutputDir_flux)
writeHist(flux_PPFX_ME_yesNuEConstraint_analysisBinning,histOutputDir_flux)
writeHist(flux_PPFX_ME_originalBinning_notBinWidthNormalized,histOutputDir_flux)
writeHist(flux_PPFX_ME_analysisBinning_notBinWidthNormalized,histOutputDir_flux)

## These need to get sorted out. They're the LE fluxes derived using alternate binnings (I think?)
#writeHist(flux_PPFX_LE_newOriginalBinning,histOutputDir_flux) ## See above RDF 2024-03-26
#writeHist(flux_PPFX_LE,histOutputDir_flux) ## See above RDF 2024-03-26

#############################################################################################################
### Common low-nu extraction infrastructure  ################################################################
#############################################################################################################

## For NOMAD normalization
NOMAD_DataPoint = (0.699*10**-38)*10 # per-GeV integrated cross section in 12-22 GeV bin, multiplied by its bin-width; in units of m^2, not cm^2

## For LE low-nu calculation
binCenterVals_list = [0.5e3,1.5e3,2.5e3,3.5e3,4.5e3,5.5e3,6.5e3,7.5e3,8.5e3,9.5e3,1.1e4,1.3e4,
                      1.5e4,1.7e4,1.9e4,2.1e4,2.5e4,3.1e4,3.7e4,4.5e4,5.5e4,7.0e4,9.0e4,1.1e5]
binCenterVals_LE = PlotUtils.MnvH1D( 'h_binCenterVals_LE' , 'h_binCenterVals_LE' , nBins_nuE_LE , array('d',bins_nuE_LE)) ## Preserve ability to use distinct binning for LE/ME in the future
binCenterVals_ME = PlotUtils.MnvH1D( 'h_binCenterVals_ME' , 'h_binCenterVals_ME' , nBins_nuE_ME, array('d',bins_nuE_ME)) ## Preserve ability to use distinct binning for LE/ME in the future

for i in range(nBins_nuE_LE):
  local_binCenterVal = binCenterVals_list[i]
  binCenterVals_LE.SetBinContent(i+1,local_binCenterVal)
  binCenterVals_ME.SetBinContent(i+1,local_binCenterVal)
binCenterVals_LE.AddMissingErrorBandsAndFillWithCV(effNumeratorHist_highNu_ME)
binCenterVals_ME.AddMissingErrorBandsAndFillWithCV(effNumeratorHist_highNu_ME)

writeHist(binCenterVals_LE,histOutputDir_xSections)
writeHist(binCenterVals_ME,histOutputDir_xSections)

## Create a transformation vector that will zero-out all entries in each lowNu xSection apart from the bins that we'll use that particular cut for
for LEMEString in ["LE","ME"]:
  for nuCut in range(1,5):
    exec("transformVec_{0}_nuCut_{1} = [0]*(nBins_nuE_lowNu_{0}+1)".format(LEMEString,nuCut))
  exec("for binNum in range(0,3): transformVec_{0}_nuCut_1[binNum] = 1.".format(LEMEString)) 
  if LEMEString == "LE":
    #for binNum in range(3,6):   transformVec_LE_nuCut_2[binNum] = 1. 
    #for binNum in range(6,8):   transformVec_LE_nuCut_3[binNum] = 1. 
    #for binNum in range(8,17):  transformVec_LE_nuCut_4[binNum] = 1. 
    #for binNum in range(3,7):   transformVec_LE_nuCut_2[binNum] = 1. 
    #for binNum in range(7,11):   transformVec_LE_nuCut_3[binNum] = 1. 
    #for binNum in range(11,21):  transformVec_LE_nuCut_4[binNum] = 1. 
    for binNum in range(3,7):   transformVec_LE_nuCut_2[binNum] = 1. 
    for binNum in range(7,11):   transformVec_LE_nuCut_3[binNum] = 1. 
    for binNum in range(11,25):  transformVec_LE_nuCut_4[binNum] = 1. 
    #for binNum in range(3,7):   transformVec_LE_nuCut_2[binNum] = 1. 
    #for binNum in range(7,12):   transformVec_LE_nuCut_3[binNum] = 1. 
    #for binNum in range(12,27):  transformVec_LE_nuCut_4[binNum] = 1. 
  else:
    #for binNum in range(3,7):   transformVec_ME_nuCut_2[binNum] = 1. 
    #for binNum in range(7,11):  transformVec_ME_nuCut_3[binNum] = 1. 
    #for binNum in range(11,21): transformVec_ME_nuCut_4[binNum] = 1. 
    for binNum in range(3,7):   transformVec_ME_nuCut_2[binNum] = 1. 
    for binNum in range(7,11):   transformVec_ME_nuCut_3[binNum] = 1. 
    for binNum in range(11,25):  transformVec_ME_nuCut_4[binNum] = 1. 
  for nuCut in range(1,5):
    # Create TH1D out of transformation vector defined above
    exec("transformVecHist_{0}_nuCut_{1} = ROOT.TH1D( 'transformVecHist_{0}_nuCut_{1}' , 'transformVecHist_{0}_nuCut_{1}' , nBins_nuE_lowNu_{0} , array('d',bins_nuE_lowNu_{0}))".format(LEMEString,nuCut))
    exec("nBins = nBins_nuE_lowNu_{0}".format(LEMEString))
    for nBin in range(nBins):
      exec("transformVecHist_{0}_nuCut_{1}.SetBinContent(nBin+1,transformVec_{0}_nuCut_{1}[nBin])".format(LEMEString,nuCut))
      exec("transformVecHist_{0}_nuCut_{1}.SetBinError(nBin+1,0.0)".format(LEMEString,nuCut))

#############################################################################################################
### Low-nu Extraction in LE and ME  #########################################################################
#############################################################################################################

#for LEMEString,isME in zip(['LE','ME'],[False,True]):
if True: 
  LEMEString = "ME"
  isME = True
### Swap out above when running with LE

  # Extract data rate and efficiency for each nu cut separately for LE and ME
  for nuCut in range(1,5):
    for sigDef in SIGNAL_DEFINITIONS:
      for component in FLUX_COMPONENTS:
        exec('{0}Hist_{3}_{1}_nuCut_{2} = {0}Hist2D_inclusive_{1}.ProjectionX("",1,{2})'.format(component,LEMEString,nuCut,sigDef))
        exec("{0}Hist_{3}_{1}_nuCut_{2}.SetName('{0}Hist_{3}_{1}_nuCut_{2}')".format(component,LEMEString,nuCut,sigDef))
      exec("eff_{2}_{0}_nuCut_{1} = effNumeratorHist_{2}_{0}_nuCut_{1}.Clone('eff_{2}_{0}_nuCut_{1}')".format(LEMEString,nuCut,sigDef))
      exec('eff_{2}_{0}_nuCut_{1}.Divide(effNumeratorHist_{2}_{0}_nuCut_{1},effDenominatorHist_{2}_{0}_nuCut_{1})'.format(LEMEString,nuCut,sigDef))
      # Write data rate and efficiency components to output file for each nu cut
      exec('writeHist(dataRateHist_{2}_{0}_nuCut_{1},histOutputDir_all{0})'.format(LEMEString,nuCut,sigDef))
      exec('writeHist(effNumeratorHist_{2}_{0}_nuCut_{1},histOutputDir_all{0})'.format(LEMEString,nuCut,sigDef))
      exec('writeHist(effDenominatorHist_{2}_{0}_nuCut_{1},histOutputDir_all{0})'.format(LEMEString,nuCut,sigDef))
      exec('writeHist(eff_{2}_{0}_nuCut_{1},histOutputDir_all{0})'.format(LEMEString,nuCut,sigDef))

  # Create GENIE xsection object, which will be used in the denominator of the low-nu xsection extraction
  exec("xSection_inclusive_{0}_GENIE_2D_local = effDenominatorHist2D_inclusive_{0}.Clone()".format(LEMEString))

  for nuCut in range(1,5):
    exec('xSection_lowNu_{0}_GENIE_nuCut_{1} = xSection_inclusive_{0}_GENIE_2D_local.ProjectionX("",1,{1})'.format(LEMEString,nuCut))

    # Duplicate, then write out an unaltered ("RAW") copy of the GENIE event rate
    exec('xSection_lowNu_{0}_GENIE_RAW_nuCut_{1} = xSection_lowNu_{0}_GENIE_nuCut_{1}.Clone("xSection_lowNu_{0}_GENIE_RAW_nuCut_{1}")'.format(LEMEString,nuCut))
    exec('writeHist(xSection_lowNu_{0}_GENIE_RAW_nuCut_{1},histOutputDir_xSections)'.format(LEMEString,nuCut))
 
    exec("xSection_lowNu_{0}_GENIE_nuCut_{1}.SetName('xSection_lowNu_{0}_GENIE_nuCut_{1}')".format(LEMEString,nuCut))

    ## If ME, we need to use the PPFX flux that has the ME-lowNu-specific binning
    if LEMEString == "ME":
      exec('xSection_lowNu_{0}_GENIE_nuCut_{1}.Divide(xSection_lowNu_{0}_GENIE_nuCut_{1},flux_PPFX_{0}_analysisBinning_notBinWidthNormalized)'.format(LEMEString,nuCut))
    else:
      exec('xSection_lowNu_{0}_GENIE_nuCut_{1}.Divide(xSection_lowNu_{0}_GENIE_nuCut_{1},flux_PPFX_{0}_notBinWidthNormalized)'.format(LEMEString,nuCut))
    exec('xSection_lowNu_{0}_GENIE_nuCut_{1}.Scale(1./(nTargets*totalDataPOT_{0}))'.format(LEMEString,nuCut)) 
    exec('xSection_lowNu_{0}_GENIE_nuCut_{1}.Scale(10**4)'.format(LEMEString,nuCut)) # Change units of flux from m^{-2} to cm^{-2}
  
    # Calculate the flux expliclity, though we don't actually use this to calculate the inclusive cross section
    exec("flux_lowNu_{0}_nuCut_{1} = dataRateHist_lowNu_{0}_nuCut_{1}.Clone('flux_lowNu_{0}_nuCut_{1}')".format(LEMEString,nuCut))
    exec('flux_lowNu_{0}_nuCut_{1}.Divide(flux_lowNu_{0}_nuCut_{1},eff_lowNu_{0}_nuCut_{1})'.format(LEMEString,nuCut))
    exec('flux_lowNu_{0}_nuCut_{1}.Divide(flux_lowNu_{0}_nuCut_{1},xSection_lowNu_{0}_GENIE_nuCut_{1})'.format(LEMEString,nuCut))
    exec('flux_lowNu_{0}_nuCut_{1}.Scale(1./(nTargets*totalDataPOT_{0}))'.format(LEMEString,nuCut))
    exec('flux_lowNu_{0}_nuCut_{1}.Scale(10**4)'.format(LEMEString,nuCut)) # Change units of flux from cm^{-2} to m^{-2}, via the xSection, hence the inverted factor
    # Write flux out
    exec('writeHist(flux_lowNu_{0}_nuCut_{1},histOutputDir_flux)'.format(LEMEString,nuCut))
 
    # Extract inclusive cross section for each nu cut, which we do with a shortcut that side-steps calculating the flux explicitly
    exec("xSection_inclusive_{0}_nuCut_{1} = dataRateHist_inclusive_{0}.Clone('xSection_inclusive_{0}_nuCut_{1}')".format(LEMEString,nuCut))
    exec('xSection_inclusive_{0}_nuCut_{1}.Divide(dataRateHist_inclusive_{0},dataRateHist_lowNu_{0}_nuCut_{1})'.format(LEMEString,nuCut))
    exec('xSection_inclusive_{0}_nuCut_{1}.Multiply(xSection_inclusive_{0}_nuCut_{1},eff_lowNu_{0}_nuCut_{1})'.format(LEMEString,nuCut))
    exec('xSection_inclusive_{0}_nuCut_{1}.Divide(xSection_inclusive_{0}_nuCut_{1},eff_inclusive_{0})'.format(LEMEString,nuCut))
  
    exec('xSection_inclusive_{0}_nuCut_{1}.Multiply(xSection_inclusive_{0}_nuCut_{1},xSection_lowNu_{0}_GENIE_nuCut_{1})'.format(LEMEString,nuCut))
  
    exec("xSectionPerE_inclusive_{0}_nuCut_{1} = xSection_inclusive_{0}_nuCut_{1}.Clone('xSectionPerE_inclusive_{0}_nuCut_{1}')".format(LEMEString,nuCut))
    ## DEBUG
    exec("testA = xSectionPerE_inclusive_{0}_nuCut_{1}.Clone('testA')".format(LEMEString,nuCut))
    writeHist(testA,histOutputDir_xSections)
    exec("xSectionPerE_inclusive_{0}_nuCut_{1}.Divide(xSectionPerE_inclusive_{0}_nuCut_{1},binCenterVals_{0})".format(LEMEString,nuCut))
    ## DEBUG
    exec("testB = xSectionPerE_inclusive_{0}_nuCut_{1}.Clone('testB')".format(LEMEString,nuCut))
    writeHist(testB,histOutputDir_xSections)
  
    # Get list of systematic universes to loop through 
    exec('errorBandNames = xSectionPerE_inclusive_{0}_nuCut_{1}.GetVertErrorBandNames()'.format(LEMEString,nuCut))
  
    # Make MnvH1D to store scale factors, by cloning one of the xSection MnvH1Ds and dividing it by itself, so this new MnvH1D has the correct set of SUs, each filled with 1 in all its bins
    exec('externalNormalizationScaleFactors_{0}_nuCut_{1} = xSectionPerE_inclusive_{0}_nuCut_{1}.Clone()'.format(LEMEString,nuCut))
    exec('externalNormalizationScaleFactors_{0}_nuCut_{1}.Divide(externalNormalizationScaleFactors_{0}_nuCut_{1},xSectionPerE_inclusive_{0}_nuCut_{1})'.format(LEMEString,nuCut))
  
    # Construct CV scale factor and add into MnvH1D of scale factors 
    #integralLowerBound = 9 if LEMEString == "LE" else 12
    #integralUpperBound = 11 if LEMEString == "LE" else 15
    #integralLowerBound = 12 if LEMEString == "LE" else 12
    #integralUpperBound = 15 if LEMEString == "LE" else 15
    #integralLowerBound = 13 if LEMEString == "LE" else 12
    #integralUpperBound = 21 if LEMEString == "LE" else 15
    integralLowerBound = 12
    integralUpperBound = 16
    exec('CV_xSectionIntegral_{0} = xSectionPerE_inclusive_{0}_nuCut_{1}.Integral(integralLowerBound,integralUpperBound,"width")'.format(LEMEString,nuCut))
    exec('CV_scaleFactor = NOMAD_DataPoint/CV_xSectionIntegral_{0}'.format(LEMEString))
    exec('externalNormalizationScaleFactors_{0}_nuCut_{1}.Scale({2},\"\",False)'.format(LEMEString,nuCut,CV_scaleFactor)) #False here tells Scale to not scale the universes, rather only the CV
 
    # Loop over SUs
    for errorBandName in errorBandNames:
      exec('errorBand = xSectionPerE_inclusive_{0}_nuCut_{1}.GetVertErrorBand(\"{2}\")'.format(LEMEString,nuCut,errorBandName))
      for iHist in range(errorBand.GetNHists()):
        systematicUniverse = errorBand.GetHist(iHist)
        xSectionIntegral = systematicUniverse.Integral(integralLowerBound,integralUpperBound,"width")
        scaleFactor = NOMAD_DataPoint/xSectionIntegral
        exec('externalNormalizationScaleFactors_{0}_nuCut_{1}.GetVertErrorBand(\"{2}\").GetHist({3}).Scale({4})'.format(LEMEString,nuCut,errorBandName,iHist,scaleFactor))
    
    for sigDef in ['xSection_lowNu_{0}_GENIE'.format(LEMEString),'xSection_inclusive_{0}'.format(LEMEString),'xSectionPerE_inclusive_{0}'.format(LEMEString)]:
      # Write pre-normalization lowNu and inclusive xSections to output file
      exec('writeHist({0}_nuCut_{1},histOutputDir_xSections)'.format(sigDef,nuCut))
      # Normalize xSection MnvH1D
      exec('{0}_nuCut_{1}.Multiply({0}_nuCut_{1},externalNormalizationScaleFactors_{2}_nuCut_{1})'.format(sigDef,nuCut,LEMEString))
      # Rename to clarify that this is the normalized distribution
      exec("{0}_nuCut_{1}.SetName('{0}_nuCut_{1}_normalized')".format(sigDef,nuCut))
      # Write normalized lowNu and inclusive xSections to output file
      exec('writeHist({0}_nuCut_{1},histOutputDir_xSections)'.format(sigDef,nuCut))

  ## Scale the low-nu fluxes
  for nuCut in range(1,5):
    exec("flux_lowNu_normalized_{0}_nuCut_{1} = flux_lowNu_{0}_nuCut_{1}.Clone('flux_lowNu_normalized_{0}_nuCut_{1}')".format(LEMEString,nuCut))
    #exec("flux_lowNu_normalized_{0}_nuCut_{1}.Multiply(flux_lowNu_normalized_{0}_nuCut_{1},externalNormalizationScaleFactors_{0}_nuCut_{1})".format(LEMEString,nuCut))
    exec("flux_lowNu_normalized_{0}_nuCut_{1}.Divide(flux_lowNu_normalized_{0}_nuCut_{1},externalNormalizationScaleFactors_{0}_nuCut_{1})".format(LEMEString,nuCut))
    exec("writeHist(flux_lowNu_normalized_{0}_nuCut_{1},histOutputDir_flux)".format(LEMEString,nuCut))

  ## Now take a step back to the low-nu cross section and assemble the final version that we'll use
  for nuCut in range(1,5):
    # Create local copy of the lowNu xSection corresponding to each nu cut
    exec("local_xSection_lowNu_{0}_nuCut_{1} = xSection_lowNu_{0}_GENIE_nuCut_{1}.Clone('local_xSection_lowNu_{0}_nuCut_{1}')".format(LEMEString,nuCut))
    # Scale local copy by the transformation vector
    exec("local_xSection_lowNu_{0}_nuCut_{1}.MultiplySingle(local_xSection_lowNu_{0}_nuCut_{1},transformVecHist_{0}_nuCut_{1})".format(LEMEString,nuCut))
    # Write this local hist to output
    exec("writeHist(local_xSection_lowNu_{0}_nuCut_{1},histOutputDir_xSections)".format(LEMEString,nuCut))
 
  # Add the above hists together and this is the cumulative lowNu XS
  exec("xSection_lowNu_{0} = local_xSection_lowNu_{0}_nuCut_1.Clone('xSection_lowNu_{0}')".format(LEMEString))
  for nuCut in range(2,5):
    exec("xSection_lowNu_{0}.Add(local_xSection_lowNu_{0}_nuCut_{1})".format(LEMEString,nuCut))
  
  # Give the LE low-nu xSection the actual NOMAD normalization uncertainty
  exec("AddNOMADErrorBand(xSection_lowNu_{0},True)".format(LEMEString))
  
  # Write this xSection to output file
  exec("writeHist(xSection_lowNu_{0},histOutputDir_xSections)".format(LEMEString))

  ## Propagate this lowNu XS to LE lowNu flux and inclusive XS
 
  # Extract LE lowNu flux using the lowNu XS constructed above
  exec("flux_lowNu_{0} = dataRateHist_lowNu_{0}.Clone('flux_lowNu_{0}')".format(LEMEString))
  exec("flux_lowNu_{0}.Divide(flux_lowNu_{0},eff_lowNu_{0})".format(LEMEString))
  exec("flux_lowNu_{0}.Divide(flux_lowNu_{0},xSection_lowNu_{0})".format(LEMEString))
  exec("flux_lowNu_{0}.Scale(1./(nTargets*totalDataPOT_{0}))".format(LEMEString))
  exec("flux_lowNu_{0}.Scale(10**4)".format(LEMEString)) # Change units of flux from cm^{-2} to m^{-2}, via the xSection, hence the inverted factor
  
  exec("writeHist(flux_lowNu_{0},histOutputDir_flux)".format(LEMEString))
  
  # Extract LE xSections
  exec("xSection_inclusive_{0} = dataRateHist_inclusive_{0}.Clone('xSection_inclusive_{0}')".format(LEMEString))
  exec("xSection_inclusive_{0}.Divide(dataRateHist_inclusive_{0},dataRateHist_lowNu_{0})".format(LEMEString))
  exec("xSection_inclusive_{0}.Multiply(xSection_inclusive_{0},eff_lowNu_{0})".format(LEMEString))
  exec("xSection_inclusive_{0}.Divide(xSection_inclusive_{0},eff_inclusive_{0})".format(LEMEString))
  
  # Propagate nominal LE/ME low nu cross section to inclusive measurements
  exec("xSection_inclusive_{0}.Multiply(xSection_inclusive_{0},xSection_lowNu_{0})".format(LEMEString))
  
  exec("xSection_highNu_{0} = xSection_inclusive_{0}.Clone('xSection_highNu_{0}')".format(LEMEString))
  exec("xSection_highNu_{0}.Add(xSection_lowNu_{0},-1.)".format(LEMEString))
  
  exec("writeHist(xSection_inclusive_{0},histOutputDir_xSections)".format(LEMEString))
  exec("writeHist(xSection_highNu_{0},histOutputDir_xSections)".format(LEMEString))

os.sys.exit(1)

#############################################################################################################
### Extract Lu-style Deliverables ###########################################################################
#############################################################################################################

# Also extract LE lowNu flux using Lu's lowNu XS
flux_lowNu_LE_Lu = dataRateHist_lowNu_LE.Clone('flux_lowNu_LE_Lu')
flux_lowNu_LE_Lu.Divide(flux_lowNu_LE_Lu,eff_lowNu_LE)
flux_lowNu_LE_Lu.Divide(flux_lowNu_LE_Lu,xSection_lowNu_LE_Lu)
flux_lowNu_LE_Lu.Scale(1./(nTargets*totalDataPOT_LE))
flux_lowNu_LE_Lu.Scale(10**4) # Change units of flux from cm^{-2} to m^{-2}, via the xSection, hence the inverted factor

# Propagate Lu LE low nu cross section to inclusive measurement
xSection_inclusive_LE_Lu = xSection_inclusive_LE.Clone('xSection_inclusive_LE_Lu')
xSection_inclusive_LE_Lu.Multiply(xSection_inclusive_LE_Lu,xSection_lowNu_LE_Lu)

writeHist(flux_lowNu_LE_Lu,histOutputDir_flux)
writeHist(xSection_inclusive_LE_Lu,histOutputDir_xSections)

#############################################################################################################
### Extract ME high-nu Flux #################################################################################
#############################################################################################################

flux_highNu_ME = dataRateHist_highNu_ME.Clone('flux_highNu_ME')
flux_highNu_ME.Divide(flux_highNu_ME,eff_highNu_ME)

## Temp until the LE hists are updated with the full suite of error bands
xSection_highNu_LE.AddMissingErrorBandsAndFillWithCV(flux_highNu_ME)
##

flux_highNu_ME.Divide(flux_highNu_ME,xSection_highNu_LE)
flux_highNu_ME.Scale(1./(nTargets*totalDataPOT_ME))
writeHist(flux_highNu_ME,histOutputDir_flux)

#############################################################################################################
### Close output file; other business #######################################################################
#############################################################################################################
histOutputFile.Close()

