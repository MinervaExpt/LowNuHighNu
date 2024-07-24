

#### Neutrino Energy Binning
#######################################
#bins_nuE_MeV = [0.,1.e3,2.e3,3.e3,4.e3,5.e3,6.e3,7.e3,8.e3,9.e3,1.e4,1.2e4,1.4e4,1.6e4,
#                1.8e4,2.0e4,2.2e4,2.8e4,3.4e4,4.0e4,5.0e4,6.0e4,8.0e4,1.0e5,1.2e5]
#bins_nuE_GeV = [0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,28,34,40,50,60,80,100,120]
bins_nuE_MeV = [0.,1.e3,2.e3,3.e3,4.e3,5.e3,6.e3,7.e3,8.e3,9.e3,1.e4,1.2e4,1.4e4,1.6e4,
                1.8e4,2.0e4,3.0e4,4.0e4,5.0e4,6.0e4,8.0e4,1.0e5,1.2e5]
bins_nuE_GeV = [0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,30,40,50,60,80,100,120]

bins_nuE = bins_nuE_GeV
nBins_nuE = len(bins_nuE)-1

bins_nuE_LE = bins_nuE
nBins_nuE_LE = nBins_nuE

## All signal definitions in LE get the default LE binning
for sigDef in ["lowNu","highNu","inclusive"]:
  exec("bins_nuE_{0}_LE = bins_nuE".format(sigDef))
  exec("nBins_nuE_{0}_LE = nBins_nuE".format(sigDef))

## The ME high-nu analysis uses the LE binning
bins_nuE_highNu_ME = bins_nuE
nBins_nuE_highNu_ME = nBins_nuE

## Signal definitions in ME used for the ME low-nu analysis may utilize different binning
#bins_nuE_ME = [0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,28,34,40,50,60,80,100,120]
#nBins_nuE_ME = len(bins_nuE_ME)-1
bins_nuE_ME = bins_nuE
nBins_nuE_ME = nBins_nuE

for sigDef in ["lowNu","inclusive"]:
  exec("bins_nuE_{0}_ME = bins_nuE_ME".format(sigDef))
  exec("nBins_nuE_{0}_ME = nBins_nuE_ME".format(sigDef))

#### Alternate Neutrino Energy Binnings
#######################################
bins_nuE_Amit = [0,1.5,3,4.5,6,7.5,9,10.5,12,13.5,15,20,25,30,50,75,100,120] # Amit's binning
bins_nuE_uniform = range(51) # Another alternate binning


#### Nu (Recoil Energy) Binning
#######################################
bins_nu = [0,0.3,0.5,1,2,5]
nBins_nu = len(bins_nu)-1

### PPFX flux full binning (borrowed fr
#om PlotUtils/scripts/fluxRearranger.py
#######################################

binMap = [#[min,max,size] (all in units of GeV)
          [0,1.5,0.5],
          [1.5,10,0.1],
          [10,20,0.5],
          [20,30,1],
          [30,60,5],
          [60,100,10] # The flux generation only goes up to 100 GeV
]
bins_nuE_ME_fineBinning = []
for minBin,maxBin,binSize in binMap:
  nBins = (maxBin-minBin)/binSize
  for iBin in range(int(nBins)):
    binBoundary = minBin + iBin*binSize
    bins_nuE_ME_fineBinning.append(binBoundary)
bins_nuE_ME_fineBinning.append(100.) # The above prescription doesn't add the last bin
nBins_nuE_ME_fineBinning = len(bins_nuE_ME_fineBinning)-1

binCenters_nuE_ME_fineBinning = []
for i in range(nBins_nuE_ME_fineBinning):
  binCenter = (bins_nuE_ME_fineBinning[i]+bins_nuE_ME_fineBinning[i+1])/2.
  binCenters_nuE_ME_fineBinning.append(binCenter)

