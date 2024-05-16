

#### Neutrino Energy Binning
#######################################
#bins_nuE_LE = [0,1,2,3,4,5,7,9,12,15,18,22,36,50,75,100,120]
#bins_nuE_LE = [0,1,2,3,4,5,6,7,8,9,10.5,12,13.5,15,18,22,36,50,75,100,120]
#bins_nuE_LE = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,36,50,75,100,120]
bins_nuE_LE = [0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,28,34,40,50,60,80,100,120]
nBins_nuE_LE = len(bins_nuE_LE)-1

## All signal definitions in LE get the default LE binning
for sigDef in ["lowNu","highNu","inclusive"]:
  exec("bins_nuE_{0}_LE = bins_nuE_LE".format(sigDef))
  exec("nBins_nuE_{0}_LE = nBins_nuE_LE".format(sigDef))

## The ME high-nu analysis uses the LE binning
bins_nuE_highNu_ME = bins_nuE_LE
nBins_nuE_highNu_ME = nBins_nuE_LE

## Signal definitions in ME used for the ME low-nu analysis get different binning
#bins_nuE_ME = [0,1,2,3,4,5,6,7,8,9,10.5,12,13.5,15,18,22,36,50,75,100,120]
bins_nuE_ME = [0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,28,34,40,50,60,80,100,120]
binCenters_nuE_ME = [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,11,13,15,17,19,21,25,31,37,45,55,70,90,110]
#bins_nuE_ME = [0,1.5,3,4.5,6,7.5,9,10.5,12,13.5,15,20,25,30,50,75,100,120] # Amit's binning
nBins_nuE_ME = len(bins_nuE_ME)-1

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

