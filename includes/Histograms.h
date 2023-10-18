//==============================================================================
// Container class that holds a lot of histograms for a specific variable
// Also knows how to initialize those histograms
//==============================================================================
#ifndef Histograms_h
#define Histograms_h

#include "Binning.h"  // MakeUniformBinArray
#include "CVUniverse.h"
#include "Constants.h"  // typedefs MH1D, CVHW, CVH2DW
#include "StackedHistogram.h"
#include "TArrayD.h"
#include "TFile.h"
#include "TruthMatching.h"
#include "utilities.h"  // uniq

class Histograms {
 public:
  //==========================================================================
  // Constructors
  //==========================================================================
  Histograms();

  Histograms(const std::string label, const std::string xlabel, const int nbins,
             const double xmin, const double xmax);

  Histograms(const std::string label, const std::string xlabel,
             const TArrayD& bins_array);

  //  Histograms(const Histograms&);

  //==========================================================================
  // Data Members
  //==========================================================================
  // Basic Data Members
  std::string m_label;
  std::string m_xlabel;
  TArrayD m_bins_array;

  // Cross Section Pipeline
  CVHW m_bg;
  CVHW m_bg_hiW;
  CVHW m_bg_loW;
  CVHW m_bg_midW;
  MH1D* m_bg_subbed_data;
  MH1D* m_cross_section;
  CVHW m_effden;
  MH1D* m_efficiency;
  CVHW m_effnum;
  CVH2DW m_migration;
  MH1D* m_selection_data;
  CVHW m_selection_mc;
  MH1D* m_tuned_bg;
  MH1D* m_unfolded;
  MH1D* m_wsideband_data;
  MH1D* m_wsidebandfit_data;
  CVHW m_wsidebandfit_hiW;
  CVHW m_wsidebandfit_loW;
  CVHW m_wsidebandfit_midW;
  CVHW m_wsidebandfit_sig;

  // Stacked Histograms
  StackedHistogram<ChannelType> m_stacked_channel;
  StackedHistogram<FSParticleType> m_stacked_fspart;
  StackedHistogram<MesonBackgroundType> m_stacked_mesonbg;
  StackedHistogram<SignalBackgroundType> m_stacked_sigbg;
  StackedHistogram<WBackgroundType> m_stacked_wbg;
  StackedHistogram<WSidebandType> m_stacked_wsideband;

  //==========================================================================
  // Functions
  //==========================================================================
  int NBins() const { return m_bins_array.GetSize() - 1; }
  double XMin() const { return m_bins_array[0]; }
  double XMax() const { return m_bins_array[NBins()]; }
  void PrintBinning() const {
    for (int i = 0; i <= NBins(); ++i) std::cout << m_bins_array[i] << " ";
    std::cout << "\n";
  }

  // Histogram Initialization
  template <typename T>
  void InitializeAllHists(T systematic_univs, T systematic_univs_truth);
  template <typename T>
  void InitializeSelectionHists(T systematic_univs, T systematic_univs_truth);
  template <typename T>
  void InitializeSidebandHists(T systematic_univs);
  template <typename T>
  void InitializeMigrationHists(T systematic_univs);
  void InitializeDataHists();
  void InitializeStackedHists();

  // Stack Map Access
  std::map<SignalBackgroundType, MH1D*> GetStackMap(
      SignalBackgroundType type) const;
  std::map<WBackgroundType, MH1D*> GetStackMap(WBackgroundType type) const;
  std::map<MesonBackgroundType, MH1D*> GetStackMap(
      MesonBackgroundType type) const;
  std::map<FSParticleType, MH1D*> GetStackMap(FSParticleType type) const;
  std::map<ChannelType, MH1D*> GetStackMap(ChannelType type) const;
  std::map<WSidebandType, MH1D*> GetStackMap(WSidebandType type) const;

  // Load MC hists from file
  void LoadDataHistsFromFile(TFile& fin);
  void LoadMCHistsFromFile(TFile& fin, UniverseMap& error_bands, bool is_true);
  CVHW LoadHWFromFile(TFile& fin, UniverseMap& error_bands, std::string name);
  CVH2DW LoadH2DWFromFile(TFile& fin, UniverseMap& error_bands,
                          std::string name);
};

// Template member functions need to be available in the header.
#include "Histograms.cxx"

#endif  // Histograms_h
