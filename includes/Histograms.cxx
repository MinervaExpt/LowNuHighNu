#ifndef Histograms_cxx
#define Histograms_cxx

#include "Histograms.h"

#include <algorithm>

// CTOR -- default
Histograms::Histograms()
    : m_label(),
      m_xlabel(),
      m_bins_array(0),
      m_bg(),
      m_bg_subbed_data(),
      m_cross_section(),
      m_effden(),
      m_efficiency(),
      m_effnum(),
      m_migration(),
      m_selection_data(),
      m_selection_mc(),
      m_tuned_bg(),
      m_unfolded(),
      m_stacked_channel(),
      m_stacked_fspart(),
      m_stacked_mesonbg(),
      m_stacked_sigbg(),
      m_stacked_wbg() {}

// CTOR -- uniform binning
Histograms::Histograms(const std::string label, const std::string xlabel,
                       const int nbins, const double xmin, const double xmax)
    : m_label(label),
      m_xlabel(xlabel),
      m_bins_array(MakeUniformBinArray(nbins, xmin, xmax)),
      m_bg(),
      m_bg_subbed_data(),
      m_cross_section(),
      m_effden(),
      m_efficiency(),
      m_effnum(),
      m_migration(),
      m_selection_data(),
      m_selection_mc(),
      m_tuned_bg(),
      m_unfolded(),
      m_stacked_channel(),
      m_stacked_fspart(),
      m_stacked_mesonbg(),
      m_stacked_sigbg(),
      m_stacked_wbg() {}

// CTOR -- variable binning
Histograms::Histograms(const std::string label, const std::string xlabel,
                       const TArrayD& bins_array)
    : m_label(label),
      m_xlabel(xlabel),
      m_bins_array(GetSortedArray(bins_array)),
      m_bg(),
      m_bg_subbed_data(),
      m_cross_section(),
      m_effden(),
      m_efficiency(),
      m_effnum(),
      m_migration(),
      m_selection_data(),
      m_selection_mc(),
      m_tuned_bg(),
      m_unfolded(),
      m_stacked_channel(),
      m_stacked_fspart(),
      m_stacked_mesonbg(),
      m_stacked_sigbg(),
      m_stacked_wbg() {}

// Load hists from file
void Histograms::LoadMCHistsFromFile(TFile& fin, UniverseMap& error_bands,
                                     bool is_true) {
  m_selection_mc = LoadHWFromFile(fin, error_bands, "selection_mc");

  m_bg = LoadHWFromFile(fin, error_bands, "bg");
  m_tuned_bg =
      (PlotUtils::MnvH1D*)fin.Get(Form("tuned_bg_%s", m_label.c_str()));
  m_effnum = LoadHWFromFile(fin, error_bands, "effnum");
  m_effden = LoadHWFromFile(fin, error_bands, "effden");

  if (!is_true) {
    m_migration = LoadH2DWFromFile(fin, error_bands, "migration");
    assert(m_migration.hist);
  }
}

CVHW Histograms::LoadHWFromFile(TFile& fin, UniverseMap& error_bands,
                                std::string name) {
  const bool do_erase_bands = false;
  PlotUtils::MnvH1D* hist =
      (PlotUtils::MnvH1D*)fin.Get(Form("%s_%s", name.c_str(), m_label.c_str()));

  if (hist == 0) {
    std::cout << name << " hist doesn't exist. skipping...\n";
    return CVHW();
  } else {
    TArrayD bins_array = *(hist->GetXaxis()->GetXbins());
    // Source histo has uniform binning
    if (bins_array.GetSize() == 0) {
      bins_array.Reset();
      bins_array = MakeUniformBinArray(hist->GetXaxis()->GetNbins(),
                                       hist->GetXaxis()->GetXmin(),
                                       hist->GetXaxis()->GetXmax());
      hist = dynamic_cast<PlotUtils::MnvH1D*>(hist->Rebin(
          bins_array.GetSize() - 1, hist->GetName(), bins_array.GetArray()));
    }

    // Compare source and destination binnings
    for (int i = 0; i < NBins(); ++i) {
      if (m_bins_array[i] != bins_array[i]) {
        std::cout << "WARNING! Binning mismatch for " << m_label << "\n";
        std::cout << "Aligning output binning to match source binning\n";
        m_bins_array = bins_array;
        break;
      }
    }

    return CVHW(hist, error_bands, do_erase_bands);
  }
}

CVH2DW Histograms::LoadH2DWFromFile(TFile& fin, UniverseMap& error_bands,
                                    std::string name) {
  const bool do_erase_bands = false;
  PlotUtils::MnvH2D* hist =
      (PlotUtils::MnvH2D*)fin.Get(Form("%s_%s", name.c_str(), m_label.c_str()));
  assert(hist);
  return CVH2DW(hist, error_bands, do_erase_bands);
}

void Histograms::LoadDataHistsFromFile(TFile& fin) {
  m_selection_data =
      (PlotUtils::MnvH1D*)fin.Get(Form("selection_data_%s", m_label.c_str()));
  m_bg_subbed_data =
      (PlotUtils::MnvH1D*)fin.Get(Form("bg_subbed_data_%s", m_label.c_str()));
  m_unfolded =
      (PlotUtils::MnvH1D*)fin.Get(Form("unfolded_%s", m_label.c_str()));
  m_cross_section =
      (PlotUtils::MnvH1D*)fin.Get(Form("cross_section_%s", m_label.c_str()));
}

// Initialize Hists
template <typename T>
void Histograms::InitializeAllHists(T systematic_univs,
                                    T systematic_univs_truth) {
  // Event Section Analysis
  InitializeSelectionHists(systematic_univs, systematic_univs_truth);

  // Migration Matrix
  InitializeMigrationHists(systematic_univs);

  // Data
  InitializeDataHists();

  // Event Selection Stacked
  InitializeStackedHists();
}

template <typename T>
void Histograms::InitializeSelectionHists(T systematic_univs,
                                          T systematic_univs_truth) {
  const Double_t* bins = m_bins_array.GetArray();
  const char* label = m_label.c_str();

  MH1D* selection_mc =
      new MH1D(Form("selection_mc_%s", label), label, NBins(), bins);

  MH1D* bg = new MH1D(Form("bg_%s", label), label, NBins(), bins);

  MH1D* effnum = new MH1D(Form("effnum_%s", label), label, NBins(), bins);

  MH1D* effden = new MH1D(Form("effden_%s", label), label, NBins(), bins);

  const bool clear_bands = true;
  m_selection_mc = CVHW(selection_mc, systematic_univs, clear_bands);
  m_bg = CVHW(bg, systematic_univs, clear_bands);
  m_effnum = CVHW(effnum, systematic_univs, clear_bands);
  m_effden = CVHW(effden, systematic_univs_truth, clear_bands);

  delete selection_mc;
  delete bg;
  delete effnum;
  delete effden;
}

void Histograms::InitializeStackedHists() {
  // Event Selection Stacked
  m_stacked_sigbg = StackedHistogram<SignalBackgroundType>(
      m_label, m_xlabel, m_bins_array, int(kNSignalBackgroundTypes));

  m_stacked_wbg = StackedHistogram<WBackgroundType>(
      m_label, m_xlabel, m_bins_array, int(kNWBackgroundTypes), 4);

  m_stacked_mesonbg = StackedHistogram<MesonBackgroundType>(
      m_label, m_xlabel, m_bins_array, int(kNMesonBackgroundTypes));

  m_stacked_fspart = StackedHistogram<FSParticleType>(
      m_label, m_xlabel, m_bins_array, int(kNFSParticleTypes));

  m_stacked_channel = StackedHistogram<ChannelType>(
      m_label, m_xlabel, m_bins_array, int(kNChannelTypes));
}

void Histograms::InitializeDataHists() {
  const Double_t* bins = m_bins_array.GetArray();
  const char* label = m_label.c_str();
  m_selection_data =
      new MH1D(Form("selection_data_%s", label), label, NBins(), bins);

  m_bg_subbed_data =
      new MH1D(Form("bg_subbed_data_%s", label), label, NBins(), bins);

  m_unfolded = new MH1D(Form("unfolded_%s", label), label, NBins(), bins);

  m_cross_section =
      new MH1D(Form("cross_section_%s", label), label, NBins(), bins);
}

template <typename T>
void Histograms::InitializeMigrationHists(T systematic_univs) {
  const Double_t* bins = m_bins_array.GetArray();
  const char* label = m_label.c_str();
  PlotUtils::MnvH2D* migration = new PlotUtils::MnvH2D(
      Form("migration_%s", label), label, NBins(), bins, NBins(), bins);

  const bool clear_bands = true;
  m_migration = CVH2DW(migration, systematic_univs, clear_bands);

  delete migration;
}

// Accessor functions
std::map<SignalBackgroundType, MH1D*> Histograms::GetStackMap(
    SignalBackgroundType type) const {
  return m_stacked_sigbg.m_hist_map;
}

std::map<WBackgroundType, MH1D*> Histograms::GetStackMap(
    WBackgroundType type) const {
  return m_stacked_wbg.m_hist_map;
}

std::map<MesonBackgroundType, MH1D*> Histograms::GetStackMap(
    MesonBackgroundType type) const {
  return m_stacked_mesonbg.m_hist_map;
}

std::map<FSParticleType, MH1D*> Histograms::GetStackMap(
    FSParticleType type) const {
  return m_stacked_fspart.m_hist_map;
}

std::map<ChannelType, MH1D*> Histograms::GetStackMap(ChannelType type) const {
  return m_stacked_channel.m_hist_map;
}

#endif  // Histograms_cxx
