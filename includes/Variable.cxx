#ifndef Variable_cxx
#define Variable_cxx

#include "Variable.h"

#include <stdlib.h>  // exit

#include "TruthMatching.h"

// CTOR -- default
Variable::Variable()
    : m_label(),
      m_units(),
      m_pointer_to_GetValue(&CVUniverse::GetDummyVar),
      m_hists(),
      m_is_true(false) {}

// CTOR -- uniform binning
Variable::Variable(const std::string label, const std::string xlabel,
                   const std::string units, const int nbins, const double xmin,
                   const double xmax, PointerToCVUniverseFunction p,
                   const bool is_true)
    : m_label(label),
      m_units(units),
      m_pointer_to_GetValue(p),
      m_hists(m_label, xlabel, nbins, xmin, xmax),
      m_is_true(is_true) {}

// CTOR -- variable binning
Variable::Variable(const std::string label, const std::string xlabel,
                   const std::string units, const TArrayD& bins_array,
                   PointerToCVUniverseFunction p, const bool is_true)
    : m_label(label),
      m_units(units),
      m_pointer_to_GetValue(p),
      m_hists(m_label, xlabel, bins_array),
      m_is_true(is_true) {}

// GetValue defines this variable
double Variable::GetValue(const CVUniverse& universe) const {
  return m_pointer_to_GetValue(universe);
}

// Histogram Initialization
template <typename T>
void Variable::InitializeAllHists(T systematic_univs,
                                  T systematic_univs_truth) {
  m_hists.InitializeAllHists(systematic_univs, systematic_univs_truth);
}

void Variable::InitializeStackedHists() { m_hists.InitializeStackedHists(); }

void Variable::InitializeDataHists() { m_hists.InitializeDataHists(); }

// Histogram Access
template <typename T>
PlotUtils::MnvH1D* Variable::GetStackComponentHist(T type) const {
  std::map<T, PlotUtils::MnvH1D*> stack_map = m_hists.GetStackMap(type);
  return stack_map[type];
}

// Save with the object names that hists were initialized with
void Variable::WriteAllHistogramsToFile(TFile& fout) const {
  fout.cd();
  m_hists.m_selection_data->Write();
  m_hists.m_selection_mc.hist->Write();
  m_hists.m_bg.hist->Write();
  m_hists.m_effnum.hist->Write();
  m_hists.m_effden.hist->Write();
  if (!m_is_true)
    m_hists.m_migration.hist->Write();
}

void Variable::LoadDataHistsFromFile(TFile& fin) {
  m_hists.LoadDataHistsFromFile(fin);
}

void Variable::LoadMCHistsFromFile(TFile& fin, UniverseMap& error_bands) {
  m_hists.LoadMCHistsFromFile(fin, error_bands, m_is_true);
}

template <typename T>
TObjArray Variable::GetStackArray(T type) const {
  if (std::is_same<T, SignalBackgroundType>::value)
    return m_hists.m_stacked_sigbg.m_hist_array;
  else if (std::is_same<T, WBackgroundType>::value)
    return m_hists.m_stacked_wbg.m_hist_array;
  else if (std::is_same<T, MesonBackgroundType>::value)
    return m_hists.m_stacked_mesonbg.m_hist_array;
  else if (std::is_same<T, FSParticleType>::value)
    return m_hists.m_stacked_fspart.m_hist_array;
  else if (std::is_same<T, ChannelType>::value)
    return m_hists.m_stacked_channel.m_hist_array;
  else {
    std::cerr << "GetStackArray: Unknown truth type.\n";
    std::exit(1);
  }
}

#endif  // Variable_cxx
