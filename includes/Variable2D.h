#ifndef VARIABLE_H
#define VARIABLE_H

#include "CVUniverse.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/MnvH2D.h"

#ifndef __CINT__  // CINT doesn't know about std::function
#include "PlotUtils/Variable2DBase.h"
#include "PlotUtils/VariableBase.h"
#endif  // __CINT__

class VariableMAT : public PlotUtils::VariableBase<CVUniverse> {
 private:
  typedef PlotUtils::HistWrapper<CVUniverse> HW;
  typedef PlotUtils::MnvH1D MH1D;
  typedef PlotUtils::Hist2DWrapper<CVUniverse> HW2D;
  typedef PlotUtils::MnvH2D MH2D;

 public:
  //=======================================================================================
  // CTOR
  //=======================================================================================
  template <class... ARGS>
  VariableMAT(ARGS... args) : PlotUtils::VariableBase<CVUniverse>(args...) {}

  //=======================================================================================
  // DECLARE NEW HISTOGRAMS
  //=======================================================================================
  // HISTWRAPPER
  HW2D m_migration;

  //=======================================================================================
  // INITIALIZE ALL HISTOGRAMS
  //=======================================================================================
  template <typename T>
  void InitializeAllHists(T univs, T univs_truth) {

    const bool clear_bands = true;  // we want empty histograms

    MH2D* temp_migration = new MH2D(
        Form("migration_%s", GetName().c_str()), GetName().c_str(),
        GetNBins(), GetBinVec().data(), GetNBins(), GetBinVec().data());
    m_migration = HW2D(temp_migration, univs, clear_bands);

    delete temp_migration;

  }

  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
  void WriteAllHistogramsToFile(TFile& f, bool isMC) const {
    f.cd();

    m_migration.hist->Write();
  }
};

class Variable2D : public PlotUtils::Variable2DBase<CVUniverse> {
 private:
  typedef PlotUtils::Hist2DWrapper<CVUniverse> HW2D;
  typedef PlotUtils::MnvH2D MH2D;

 public:
  //=======================================================================================
  // CTOR
  //=======================================================================================
  template <class... ARGS>
  Variable2D(ARGS... args) : PlotUtils::Variable2DBase<CVUniverse>(args...) {}

  //==========================================================================
  // Data members
  //==========================================================================
  std::string m_label;
  bool m_is_true;

  //==========================================================================
  // Functions
  //==========================================================================
  std::string Name() const { return m_label; }

  //=======================================================================================
  // DECLARE NEW HISTOGRAMS
  //=======================================================================================
  // HISTWRAPPER
  HW2D m_selection_data_inclusive;
  HW2D m_selection_data_lowNu;
  HW2D m_selection_data_highNu;
  HW2D m_selection_mc_inclusive;
  HW2D m_selection_mc_lowNu;
  HW2D m_selection_mc_highNu;
  HW2D m_bg_inclusive;
  HW2D m_bg_lowNu;
  HW2D m_bg_highNu;
  HW2D m_effnum_inclusive;
  HW2D m_effnum_lowNu;
  HW2D m_effnum_highNu;
  HW2D m_effdenom_inclusive;
  HW2D m_effdenom_lowNu;
  HW2D m_effdenom_highNu;

  //=======================================================================================
  // INITIALIZE ALL HISTOGRAMS
  //=======================================================================================
  template <typename T>
  void InitializeAllHists(T univs, T univs_truth) {
    const bool clear_bands = true;  // we want empty histograms

    // Data selected sample
    MH2D* temp_selection_data_inclusive = new MH2D(
        Form("selection_data_%s_inclusive", GetName().c_str()), GetName().c_str(),
        GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selection_data_inclusive = HW2D(temp_selection_data_inclusive, univs, clear_bands);

    MH2D* temp_selection_data_lowNu = new MH2D(
        Form("selection_data_%s_lowNu", GetName().c_str()), GetName().c_str(),
        GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selection_data_lowNu = HW2D(temp_selection_data_lowNu, univs, clear_bands);

    MH2D* temp_selection_data_highNu = new MH2D(
        Form("selection_data_%s_highNu", GetName().c_str()), GetName().c_str(),
        GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selection_data_highNu = HW2D(temp_selection_data_highNu, univs, clear_bands);

    // MC selected sample
    MH2D* temp_selection_mc_inclusive = new MH2D(
        Form("selection_mc_%s_inclusive", GetName().c_str()), GetName().c_str(),
        GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selection_mc_inclusive = HW2D(temp_selection_mc_inclusive, univs, clear_bands);

    MH2D* temp_selection_mc_lowNu = new MH2D(
        Form("selection_mc_%s_lowNu", GetName().c_str()), GetName().c_str(),
        GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selection_mc_lowNu = HW2D(temp_selection_mc_lowNu, univs, clear_bands);

    MH2D* temp_selection_mc_highNu = new MH2D(
        Form("selection_mc_%s_highNu", GetName().c_str()), GetName().c_str(),
        GetNBinsX(), GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selection_mc_highNu = HW2D(temp_selection_mc_highNu, univs, clear_bands);

    // MC background
    MH2D* temp_bg_inclusive = new MH2D(
        Form("bg_%s_inclusive", GetName().c_str()), GetName().c_str(), GetNBinsX(),
        GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_bg_inclusive = HW2D(temp_bg_inclusive, univs, clear_bands);

    MH2D* temp_bg_lowNu = new MH2D(
        Form("bg_%s_lowNu", GetName().c_str()), GetName().c_str(), GetNBinsX(),
        GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_bg_lowNu = HW2D(temp_bg_lowNu, univs, clear_bands);

    MH2D* temp_bg_highNu = new MH2D(
        Form("bg_%s_highNu", GetName().c_str()), GetName().c_str(), GetNBinsX(),
        GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_bg_highNu = HW2D(temp_bg_highNu, univs, clear_bands);

    // MC efficiency numerator
    MH2D* temp_effnum_inclusive = new MH2D(
        Form("effnum_%s_inclusive", GetName().c_str()), GetName().c_str(), GetNBinsX(),
        GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_effnum_inclusive = HW2D(temp_effnum_inclusive, univs, clear_bands);

    MH2D* temp_effnum_lowNu = new MH2D(
        Form("effnum_%s_lowNu", GetName().c_str()), GetName().c_str(), GetNBinsX(),
        GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_effnum_lowNu = HW2D(temp_effnum_lowNu, univs, clear_bands);

    MH2D* temp_effnum_highNu = new MH2D(
        Form("effnum_%s_highNu", GetName().c_str()), GetName().c_str(), GetNBinsX(),
        GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_effnum_highNu = HW2D(temp_effnum_highNu, univs, clear_bands);

    // MC efficiency denominator
    MH2D* temp_effdenom_inclusive = new MH2D(
        Form("effdenom_%s_inclusive", GetName().c_str()), GetName().c_str(), GetNBinsX(),
        GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_effdenom_inclusive = HW2D(temp_effdenom_inclusive, univs_truth, clear_bands);

    MH2D* temp_effdenom_lowNu = new MH2D(
        Form("effdenom_%s_lowNu", GetName().c_str()), GetName().c_str(), GetNBinsX(),
        GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_effdenom_lowNu = HW2D(temp_effdenom_lowNu, univs_truth, clear_bands);

    MH2D* temp_effdenom_highNu = new MH2D(
        Form("effdenom_%s_highNu", GetName().c_str()), GetName().c_str(), GetNBinsX(),
        GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_effdenom_highNu = HW2D(temp_effdenom_highNu, univs_truth, clear_bands);

    delete temp_selection_data_inclusive;
    delete temp_selection_data_lowNu;
    delete temp_selection_data_highNu;
    delete temp_selection_mc_inclusive;
    delete temp_selection_mc_lowNu;
    delete temp_selection_mc_highNu;
    delete temp_bg_inclusive;
    delete temp_bg_lowNu;
    delete temp_bg_highNu;
    delete temp_effnum_inclusive;
    delete temp_effnum_lowNu;
    delete temp_effnum_highNu;
    delete temp_effdenom_inclusive;
    delete temp_effdenom_lowNu;
    delete temp_effdenom_highNu;
  }

  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
  void WriteAllHistogramsToFile(TFile& f, bool isMC) const {
    f.cd();

    m_selection_data_inclusive.hist->Write();
    m_selection_data_lowNu.hist->Write();
    m_selection_data_highNu.hist->Write();
    m_selection_mc_inclusive.hist->Write();
    m_selection_mc_lowNu.hist->Write();
    m_selection_mc_highNu.hist->Write();
    m_bg_inclusive.hist->Write();
    m_bg_lowNu.hist->Write();
    m_bg_highNu.hist->Write();
    m_effnum_inclusive.hist->Write();
    m_effnum_lowNu.hist->Write();
    m_effnum_highNu.hist->Write();
    m_effdenom_inclusive.hist->Write();
    m_effdenom_lowNu.hist->Write();
    m_effdenom_highNu.hist->Write();
  }
};

#endif  // VARIABLE_H
