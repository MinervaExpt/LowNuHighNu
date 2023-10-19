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

 public:
  //=======================================================================================
  // CTOR
  //=======================================================================================
  template <class... ARGS>
  VariableMAT(ARGS... args)
      : PlotUtils::VariableBase<CVUniverse>(args...) {}

  //=======================================================================================
  // INITIALIZE ALL HISTOGRAMS
  //=======================================================================================
  template <typename T>
  void InitializeAllHists(T univs) {
    return;
  }

  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
  void WriteAllHistogramsToFile(TFile& f, bool isMC) const {
    return;
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
  Variable2D(ARGS... args)
      : PlotUtils::Variable2DBase<CVUniverse>(args...) {}

  //=======================================================================================
  // DECLARE NEW HISTOGRAMS
  //=======================================================================================
  // HISTWRAPPER
  // selected mc reco histwrapper
  HW2D m_selection_mc;
  HW2D m_selection_data;
  HW2D m_effnum;
  HW2D m_bg;

  //=======================================================================================
  // INITIALIZE ALL HISTOGRAMS
  //=======================================================================================
  template <typename T>
  void InitializeAllHists(T univs) {
    const bool clear_bands = true;  // we want empty histograms

    MH2D* temp_selection_mc =
        new MH2D(Form("selection_mc_%s", GetName().c_str()), GetName().c_str(), GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selection_mc = HW2D(temp_selection_mc, univs, clear_bands);

    MH2D* temp_selection_data =
        new MH2D(Form("selection_data_%s", GetName().c_str()), GetName().c_str(), GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selection_data = HW2D(temp_selection_data, univs, clear_bands);

    MH2D* temp_effnum =
        new MH2D(Form("effnum_%s", GetName().c_str()), GetName().c_str(), GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_effnum = HW2D(temp_effnum, univs, clear_bands);

    MH2D* temp_bg =
        new MH2D(Form("bg_%s", GetName().c_str()), GetName().c_str(), GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_bg = HW2D(temp_bg, univs, clear_bands);

    delete temp_selection_mc;
    delete temp_selection_data;
    delete temp_effnum;
    delete temp_bg;
  }

  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
  void WriteAllHistogramsToFile(TFile& f, bool isMC) const {
    f.cd();

    m_selection_mc.hist->Write();
    m_selection_data.hist->Write();
    m_effnum.hist->Write();
    m_bg.hist->Write();

  }
};

#endif  // VARIABLE_H
