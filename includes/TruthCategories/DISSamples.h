#ifndef Sidebands_H
#define Sidebands_H

#include <algorithm>

#include "../CVUniverse.h"
#include "../SignalDefinition.h"  // IsSignal

namespace sidebands {
const int kNWFitCategories = 3;
const std::string kFitVarString("wexp_fit");
const double kSidebandCutVal = 1500;
const int kWSideband_ColorScheme = 3;
}  // namespace sidebands

enum SampleType {
  kDIS_Signal,
  kSideband_LowWHighQ2,
  kSideband_LowQ2HighW,
  kOther
};

// n_categories: should we chop up the sideband into 2 or 3 categories?
// n_categories == 2: WexpTrue [0, 1.8],[1.8, up]
// n_categories == 3: WexpTrue [0, 1.4],[1.4, 1.8],[1.8, up]
// This is separate from the question of whether we float or fix each category.
// For that question, see the implementation of the fitter.
SampleType GetTrueSampleType(const CVUniverse& universe,
                             const SignalDefinition signal_definition,
                             const int n_categories = 4) {
  //  if (n_categories != 2 && n_categories != 3)
  //    throw std::invalid_argument("GetWSidebandType: n_categories is 2 or 3");

  double Wexp_true = universe.GetWexpTrue();

  if (IsSignal(universe, signal_definition)) {
    if (Wexp_true > 2000 && Q2exp_true > 1000000) return kDIS_Signal;
    if (Wexp_true > 2000 && Q2exp_true < 9000000) return kSideband_LowQ2HighW;
    if (Wexp_true < 1900 && Q2exp_true > 1000000)
      return kSideband_LowWHighQ2;
    else
      return kOther;
  }
}

std::string GetTruthClassification_LegendLabel(const SampleType sample,
                                               const int n_categories = 4) {
  //  if (n_categories != 2 && n_categories != 3)
  //    throw std::invalid_argument("GetWSidebandType: n_categories is 2 or 3");

  switch (sample) {
    case kDIS_Signal:
      return "DIS";
    case kSideband_LowWHighQ2:
      return "Low W High Q2";
    case kSideband_LowQ2HighW:
      return "Low Q2 High W";
    default:
      return "ERROR";
  }
}

std::string GetTruthClassification_Name(SampleType sample) {
  switch (sample) {
    case kDIS_Signal:
      return "DIS";
    case kSideband_LowWHighQ2:
      return "Low W High Q2";
    case kSideband_LowQ2HighW:
      return "Low Q2 High W";
    default:
      return "ERROR";
  }
}

//==============================================================================

#endif
