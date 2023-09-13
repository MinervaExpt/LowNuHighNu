#ifndef Constants_h
#define Constants_h

#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/MnvH1D.h"
#include "TObjArray.h"

// example code for timing
#include <chrono>
std::chrono::duration<double> t01; // seconds, global
std::chrono::duration<double> t12; // seconds, global
std::chrono::duration<double> t23; // seconds, global
//auto td = std::chrono::steady_clock::now();
//<...>
//auto td = std::chrono::steady_clock::now();
//tcd += std::chrono::duration_cast<std::chrono::microseconds>(td - tc);

class CVUniverse;

//==============================================================================
// Typedefs, enums, colors, constants
//==============================================================================
struct PassesCutsInfo {
  bool passes_all_cuts;
  bool GetAll() {
    return passes_all_cuts;
  };
};

struct PassesCutsInfoInclusive {
  bool passes_incl_cuts;
  bool is_dis_signal;
  bool is_lowwhighq2_sideband;
  bool is_lowq2highw_sideband;
  std::tuple<bool,bool,bool,bool> GetAll(){
    return {passes_incl_cuts, is_dis_signal, is_lowwhighq2_sideband, is_lowq2highw_sideband};
  };
};

enum ECuts {
  kNoCuts,
  kGoodObjects,
  kGoodVertex,
  kFiducialVolume,
  kMinosActivity,
  kPrecuts,
  kVtx,
  kMinosMatch,
  kMinosCharge,
  kMinosMuon,
  kPmu,
  kAllCuts,
};

enum EDataMCTruth { kData, kMC, kTruth, kNDataMCTruthTypes };

typedef std::map<ECuts, double> EventCount;

typedef PlotUtils::MnvH1D MH1D;
typedef PlotUtils::HistWrapper<CVUniverse> CVHW;
typedef PlotUtils::Hist2DWrapper<CVUniverse> CVH2DW;

typedef std::map<std::string, std::vector<CVUniverse*>> UniverseMap;

//==============================================================================

namespace CCNuPionIncConsts {
const int N_LAT_UNIVERSES = 60;
const double CHARGED_PION_MASS = 139.569;
const double CHARGED_KAON_MASS = 493.677;
const double MUON_MASS = 105.659;
const double PROTON_MASS = 938.3;
const double PI = 3.14159265358979323846;

const double bindE = 25.;
const double numi_beam_angle_rad = -0.05887;

const unsigned int MAX_N_REC_HADRONS = 10;

const int kIsVertexPion = -1;
const int kEmptyPionCandidateVector = -2;

const int kIsoProngCutVal = 2;
const double kThetamuMaxCutVal = 0.3491;  // rad (20 deg)
const double kPmuMinCutVal = 1500.;       // MeV/c
const double kPmuMaxCutVal = 20000.;      // MeV/c
const double kZVtxMinCutVal = 5990.;      // cm
const double kZVtxMaxCutVal = 8340.;      // cm
const double kApothemCutVal = 850.;       // cm
const double kTpiLoCutVal = 35.;          // MeV
const double kTpiHiCutVal = 350.;         // MeV

const bool kUseNueConstraint = true;
const int kAnaNuPDG = 14;
const bool kUseNonResPiWgt = true;
const bool kUseDeuteriumGeniePiTune = true;

const int kNFluxUniverses = 100;
}  // namespace CCNuPionIncConsts

//==============================================================================

namespace CCNuPionIncShifts {
const std::string kLowQ2PiChannel("MENU1PI");  // MENU1PI is the weight
                                               // for GENIE v4.3

const double muon_angle_res = 0.002;                 // radians
const double pion_angle_res = 5.0 * muon_angle_res;  // radians
const double vtx_smear_correction =
    0.9995;  // difference between nominal and smeared fiducial volume
             // selections

// hadron response systematics
// TODO subdivide by kinematic ranges
const double frac_proton_resp = 0.04;
const double frac_pion_resp = 0.05;
const double frac_kaon_resp = 0.05;
const double frac_em_resp = 0.03;
const double frac_neutron_resp = 0.2;
const double frac_xtalk_resp = 0.2;
const double frac_other_resp = 0.1;

// testbeam hadron efficiency constraints
const double hadron_tracking_frac_err[3] = {
    0.056, 0.015, 0.006};  // < 8 planes, 8 planes, and > 8 planes
const double proton_pid_frac_err = 0.32;
}  // namespace CCNuPionIncShifts

//==============================================================================

namespace CCNuPionIncPlotting {
// const int N_COLORS=10;
// const int colors[N_COLORS] = {kRed+1, kOrange+8, kGray+1, kOrange-7,
// kGreen+2, kYellow+2, kCyan-3, kBlue-7, kMagenta-6, kMagenta+2};

// const int colors[N_COLORS] = {kRed+3, kOrange+9, kGray+2, kOrange+3,
// kGreen+3, kYellow+3, kCyan-2, kBlue-7, kMagenta-6, kMagenta+2}; const int
// fill_colors[N_COLORS] = {kRed+1, kOrange+6, kGray+1, kOrange-5, kGreen-2,
// kYellow+2, kCyan-8, kBlue-9, kMagenta-9, kMagenta-2};

const int N_COLORS = 8;
const int colors[N_COLORS] = {kRed + 3,     kOrange + 9, kGray + 2,
                              kGreen + 3,   kCyan - 2,   kBlue,
                              kMagenta - 6, kMagenta + 2};
const int fill_colors[N_COLORS] = {kRed + 1,     kOrange + 6, kGray + 1,
                                   kGreen - 2,   kCyan - 8,   kBlue - 6,
                                   kMagenta - 9, kMagenta - 2};

const int nColors2 = 12;
const int colors2[nColors2] = {
    kMagenta + 1, kOrange + 2, kTeal - 7,   kRed - 6,   kSpring - 1, kBlue,
    kAzure + 7,   kCyan + 2,   kViolet - 4, kGreen + 3, kYellow + 2, kGray + 2,
};
const int fill_colors2[nColors2] = {
    kMagenta - 1, kOrange, kTeal - 9,   kRed - 8,   kSpring - 3, kBlue - 2,
    kAzure + 5,   kCyan,   kViolet - 6, kGreen + 1, kYellow,     kGray,
};

const int nColors3 = 7;
const int colors3[nColors3] = {kGreen + 3, kYellow + 2, kOrange + 6, kBlue,
                               kCyan + 2,  kRed + 3,    kAzure + 7};
const int fill_colors3[nColors3] = {
    kGreen + 1, kYellow, kOrange + 5, kAzure - 2, kCyan, kRed + 1, kAzure + 5};

const int nColors4 = 6;
const int colors4[nColors4] = {kRed + 3,    kBlue - 6,  kMagenta - 6,
                               kOrange + 6, kGreen - 2, kAzure + 8};
const int fill_colors4[nColors4] = {kRed + 1,    kBlue - 2,  kMagenta - 3,
                                    kOrange + 5, kGreen - 6, kAzure + 7};

const int nColors5 = 7;
const int colors5[nColors5] = {kGreen + 3, kGray + 2, kBlue,     kYellow + 2,
                               kCyan + 2,  kRed + 3,  kAzure + 7};
const int fill_colors5[nColors5] = {kGreen + 1, kGray,    kBlue - 2, kYellow,
                                    kCyan,      kRed + 1, kAzure + 5};

const int coherent_fill_color = kRed - 7;
const int coherent_color = kRed - 3;
}  // namespace CCNuPionIncPlotting

void SetHistColorScheme(PlotUtils::MnvH1D* h, const int type,
                        const int color_scheme) {
  switch (color_scheme) {
    case 2:
      h->SetFillColor(CCNuPionIncPlotting::fill_colors2[type]);
      h->SetLineColor(CCNuPionIncPlotting::colors2[type]);
      break;
    case 3:
      h->SetFillColor(CCNuPionIncPlotting::fill_colors3[type]);
      h->SetLineColor(CCNuPionIncPlotting::colors3[type]);
      break;
    case 4:
      h->SetFillColor(CCNuPionIncPlotting::fill_colors4[type]);
      h->SetLineColor(CCNuPionIncPlotting::colors4[type]);
      break;
    default:
      h->SetFillColor(CCNuPionIncPlotting::fill_colors[type]);
      h->SetLineColor(CCNuPionIncPlotting::colors[type]);
  }
}

#endif  // Constants_h
