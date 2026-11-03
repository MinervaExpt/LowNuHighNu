#ifndef makeXsecMCInputs_C
#define makeXsecMCInputs_C

#include <cassert>
#include <ctime>

#include "TNamed.h"
#include "TString.h"
#include "TSystem.h"

#include "includes/Binning.h"
#include "includes/CVUniverse.h"
#include "includes/Constants.h"
#include "includes/Cuts.h"
#include "includes/LowNuHighNuEvent.h"
#include "includes/MacroUtil.h"
#include "includes/SignalDefinition.h"
#include "includes/TruthCategories/Sidebands.h"  // sidebands::kFitVarString, IsWSideband
#include "includes/Variable.h"
#include "includes/Variable2D.h"
#include "includes/common_functions.h"  // GetVar, WritePOT

class Variable;

//==============================================================================
// Helper Functions
//==============================================================================
namespace make_xsec_mc_inputs {
typedef Variable Var;
typedef VariableMAT VarMAT;
typedef Variable2D Var2D;

// Temporary solution to the fact that Ben's Variable class accepts a TArrayD
// binning but the MAT VariableBase class only accepts a std::vector<double>
// binning
std::vector<double> ConvertTArrayDToStdVector(const TArrayD& array) {
  std::vector<double> result;
  result.reserve(array.GetSize());  // Optional for optimization

  for (int i = 0; i < array.GetSize(); i++) {
    result.push_back(array[i]);
  }

  return result;
}

std::vector<Variable*> GetLowNuHighNuVariables(bool include_truth_vars = true) {
  Var* pmu = new Var("pmu", "p_{#mu}", "MeV", CCPi::GetBinning("pmu"),
                     &CVUniverse::GetPmu);

  Var* thetamu_deg =
      new Var("thetamu_deg", "#theta_{#mu}", "deg",
              CCPi::GetBinning("thetamu_deg"), &CVUniverse::GetThetamuDeg);

  Var* enu = new Var("enu", "E_{#nu}", "MeV", CCPi::GetBinning("enu"),
                     &CVUniverse::GetEnu_GeV);

  Var* ehad = new Var("ehad", "ehad", "MeV", CCPi::GetBinning("ehad"),
                      &CVUniverse::GetEhad_GeV);

  // Boundary-region diagnostic: signed distance from the low-nu staircase.
  // Negative => classified low-nu, positive => high-nu.  See LowNuBoundary.h.
  Var* delta_nu =
      new Var("delta_nu", "#nu - #nu_{cut}(E_{#nu})", "GeV",
              CCPi::GetBinning("delta_nu"), &CVUniverse::GetDeltaNu_GeV);

  // True Variables
  bool is_true = true;
  Var* pmu_true =
      new Var("pmu_true", "p_{#mu} True", pmu->m_units,
              pmu->m_hists.m_bins_array, &CVUniverse::GetPmuTrue, is_true);

  Var* thetamu_deg_true =
      new Var("thetamu_deg_true", "#theta_{#mu} True", thetamu_deg->m_units,
              thetamu_deg->m_hists.m_bins_array, &CVUniverse::GetThetamuTrueDeg,
              is_true);

  Var* enu_true =
      new Var("enu_true", "E_{#nu} True", enu->m_units,
              enu->m_hists.m_bins_array, &CVUniverse::GetEnuTrue_GeV, is_true);

  Var* ehad_true =
      new Var("ehad_true", "ehad True", "MeV", ehad->m_hists.m_bins_array,
              &CVUniverse::GetEhadTrue_GeV, is_true);

  Var* delta_nu_true =
      new Var("delta_nu_true", "#nu - #nu_{cut}(E_{#nu}) True",
              delta_nu->m_units, delta_nu->m_hists.m_bins_array,
              &CVUniverse::GetDeltaNuTrue_GeV, is_true);

  std::vector<Var*> variables = {pmu, thetamu_deg, enu, ehad, delta_nu};

  if (include_truth_vars) {
    variables.push_back(pmu_true);
    variables.push_back(thetamu_deg_true);
    variables.push_back(enu_true);
    variables.push_back(ehad_true);
    variables.push_back(delta_nu_true);
  }

  return variables;
}

std::vector<VariableMAT*> GetLowNuHighNuMATVariables(bool include_truth_vars = true) {
  VarMAT* vmat_ehad = new VarMAT(
    "vmat_ehad", "vmat_ehad", ConvertTArrayDToStdVector(CCPi::GetBinning("ehad_fine")),
    &CVUniverse::GetEhad_GeV, &CVUniverse::GetEhadTrue_GeV);
  return std::vector<VarMAT*>{vmat_ehad};
}

std::vector<Variable2D*> GetLowNuHighNu2DVariables(
    bool include_truth_vars = true) {
  VarMAT* vmat_enu = new VarMAT(
      "enu", "enu", ConvertTArrayDToStdVector(CCPi::GetBinning("enu")),
      &CVUniverse::GetEnu_GeV, &CVUniverse::GetEnuTrue_GeV);
  VarMAT* vmat_ehad = new VarMAT(
      "ehad", "ehad", ConvertTArrayDToStdVector(CCPi::GetBinning("ehad")),
      &CVUniverse::GetEhad_GeV, &CVUniverse::GetEhadTrue_GeV);
  VarMAT* vmat_ehad_fine = new VarMAT(
      "ehad_fine", "ehad_fine", ConvertTArrayDToStdVector(CCPi::GetBinning("ehad_fine")),
      &CVUniverse::GetEhad_GeV, &CVUniverse::GetEhadTrue_GeV);
  Var2D* enu_ehad = new Var2D(*vmat_enu, *vmat_ehad);
  Var2D* enu_ehad_fine = new Var2D(*vmat_enu, *vmat_ehad_fine);

  // Boundary-region 2D map: fine-binned (E_nu, nu) for visualizing the
  // staircase cut directly with events overlaid.  The standard `enu` binning
  // doesn't resolve the corners at 3, 7, 12 GeV cleanly; `enu_boundary`
  // and `ehad_boundary` are uniform fine binnings defined in Binning.h.
  VarMAT* vmat_enu_boundary = new VarMAT(
      "enu_boundary", "enu_boundary",
      ConvertTArrayDToStdVector(CCPi::GetBinning("enu_boundary")),
      &CVUniverse::GetEnu_GeV, &CVUniverse::GetEnuTrue_GeV);
  VarMAT* vmat_ehad_boundary = new VarMAT(
      "ehad_boundary", "ehad_boundary",
      ConvertTArrayDToStdVector(CCPi::GetBinning("ehad_boundary")),
      &CVUniverse::GetEhad_GeV, &CVUniverse::GetEhadTrue_GeV);
  // Explicit name: the auto-joined default would be
  // "enu_boundary_ehad_boundary" (component names joined with "_").
  Var2D* enu_ehad_boundary = new Var2D(std::string("enu_ehad_boundary"),
                                       *vmat_enu_boundary, *vmat_ehad_boundary);

  // Native (E_nu, Delta-nu) map.  x = nominal analysis E_nu binning (edges at
  // the 3/7/12 GeV staircase corners and at 20 GeV, the top of the NOMAD
  // normalization integral); y = fine Delta-nu on [-2, 3] GeV.  Delta-nu is
  // filled directly per universe, so the plotter no longer has to rebuild it
  // from enu_ehad_boundary by a column shift (that map stops at nu = 2.5 GeV
  // and E_nu = 15 GeV).  Auto-joined name: "enu_delta_nu".
  VarMAT* vmat_delta_nu_2d = new VarMAT(
      "delta_nu", "delta_nu",
      ConvertTArrayDToStdVector(CCPi::GetBinning("delta_nu_2d")),
      &CVUniverse::GetDeltaNu_GeV, &CVUniverse::GetDeltaNuTrue_GeV);
  Var2D* enu_delta_nu = new Var2D(*vmat_enu, *vmat_delta_nu_2d);

  return std::vector<Var2D*>{enu_ehad, enu_ehad_fine, enu_ehad_boundary,
                             enu_delta_nu};
}

}  // namespace make_xsec_mc_inputs

std::vector<Variable*> GetAnalysisVariables(SignalDefinition signal_definition,
                                            bool include_truth_vars = false) {
  std::vector<Variable*> variables;
  std::cout << "My signal is " << GetSignalName(signal_definition) << std::endl;
  switch (signal_definition) {
    case kLowNuHighNu:
      variables =
          make_xsec_mc_inputs::GetLowNuHighNuVariables(include_truth_vars);
      break;
    default:
      std::cerr << "Variables for other SDs not yet implemented.\n";
      std::exit(1);
  }
  return variables;
}

void SyncAllHists(Variable& v) {
  v.m_hists.m_selection_mc.SyncCVHistos();
  v.m_hists.m_bg.SyncCVHistos();
  v.m_hists.m_bg_loW.SyncCVHistos();
  v.m_hists.m_bg_midW.SyncCVHistos();
  v.m_hists.m_bg_hiW.SyncCVHistos();
  v.m_hists.m_wsidebandfit_sig.SyncCVHistos();
  v.m_hists.m_wsidebandfit_loW.SyncCVHistos();
  v.m_hists.m_wsidebandfit_midW.SyncCVHistos();
  v.m_hists.m_wsidebandfit_hiW.SyncCVHistos();
  v.m_hists.m_effnum.SyncCVHistos();
  v.m_hists.m_effden.SyncCVHistos();
}

// 2D analogue of SyncAllHists(Variable&): push each filled universe into the
// MnvH2D vertical error bands before writing.  Only the MC wrappers are synced;
// the data wrappers are filled on the CV hist directly (MnvH2D propagates a
// plain Fill to its bands) and carry only statistical uncertainty, matching the
// 1D data treatment.  This is the "// to-do follow up on this" from the write
// loop below -- needed for systematic bands on the (E_nu, nu) map + E_nu slices.
void SyncAllHists(Variable2D& v) {
  v.m_selection_mc_inclusive.SyncCVHistos();
  v.m_selection_mc_lowNu.SyncCVHistos();
  v.m_selection_mc_highNu.SyncCVHistos();
  v.m_bg_inclusive.SyncCVHistos();
  v.m_bg_lowNu.SyncCVHistos();
  v.m_bg_highNu.SyncCVHistos();
  v.m_effnum_inclusive.SyncCVHistos();
  v.m_effnum_lowNu.SyncCVHistos();
  v.m_effnum_highNu.SyncCVHistos();
  v.m_effdenom_inclusive.SyncCVHistos();
  v.m_effdenom_lowNu.SyncCVHistos();
  v.m_effdenom_highNu.SyncCVHistos();
}

//==============================================================================
// Low-nu misclassification of selected signal MC (reco call vs truth call)
//==============================================================================
// x = E_nu (nominal analysis binning, GeV); y = category, one unit-wide bin each:
//   y bin 1 [0,1) LL : reco low-nu,  truth low-nu
//   y bin 2 [1,2) LH : reco low-nu,  truth high-nu
//   y bin 3 [2,3) HL : reco high-nu, truth low-nu
//   y bin 4 [3,4) HH : reco high-nu, truth high-nu
// Same population and weights as migration_delta_nu (selected signal MC, reco
// loop, every universe), so summed over E_nu (incl. under/overflow) the four
// categories reproduce its sign quadrants.  Two versions:
//   misclass_lownu_enu_reco -- x = reco E_nu: composition of each analysis bin
//                              (what the reco-binned efficiency numerator sees)
//   misclass_lownu_enu_true -- x = true E_nu (what the truth-binned denominator sees)
// No axis labels on purpose: labelled axes change how ROOT adds/merges hists.
// (A typedef rather than a struct: ACLiC would try to build a dictionary for a
// class defined in this file, and Hist2DWrapper<CVUniverse> has none.)
typedef PlotUtils::Hist2DWrapper<CVUniverse> MisclassHist;

MisclassHist MakeMisclassHist(const char* name, const char* title,
                              UniverseMap& error_bands) {
  const TArrayD enu_bins = CCPi::GetBinning("enu");
  const double cat_bins[5] = {0., 1., 2., 3., 4.};
  PlotUtils::MnvH2D* temp = new PlotUtils::MnvH2D(
      name, title, enu_bins.GetSize() - 1, enu_bins.GetArray(), 4, cat_bins);
  const bool clear_bands = true;  // we want empty histograms
  MisclassHist hw(temp, error_bands, clear_bands);
  delete temp;  // the wrapper holds its own copy (as in VariableMAT)
  return hw;
}

// y-bin centre for the category of one selected signal MC event
double MisclassCategory(const bool reco_low, const bool true_low) {
  return (reco_low ? 0. : 2.) + (true_low ? 0. : 1.) + 0.5;
}

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillData(const CCPi::MacroUtil& util,
                     std::vector<Variable*>& variables,
                     std::vector<VariableMAT*>& variables_MAT,
                     std::vector<Variable2D*>& variables2D,
                     bool truncate_run) {
  // Fill data distributions.
  const bool is_mc = false;
  const bool is_truth = false;
  std::cout << "*** Starting Data Loop ***" << std::endl;
  Long64_t n_entries = util.GetDataEntries();
  Long64_t n_entries_to_loop_over = truncate_run ? 1000 : n_entries;
  for (Long64_t i_event = 0; i_event < n_entries_to_loop_over; ++i_event) {
    if (i_event % (n_entries_to_loop_over / 10) == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;
    util.m_data_universe->SetEntry(i_event);

    LowNuHighNuEvent event(is_mc, is_truth, util.m_signal_definition,
                    util.m_data_universe);

    // Check cuts
    // And extract whether this is w sideband and get candidate pion indices
    PassesCutsInfo cuts_info = PassesCuts(event);

    // Set what we've learned to the event
    event.m_passes_cuts = cuts_info.GetAll();

    lownuhighnu_event::FillRecoEvent(event, variables, variables_MAT, variables2D);
  }
  std::cout << "*** Done Data Loop ***\n\n";
}

void LoopAndFillMCXSecInputs(const CCPi::MacroUtil& util,
                             const EDataMCTruth& type,
                             std::vector<Variable*>& variables,
                             std::vector<VariableMAT*>& variables_MAT,
                             std::vector<Variable2D*>& variables2D,
                             bool truncate_run,
                             MisclassHist* misclass_enu_reco = nullptr,
                             MisclassHist* misclass_enu_true = nullptr) {
  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);
  const UniverseMap error_bands =
      is_truth ? util.m_error_bands_truth : util.m_error_bands;
  Long64_t n_entries_to_loop_over = truncate_run ? 1000 : n_entries;
  for (Long64_t i_event = 0; i_event < n_entries_to_loop_over; ++i_event) {
    if (i_event % (n_entries_to_loop_over / 10) == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;

    // Variables that hold info about whether the CVU passes cuts
    PassesCutsInfo cv_cuts_info;
    bool checked_cv = false;

    // Loop universes, make cuts, and fill
    for (auto error_band : error_bands) {
      std::vector<CVUniverse*> universes = error_band.second;
      for (auto universe : universes) {
        universe->SetEntry(i_event);
        // std::cout << universe->ShortName() << "\n";
        // if (universe->GetDouble("mc_incoming") == 12 &&
        //    universe->ShortName() == "cv")
        //  universe->PrintArachneLink();

        // MC weight is set below, once, only where it is needed
        LowNuHighNuEvent event(is_mc, is_truth, util.m_signal_definition,
                               universe);

        //===============
        // FILL TRUTH
        //===============
        if (type == kTruth) {
          // Truth fills use the weight only for signal events
          if (event.m_is_signal) event.m_weight = universe->GetWeight();
          lownuhighnu_event::FillTruthEvent(event, variables, variables_MAT, variables2D);
          continue;
        }

        //===============
        // CHECK CUTS
        //===============
        // Universe only affects weights
        if (universe->IsVerticalOnly()) {
          // Only check vertical-only universes once.
          if (!checked_cv) {
            // Check cuts
            cv_cuts_info = PassesCuts(event);
            checked_cv = true;
          }

          // Already checked a vertical-only universe
          if (checked_cv) {
            event.m_passes_cuts = cv_cuts_info.GetAll();
          }
          // Universe shifts something laterally
        } else {
          PassesCutsInfo cuts_info = PassesCuts(event);
          event.m_passes_cuts = cuts_info.GetAll();
        }

        // Weight: computed once per (event, universe), and only for events
        // that pass the cuts -- every reco fill below requires m_passes_cuts.
        // (This used to be a second GetWeight() call on top of the one in the
        // LowNuHighNuEvent constructor, for the CCPi analysis's pion-dependent
        // node-cut systematic; not used here.)
        if (event.m_passes_cuts) event.m_weight = universe->GetWeight();

        //===============
        // FILL RECO
        //===============
        lownuhighnu_event::FillRecoEvent(event, variables, variables_MAT, variables2D);

        //===============
        // FILL LOW-NU MISCLASSIFICATION (selected signal MC; same population
        // as migration_delta_nu)
        //===============
        if (misclass_enu_reco && misclass_enu_true && event.m_is_mc &&
            event.m_is_signal && event.m_passes_cuts) {
          // m_is_lowNu is the reco call for THIS universe (lateral shifts move it)
          const double cat =
              MisclassCategory(event.m_is_lowNu, IsLowNuTruth(*universe));
          misclass_enu_reco->FillUniverse(*universe, universe->GetEnu_GeV(),
                                          cat, event.m_weight);
          misclass_enu_true->FillUniverse(*universe, universe->GetEnuTrue_GeV(),
                                          cat, event.m_weight);
        }
      }  // universes
    }    // error bands
  }      // events
  std::cout << "*** Done MC Loop ***\n\n";
}

//==============================================================================
// Main
//==============================================================================
void makeCrossSectionMCInputs(int signal_definition_int = 1,
                              std::string plist = "ME1A",
                              bool do_systematics = false,
                              bool do_truth = false, bool is_grid = false,
                              std::string mc_input_file = "",
                              std::string data_input_file = "",
                              bool test_run = false,
                              std::string out_tag = "") {
  // Only signal definition 1 (kLowNuHighNu) is supported; it writes the
  // inclusive, lowNu and highNu samples.  Definition 0 (kInclusive: Anne's
  // DIS/sideband variant, InclusiveEvent) was removed 2026-11-03.
  if (signal_definition_int != kLowNuHighNu) {
    std::cerr << "ERROR: only signal definition 1 (kLowNuHighNu) is supported.\n";
    std::exit(1);
  }

  // INPUT TUPLES
  const bool is_mc = true;
  // Playlist files must be given explicitly (run_playlists.sh does this).  The
  // old fallback, GetPlaylistFile() in playlist_methods.h, points at
  // pre-open-data MAD playlists that no longer exist.
  if (mc_input_file.empty() || data_input_file.empty()) {
    std::cerr << "ERROR: pass both mc_input_file and data_input_file "
                 "(full paths to the playlist .txt files).\n";
    std::exit(1);
  }
  const std::string mc_file_list = mc_input_file;
  const std::string data_file_list = data_input_file;

  // INIT MACRO UTILITY
  const std::string macro("XSecInputs");
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list, 
                       plist, do_truth, is_grid, do_systematics);
  util.PrintMacroConfiguration(macro);

  // INIT OUTPUT
  auto time = std::time(nullptr);
  char tchar[100];
  std::strftime(tchar, sizeof(tchar), "%F", std::gmtime(&time));  // YYYY-MM-dd
  const std::string tag = tchar;
  // Optional out_tag (e.g. "wk45" or "wk45_p1of2" from run_playlists.sh) is
  // appended so tagged or split runs of the same playlist never collide.  The
  // playlist name itself must stay real: MAT uses it for playlist-dependent
  // weights.
  const std::string out_sfx = out_tag.empty() ? "" : "_" + out_tag;
  std::string outfile_name(Form("%s_%d%d%d%d_%s_%s%s.root", macro.c_str(),
                                signal_definition_int, int(do_systematics),
                                int(do_truth), int(is_grid), plist.c_str(),
                                tag.c_str(), out_sfx.c_str()));
  std::cout << "Saving output to " << outfile_name << "\n\n";
  TFile fout(outfile_name.c_str(), "RECREATE");

  // PROVENANCE: record which code built this file, read at launch (the
  // checkout may change while a long run is going).  git_describe names the
  // commit and ends in "-dirty" if there were uncommitted edits; git_diff_md5
  // fingerprints those edits (the md5 of an empty diff, d41d8cd98f00, means a
  // clean tree).  --no-optional-locks: parallel jobs must not race on
  // .git/index.lock.
  {
    const TString repo =
        TString::Format("%s/LowNuHighNu", gSystem->Getenv("TOPDIR"));
    const TString describe = gSystem->GetFromPipe(TString::Format(
        "git -C %s --no-optional-locks describe --always --dirty 2>&1",
        repo.Data()));
    const TString diff_md5 = gSystem->GetFromPipe(TString::Format(
        "git -C %s --no-optional-locks diff HEAD 2>/dev/null | md5sum | cut -c1-12",
        repo.Data()));
    std::cout << "Code version: " << describe << "  (diff md5 " << diff_md5
              << ")\n\n";
    fout.cd();
    TNamed("git_describe", describe.Data()).Write();
    TNamed("git_diff_md5", diff_md5.Data()).Write();
  }

  // INIT VARS, HISTOS, AND EVENT COUNTERS
  const bool do_truth_vars = true;
  std::vector<Variable*> variables =
      GetAnalysisVariables(util.m_signal_definition, do_truth_vars);
  std::vector<VariableMAT*> variables_MAT =
      make_xsec_mc_inputs::GetLowNuHighNuMATVariables(do_truth_vars);
  std::vector<Variable2D*> variables2D =
      make_xsec_mc_inputs::GetLowNuHighNu2DVariables(do_truth_vars);

  for (auto v : variables)
    v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);
  for (auto v : variables_MAT) v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);
  for (auto v : variables2D) v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);

  // Low-nu misclassification hists (filled in the MC reco loop only)
  MisclassHist misclass_enu_reco = MakeMisclassHist(
      "misclass_lownu_enu_reco",
      "low-nu reco/truth call (LL,LH,HL,HH) vs reco E_nu", util.m_error_bands);
  MisclassHist misclass_enu_true = MakeMisclassHist(
      "misclass_lownu_enu_true",
      "low-nu reco/truth call (LL,LH,HL,HH) vs true E_nu", util.m_error_bands);

  // LOOP DATA
  // Uses the full `variables` (incl. truth vars); FillSelected skips truth vars
  // per-variable on data (the `continue` in LowNuHighNuEvent.cxx), so the 1D
  // selection_data histos fill into the same objects the write loop below emits.
  // (Previously a separate truth-free `variables_lessTruth` was looped here but
  // never written, leaving selection_data_<var> empty in the output.)
  LoopAndFillData(util, variables, variables_MAT, variables2D, test_run);

  // LOOP MC RECO
  for (auto band : util.m_error_bands) {
    std::vector<CVUniverse*> universes = band.second;
    for (auto universe : universes) universe->SetTruth(false);
  }
  LoopAndFillMCXSecInputs(util, kMC, variables, variables_MAT, variables2D, test_run,
                          &misclass_enu_reco, &misclass_enu_true);

  // LOOP TRUTH
  if (util.m_do_truth) {
    // m_is_truth is static, so we turn it on now
    for (auto band : util.m_error_bands_truth) {
      std::vector<CVUniverse*> universes = band.second;
      for (auto universe : universes) universe->SetTruth(true);
    }
    LoopAndFillMCXSecInputs(util, kTruth, variables, variables_MAT, variables2D, test_run);
  }

  // WRITE TO FILE
  std::cout << "Synching and Writing\n\n";
  WritePOT(fout, true, util.m_mc_pot);
  WritePOT(fout, false, util.m_data_pot);
  fout.cd();
  for (auto v : variables) {
    SyncAllHists(*v);
    v->WriteAllHistogramsToFile(fout);
  }
  for (auto v: variables_MAT) {
    v->WriteAllHistogramsToFile(fout, true);
  }
  for (auto v : variables2D) {
    SyncAllHists(*v);
    v->WriteAllHistogramsToFile(fout, true);
  }
  fout.cd();
  misclass_enu_reco.SyncCVHistos();
  misclass_enu_reco.hist->Write();
  misclass_enu_true.SyncCVHistos();
  misclass_enu_true.hist->Write();
}

#endif  // makeXsecMCInputs_C
