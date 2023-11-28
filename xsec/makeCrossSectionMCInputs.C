#ifndef makeXsecMCInputs_C
#define makeXsecMCInputs_C

#include <cassert>
#include <ctime>

#include "includes/Binning.h"
#include "includes/CVUniverse.h"
#include "includes/Constants.h"
#include "includes/Cuts.h"
#include "includes/InclusiveEvent.h"
#include "includes/LowNuHighNuEvent.h"
#include "includes/MacroUtil.h"
#include "includes/SignalDefinition.h"
#include "includes/TruthCategories/Sidebands.h"  // sidebands::kFitVarString, IsWSideband
#include "includes/Variable.h"
#include "includes/Variable2D.h"
#include "includes/common_functions.h"  // GetVar, WritePOT
#include "playlist_methods.h"           // GetPlaylistFile

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
                     &CVUniverse::GetEnu);

  Var* ehad = new Var("ehad", "ehad", "MeV", CCPi::GetBinning("ehad"),
                      &CVUniverse::GetEhad);

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
              enu->m_hists.m_bins_array, &CVUniverse::GetEnuTrue, is_true);

  Var* ehad_true =
      new Var("ehad_true", "ehad True", "MeV", ehad->m_hists.m_bins_array,
              &CVUniverse::GetEhadTrue, is_true);

  std::vector<Var*> variables = {pmu, thetamu_deg, enu, ehad};

  if (include_truth_vars) {
    variables.push_back(pmu_true);
    variables.push_back(thetamu_deg_true);
    variables.push_back(enu_true);
    variables.push_back(ehad_true);
  }

  return variables;
}

std::vector<VariableMAT*> GetLowNuHighNuMATVariables(bool include_truth_vars = true) {
  VarMAT* vmat_ehad = new VarMAT(
    "vmat_ehad", "vmat_ehad", ConvertTArrayDToStdVector(CCPi::GetBinning("ehad_fine")),
    &CVUniverse::GetEhad, &CVUniverse::GetEhadTrue);
  return std::vector<VarMAT*>{vmat_ehad};
}

std::vector<Variable2D*> GetLowNuHighNu2DVariables(
    bool include_truth_vars = true) {
  VarMAT* vmat_enu = new VarMAT(
      "enu", "enu", ConvertTArrayDToStdVector(CCPi::GetBinning("enu")),
      &CVUniverse::GetEnu, &CVUniverse::GetEnuTrue);
  VarMAT* vmat_ehad = new VarMAT(
      "ehad", "ehad", ConvertTArrayDToStdVector(CCPi::GetBinning("ehad")),
      &CVUniverse::GetEhad, &CVUniverse::GetEhadTrue);
  Var2D* enu_ehad = new Var2D(*vmat_enu, *vmat_ehad);

  return std::vector<Var2D*>{enu_ehad};
}

std::vector<Variable*> GetInclusiveVariables(bool include_truth_vars = true) {
  Var* pmu = new Var("pmu", "p_{#mu}", "MeV", CCPi::GetBinning("pmu"),
                     &CVUniverse::GetPmu);

  Var* thetamu_deg =
      new Var("thetamu_deg", "#theta_{#mu}", "deg",
              CCPi::GetBinning("thetamu_deg"), &CVUniverse::GetThetamuDeg);

  Var* enu = new Var("enu", "E_{#nu}", "MeV", CCPi::GetBinning("enu"),
                     &CVUniverse::GetEnu);

  Var* q2 = new Var("q2", "Q^{2}", "MeV^{2}", CCPi::GetBinning("q2"),
                    &CVUniverse::GetQ2);

  Var* wexp = new Var("wexp", "W_{exp}", "MeV", CCPi::GetBinning("wexp"),
                      &CVUniverse::GetWexp);

  // Ehad variables
  Var* ehad = new Var("ehad", "ehad", "MeV", CCPi::GetBinning("ehad"),
                      &CVUniverse::GetEhad);

  bool is_true = true;
  Var* ehad_true =
      new Var("ehad_true", "ehad True", "MeV", ehad->m_hists.m_bins_array,
              &CVUniverse::GetEhadTrue);
  ehad_true->m_is_true = true;

  Var* thetamu_deg_true =
      new Var("thetamu_deg_true", "#theta_{#mu} True", thetamu_deg->m_units,
              thetamu_deg->m_hists.m_bins_array, &CVUniverse::GetThetamuTrueDeg,
              is_true);

  Var* enu_true =
      new Var("enu_true", "E_{#nu} True", enu->m_units,
              enu->m_hists.m_bins_array, &CVUniverse::GetEnuTrue, is_true);

  Var* q2_true =
      new Var("q2_true", "Q^{2} True", q2->m_units, q2->m_hists.m_bins_array,
              &CVUniverse::GetQ2True, is_true);
  Var* pmu_true =
      new Var("pmu_true", "p_{#mu} True", pmu->m_units,
              pmu->m_hists.m_bins_array, &CVUniverse::GetPmuTrue, is_true);

  Var* wexp_true =
      new Var("wexp_true", "W_{exp} True", wexp->m_units,
              wexp->m_hists.m_bins_array, &CVUniverse::GetWexpTrue, is_true);

  std::vector<Var*> variables = {pmu, thetamu_deg, enu, q2, wexp, ehad};

  if (include_truth_vars) {
    variables.push_back(pmu_true);
    variables.push_back(thetamu_deg_true);
    variables.push_back(enu_true);
    variables.push_back(q2_true);
    variables.push_back(wexp_true);
    variables.push_back(ehad_true);
  }

  return variables;
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
    case kInclusive:
      variables =
          make_xsec_mc_inputs::GetInclusiveVariables(include_truth_vars);
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

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillMCXSecInputs(const CCPi::MacroUtil& util,
                             const EDataMCTruth& type,
                             std::vector<Variable*>& variables,
                             std::vector<VariableMAT*>& variables_MAT,
                             std::vector<Variable2D*>& variables2D,
                             bool truncate_run) {
  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);
  const UniverseMap error_bands =
      is_truth ? util.m_error_bands_truth : util.m_error_bands;
  Long64_t n_entries_to_loop_over = truncate_run ? 10000 : n_entries;
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

        // calls GetWeight
        // LowNuHighNuEvent event(is_mc, is_truth, util.m_signal_definition,
        // universe);
        InclusiveEvent event(is_mc, is_truth, util.m_signal_definition,
                             universe);

        //===============
        // FILL TRUTH
        //===============
        if (type == kTruth) {
          // lownuhighnu_event::FillTruthEvent(event, variables);
          inclusive_event::FillTruthEvent(event, variables);
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

        // RDF: Still need this in its current location? Or can be moved
        // upstream? Need to re-call this because the node cut efficiency
        // systematic needs a pion candidate to calculate its weight.
        event.m_weight = universe->GetWeight();

        //===============
        // FILL RECO
        //===============
        // lownuhighnu_event::FillRecoEvent(event, variables, variables_MAT, variables2D);
        inclusive_event::FillRecoEvent(event, variables);
      }  // universes
    }    // error bands
  }      // events
  std::cout << "*** Done ***\n\n";
}

//==============================================================================
// Main
//==============================================================================
void makeCrossSectionMCInputs(int signal_definition_int = 0,
                              std::string plist = "ME1A",
                              bool do_systematics = false,
                              bool do_truth = false, bool is_grid = false,
                              std::string input_file = "", int run = 0,
                              bool test_run = false) {
  // INPUT TUPLES
  const bool is_mc = true;
  std::string mc_file_list;
  assert(!(is_grid && input_file.empty()) &&
         "On the grid, infile must be specified.");
  // const bool use_xrootd = false;
  mc_file_list = input_file.empty()
                     ? GetPlaylistFile(plist, is_mc /*, use_xrootd*/)
                     : input_file;

  // INIT MACRO UTILITY
  const std::string macro("MCXSecInputs");
  CCPi::MacroUtil util(signal_definition_int, mc_file_list, plist, do_truth,
                       is_grid, do_systematics);
  util.PrintMacroConfiguration(macro);

  // INIT OUTPUT
  auto time = std::time(nullptr);
  char tchar[100];
  std::strftime(tchar, sizeof(tchar), "%F", std::gmtime(&time));  // YYYY-MM-dd
  const std::string tag = tchar;
  std::string outfile_name(Form("%s_%d%d%d%d_%s_%d_%s.root", macro.c_str(),
                                signal_definition_int, int(do_systematics),
                                int(do_truth), int(is_grid), plist.c_str(), run,
                                tag.c_str()));
  std::cout << "Saving output to " << outfile_name << "\n\n";
  TFile fout(outfile_name.c_str(), "RECREATE");

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

  for (auto v : variables_MAT) v->InitializeAllHists(util.m_error_bands);
  for (auto v : variables2D) v->InitializeAllHists(util.m_error_bands);

  // LOOP MC RECO
  for (auto band : util.m_error_bands) {
    std::vector<CVUniverse*> universes = band.second;
    for (auto universe : universes) universe->SetTruth(false);
  }
  LoopAndFillMCXSecInputs(util, kMC, variables, variables_MAT, variables2D, test_run);

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
  fout.cd();
  for (auto v : variables) {
    SyncAllHists(*v);
    v->WriteMCHists(fout);
  }
  for (auto v: variables_MAT) {
    v->WriteAllHistogramsToFile(fout, true);
  }
  for (auto v : variables2D) {
    // SyncAllHists(*v); // to-do follow up on this
    v->WriteAllHistogramsToFile(fout, true);
  }
}

#endif  // makeXsecMCInputs_C
