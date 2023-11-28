//==============================================================================
// Compare various options for calorimetry that factors into calculation of Ehad.
//==============================================================================
#include <iostream>
#include <vector>

#include "includes/Binning.h"
#include "includes/CVUniverse.h"
#include "includes/Cuts.h"
#include "includes/LowNuHighNuEvent.h"
#include "includes/MacroUtil.h"
#include "includes/Variable.h"
#include "includes/Variable2D.h"
#include "includes/common_functions.h"
#include "playlist_methods.h"  // GetPlaylistFile

// Forward declare my variables because we're hiding the header.
class Variable;

namespace run_recoil_study {
typedef Variable Var;
typedef VariableMAT VarMAT;

//==============================================================================
// Helper Functions
//==============================================================================
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

//==============================================================================
// Do some event processing (e.g. make cuts, get best pion) and fill hists
//==============================================================================
void FillVars(LowNuHighNuEvent& event, const std::vector<MnvH1D*>& study_hists,
              const std::vector<Variable*>& variables,
              const std::vector<VariableMAT*>& variables_MAT) {
  const CVUniverse* universe = event.m_universe;
  const double wgt = event.m_weight;
  const bool is_mc = event.m_is_mc;
  const SignalDefinition sd = event.m_signal_definition;

  // No systematics considered here
  if (universe->ShortName() != "cv") return;

  PassesCutsInfo cuts_info = PassesCuts(event);
  event.m_passes_cuts = cuts_info.GetAll();
  if(!event.m_passes_cuts) return;

  double ehad_reco = event.m_universe->GetEhad();
  double ehad_true = event.m_universe->GetEhadTrue();

  if(ehad_true>200&&ehad_true<400){ 
    study_hists[0]->Fill(ehad_reco,event.m_weight);
  }
  else if(ehad_true>400&&ehad_true<600){
    study_hists[1]->Fill(ehad_reco,event.m_weight);
  }

  for (auto var: variables_MAT) {
    // Get fill value
    double fill_val_reco = var->GetRecoValue(*event.m_universe);
    double fill_val_truth = var->GetTrueValue(*event.m_universe);

    // all MC events
    if (event.m_is_mc) {
      var->m_migration.FillUniverse(*event.m_universe, fill_val_reco,
                                    fill_val_truth, event.m_weight);
    }
  } //end variables_MAT
}

//==============================================================================
// Get Variables
//==============================================================================
std::vector<Variable*> GetVariables() {
  Var* ehad = new Var("ehad", "ehad", "MeV", CCPi::GetBinning("ehad"),
                     &CVUniverse::GetEhad);
  return std::vector<Var*>{ehad};
}

std::vector<VariableMAT*> GetMATVariables() {
  VarMAT* vmat_ehad = new VarMAT("vmat_ehad", "vmat_ehad", ConvertTArrayDToStdVector(CCPi::GetBinning("ehad_fine")),
                     &CVUniverse::GetEhad, &CVUniverse::GetEhadTrue);
  return std::vector<VarMAT*>{vmat_ehad};
}

//==============================================================================
// Get Study-Specific Histograms 
//==============================================================================
std::vector<MnvH1D*> GetStudyHists() {
  // reco-true/reco, etc.

  int n_bins = CCPi::GetBinning("ehad_fine").GetSize()-1;

  MnvH1D* ehad_reco_residual_03 = new MnvH1D("ehad_reco_residual_03", "ehad_reco_residual_03", n_bins, CCPi::GetBinning("ehad_fine").GetArray());
  MnvH1D* ehad_reco_residual_05 = new MnvH1D("ehad_reco_residual_05", "ehad_reco_residual_05", n_bins, CCPi::GetBinning("ehad_fine").GetArray());

  return std::vector<MnvH1D*>{ehad_reco_residual_03,ehad_reco_residual_05};
}

}  // namespace run_recoil_study

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFill(const CCPi::MacroUtil& util, CVUniverse* universe,
                 const EDataMCTruth& type, std::vector<MnvH1D*>& study_hists,
                 std::vector<Variable*>& variables,
                 std::vector<VariableMAT*>& variables_MAT) {
  std::cout << "Loop and Fill CutVars\n";
  bool is_mc, is_truth;
  Long64_t n_entries;
  SetupLoop(type, util, is_mc, is_truth, n_entries);

  for (Long64_t i_event = 0; i_event < n_entries; ++i_event) {
    if (i_event % 500000 == 0)
      std::cout << (i_event / 1000) << "k " << std::endl;
    universe->SetEntry(i_event);

    // For mc, get weight, check signal, and sideband
    LowNuHighNuEvent event(is_mc, is_truth, util.m_signal_definition, universe);

    // WRITE THE FILL FUNCTION
    run_recoil_study::FillVars(event, study_hists, variables, variables_MAT);
  }  // events
  std::cout << "*** Done ***\n\n";
}

//==============================================================================
// Main
//==============================================================================
void runRecoilEnergyStudy(std::string plist = "ME1L", std::string input_file = "") {

  //=========================================
  // Input tuples
  //=========================================
  bool is_mc = true;
  std::string mc_file_list;
  // const bool use_xrootd = false;
  mc_file_list = input_file.empty()
                     ? GetPlaylistFile(plist, is_mc)
                     : input_file;

  //=========================================
  // Init macro utility
  //=========================================
  const int signal_definition_int = 0;
  const std::string macro("runStudyTemplate");
  const bool is_grid = false;
  const bool do_truth = false;
  const bool do_systematics = false;

  CCPi::MacroUtil util(signal_definition_int, mc_file_list, 
                       plist, do_truth, is_grid, do_systematics);
  util.PrintMacroConfiguration(macro);

  //=========================================
  // Init output
  //=========================================
  auto time = std::time(nullptr);
  char tchar[100];
  std::strftime(tchar, sizeof(tchar), "%F", std::gmtime(&time));  // YYYY-MM-dd
  const std::string tag = tchar;
  std::string outfile_name(Form("%s_%d%d%d%d_%s_%s.root", macro.c_str(),
                                signal_definition_int, int(do_systematics),
                                int(do_truth), int(is_grid), plist.c_str(),
                                tag.c_str()));
  std::cout << "Saving output to " << outfile_name << "\n\n";
  TFile fout(outfile_name.c_str(), "RECREATE");

  //=========================================
  // Get variables and initialize their hists
  //=========================================
  std::vector<MnvH1D*> study_hists = run_recoil_study::GetStudyHists();
  std::vector<Variable*> variables = run_recoil_study::GetVariables();
  std::vector<VariableMAT*> variables_MAT = run_recoil_study::GetMATVariables();
  for (auto v : variables)
    v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);
  for (auto v : variables_MAT) v->InitializeAllHists(util.m_error_bands);

  //=========================================
  // Loop and Fill
  //=========================================
  LoopAndFill(util, util.m_data_universe, kData, study_hists, variables, variables_MAT);
  LoopAndFill(util, util.m_error_bands.at("cv").at(0), kMC, study_hists, variables, variables_MAT);

  //=========================================
  // Write to file
  //=========================================
  std::cout << "Synching and Writing\n\n";
  WritePOT(fout, true, util.m_mc_pot);
  fout.cd();
  //for (auto v : variables) {
  //  v->WriteMCHists(fout);
  //}
  for (auto v: variables_MAT) {
    v->WriteAllHistogramsToFile(fout, true);
  }

  for (auto h: study_hists) {
    h->Write();
  }

  std::cout << "Success" << std::endl;

}
