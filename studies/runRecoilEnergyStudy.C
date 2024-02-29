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
#include "PlotUtils/LowRecoilPionReco.h" // Cluster struct
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

  double ehad_true = event.m_universe->GetEhadTrue();
  double ehad_reco = event.m_universe->GetEhad();
  double ehad_residual = (ehad_reco-ehad_true)/ehad_true;
  double ehad_reco_lessOD = event.m_universe->GetEhadLessOD();
  double ehad_residual_lessOD = (ehad_reco_lessOD-ehad_true)/ehad_true;

  if(ehad_true>200&&ehad_true<400){ 
    study_hists[0]->Fill(ehad_reco,event.m_weight);
    study_hists[1]->Fill(ehad_reco_lessOD,event.m_weight);
    study_hists[2]->Fill(ehad_residual,event.m_weight);
    study_hists[3]->Fill(ehad_residual/100.,event.m_weight);
    study_hists[4]->Fill(ehad_residual_lessOD,event.m_weight);
    study_hists[5]->Fill(ehad_residual_lessOD/100.,event.m_weight);
  }
  else if(ehad_true>400&&ehad_true<600){
    study_hists[6]->Fill(ehad_reco,event.m_weight);
    study_hists[7]->Fill(ehad_reco_lessOD,event.m_weight);
    study_hists[8]->Fill(ehad_residual,event.m_weight);
    study_hists[9]->Fill(ehad_residual/100.,event.m_weight);
    study_hists[10]->Fill(ehad_residual_lessOD,event.m_weight);
    study_hists[11]->Fill(ehad_residual_lessOD/100.,event.m_weight);
  }
  else if(ehad_true>900&&ehad_true<1100){
    study_hists[12]->Fill(ehad_reco,event.m_weight);
    study_hists[13]->Fill(ehad_reco_lessOD,event.m_weight);
    study_hists[14]->Fill(ehad_residual,event.m_weight);
    study_hists[15]->Fill(ehad_residual/100.,event.m_weight);
    study_hists[16]->Fill(ehad_residual_lessOD,event.m_weight);
    study_hists[17]->Fill(ehad_residual_lessOD/100.,event.m_weight);
  }
  else if(ehad_true>1900&&ehad_true<2100){
    study_hists[18]->Fill(ehad_reco,event.m_weight);
    study_hists[19]->Fill(ehad_reco_lessOD,event.m_weight);
    study_hists[10]->Fill(ehad_residual,event.m_weight);
    study_hists[21]->Fill(ehad_residual/100.,event.m_weight);
    study_hists[22]->Fill(ehad_residual_lessOD,event.m_weight);
    study_hists[23]->Fill(ehad_residual_lessOD/100.,event.m_weight);
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

  // Loop over clusters
  int nclusters = universe->GetInt("cluster_view_sz");
  double interaction_time = universe->GetVecElem("vtx",3);///pow(10,3);
  for (int i = 0; i < nclusters; i++) {

    // Grab ith cluster
    LowRecoilPion::Cluster clus(*universe,i);

    // Some of the clusters have energy = nan (don't know why). Skip these clusters 
    if(clus.energy!=clus.energy){
      continue;
    }

    study_hists[24]->Fill(clus.time*pow(10,3)-interaction_time,clus.energy); // weight by cluster energy because low-energy clusters have bad timing reco
    study_hists[25]->Fill(clus.energy,event.m_weight);

    if(clus.subdet==1){ // Nuclear Targets
      study_hists[26]->Fill(clus.time*pow(10,3)-interaction_time,clus.energy);
      study_hists[27]->Fill(clus.energy,event.m_weight);
    }
    else if(clus.subdet==2){ // Tracker
      study_hists[28]->Fill(clus.time*pow(10,3)-interaction_time,clus.energy);
      study_hists[29]->Fill(clus.energy,event.m_weight);
    } 
    else if(clus.subdet==3){ // ECAL
      study_hists[30]->Fill(clus.time*pow(10,3)-interaction_time,clus.energy);
      study_hists[31]->Fill(clus.energy,event.m_weight);
    } 
    else if(clus.subdet==4){ // HCAL 
      study_hists[32]->Fill(clus.time*pow(10,3)-interaction_time,clus.energy);
      study_hists[33]->Fill(clus.energy,event.m_weight);
    }
    else{
      std::cout << "Maybe I don't understand clus.subdet, because it's value here is: " << clus.subdet << std::endl;
    }

  }
}

void FillVarsData(LowNuHighNuEvent& event, const std::vector<MnvH1D*>& study_hists) {
  const CVUniverse* universe = event.m_universe;

  // No systematics considered here
  if (universe->ShortName() != "cv") return;

  PassesCutsInfo cuts_info = PassesCuts(event);
  event.m_passes_cuts = cuts_info.GetAll();
  if(!event.m_passes_cuts) return;

  // Loop over clusters
  int nclusters = universe->GetInt("cluster_view_sz");
  double interaction_time = universe->GetVecElem("vtx",3);///pow(10,3);
  for (int i = 0; i < nclusters; i++) {

    // Grab ith cluster
    LowRecoilPion::Cluster clus(*universe,i);

    // Some of the clusters have energy = nan (don't know why). Skip these clusters 
    if(clus.energy!=clus.energy){
      continue;
    }

    study_hists[0]->Fill(clus.time*pow(10,3)-interaction_time,clus.energy); // weight by cluster energy because low-energy clusters have bad timing reco
    study_hists[1]->Fill(clus.energy,event.m_weight);

    if(clus.subdet==1){ // Nuclear Targets
      study_hists[2]->Fill(clus.time*pow(10,3)-interaction_time,clus.energy);
      study_hists[3]->Fill(clus.energy,event.m_weight);
    }
    else if(clus.subdet==2){ // Tracker
      study_hists[4]->Fill(clus.time*pow(10,3)-interaction_time,clus.energy);
      study_hists[5]->Fill(clus.energy,event.m_weight);
    } 
    else if(clus.subdet==3){ // ECAL
      study_hists[6]->Fill(clus.time*pow(10,3)-interaction_time,clus.energy);
      study_hists[7]->Fill(clus.energy,event.m_weight);
    } 
    else if(clus.subdet==4){ // HCAL 
      study_hists[8]->Fill(clus.time*pow(10,3)-interaction_time,clus.energy);
      study_hists[9]->Fill(clus.energy,event.m_weight);
    }
    else{
      std::cout << "Maybe I don't understand clus.subdet, because it's value here is: " << clus.subdet << std::endl;
    }

  }
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
  VarMAT* vmat_ehad_lessOD = new VarMAT("vmat_ehad_lessOD", "vmat_ehad_lessOD", ConvertTArrayDToStdVector(CCPi::GetBinning("ehad_fine")),
                             &CVUniverse::GetEhadLessOD, &CVUniverse::GetEhadTrue);
  return std::vector<VarMAT*>{vmat_ehad,vmat_ehad_lessOD};
}

//==============================================================================
// Get Study-Specific Histograms 
//==============================================================================
std::vector<MnvH1D*> GetStudyHists() {
  int n_bins_ehad_reco = CCPi::GetBinning("ehad_fine").GetSize()-1;
  int n_bins_ehad_res = CCPi::GetBinning("ehad_res").GetSize()-1;
  int n_bins_ehad_res_zoom = CCPi::GetBinning("ehad_res_zoom").GetSize()-1;
  int n_bins_cluster_timing = CCPi::GetBinning("cluster_timing").GetSize()-1;
  int n_bins_cluster_energy = CCPi::GetBinning("cluster_energy").GetSize()-1;

  MnvH1D* ehad_reco_03 = new MnvH1D("ehad_reco_03", "ehad_reco_03", n_bins_ehad_reco, CCPi::GetBinning("ehad_fine").GetArray());
  MnvH1D* ehad_reco_lessOD_03 = new MnvH1D("ehad_reco_lessOD_03", "ehad_reco_lessOD_03", n_bins_ehad_reco, CCPi::GetBinning("ehad_fine").GetArray());
  MnvH1D* ehad_residual_03 = new MnvH1D("ehad_residual_03", "ehad_residual_03", n_bins_ehad_res, CCPi::GetBinning("ehad_res").GetArray());
  MnvH1D* ehad_residual_03_zoom = new MnvH1D("ehad_residual_zoom_03", "ehad_residual_zoom_03", n_bins_ehad_res_zoom, CCPi::GetBinning("ehad_res_zoom").GetArray());
  MnvH1D* ehad_residual_lessOD_03 = new MnvH1D("ehad_residual_lessOD_03", "ehad_residual_lessOD_03", n_bins_ehad_res, CCPi::GetBinning("ehad_res").GetArray());
  MnvH1D* ehad_residual_lessOD_03_zoom = new MnvH1D("ehad_residual_lessOD_zoom_03", "ehad_residual_lessOD_zoom_03", n_bins_ehad_res_zoom, CCPi::GetBinning("ehad_res_zoom").GetArray());

  MnvH1D* ehad_reco_05 = new MnvH1D("ehad_reco_05", "ehad_reco_05", n_bins_ehad_reco, CCPi::GetBinning("ehad_fine").GetArray());
  MnvH1D* ehad_reco_lessOD_05 = new MnvH1D("ehad_reco_lessOD_05", "ehad_reco_lessOD_05", n_bins_ehad_reco, CCPi::GetBinning("ehad_fine").GetArray());
  MnvH1D* ehad_residual_05 = new MnvH1D("ehad_residual_05", "ehad_residual_05", n_bins_ehad_res, CCPi::GetBinning("ehad_res").GetArray());
  MnvH1D* ehad_residual_05_zoom = new MnvH1D("ehad_residual_zoom_05", "ehad_residual_zoom_05", n_bins_ehad_res_zoom, CCPi::GetBinning("ehad_res_zoom").GetArray());
  MnvH1D* ehad_residual_lessOD_05 = new MnvH1D("ehad_residual_lessOD_05", "ehad_residual_lessOD_05", n_bins_ehad_res, CCPi::GetBinning("ehad_res").GetArray());
  MnvH1D* ehad_residual_lessOD_05_zoom = new MnvH1D("ehad_residual_lessOD_zoom_05", "ehad_residual_lessOD_zoom_05", n_bins_ehad_res_zoom, CCPi::GetBinning("ehad_res_zoom").GetArray());

  MnvH1D* ehad_reco_10 = new MnvH1D("ehad_reco_10", "ehad_reco_10", n_bins_ehad_reco, CCPi::GetBinning("ehad_fine").GetArray());
  MnvH1D* ehad_reco_lessOD_10 = new MnvH1D("ehad_reco_lessOD_10", "ehad_reco_lessOD_10", n_bins_ehad_reco, CCPi::GetBinning("ehad_fine").GetArray());
  MnvH1D* ehad_residual_10 = new MnvH1D("ehad_residual_10", "ehad_residual_10", n_bins_ehad_res, CCPi::GetBinning("ehad_res").GetArray());
  MnvH1D* ehad_residual_10_zoom = new MnvH1D("ehad_residual_zoom_10", "ehad_residual_zoom_10", n_bins_ehad_res_zoom, CCPi::GetBinning("ehad_res_zoom").GetArray());
  MnvH1D* ehad_residual_lessOD_10 = new MnvH1D("ehad_residual_lessOD_10", "ehad_residual_lessOD_10", n_bins_ehad_res, CCPi::GetBinning("ehad_res").GetArray());
  MnvH1D* ehad_residual_lessOD_10_zoom = new MnvH1D("ehad_residual_lessOD_zoom_10", "ehad_residual_lessOD_zoom_10", n_bins_ehad_res_zoom, CCPi::GetBinning("ehad_res_zoom").GetArray());

  MnvH1D* ehad_reco_20 = new MnvH1D("ehad_reco_20", "ehad_reco_20", n_bins_ehad_reco, CCPi::GetBinning("ehad_fine").GetArray());
  MnvH1D* ehad_reco_lessOD_20 = new MnvH1D("ehad_reco_lessOD_20", "ehad_reco_lessOD_20", n_bins_ehad_reco, CCPi::GetBinning("ehad_fine").GetArray());
  MnvH1D* ehad_residual_20 = new MnvH1D("ehad_residual_20", "ehad_residual_20", n_bins_ehad_res, CCPi::GetBinning("ehad_res").GetArray());
  MnvH1D* ehad_residual_20_zoom = new MnvH1D("ehad_residual_zoom_20", "ehad_residual_zoom_20", n_bins_ehad_res_zoom, CCPi::GetBinning("ehad_res_zoom").GetArray());
  MnvH1D* ehad_residual_lessOD_20 = new MnvH1D("ehad_residual_lessOD_20", "ehad_residual_lessOD_20", n_bins_ehad_res, CCPi::GetBinning("ehad_res").GetArray());
  MnvH1D* ehad_residual_lessOD_20_zoom = new MnvH1D("ehad_residual_lessOD_zoom_20", "ehad_residual_lessOD_zoom_20", n_bins_ehad_res_zoom, CCPi::GetBinning("ehad_res_zoom").GetArray());

  MnvH1D* cluster_timing = new MnvH1D("cluster_timing", "cluster_timing", n_bins_cluster_timing, CCPi::GetBinning("cluster_timing").GetArray());
  MnvH1D* cluster_energy = new MnvH1D("cluster_energy", "cluster_energy", n_bins_cluster_energy, CCPi::GetBinning("cluster_energy").GetArray());
  MnvH1D* cluster_timing_nuclear_targets = new MnvH1D("cluster_timing_nuclear_targets", "cluster_timing_nuclear_targets", n_bins_cluster_timing, CCPi::GetBinning("cluster_timing").GetArray());
  MnvH1D* cluster_energy_nuclear_targets = new MnvH1D("cluster_energy_nuclear_targets", "cluster_energy_nuclear_targets", n_bins_cluster_energy, CCPi::GetBinning("cluster_energy").GetArray());
  MnvH1D* cluster_timing_tracker = new MnvH1D("cluster_timing_tracker", "cluster_timing_tracker", n_bins_cluster_timing, CCPi::GetBinning("cluster_timing").GetArray());
  MnvH1D* cluster_energy_tracker = new MnvH1D("cluster_energy_tracker", "cluster_energy_tracker", n_bins_cluster_energy, CCPi::GetBinning("cluster_energy").GetArray());
  MnvH1D* cluster_timing_ECAL = new MnvH1D("cluster_timing_ECAL", "cluster_timing_ECAL", n_bins_cluster_timing, CCPi::GetBinning("cluster_timing").GetArray());
  MnvH1D* cluster_energy_ECAL = new MnvH1D("cluster_energy_ECAL", "cluster_energy_ECAL", n_bins_cluster_energy, CCPi::GetBinning("cluster_energy").GetArray());
  MnvH1D* cluster_timing_HCAL = new MnvH1D("cluster_timing_HCAL", "cluster_timing_HCAL", n_bins_cluster_timing, CCPi::GetBinning("cluster_timing").GetArray());
  MnvH1D* cluster_energy_HCAL = new MnvH1D("cluster_energy_HCAL", "cluster_energy_HCAL", n_bins_cluster_energy, CCPi::GetBinning("cluster_energy").GetArray());
  
  return std::vector<MnvH1D*>{
    ehad_reco_03,
    ehad_reco_lessOD_03,
    ehad_residual_03,
    ehad_residual_03_zoom,
    ehad_residual_lessOD_03,
    ehad_residual_lessOD_03_zoom,
    ehad_reco_05,
    ehad_reco_lessOD_05,
    ehad_residual_05,
    ehad_residual_05_zoom,
    ehad_residual_lessOD_05,
    ehad_residual_lessOD_05_zoom,
    ehad_reco_10,
    ehad_reco_lessOD_10,
    ehad_residual_10,
    ehad_residual_10_zoom,
    ehad_residual_lessOD_10,
    ehad_residual_lessOD_10_zoom,
    ehad_reco_20,
    ehad_reco_lessOD_20,
    ehad_residual_20,
    ehad_residual_20_zoom,
    ehad_residual_lessOD_20,
    ehad_residual_lessOD_20_zoom,
    cluster_timing,
    cluster_energy,
    cluster_timing_nuclear_targets,
    cluster_energy_nuclear_targets,
    cluster_timing_tracker,
    cluster_energy_tracker,
    cluster_timing_ECAL,
    cluster_energy_ECAL,
    cluster_timing_HCAL,
    cluster_energy_HCAL
  };
}

std::vector<MnvH1D*> GetStudyHistsData() {
  int n_bins_cluster_timing_data = CCPi::GetBinning("cluster_timing").GetSize()-1;
  int n_bins_cluster_energy_data = CCPi::GetBinning("cluster_energy").GetSize()-1;

  MnvH1D* cluster_timing_data = new MnvH1D("cluster_timing_data", "cluster_timing_data", n_bins_cluster_timing_data, CCPi::GetBinning("cluster_timing").GetArray());
  MnvH1D* cluster_energy_data = new MnvH1D("cluster_energy_data", "cluster_energy_data", n_bins_cluster_energy_data, CCPi::GetBinning("cluster_energy").GetArray());
  MnvH1D* cluster_timing_nuclear_targets_data = new MnvH1D("cluster_timing_nuclear_targets_data", "cluster_timing_nuclear_targets_data", n_bins_cluster_timing_data, CCPi::GetBinning("cluster_timing").GetArray());
  MnvH1D* cluster_energy_nuclear_targets_data = new MnvH1D("cluster_energy_nuclear_targets_data", "cluster_energy_nuclear_targets_data", n_bins_cluster_energy_data, CCPi::GetBinning("cluster_energy").GetArray());
  MnvH1D* cluster_timing_tracker_data = new MnvH1D("cluster_timing_tracker", "cluster_timing_tracker", n_bins_cluster_timing_data, CCPi::GetBinning("cluster_timing").GetArray());
  MnvH1D* cluster_energy_tracker_data = new MnvH1D("cluster_energy_tracker", "cluster_energy_tracker", n_bins_cluster_energy_data, CCPi::GetBinning("cluster_energy").GetArray());
  MnvH1D* cluster_timing_ECAL_data = new MnvH1D("cluster_timing_ECAL_data", "cluster_timing_ECAL_data", n_bins_cluster_timing_data, CCPi::GetBinning("cluster_timing").GetArray());
  MnvH1D* cluster_energy_ECAL_data = new MnvH1D("cluster_energy_ECAL_data", "cluster_energy_ECAL_data", n_bins_cluster_energy_data, CCPi::GetBinning("cluster_energy").GetArray());
  MnvH1D* cluster_timing_HCAL_data = new MnvH1D("cluster_timing_HCAL_data", "cluster_timing_HCAL_data", n_bins_cluster_timing_data, CCPi::GetBinning("cluster_timing").GetArray());
  MnvH1D* cluster_energy_HCAL_data = new MnvH1D("cluster_energy_HCAL_data", "cluster_energy_HCAL_data", n_bins_cluster_energy_data, CCPi::GetBinning("cluster_energy").GetArray());
  
  return std::vector<MnvH1D*>{
    cluster_timing_data,
    cluster_energy_data,
    cluster_timing_nuclear_targets_data,
    cluster_energy_nuclear_targets_data,
    cluster_timing_tracker_data,
    cluster_energy_tracker_data,
    cluster_timing_ECAL_data,
    cluster_energy_ECAL_data,
    cluster_timing_HCAL_data,
    cluster_energy_HCAL_data
  };
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
    if (is_mc){
      run_recoil_study::FillVars(event, study_hists, variables, variables_MAT);
    }
    else{
      run_recoil_study::FillVarsData(event, study_hists);
    }
  }  // events
  std::cout << "*** Done ***\n\n";
}

//==============================================================================
// Main
//==============================================================================
void runRecoilEnergyStudy(std::string plist = "ME1L", std::string input_file = "", std::string input_file_data = "") {

  //=========================================
  // Input tuples
  //=========================================
  bool is_mc = true;
  std::string mc_file_list;
  // const bool use_xrootd = false;
  mc_file_list = input_file.empty()
                     ? GetPlaylistFile(plist, is_mc)
                     : input_file;
  std::string data_file_list;
  data_file_list = input_file_data.empty()
                     ? GetPlaylistFile(plist,false)
                     : input_file_data;

  //=========================================
  // Init macro utility
  //=========================================
  const int signal_definition_int = 0;
  const std::string macro("RecoilEnergyStudy");
  const bool is_grid = false;
  const bool do_truth = false;
  const bool do_systematics = false;

  CCPi::MacroUtil util(signal_definition_int, mc_file_list, data_file_list,
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
  std::vector<MnvH1D*> study_hists_data = run_recoil_study::GetStudyHistsData();
  std::vector<Variable*> variables = run_recoil_study::GetVariables();
  std::vector<VariableMAT*> variables_MAT = run_recoil_study::GetMATVariables();
  for (auto v : variables)
    v->InitializeAllHists(util.m_error_bands, util.m_error_bands_truth);
  for (auto v : variables_MAT) v->InitializeAllHists(util.m_error_bands);

  //=========================================
  // Loop and Fill
  //=========================================
  LoopAndFill(util, util.m_data_universe, kData, study_hists_data, variables, variables_MAT);
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

  for (auto h: study_hists_data) {
    h->Write();
  }

  std::cout << "Success" << std::endl;

}
