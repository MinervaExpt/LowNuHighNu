#ifndef InclusiveEvent_cxx
#define InclusiveEvent_cxx

#include "InclusiveEvent.h"

#include "Cuts.h"    // kCutsVector
#include "Michel.h"  // class endpoint::Michel, typdef endpoint::MichelMap, endpoint::GetQualityMichels
#include "MichelTrackless.h"
#include "common_functions.h"  // GetVar, HasVar

//==============================================================================
// CTOR
//==============================================================================
InclusiveEvent::InclusiveEvent(const bool is_mc, const bool is_truth,
                     const SignalDefinition signal_definition,
                     CVUniverse* universe)
    : m_is_mc(is_mc),
      m_is_truth(is_truth),
      m_signal_definition(signal_definition),
      m_universe(universe),
      m_reco_pion_candidate_idxs(),
      m_highest_energy_pion_idx(-300) {
  m_is_signal = is_mc ? IsSignal(*universe, signal_definition) : false;
  m_weight = is_mc ? universe->GetWeight() : 1.;
  m_w_type = is_mc ? GetWSidebandType(*universe, signal_definition,
                                      sidebands::kNWFitCategories)
                   : kNWSidebandTypes;
}

//==============================================================================
// Helper Functions
//==============================================================================
// PassesCutsInfo {passes_all_cuts, is_w_sideband, passes_all_except_w, pion_candidate_idxs}
PassesCutsInfo PassesCuts(const InclusiveEvent& e) {
  return PassesCuts(*e.m_universe, e.m_is_mc, e.m_signal_definition);
}

SignalBackgroundType GetSignalBackgroundType(const InclusiveEvent& e) {
  return GetSignalBackgroundType(*e.m_universe, e.m_signal_definition);
}

RecoPionIdx GetHighestEnergyPionCandidateIndex(const InclusiveEvent& e) {
  return e.m_universe->GetHighestEnergyPionCandidateIndex(
      e.m_reco_pion_candidate_idxs);
}

TruePionIdx GetHighestEnergyTruePionIndex(const InclusiveEvent& e) {
  return e.m_universe->GetHighestEnergyTruePionIndex();
}

//==============================================================================
// Fill all histos for an entire event -- call other specialized fill functions
//==============================================================================
void inclusive_event::FillRecoEvent(const InclusiveEvent& event,
                               const std::vector<Variable*>& variables) {
  // Fill selection -- total, signal-only, and bg-only
  if (event.m_passes_cuts) {
    inclusive_event::FillSelected(event, variables);
  }
  // // Fill W Sideband
  // if (event.m_is_w_sideband) {
  //   inclusive_event::FillWSideband(event, variables);
  // }

  // // Fill W Sideband Study
  // if (event.m_passes_all_cuts_except_w && event.m_universe->ShortName() == "cv") {
  //   inclusive_event::FillWSideband_Study(event, variables);
  // }

  // Fill Migration
  if (event.m_is_mc && event.m_is_signal && event.m_passes_cuts) {
    if (HasVar(variables, "tpi") && HasVar(variables, "tpi_true"))
      FillMigration(event, variables, std::string("tpi"));
    if (HasVar(variables, "thetapi_deg") &&
        HasVar(variables, "thetapi_deg_true"))
      FillMigration(event, variables, std::string("thetapi_deg"));
    if (HasVar(variables, "pmu") && HasVar(variables, "pmu_true"))
      FillMigration(event, variables, std::string("pmu"));
    if (HasVar(variables, "pzmu") && HasVar(variables, "pzmu_true"))
      FillMigration(event, variables, std::string("pzmu"));
    if (HasVar(variables, "ptmu") && HasVar(variables, "ptmu_true"))
      FillMigration(event, variables, std::string("ptmu"));
    if (HasVar(variables, "thetamu_deg") &&
        HasVar(variables, "thetamu_deg_true"))
      FillMigration(event, variables, std::string("thetamu_deg"));
    if (HasVar(variables, "q2") && HasVar(variables, "q2_true"))
      FillMigration(event, variables, std::string("q2"));
    if (HasVar(variables, "enu") && HasVar(variables, "enu_true"))
      FillMigration(event, variables, std::string("enu"));
    if (HasVar(variables, "wexp") && HasVar(variables, "wexp_true"))
      FillMigration(event, variables, std::string("wexp"));
    if (HasVar(variables, "ehad") && HasVar(variables, "ehad_true"))
      FillMigration(event, variables, std::string("ehad"));
    if (HasVar(variables, "cosadtheta") && HasVar(variables, "cosadtheta_true"))
      FillMigration(event, variables, std::string("cosadtheta"));
    if (HasVar(variables, "adphi") && HasVar(variables, "adphi_true"))
      FillMigration(event, variables, std::string("adphi"));
    if (HasVar(variables, "pimuAngle") && HasVar(variables, "pimuAngle_true"))
      FillMigration(event, variables, std::string("pimuAngle"));
    if (HasVar(variables, "PT") && HasVar(variables, "PT_true"))
      FillMigration(event, variables, std::string("PT"));
  }
}

void inclusive_event::FillTruthEvent(const InclusiveEvent& event,
                                const std::vector<Variable*>& variables) {
  // Fill Efficiency Denominator
  if (event.m_is_signal)
    inclusive_event::FillEfficiencyDenominator(event, variables);
}

//==============================================================================
// Specialized fill functions -- for xsec calculation
//==============================================================================
// Fill histos with selected (i.e. passes_cuts) events:
// ** sig + bg (true and reco vars, data and mc)
// ** signal only (true vars, for eff num & closure)
// ** bg only (reco and true vars)
void inclusive_event::FillSelected(const InclusiveEvent& event,
                              const std::vector<Variable*>& variables) {
  for (auto var : variables) {
    // Sanity Checks
    if (var->m_is_true && !event.m_is_mc) return;  // truth, but not MC?

    // Get fill value
    double fill_val = var->GetValue(*event.m_universe);

    // total = signal & background, together
    if (event.m_is_mc) {
      var->m_hists.m_selection_mc.FillUniverse(*event.m_universe, fill_val,
                                               event.m_weight);
    } else {
      var->m_hists.m_selection_data->Fill(fill_val);
    }

    // done with data
    if (!event.m_is_mc) continue;

    // signal and background individually
    if (event.m_is_signal) {
      var->m_hists.m_effnum.FillUniverse(*event.m_universe, fill_val,
                                         event.m_weight);
    } else {
      var->m_hists.m_bg.FillUniverse(*event.m_universe, fill_val,
                                     event.m_weight);

      // Fill bg by W sideband category
      switch (event.m_w_type) {
        case kWSideband_Signal:
          break;
        case kWSideband_Low:
          var->m_hists.m_bg_loW.FillUniverse(*event.m_universe, fill_val,
                                             event.m_weight);
          break;
        case kWSideband_Mid:
          var->m_hists.m_bg_midW.FillUniverse(*event.m_universe, fill_val,
                                              event.m_weight);
          break;
        case kWSideband_High:
          var->m_hists.m_bg_hiW.FillUniverse(*event.m_universe, fill_val,
                                             event.m_weight);
          break;
        default:
          std::cerr << "FillBackgrounds: no such W category\n";
          std::exit(2);
      }
    }
  }  // end variables
}

// Fill histograms of all variables with events in the sideband region
void inclusive_event::FillWSideband(const InclusiveEvent& event,
                               const std::vector<Variable*>& variables) {
  if (!event.m_is_w_sideband) {
    std::cerr << "FillWSideband Warning: This event is not in the wsideband "
                 "region, are you sure you want to be filling?\n";
  }
  if (!HasVar(variables, sidebands::kFitVarString)) {
    std::cerr << "FillWSideband: variables container is missing fit var\n";
    std::exit(1);
  }
  if (event.m_reco_pion_candidate_idxs.empty()) {
    std::cerr << "FillWSideband: member pion idxs is empty\n";
    std::exit(1);
  }

  const RecoPionIdx idx = event.m_highest_energy_pion_idx;

  for (auto var : variables) {
    // if (var->m_is_true && !event.m_is_mc) continue; // truth, but not MC?
    if (var->m_is_true) continue;  // truth pion variables don't generally work
    const double fill_val = var->GetValue(*event.m_universe, idx);

    if (event.m_is_mc) {
      switch (event.m_w_type) {
        case kWSideband_Signal:
          var->m_hists.m_wsidebandfit_sig.FillUniverse(
              *event.m_universe, fill_val, event.m_weight);
          break;
        case kWSideband_Low:
          var->m_hists.m_wsidebandfit_loW.FillUniverse(
              *event.m_universe, fill_val, event.m_weight);
          break;
        case kWSideband_Mid:
          var->m_hists.m_wsidebandfit_midW.FillUniverse(
              *event.m_universe, fill_val, event.m_weight);
          break;
        case kWSideband_High:
          var->m_hists.m_wsidebandfit_hiW.FillUniverse(
              *event.m_universe, fill_val, event.m_weight);
          break;
        default:
          std::cerr << "FillWSideband: invalid W category\n";
          std::exit(2);
      }
    } else {
      var->m_hists.m_wsidebandfit_data->Fill(fill_val);
    }
  }  // end variables
}

void inclusive_event::FillMigration(const InclusiveEvent& event,
                               const vector<Variable*>& variables,
                               std::string name) {
  Variable* reco_var = GetVar(variables, name);
  Variable* true_var = GetVar(variables, name + string("_true"));
  if (true_var == 0) return;
  RecoPionIdx reco_idx = event.m_highest_energy_pion_idx;
  TruePionIdx true_idx = GetHighestEnergyTruePionIndex(event);
  double reco_fill_val = reco_var->GetValue(*event.m_universe, reco_idx);
  double true_fill_val = true_var->GetValue(*event.m_universe, true_idx);
  reco_var->m_hists.m_migration.FillUniverse(*event.m_universe, reco_fill_val,
                                             true_fill_val, event.m_weight);
}

// Only for true variables
void inclusive_event::FillEfficiencyDenominator(
    const InclusiveEvent& event, const std::vector<Variable*>& variables) {
  for (auto var : variables) {
    if (!var->m_is_true) continue;
    double fill_val = var->GetValue(*event.m_universe);
    try {
      var->m_hists.m_effden.FillUniverse(*event.m_universe, fill_val,
                                         event.m_weight);
    } catch (...) {
      std::cerr << "From inclusive_event::FillEfficiencyDenominator\n";
      std::cerr << "Variable is " << var->Name() << "\n";
      throw;
    }
  }
}

//==============================================================================
// Specialized fill functions -- for studies
//==============================================================================
// Fill Stacked components of the wexp_fit variable without the W cut. For
// visualizing the sideband sample. These hists are filled for study purposes
// only. Other hists owned by this variable are used to perform the fit (those
// are filled in FillWSideband.)
void inclusive_event::FillWSideband_Study(const InclusiveEvent& event,
                                     std::vector<Variable*> variables) {
  if (event.m_universe->ShortName() != "cv") {
    std::cerr << "FillWSideband_Study Warning: you're filling the wexp_fit "
                 "variable w/o the W-cut for a universe other than the CV\n";
  }

  if (!event.m_passes_all_cuts_except_w) {
    std::cerr << "FillWSideband_Study Warning: This event does not pass "
                 "correct cuts, are you sure you want to be filling?\n";
  }

  const RecoPionIdx pion_idx = event.m_highest_energy_pion_idx;

  Variable* var = GetVar(variables, sidebands::kFitVarString);
  double fill_val = var->GetValue(*event.m_universe, pion_idx);
  if (event.m_is_mc) {
    var->GetStackComponentHist(event.m_w_type)
        ->Fill(fill_val, event.m_weight);
  } else {
    var->m_hists.m_wsideband_data->Fill(fill_val);
  }
}

// Like FillCutVars, this function loops through cuts and calls PassesCut.
// Michel containers updated as we go, but thrown away at the end.
void inclusive_event::FillCounters(
    const InclusiveEvent& event,
    const std::pair<EventCount*, EventCount*>& counters) {
  EventCount* signal = counters.first;
  EventCount* bg = event.m_is_mc ? counters.second : nullptr;
  endpoint::MichelMap dummy1;
  trackless::MichelEvent<CVUniverse> dummy2;
  bool pass = true;
  // Purity and efficiency
  for (auto i_cut : kCutsVector) {
    if (event.m_is_truth != IsPrecut(i_cut))
      continue;  // truth loop does precuts

    bool passes_this_cut = true;
    std::tie(passes_this_cut, dummy1, dummy2) =
        PassesCut(*event.m_universe, i_cut, event.m_is_mc,
                  event.m_signal_definition, dummy1, dummy2);

    pass = pass && passes_this_cut;

    if (pass) {
      if (!event.m_is_mc) {
        (*signal)[i_cut] += event.m_weight;
        continue;
      }
      if (event.m_is_signal) {
        (*signal)[i_cut] += event.m_weight;  // eff/pur numerator
      } else {
        (*bg)[i_cut] += event.m_weight;
      }
    }
  }  // cuts
}

std::pair<EventCount, EventCount> inclusive_event::FillCounters(
    const InclusiveEvent& event, const EventCount& s, const EventCount& b) {
  EventCount signal = s;
  EventCount bg = b;

  endpoint::MichelMap endpoint_michels;
  trackless::MichelEvent<CVUniverse> vtx_michels;
  bool pass = true;
  for (auto i_cut : kCutsVector) {
    if (event.m_is_truth != IsPrecut(i_cut)) continue;

    bool passes_this_cut = true;
    std::tie(passes_this_cut, endpoint_michels, vtx_michels) =
        PassesCut(*event.m_universe, i_cut, event.m_is_mc,
                  event.m_signal_definition, endpoint_michels, vtx_michels);

    event.m_universe->SetPionCandidates(
        GetHadIdxsFromMichels(endpoint_michels, vtx_michels));

    pass = pass && passes_this_cut;

    if (!pass) break;

    if (!event.m_is_mc) {
      signal[i_cut] += event.m_weight;  // selected data
    } else {
      if (event.m_is_signal)
        signal[i_cut] += event.m_weight;  // selected mc signal
      else
        bg[i_cut] += event.m_weight;  // selected mc bg
    }
  }  // cuts loop
  return {signal, bg};
}

void inclusive_event::FillCutVars(InclusiveEvent& event,
                             const std::vector<Variable*>& variables) {
  const CVUniverse* universe = event.m_universe;
  const double wgt = event.m_weight;
  const bool is_mc = event.m_is_mc;
  const SignalDefinition sd = event.m_signal_definition;

  if (universe->ShortName() != "cv") return;

  endpoint::MichelMap endpoint_michels;
  endpoint_michels.clear();

  trackless::MichelEvent<CVUniverse> vtx_michels;
  // vtx_michels.clear();

  // loop cuts
  bool pass = true;
  for (unsigned int i = 0; i < kCutsVector.size(); ++i) {
    ECuts cut = (ECuts)kCutsVector[i];
    ECuts next_cut;
    try {
      next_cut = (ECuts)(kCutsVector[i + 1]);
    } catch (const std::out_of_range& e) {
      next_cut = (ECuts)(-1);
    }
    event.m_reco_pion_candidate_idxs.clear();

    bool passes_this_cut = true;
    std::tie(passes_this_cut, endpoint_michels, vtx_michels) =
        PassesCut(*universe, cut, is_mc, sd, endpoint_michels, vtx_michels);

    pass = pass && passes_this_cut;
    if (!pass) continue;

    // fill container of pion candidate idxs
    for (auto m : endpoint_michels)
      event.m_reco_pion_candidate_idxs.push_back(m.second.had_idx);

    // Get the highest energy pion candidate
    // This quantity is only well-defined after you've made the
    // AtLeastOneMichel cut. This cut identifies our pion candidates and their
    // associated idxs.
    int pion_idx = -200;
    if (cut == kAtLeastOneMichel || cut == kLLR || cut == kNode ||
        cut == kIsoProngs || cut == kPionMult) {
      pion_idx = GetHighestEnergyPionCandidateIndex(event);
      event.m_highest_energy_pion_idx = pion_idx;
    }

    // Fill Wexp for each cut
    if (HasVar(variables, Form("wexp%d", i)))
      FillStackedHists(event, GetVar(variables, Form("wexp%d", i)));

    // N Hadron Tracks
    if (next_cut == kAtLeastOnePionCandidateTrack &&
        HasVar(variables, "n_had_tracks")) {
      int fill_val = universe->GetInt("MasterAnaDev_hadron_number");
      FillStackedHists(event, GetVar(variables, "n_had_tracks"), fill_val);
    }
    // Wexp
    if (next_cut == kWexp && HasVar(variables, "wexp_cut")) {
      FillStackedHists(event, GetVar(variables, "wexp_cut"));
    }
    // N michels
    if (next_cut == kAtLeastOneMichel && HasVar(variables, "michel_count")) {
      double fill_val = endpoint::GetQualityMichels(*universe).size();
      FillStackedHists(event, GetVar(variables, "michel_count"), fill_val);
      // if (fill_val == 0 && event.m_is_signal)
      //  universe->PrintArachneLink();
    }
    // New Tracking Variables -- check before LLR cut
    if (next_cut == kLLR) {
      if (HasVar(variables, "n_short_tracks"))
        FillStackedHists(event, GetVar(variables, "n_short_tracks"));
      if (HasVar(variables, "n_long_tracks"))
        FillStackedHists(event, GetVar(variables, "n_long_tracks"));
      if (HasVar(variables, "fit_vtx_x"))
        FillStackedHists(event, GetVar(variables, "fit_vtx_x"));
      if (HasVar(variables, "fit_vtx_y"))
        FillStackedHists(event, GetVar(variables, "fit_vtx_y"));
      if (HasVar(variables, "fit_vtx_z"))
        FillStackedHists(event, GetVar(variables, "fit_vtx_z"));
      if (HasVar(variables, "track_reco_meth"))
        FillStackedHists(event, GetVar(variables, "track_reco_meth"));
      if (HasVar(variables, "n_nodes"))
        FillStackedHists(event, GetVar(variables, "n_nodes"));
    }
    // LLR
    if (next_cut == kLLR && HasVar(variables, "llr")) {
      FillStackedHists(event, GetVar(variables, "llr"));
    }
    // Node
    if (next_cut == kNode && HasVar(variables, "enode01")) {
      FillStackedHists(event, GetVar(variables, "enode01"));
      FillStackedHists(event, GetVar(variables, "enode2"));
      FillStackedHists(event, GetVar(variables, "enode3"));
      FillStackedHists(event, GetVar(variables, "enode4"));
      FillStackedHists(event, GetVar(variables, "enode5"));
    }
    // N Isolated Prongs
    if (next_cut == kIsoProngs && HasVar(variables, "n_iso_prongs")) {
      FillStackedHists(event, GetVar(variables, "n_iso_prongs"));
    }
    // Pion Multiplicity
    if (next_cut == kPionMult && HasVar(variables, "n_pions")) {
      int fill_val = endpoint_michels.size();
      FillStackedHists(event, GetVar(variables, "n_pions"), fill_val);
    }

    if (i == kCutsVector.size() - 1) {
      if (HasVar(variables, "wexp"))
        FillStackedHists(event, GetVar(variables, "wexp"));
      if (HasVar(variables, "pmu"))
        FillStackedHists(event, GetVar(variables, "pmu"));
      if (HasVar(variables, "ptmu"))
        FillStackedHists(event, GetVar(variables, "ptmu"));
      if (HasVar(variables, "pzmu"))
        FillStackedHists(event, GetVar(variables, "pzmu"));
      if (HasVar(variables, "tpi"))
        FillStackedHists(event, GetVar(variables, "tpi"));
      if (HasVar(variables, "tpi_mbr"))
        FillStackedHists(event, GetVar(variables, "tpi_mbr"));
      if (HasVar(variables, "enu"))
        FillStackedHists(event, GetVar(variables, "enu"));
      if (HasVar(variables, "enu"))
        FillStackedHists(event, GetVar(variables, "enu"));
      if (HasVar(variables, "q2"))
        FillStackedHists(event, GetVar(variables, "q2"));
      if (HasVar(variables, "thetamu_deg"))
        FillStackedHists(event, GetVar(variables, "thetamu_deg"));
      if (HasVar(variables, "thetapi_deg"))
        FillStackedHists(event, GetVar(variables, "thetapi_deg"));
      if (HasVar(variables, "ehad"))
        FillStackedHists(event, GetVar(variables, "ehad"));
      if (HasVar(variables, "cosadtheta"))
        FillStackedHists(event, GetVar(variables, "cosadtheta"));
      if (HasVar(variables, "adphi"))
        FillStackedHists(event, GetVar(variables, "adphi"));
      if (HasVar(variables, "pimuAngle"))
        FillStackedHists(event, GetVar(variables, "pimuAngle"));
      if (HasVar(variables, "PT"))
        FillStackedHists(event, GetVar(variables, "PT"));
    }
  }  // end cuts loop
}

void inclusive_event::FillStackedHists(const InclusiveEvent& event,
                                  const std::vector<Variable*>& variables) {
  for (auto var : variables) FillStackedHists(event, var);
}

void inclusive_event::FillStackedHists(const InclusiveEvent& event, Variable* v,
                                  double fill_val) {
  if (!event.m_is_mc && v->m_is_true) return;

  const RecoPionIdx pion_idx = event.m_highest_energy_pion_idx;
  if (fill_val == -999.) fill_val = v->GetValue(*event.m_universe, pion_idx);

  if (!event.m_is_mc) {
    v->m_hists.m_selection_data->Fill(fill_val);
    return;
  }

  v->GetStackComponentHist(GetFSParticleType(*event.m_universe))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(GetChannelType(*event.m_universe))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(GetHadronType(*event.m_universe, pion_idx))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(GetNPionsType(*event.m_universe))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(GetNPi0Type(*event.m_universe))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(GetNPipType(*event.m_universe))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(
       GetSignalBackgroundType(*event.m_universe, event.m_signal_definition))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(
       GetWSidebandType(*event.m_universe, event.m_signal_definition))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(
       GetMesonBackgroundType(*event.m_universe, event.m_signal_definition))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(
       GetWBackgroundType(*event.m_universe, event.m_signal_definition))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(
       GetTruthWType(*event.m_universe, event.m_signal_definition))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(
       GetCoherentType(*event.m_universe, event.m_signal_definition))
      ->Fill(fill_val, event.m_weight);
}

#endif  // InclusiveEvent_cxx
