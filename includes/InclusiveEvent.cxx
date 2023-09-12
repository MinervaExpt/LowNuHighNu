#ifndef InclusiveEvent_cxx
#define InclusiveEvent_cxx

#include "InclusiveEvent.h"

#include "Cuts.h"    // kCutsVector
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
      m_universe(universe) {
  m_is_signal = is_mc ? IsSignal(*universe, signal_definition) : false;
  m_weight = is_mc ? universe->GetWeight() : 1.;
  m_w_type = is_mc ? GetWSidebandType(*universe, signal_definition,
                                      sidebands::kNWFitCategories)
                   : kNWSidebandTypes;
}

//==============================================================================
// Helper Functions
//==============================================================================
// PassesCutsInfo {passes_all_cuts, is_w_sideband, passes_all_except_w}
PassesCutsInfo PassesCuts(const InclusiveEvent& e) {
  return PassesCuts(*e.m_universe, e.m_is_mc, e.m_signal_definition);
}

SignalBackgroundType GetSignalBackgroundType(const InclusiveEvent& e) {
  return GetSignalBackgroundType(*e.m_universe, e.m_signal_definition);
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

  for (auto var : variables) {
    // if (var->m_is_true && !event.m_is_mc) continue; // truth, but not MC?
    if (var->m_is_true) continue;  // truth pion variables don't generally work
    const double fill_val = var->GetValue(*event.m_universe);

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
  double reco_fill_val = reco_var->GetValue(*event.m_universe);
  double true_fill_val = true_var->GetValue(*event.m_universe);
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

  Variable* var = GetVar(variables, sidebands::kFitVarString);
  double fill_val = var->GetValue(*event.m_universe);
  if (event.m_is_mc) {
    var->GetStackComponentHist(event.m_w_type)
        ->Fill(fill_val, event.m_weight);
  } else {
    var->m_hists.m_wsideband_data->Fill(fill_val);
  }
}

// Like FillCutVars, this function loops through cuts and calls PassesCut.
void inclusive_event::FillCounters(
    const InclusiveEvent& event,
    const std::pair<EventCount*, EventCount*>& counters) {
  EventCount* signal = counters.first;
  EventCount* bg = event.m_is_mc ? counters.second : nullptr;
  bool pass = true;
  // Purity and efficiency
  for (auto i_cut : kCutsVector) {
    if (event.m_is_truth != IsPrecut(i_cut))
      continue;  // truth loop does precuts

    bool passes_this_cut = true;
    passes_this_cut =
        PassesCut(*event.m_universe, i_cut, event.m_is_mc,
                  event.m_signal_definition);

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

  bool pass = true;
  for (auto i_cut : kCutsVector) {
    if (event.m_is_truth != IsPrecut(i_cut)) continue;

    bool passes_this_cut = true;
    passes_this_cut =
        PassesCut(*event.m_universe, i_cut, event.m_is_mc,
                  event.m_signal_definition);

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

    bool passes_this_cut = true;
    passes_this_cut =
        PassesCut(*universe, cut, is_mc, sd);

    pass = pass && passes_this_cut;
    if (!pass) continue;

    if (i == kCutsVector.size() - 1) {
      if (HasVar(variables, "wexp"))
        FillStackedHists(event, GetVar(variables, "wexp"));
      if (HasVar(variables, "pmu"))
        FillStackedHists(event, GetVar(variables, "pmu"));
      if (HasVar(variables, "ptmu"))
        FillStackedHists(event, GetVar(variables, "ptmu"));
      if (HasVar(variables, "pzmu"))
        FillStackedHists(event, GetVar(variables, "pzmu"));
      if (HasVar(variables, "enu"))
        FillStackedHists(event, GetVar(variables, "enu"));
      if (HasVar(variables, "enu"))
        FillStackedHists(event, GetVar(variables, "enu"));
      if (HasVar(variables, "q2"))
        FillStackedHists(event, GetVar(variables, "q2"));
      if (HasVar(variables, "thetamu_deg"))
        FillStackedHists(event, GetVar(variables, "thetamu_deg"));
      if (HasVar(variables, "ehad"))
        FillStackedHists(event, GetVar(variables, "ehad"));
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


  v->GetStackComponentHist(
       GetCoherentType(*event.m_universe, event.m_signal_definition))
      ->Fill(fill_val, event.m_weight);
}

#endif  // InclusiveEvent_cxx
