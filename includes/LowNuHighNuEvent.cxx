#ifndef LowNuHighNuEvent_cxx
#define LowNuHighNuEvent_cxx

#include "LowNuHighNuEvent.h"

#include "Cuts.h"              // kCutsVector
#include "common_functions.h"  // GetVar, HasVar

//==============================================================================
// CTOR
//==============================================================================
LowNuHighNuEvent::LowNuHighNuEvent(const bool is_mc, const bool is_truth,
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
PassesCutsInfo PassesCuts(const LowNuHighNuEvent& e) {
  return PassesCuts(*e.m_universe, e.m_is_mc, e.m_signal_definition);
}

SignalBackgroundType GetSignalBackgroundType(const LowNuHighNuEvent& e) {
  return GetSignalBackgroundType(*e.m_universe, e.m_signal_definition);
}

//==============================================================================
// Fill all histos for an entire event -- call other specialized fill functions
//==============================================================================
void lownuhighnu_event::FillRecoEvent(
    const LowNuHighNuEvent& event, const std::vector<Variable*>& variables,
    const std::vector<VariableMAT*>& variables_MAT,
    const std::vector<Variable2D*>& variables2D) {
  // Fill selection -- total, signal-only, and bg-only
  if (event.m_passes_cuts) {
    lownuhighnu_event::FillSelected(event, variables, variables_MAT, variables2D);
  }
  // // Fill W Sideband
  // if (event.m_is_w_sideband) {
  //   lownuhighnu_event::FillWSideband(event, variables);
  // }

  // Fill Migration
  if (event.m_is_mc && event.m_is_signal && event.m_passes_cuts) {
    FillMigration(event, variables);
  }
}

void lownuhighnu_event::FillTruthEvent(
    const LowNuHighNuEvent& event, const std::vector<Variable*>& variables) {
  // Fill Efficiency Denominator
  if (event.m_is_signal)
    lownuhighnu_event::FillEfficiencyDenominator(event, variables);
}

//==============================================================================
// Specialized fill functions -- for xsec calculation
//==============================================================================
// Fill histos with selected (i.e. passes_cuts) events:
// ** sig + bg (true and reco vars, data and mc)
// ** signal only (true vars, for eff num & closure)
// ** bg only (reco and true vars)
void lownuhighnu_event::FillSelected(
    const LowNuHighNuEvent& event, const std::vector<Variable*>& variables,
    const std::vector<VariableMAT*>& variables_MAT,
    const std::vector<Variable2D*>& variables2D) {
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

  for (auto var : variables2D) {
    // Get fill value
    double fill_val_x = var->GetRecoValueX(*event.m_universe);
    double fill_val_y = var->GetRecoValueY(*event.m_universe);

    // total = signal & background, together
    if (event.m_is_mc) {
      var->m_selection_mc.FillUniverse(*event.m_universe, fill_val_x,
                                       fill_val_y, event.m_weight);
    } else {
      var->m_selection_data.hist->Fill(fill_val_x, fill_val_y);
    }

    // done with data
    if (!event.m_is_mc) continue;

    // signal and background individually
    if (event.m_is_signal) {
      var->m_effnum.FillUniverse(*event.m_universe, fill_val_x, fill_val_y,
                                 event.m_weight);
    } else {
      var->m_bg.FillUniverse(*event.m_universe, fill_val_x, fill_val_y,
                             event.m_weight);
    }
  }  // end variables2D
}

// Fill histograms of all variables with events in the sideband region
void lownuhighnu_event::FillWSideband(const LowNuHighNuEvent& event,
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

void lownuhighnu_event::FillMigration(const LowNuHighNuEvent& event,
                                      const vector<Variable*>& variables) {
  for (auto var : variables) {
    if (var->m_is_true) continue;
    std::string var_name = var->Name();
    if (!HasVar(variables, var_name + string("_true"))) continue;
    Variable* reco_var = GetVar(variables, var_name);
    Variable* true_var = GetVar(variables, var_name + string("_true"));
    double reco_fill_val = reco_var->GetValue(*event.m_universe);
    double true_fill_val = true_var->GetValue(*event.m_universe);
    reco_var->m_hists.m_migration.FillUniverse(*event.m_universe, reco_fill_val,
                                               true_fill_val, event.m_weight);
  }
}

// Only for true variables
void lownuhighnu_event::FillEfficiencyDenominator(
    const LowNuHighNuEvent& event, const std::vector<Variable*>& variables) {
  for (auto var : variables) {
    if (!var->m_is_true) continue;
    double fill_val = var->GetValue(*event.m_universe);
    try {
      var->m_hists.m_effden.FillUniverse(*event.m_universe, fill_val,
                                         event.m_weight);
    } catch (...) {
      std::cerr << "From lownuhighnu_event::FillEfficiencyDenominator\n";
      std::cerr << "Variable is " << var->Name() << "\n";
      throw;
    }
  }
}

//==============================================================================
// Specialized fill functions -- for studies
//==============================================================================
// Like FillCutVars, this function loops through cuts and calls PassesCut.
void lownuhighnu_event::FillCounters(
    const LowNuHighNuEvent& event,
    const std::pair<EventCount*, EventCount*>& counters) {
  EventCount* signal = counters.first;
  EventCount* bg = event.m_is_mc ? counters.second : nullptr;
  bool pass = true;
  // Purity and efficiency
  for (auto i_cut : kCutsVector) {
    if (event.m_is_truth != IsPrecut(i_cut))
      continue;  // truth loop does precuts

    bool passes_this_cut = true;
    passes_this_cut = PassesCut(*event.m_universe, i_cut, event.m_is_mc,
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

std::pair<EventCount, EventCount> lownuhighnu_event::FillCounters(
    const LowNuHighNuEvent& event, const EventCount& s, const EventCount& b) {
  EventCount signal = s;
  EventCount bg = b;

  bool pass = true;
  for (auto i_cut : kCutsVector) {
    if (event.m_is_truth != IsPrecut(i_cut)) continue;

    bool passes_this_cut = true;
    passes_this_cut = PassesCut(*event.m_universe, i_cut, event.m_is_mc,
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

void lownuhighnu_event::FillCutVars(LowNuHighNuEvent& event,
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
    passes_this_cut = PassesCut(*universe, cut, is_mc, sd);

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

void lownuhighnu_event::FillStackedHists(
    const LowNuHighNuEvent& event, const std::vector<Variable*>& variables) {
  for (auto var : variables) FillStackedHists(event, var);
}

void lownuhighnu_event::FillStackedHists(const LowNuHighNuEvent& event,
                                         Variable* v, double fill_val) {
  if (!event.m_is_mc && v->m_is_true) return;

  if (fill_val == -999.) fill_val = v->GetValue(*event.m_universe);

  if (!event.m_is_mc) {
    v->m_hists.m_selection_data->Fill(fill_val);
    return;
  }

  v->GetStackComponentHist(GetFSParticleType(*event.m_universe))
      ->Fill(fill_val, event.m_weight);

  v->GetStackComponentHist(GetChannelType(*event.m_universe))
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
}

#endif  // LowNuHighNuEvent_cxx
