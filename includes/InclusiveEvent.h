#ifndef InclusiveEvent_h
#define InclusiveEvent_h

#include "CVUniverse.h"
#include "Constants.h"  // typedef RecoPionIdx, EventCount, PassesCutsInfo
#include "SignalDefinition.h"
#include "TruthCategories/Sidebands.h"         // WSidebandType
#include "TruthCategories/SignalBackground.h"  // SignalBackgroundType
#include "Variable.h"
class Variable;

//==============================================================================
// InclusiveEvent is a container struct holding misc info about a universe-event:
// * passes cuts
// * is signal
// * event weight
// * is mc, is truth
// * signal definition currently being used
// * whether and what kind of w sideband it is
// And it owns a pointer to the corresponding CVUniverse.
//
// This container exists to be passed around to functions that need all this
// info: functions that fill histograms, check whether the event-universe
// passed cuts, etc.
//
// InclusiveEvent has no functions of its own.
//==============================================================================

struct InclusiveEvent {
  InclusiveEvent(const bool is_mc, const bool is_truth,
            const SignalDefinition signal_definition, CVUniverse* universe);

  // Fixed by the constructor
  const bool m_is_mc;
  const bool m_is_truth;
  const SignalDefinition m_signal_definition;
  CVUniverse* m_universe;
  bool m_is_signal;
  double m_weight;
  WSidebandType m_w_type;

  // Fixed (directly) outside of constructor -- with time-intensive functions
  bool m_passes_cuts;
  bool m_is_w_sideband;
  bool m_passes_all_cuts_except_w;
};

// Helper Functions
// bool IsWSideband(InclusiveEvent&);
PassesCutsInfo PassesCuts(const InclusiveEvent&);
SignalBackgroundType GetSignalBackgroundType(const InclusiveEvent&);

// Helper Fill Histo Functions
namespace inclusive_event {
// Xsec analysis fill functions
void FillSelected(const InclusiveEvent&, const std::vector<Variable*>&);
void FillRecoEvent(const InclusiveEvent&, const std::vector<Variable*>&);
void FillWSideband(const InclusiveEvent&, const std::vector<Variable*>&);
void FillTruthEvent(const InclusiveEvent&, const std::vector<Variable*>&);
void FillEfficiencyDenominator(const InclusiveEvent&, const std::vector<Variable*>&);
void FillMigration(const InclusiveEvent&, const std::vector<Variable*>&,
                   std::string name);

// Study functions
void FillWSideband_Study(const InclusiveEvent&, std::vector<Variable*>);
void FillCounters(const InclusiveEvent&,
                  const std::pair<EventCount*, EventCount*>& counters);
std::pair<EventCount, EventCount> FillCounters(const InclusiveEvent&,
                                               const EventCount& signal,
                                               const EventCount& bg);
void FillCutVars(InclusiveEvent&, const std::vector<Variable*>&);
void FillStackedHists(const InclusiveEvent&,
                      const std::vector<Variable*>&);  // all variables
void FillStackedHists(const InclusiveEvent&, Variable*,
                      const double fill_value = -999.);  // Single variable
}  // namespace inclusive_event


#endif  // InclusiveEvent
