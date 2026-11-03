#ifndef LowNuBoundary_h
#define LowNuBoundary_h

//==============================================================================
// LowNuBoundary.h
//
// Single source of truth for the low-nu/high-nu staircase cut.
//
// GetLowNuCutValue_GeV(Enu) returns the maximum hadronic recoil nu (in GeV)
// for which an event at the given neutrino energy is classified low-nu.
// IsLowNu / IsLowNuTruth (in SignalDefinition.h) and the boundary-region
// diagnostic Variables (Delta_nu = nu - nu_cut, in CVUniverse) all read it,
// so the staircase cannot drift between the classifier and the diagnostics.
//
// Marked inline so it can be included in multiple translation units without
// triggering one-definition-rule conflicts at link time.
//==============================================================================

inline double GetLowNuCutValue_GeV(double Enu_GeV) {
  if (Enu_GeV < 3.)  return 0.3;
  if (Enu_GeV < 7.)  return 0.5;
  if (Enu_GeV < 12.) return 1.0;
  return 2.0;
}

// The four nu0 plateaus of the staircase above, in order (GeV): the per-nu0
// samples.  Must match GetLowNuCutValue_GeV's return values (and the "ehad"
// bin edges in Binning.h, which processHists.py uses for the per-nu0 samples).
const int kNLowNuPlateaus = 4;
inline double GetLowNuPlateau_GeV(int i) {
  static const double plateaus[kNLowNuPlateaus] = {0.3, 0.5, 1.0, 2.0};
  return plateaus[i];
}

#endif  // LowNuBoundary_h
