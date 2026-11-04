#ifndef CutUtils_h
#define CutUtils_h

#include "Constants.h"  // enum ECuts, CCIncConsts

// Analysis Cuts - default vector
const std::vector<ECuts> kCutsVector = {kNoCuts, kVtx, kMinosMuon, kPmu};

// The thesis selection (SELECTION=thesis), for the comparison with Rob's 2020
// thesis (THESIS-MAP.md, "The thesis selection").  Differences from
// kCutsVector:
//  * added: MINOS curvature significance eqp/qp > -0.3 on every event,
//    upstream dead discriminator pairs <= 1, theta_mu < 0.35 rad,
//    E_mu > 1.8 GeV;
//  * dropped: the 1.5 < p_mu < 20 GeV window.  The thesis had no upper cut
//    (Rob, 2026-11-04: drop it for this comparison, revisit it for the final
//    analysis); E_mu > 1.8 GeV is its threshold;
//  * kept: the 850 apothem in kVtx.  The thesis code has no explicit reco
//    apothem cut, but it required NukeCC_pass_canonical_cut, which most likely
//    held one (Rob, 2026-11-04).  That canonical cut is the one part of the
//    thesis selection this vector does not reproduce.
const std::vector<ECuts> kCutsVectorThesis = {
    kNoCuts,
    kVtx,
    kMinosMuon,
    kThesisMinosCurvature,
    kThesisDeadDiscriminators,
    kThesisThetamu,
    kThesisEmu};

// Gaudi tool cuts - only work when checking truth tuple
bool IsPrecut(ECuts c) {
  if (c == kNoCuts || c == kGoodObjects || c == kGoodVertex ||
      c == kFiducialVolume || c == kMinosActivity)
    return true;
  else
    return false;
}

// Cut Names
std::string GetCutName(ECuts cut) {
  switch (cut) {
    case kNoCuts:
      return "No Cuts";

    case kGoodObjects:
      return "Good Objects";

    case kGoodVertex:
      return "Good Vertex";

    case kFiducialVolume:
      return "Fiducial Volume";

    case kMinosActivity:
      return "MINOS Activity";

    case kVtx:
      return "vertex position Cut";

    case kMinosMatch:
      return "MINOS Muon";

    case kMinosCharge:
      return "MINOS Charge";

    case kMinosMuon:
      return "MINOS Muon";

    case kAllCuts:
      return "Total";

    case kPmu:
      return "1.5 GeV $<$ Pmu $<$ 20 GeV";

    case kThesisMinosCurvature:
      return "MINOS eqp/qp $>$ -0.3 (thesis)";

    case kThesisDeadDiscriminators:
      return "Upstream dead discriminators $\\leq$ 1 (thesis)";

    case kThesisThetamu:
      return "$\\theta_\\mu <$ 0.35 rad (thesis)";

    case kThesisEmu:
      return "E$_\\mu >$ 1.8 GeV (thesis)";

    default:
      std::cout << "ERROR: GetCutName unknown cut!" << std::endl;
      return "";
  };
}

#endif
