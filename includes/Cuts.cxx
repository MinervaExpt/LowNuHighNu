//============================================================================
/*
This file contains the definitions of individual reco (aka event selection)
cuts.

It also contains the PassesCut(s) functions to apply all of them at once to
perform the event selection.

The default cuts used for the analysis are in the kCutsVector, located in
includes/CutUtils.h
*/
//============================================================================
#ifndef Cuts_cxx
#define Cuts_cxx

#include "Cuts.h"
#include "CutUtils.h"  // IsPrecut, kCutsVector

//==============================================================================
// Passes ALL Cuts
//==============================================================================
PassesCutsInfo PassesCuts(CVUniverse& universe, const bool is_mc,
                          const SignalDefinition signal_definition,
                          std::vector<ECuts> cuts) {
  //============================================================================
  // passes all cuts but w cut
  //============================================================================
  bool passes_all_cuts = true;
  for (auto c : cuts) {
    bool passes_this_cut = false;
    passes_this_cut = PassesCut(
        universe, c, is_mc, signal_definition);
    passes_all_cuts = passes_all_cuts && passes_this_cut;
  }

  return PassesCutsInfo{passes_all_cuts};
}

//==============================================================================
// Passes INDIVIDUAL Cut
//==============================================================================
bool PassesCut(
    const CVUniverse& univ, const ECuts cut, const bool is_mc,
    const SignalDefinition signal_definition) {
  bool pass = false;

  if (IsPrecut(cut) && !is_mc) return true;

  switch (cut) {
    case kNoCuts:
      pass = true;
      break;

    // gaudi cut AKA precut (maybe still used in MAD?)
    case kGoodObjects:
      pass = univ.IsTruth() ? GoodObjectsCut(univ) : true;
      break;

    // gaudi cut AKA precut (maybe still used in MAD?)
    case kGoodVertex:
      pass = univ.IsTruth() ? GoodVertexCut(univ) : true;
      break;

    // gaudi cut AKA precut (probably not used in MAD)
    case kFiducialVolume:
      pass = univ.IsTruth() ? FiducialVolumeCut(univ) : true;
      break;

    // gaudi cut AKA precut (probably not used in MAD)
    case kMinosActivity:
      pass = univ.IsTruth() ? MinosActivityCut(univ) : true;
      break;

    case kPrecuts:
      pass = univ.IsTruth() ? GoodObjectsCut(univ) && GoodVertexCut(univ) &&
                                  FiducialVolumeCut(univ)
                            : true;
      // MinosActivityCut(univ) : true;
      break;

    case kVtx:
      pass = vtxCut(univ);
      break;

    case kMinosMatch:
      pass = MinosMatchCut(univ);
      break;

    case kMinosCharge:
      pass = MinosChargeCut(univ);
      break;

    case kMinosMuon:
      pass = MinosMatchCut(univ) && MinosChargeCut(univ);
      break;

    case kPmu:
      pass = PmuCut(univ);
      break;

    case kAllCuts:
      pass = true;
      break;

    default:
      std::cout << "PassesCut Error Unknown Cut!" << cut << "  "
                << GetCutName(cut) << "\n";
      pass = false;
  };

  return pass;
}

//==============================================================================
// Cut Definitions
//==============================================================================
// Truth precuts
bool GoodObjectsCut(const CVUniverse& univ) {
  return univ.GetBool("truth_reco_hasGoodObjects");
}
bool GoodVertexCut(const CVUniverse& univ) {
  return univ.GetBool("truth_reco_isGoodVertex");
}
bool FiducialVolumeCut(const CVUniverse& univ) {
  return univ.GetBool("truth_reco_isFidVol_smeared");
}
bool MinosActivityCut(const CVUniverse& univ) {
  return univ.GetInt("truth_reco_muon_is_minos_match");
}

// Eventwide reco cuts
bool MinosMatchCut(const CVUniverse& univ) {
  return univ.GetBool("isMinosMatchTrack");
}
// Equivalent to Brandon's, but using standard minos branches
bool MinosChargeCut(const CVUniverse& univ) {
  return univ.GetDouble("MasterAnaDev_minos_trk_qp") < 0.0;
}

// Vtx cut for detection volume
bool vtxCut(const CVUniverse& univ) {
  bool pass = true;
  pass = pass && zVertexCut(univ, CCNuPionIncConsts::kZVtxMaxCutVal,
                            CCNuPionIncConsts::kZVtxMinCutVal);
  pass = pass && XYVertexCut(univ, CCNuPionIncConsts::kApothemCutVal);
  return pass;
}

bool zVertexCut(const CVUniverse& univ, const double upZ, const double downZ) {
  double vtxZ = univ.GetVecElem("vtx", 2);
  return vtxZ > downZ && vtxZ < upZ;
}

bool XYVertexCut(const CVUniverse& univ, const double a) {
  const double x = univ.GetVecElem("vtx", 0);
  const double y = univ.GetVecElem("vtx", 1);
  return univ.IsInHexagon(x, y, a);
}

bool PmuCut(const CVUniverse& univ) {
  double pmu = univ.GetPmu();
  return CCNuPionIncConsts::kPmuMinCutVal < pmu &&
         pmu < CCNuPionIncConsts::kPmuMaxCutVal;
}

#endif  // Cuts_cxx
