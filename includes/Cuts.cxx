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
    passes_this_cut = PassesCut(universe, c, is_mc, signal_definition);
    passes_all_cuts = passes_all_cuts && passes_this_cut;
  }

  return PassesCutsInfo{passes_all_cuts};
}

PassesCutsInfoInclusive PassesCutsInclusive(
    CVUniverse& universe, const bool is_mc,
    const SignalDefinition signal_definition, std::vector<ECuts> cuts) {
  //============================================================================
  bool passes_incl_cuts = true;
  bool is_dis_signal = true;
  bool is_lowwhighq2_sideband = true;
  bool is_lowq2highw_sideband = true;

  return {passes_incl_cuts, is_dis_signal, is_lowwhighq2_sideband,
          is_lowq2highw_sideband};
}

PassesCutsInfoMerged PassesCutsMerged(CVUniverse& universe, const bool is_mc,
                                      const SignalDefinition signal_definition,
                                      std::vector<ECuts> cuts) {

  bool passes_all_cuts = true;
  bool passes_incl_cuts = true;
  bool is_dis_signal = true;
  bool is_lownu = false;
  bool is_highnu = false;
  bool is_lowwhighq2_sideband = true;
  bool is_lowq2highw_sideband = true;

  switch (signal_definition) {
    case kLowNu:
    case kHighNu:
      for (auto c : cuts) {
        bool passes_this_cut = false;
        passes_this_cut = PassesCut(universe, c, is_mc, signal_definition);
        passes_all_cuts = passes_all_cuts && passes_this_cut;
      }
      // Do something to set is_lownu and is_highnu
      break;
    case kInclusive:
      // Do something to set the Inclusive/DIS specific bools
    default:
      std::cout << "Warning: Unknown Signal Definition! PassesCutsMerged() is not doing the correct thing!!" << std::endl;
  }

  return {passes_all_cuts, passes_incl_cuts, is_dis_signal, is_lownu, is_highnu, 
          is_lowwhighq2_sideband, is_lowq2highw_sideband};
} 

//==============================================================================
// Passes INDIVIDUAL Cut
//==============================================================================
bool PassesCut(const CVUniverse& univ, const ECuts cut, const bool is_mc,
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
// Eventwide reco cuts
bool MinosMatchCut(const CVUniverse& univ) {
  return univ.GetBool("isMinosMatchTrack");
}
// Muon helicity cut
bool MinosChargeCut(const CVUniverse& univ) {
  return univ.GetDouble("MasterAnaDev_minos_trk_qp") < 0.0;
}
// Muon coil cut
bool MinosMuonCoilCut(const CVUniverse& univ) {
    const double coilXPos = 1219.0;
    const double coilYPos = 393.0;
    const double coilR = 210.;
    const double maxR = 2500.;
    const double minos_trk_end_x = univ.GetDouble("MasterAnaDev_minos_trk_end_x");
    const double minos_trk_end_y = univ.GetDouble("MasterAnaDev_minos_trk_end_y");
    const double minos_x = minos_trk_end_x + coilXPos;
    const double minos_y = minos_trk_end_y + coilYPos;
    double minosR = sqrt( pow(minos_x,2) + pow(minos_y,2) );

    return (minosR > coilR && minosR < maxR);
}
// Muon curvature cut
bool MinosMuonCurveCut(const CVUniverse& univ) {
    //Don't make a significance cut if we reconstructed by range 
    const bool minos_used_curvature = univ.GetBool("MasterAnaDev_minos_used_curvature");
    if(minos_used_curvature != 1)
        return true;
    const double minos_trk_eqp_qp = univ.GetDouble("MasterAnaDev_minos_trk_eqp_qp");
    return ( minos_trk_eqp_qp > -0.2 );
}

// To-do: Figure out if either of these cuts from Rob's historical implementation are needed now:
  //if not chain.pass_canonical_cut == 1:                               return False
  //if not chain.phys_n_dead_discr_pair_upstream_prim_track_proj <= 1:  return False

// Vtx cut for detection volume
bool vtxCut(const CVUniverse& univ) {
  bool pass = true;
  pass = pass && zVertexCut(univ, CCIncConsts::kZVtxMaxCutVal,
                            CCIncConsts::kZVtxMinCutVal);
  pass = pass && XYVertexCut(univ, CCIncConsts::kApothemCutVal);
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

bool PmuCut(const CVUniverse& univ, bool apply_upper_bound_cut /* = true */) {
  double pmu = univ.GetPmu();
  bool pass_lower_bound_cut = CCIncConsts::kPmuMinCutVal < pmu;
  bool pass_upper_bound_cut = apply_upper_bound_cut ? (pmu < CCIncConsts::kPmuMaxCutVal) : true;
  return pass_lower_bound_cut && pass_upper_bound_cut;
}

bool ThetamuCut(const CVUniverse& univ) {
  if (univ.GetThetamu() >= CCIncConsts::kThetamuMaxCutVal) return false;
  else return true;
}

#endif  // Cuts_cxx
