#ifndef SignalDefinition_H
#define SignalDefinition_H

#include "includes/CVUniverse.h"
#include "includes/Constants.h"  // namespace CCNuPionIncConsts

enum SignalDefinition { kInclusive, kLowNu, kHighNu, kNSignalDefTypes };

double GetWCutValue(SignalDefinition signal_definition) {
  switch (signal_definition) {
    case kLowNu:
      return 1000.;  // placeholder
    case kHighNu:
      return 1100.;  // placeholder
    case kInclusive:
      return 10.;
    default:
      std::cout << "ERROR GetWCutValue" << std::endl;
      return -1.;
  }
}

// ROB -- DON'T BOTHER READING THIS STUFF
// Truth topology particle counts
// From Aaron
std::map<string, int> GetParticleTopology(
    const std::vector<int>& FS_PDG, const std::vector<double>& FS_energy) {
  std::map<std::string, int> genie_n;

  // Overarching categories: nucleons, mesons
  genie_n["muons"] = 0;     // Muons, photons (do we want electrons...)
  genie_n["photons"] = 0;   // Photons are filled if there are electrons
  genie_n["pi_zeros"] = 0;  // Pions
  genie_n["piplus"] = 0;
  genie_n["piplus_range"] = 0;  // Pi+ that passes the kinematic cuts
  genie_n["piminus"] = 0;
  genie_n["pions"] = 0;
  genie_n["kaons"] = 0;  // Other mesons
  genie_n["charms"] = 0;
  genie_n["mesons"] = 0;
  genie_n["protons"] = 0;  // Baryons
  genie_n["neutrons"] = 0;
  genie_n["nucleons"] = 0;
  genie_n["heavy_baryons"] = 0;
  genie_n["nuclei"] = 0;
  genie_n["others"] = 0;  // others

  double tpi = 0;
  for (uint p = 0; p < FS_PDG.size(); p++) {
    // Get tpi
    tpi = (FS_energy[p] - MinervaUnits::M_pion);

    // Look at every particle's pdg
    switch (FS_PDG[p]) {
      case 13:
        genie_n["muons"]++;
        break;
      case 22:
        // Check the energy.  If below 10 MeV, it's a nuclear deexcitation,
        // which we ignore
        if (FS_energy[p] > 10) genie_n["photons"]++;
        break;
      case 11:
      case -11:
        genie_n["photons"]++;
        break;
      case 211:
        if (CCNuPionIncConsts::kTpiLoCutVal < tpi &&
            tpi < CCNuPionIncConsts::kTpiHiCutVal)
          genie_n["piplus_range"]++;
        genie_n["piplus"]++;
        genie_n["pions"]++;
        genie_n["mesons"]++;
        break;
      case -211:
        genie_n["piminus"]++;
        genie_n["pions"]++;
        genie_n["mesons"]++;
        break;
      case 111:
        genie_n["pi_zeros"]++;
        // genie_n["photons"]++;
        genie_n["pions"]++;
        genie_n["mesons"]++;
        break;
      case 130:
      case 310:
      case 311:
      case 321:
      case -130:
      case -310:
      case -311:
      case -321:
        genie_n["kaons"]++;
        genie_n["mesons"]++;
        break;
      case 411:
      case 421:
      case 431:
      case -411:
      case -421:
      case -431:
        genie_n["charms"]++;
        genie_n["mesons"]++;
        break;
      case 2212:
        genie_n["protons"]++;
        genie_n["nucleons"]++;
        break;
      case 2112:
        genie_n["neutrons"]++;
        genie_n["nucleons"]++;
        break;
      case 2000000101:  // bindino is binding energy placeholder
        break;
      default:
        if (FS_PDG[p] > 3000 && FS_PDG[p] < 5000)
          genie_n["heavy_baryons"]++;
        else if (FS_PDG[p] > 1000000000 && FS_PDG[p] < 1099999999)
          genie_n["nuclei"]++;
        else
          genie_n["others"]++;
    }
  }

  return genie_n;
}

int NOtherParticles(const CVUniverse& univ) {
  int n_other_particles = 0;
  n_other_particles = univ.GetInt("truth_N_chargedK") +
                      univ.GetInt("truth_N_K0") + univ.GetInt("truth_N_sigma") +
                      univ.GetInt("truth_N_lambda");
  return n_other_particles;
}

bool ZVtxIsSignal(const CVUniverse& univ) {
  double vtx_z = univ.GetVecElem("mc_vtx", 2);
  return CCNuPionIncConsts::kZVtxMinCutVal < vtx_z &&
                 vtx_z < CCNuPionIncConsts::kZVtxMaxCutVal
             ? true
             : false;
}

bool XYVtxIsSignal(const CVUniverse& univ) {
  return univ.IsInHexagon(univ.GetVecElem("mc_vtx", 0),  // x
                          univ.GetVecElem("mc_vtx", 1),  // y
                          CCNuPionIncConsts::kApothemCutVal);
}
// ROB -- DON'T BOTHER READING THIS STUFF

// ---> IMPORTANT FINAL FUNCTION  <----
bool IsSignal(const CVUniverse& univ, SignalDefinition sig_def = kInclusive) {
  const std::map<std::string, int> particles = GetParticleTopology(
      univ.GetVec<int>("mc_FSPartPDG"), univ.GetVec<double>("mc_FSPartE"));
  if (univ.GetInt("mc_current") == 1 && univ.GetBool("truth_is_fiducial") &&
      ZVtxIsSignal(univ) && XYVtxIsSignal(univ) &&
      univ.GetInt("mc_incoming") == 14 &&
      univ.GetThetalepTrue() < CCNuPionIncConsts::kThetamuMaxCutVal &&
      // 0. < univ.GetWexpTrue() && univ.GetWexpTrue() < GetWCutValue(sig_def)
      // &&
      // && NOtherParticles(univ) == 0
      // particles.at("piplus_range") == 1 && Is1PiPlus(particles) &&
      CCNuPionIncConsts::kPmuMinCutVal < univ.GetPmuTrue() &&
      univ.GetPmuTrue() < CCNuPionIncConsts::kPmuMaxCutVal) {
  } else {
    return false;
  }

  switch (sig_def) {
    case kLowNu:
      return true;
    case kHighNu:
      return true;
    case kInclusive:
      return true;

    default:
      std::cout << "IsSignal Error Unknown Signal Definition! I think that my "
                   "sig_def is:"
                << sig_def << "." << std::endl;
      return false;
  }
}
// ---> IMPORTANT FINAL FUNCTION  <----

std::string GetSignalName(SignalDefinition sig_def) {
  switch (sig_def) {
    case kLowNu:
      return "#nu_{#mu} Tracker #rightarrow #mu^{-} [Low Nu placeholder]";
    case kHighNu:
      return "#nu_{#mu} Tracker #rightarrow #mu^{-} [High Nu placeholder]";
    case kInclusive:
      return "#nu_{#mu} Tracker #rightarrow #mu^{-} X [Inclusive]";

    default:
      return "UNKNOWN SIGNAL";
  }
}

std::string GetSignalFileTag(SignalDefinition sig_def) {
  switch (sig_def) {
    case kLowNu:
      return "LowNu";
    case kHighNu:
      return "HighNu";
    case kInclusive:
      return "Inclusive";

    default:
      return "UNKNOWN SIGNAL";
  }
}

#endif  // Signal_H
