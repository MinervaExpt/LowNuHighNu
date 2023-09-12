#ifndef CVUniverse_cxx
#define CVUniverse_cxx

#include "CVUniverse.h"

#include <TVector3.h>

#include <algorithm>  // max_element
#include <cmath>      //isfinite

#include "PlotUtils/MnvTuneSystematics.h"
#include "utilities.h"  // FixAngle

//==============================================================================
// Constructor
//==============================================================================
CVUniverse::CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma)
    : PlotUtils::MinervaUniverse(chw, nsigma) {}

//==============================================================================
// Print arachne link
//==============================================================================
void CVUniverse::PrintArachneLink() const {
  int link_size = 200;
  char link[link_size];
  int run = GetInt("mc_run");
  int subrun = GetInt("mc_subrun");
  int gate = GetInt("mc_nthEvtInFile") + 1;
  int slice = GetVecElem("slice_numbers", 0);
  sprintf(link,
          "http://minerva05.fnal.gov/Arachne/"
          "arachne.html\?det=SIM_minerva&recoVer=v21r1p1&run=%d&subrun=%d&gate="
          "%d&slice=%d",
          run, subrun, gate, slice);
  // strncpy(); // FAIL
  // memcpy();  // FAIL
  std::cout << link << std::endl;
}

//==============================================================================
// Dummy access for variable constructors
//==============================================================================
double CVUniverse::GetDummyVar() const { return -999.; }
double CVUniverse::GetDummyHadVar(const int x) const { return -999.; }

//==============================================================================
// Analysis Variables
//==============================================================================
// muon
double CVUniverse::GetPTmu() const {
  return sqrt(pow(GetPXmu(), 2.0) + pow(GetPYmu(), 2.0));
}

double CVUniverse::GetPXmu() const { return GetMuon4V().Px(); }

double CVUniverse::GetPYmu() const { return GetMuon4V().Py(); }

double CVUniverse::GetPZmu() const { return GetMuon4V().Pz(); }

double CVUniverse::GetThetamuDeg() const {
  return ConvertRadToDeg(GetThetamu());
}

// event-wide
double CVUniverse::GetEhad() const {
  return GetCalRecoilEnergy();
}
double CVUniverse::GetEnu() const { return GetEmu() + GetEhad(); }

double CVUniverse::GetQ2() const {
  return CalcQ2(GetEnu(), GetEmu(), GetThetamu());
}

double CVUniverse::GetWexp() const { return CalcWexp(GetQ2(), GetEhad()); }

double CVUniverse::Getq0() const { return Calcq0(GetEnu(), GetEmu()); }

double CVUniverse::Getq3() const { return Calcq3(GetQ2(), GetEnu(), GetEmu()); }

//==============================================================================
// Truth
//==============================================================================
// Need a function that gets truth tracked energy, total E for pions, just KE
// for protons.
double CVUniverse::GetAllTrackEnergyTrue() const {
  double etracks = 0;
  // std::cout << "\n";
  return etracks;
}

double CVUniverse::GetEmuTrue() const { return GetElepTrue(); }

double CVUniverse::GetIntVtxXTrue() const { return GetVecElem("mc_vtx", 0); }

double CVUniverse::GetIntVtxYTrue() const { return GetVecElem("mc_vtx", 1); }

double CVUniverse::GetIntVtxZTrue() const { return GetVecElem("mc_vtx", 2); }

double CVUniverse::GetPTmuTrue() const {
  return GetPlepTrue() * sin(GetThetalepTrue());
}

double CVUniverse::GetPXmuTrue() const {
  return GetVecElem("mc_primFSLepton", 0);
}

double CVUniverse::GetPYmuTrue() const {
  return GetVecElem("mc_primFSLepton", 1);
}

double CVUniverse::GetPZmuTrue() const {
  return GetPlepTrue() * cos(GetThetalepTrue());
}

double CVUniverse::GetPmuTrue() const { return GetPlepTrue(); }

double CVUniverse::GetThetamuTrue() const {
  return FixAngle(GetThetalepTrue());
}

double CVUniverse::GetThetamuTrueDeg() const {
  return ConvertRadToDeg(GetThetamuTrue());
}

double CVUniverse::GetWexpTrue() const {
  return CalcWexp(GetQ2True(), GetEhadTrue());
}


//==============================
// Ehad (GetErecoil) Variables
//==============================
// Untracked recoil energy
double CVUniverse::GetCalRecoilEnergy() const {
  return GetCalRecoilEnergy_CCIncSpline(); //AEN This really only works for the tracker now 
}

// Total recoil with CCPi spline correction.
// Spline measured from: CC, 1pi+, 1mu, NBaryons,
// True W exp < 1.4 GeV, thmu_true < 20 deg, 35 < tpi < 350 MeV, Minos match
// RecoilUtils->calcRecoilEFromClusters(event, muonProng,
// "NukeCCPion_TwoTrack_Nu_Tracker");
double CVUniverse::GetCalRecoilEnergy_CCPiSpline() const {
  return GetDouble("MasterAnaDev_hadron_recoil_two_track");
}

// RecoilUtils->calcRecoilEFromClusters(event, muonProng, "Default" );
double CVUniverse::GetCalRecoilEnergy_DefaultSpline() const {
  return GetDouble("MasterAnaDev_hadron_recoil_default");
}

// RDF: This is a placeholder; required because the Response systematics expect
// a function with this name to exist
// This is what the response universe calls our tracked recoil energy
double CVUniverse::GetNonCalRecoilEnergy() const {
  return 0.; 
}

//==============================
// ehad old variables
//==============================

// Ehad CCInclusive Spline Variables
// Ehad ccinclusive splines -- doesn't account for pion
double CVUniverse::GetCalRecoilEnergy_CCIncSpline() const {
  return GetDouble("MasterAnaDev_hadron_recoil_CCInc");
}

//==============================
// ehad truth variables
//==============================

double CVUniverse::GetEhadTrue() const { return GetEnuTrue() - GetElepTrue(); }

//==============================================================================
// Cuts, Systematics, Studies
//==============================================================================
bool CVUniverse::IsInHexagon(double x, double y, double apothem) const {
  double lenOfSide = apothem * (2 / sqrt(3));
  double slope = (lenOfSide / 2.0) / apothem;
  double xp = fabs(x);
  double yp = fabs(y);

  if ((xp * xp + yp * yp) < apothem * apothem)
    return true;
  else if (xp <= apothem && yp * yp < lenOfSide / 2.0)
    return true;
  else if (xp <= apothem && yp < lenOfSide - xp * slope)
    return true;

  return false;
}

bool CVUniverse::IsInPlastic() const {
  if (!IsInHexagon(GetVecElem("mc_vtx", 0), GetVecElem("mc_vtx", 1), 1000.0))
    return false;  // This is in the calorimeters

  double mc_vtx_z = GetVecElem("mc_vtx", 2);
  if (mc_vtx_z > 8467.0) return false;  // Ditto

  int mc_nuclei = GetInt("mc_targetZ");
  // In the carbon target?  The z is gotten from NukeBinningUtils
  if (fabs(mc_vtx_z - 4945.92) <=
          PlotUtils::TargetProp::ThicknessMC::Tgt3::C / 2 &&
      mc_nuclei == 6)
    return false;

  // In the water target?
  if (5200 < mc_vtx_z && mc_vtx_z < 5420 && (mc_nuclei == 8 || mc_nuclei == 1))
    return false;

  // Finally, do you have target material?  I'm going to say lead/iron isn't a
  // big consideration elsewhere in the detector
  if (mc_nuclei == 26 || mc_nuclei == 82) return false;

  return true;
}

double CVUniverse::GetFitVtxX() const {
  return GetVecElem("MasterAnaDev_vtx", 0);
}  // cm?

double CVUniverse::GetFitVtxY() const {
  return GetVecElem("MasterAnaDev_vtx", 1);
}  // cm?

double CVUniverse::GetFitVtxZ() const {
  return GetVecElem("MasterAnaDev_vtx", 2);
}  // cm?

int CVUniverse::GetNhadrons() const {
  return GetInt("MasterAnaDev_hadron_number");
}

//==============================================================================
// Weights
//==============================================================================

double CVUniverse::GetAnisoDeltaDecayWarpWeight() const {
  return GetVecElem("truth_genie_wgt_Theta_Delta2Npi", 4);
}

// Note, this assumes you're not using the diffractive model in GENIE
// As of 03/2021, we don't really trust our diffractive model, so
// as a rough approximation, weight every coherent event
// (diffractive is coherent on hydrogen) by 1.4368.
// Coherent xsec scales by A^(1/3), and 1/(12^(1/3)) = 0.4368
double CVUniverse::GetDiffractiveWeight() const {
  if (GetInt("mc_intType") != 4) return 1.;
  // Note: diffractive should be applied only to plastic. This approximates that
  // if( PlotUtils::TargetUtils::Get().InCarbon3VolMC( GetVecElem("mc_vtx",0),
  //                                 GetVecElem("mc_vtx",1),
  //                                 GetVecElem("mc_vtx",2) ) ) return 1.;
  // if( GetInt("mc_nucleiZ") != 6 ) 1.;
  if (!IsInPlastic() && !PlotUtils::TargetUtils::Get().InWaterTargetMC(
                            GetIntVtxXTrue(), GetIntVtxYTrue(),
                            GetIntVtxZTrue(), GetInt("mc_targetZ")))
    return 1.;

  return 1.4368;
}

double CVUniverse::GetGenieWarpWeight() const {
  double wgt = GetVecElem("truth_genie_wgt_MaRES", 4);
  wgt = 1 + (wgt - 1) *
                2;  // double the size of the shift from 1 (e.g. 1.1 --> 1.2)
  return wgt;
}

double CVUniverse::GetLowQ2PiWeight(double q2, std::string channel) const {
  if (!PlotUtils::IsCCRes(*this))
    return 1.;
  else
    return PlotUtils::weight_lowq2pi().getWeight(q2, channel, +1);
}

double CVUniverse::GetWeight() const {
  // Warping strategy is to only turn on one of these at a time.
  const bool do_genie_warping = false;
  const bool do_aniso_warping = false;
  const bool do_mk_warping = false;

  double wgt_flux = 1., wgt_2p2h = 1.;
  double wgt_rpa = 1., wgt_lowq2 = 1.;
  double wgt_genie = 1., wgt_mueff = 1.;
  double wgt_anisodd = 1.;
  double wgt_diffractive = 1.;
  double wgt_mk = 1.;
  double wgt_target = 1.;
  double wgt_fsi = 1., wgt_geant = 1.,
         wgt_sbfit = 1. /* This weight depends of the sidebands, Will we applay
                           this weight?*/
      ;

  // genie
  wgt_genie = GetGenieWeight();
  if (do_genie_warping) wgt_genie = GetGenieWarpWeight();

  // flux

  wgt_flux = GetFluxAndCVWeight();

  // rpa
  wgt_rpa = GetRPAWeight();

  // MINOS efficiency
  if (!m_is_truth && GetBool("isMinosMatchTrack"))
    wgt_mueff = GetMinosEfficiencyWeight();

  // 2p2h
  wgt_2p2h = GetLowRecoil2p2hWeight();

  // low q2
  wgt_lowq2 = (GetQ2True() > 0)
                  ? GetLowQ2PiWeight(CCNuPionIncShifts::kLowQ2PiChannel)
                  : 1.;

  // aniso delta decay weight -- currently being used for warping
  if (do_aniso_warping)
    wgt_anisodd = GetVecElem("truth_genie_wgt_Theta_Delta2Npi", 4);

  // Diffractive
  wgt_diffractive = GetDiffractiveWeight();

  // MK Weight
  if (do_mk_warping) wgt_mk = GetMKWeight();

  // Target Mass
  wgt_target = GetTargetMassWeight();

  // New Weights added taking as reference Aaron's weights

  wgt_fsi = GetFSIWeight(0);

  wgt_geant = GetGeantHadronWeight();

  return wgt_genie * wgt_flux * wgt_2p2h * wgt_rpa * wgt_lowq2 * wgt_mueff *
         wgt_anisodd * wgt_diffractive * wgt_mk * wgt_target *
         wgt_fsi * wgt_geant * wgt_sbfit;
}

//==============================================================================
// Physics Calculations
//==============================================================================

double CVUniverse::CalcQ2(const double Enu, const double Emu,
                          const double Thetamu) const {
  double Q2 =
      2.0 * Enu *
          (Emu - sqrt(pow(Emu, 2.0) - pow(CCNuPionIncConsts::MUON_MASS, 2.0)) *
                     cos(Thetamu)) -
      pow(CCNuPionIncConsts::MUON_MASS, 2.0);
  if (Q2 < 0.) Q2 = 0.0;
  return Q2;
}

double CVUniverse::CalcWexp(const double Q2, const double Ehad) const {
  double W = pow(CCNuPionIncConsts::PROTON_MASS, 2.0) - Q2 +
             2.0 * (CCNuPionIncConsts::PROTON_MASS)*Ehad;
  W = W > 0 ? sqrt(W) : 0.0;
  return W;
}

double CVUniverse::Calcq0(const double Enu, const double Emu) const {
  return Enu - Emu;
}

double CVUniverse::Calcq3(const double Q2, const double Enu,
                          const double Emu) const {
  return sqrt(Q2 + pow(Enu - Emu, 2.0));
}


#endif  // CVUniverse_cxx
