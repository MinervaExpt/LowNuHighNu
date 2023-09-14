#ifndef CVUniverse_H
#define CVUniverse_H

#include <TVector3.h>

#include "Binning.h"    // CCPi::GetBinning for ehad_nopi
#include "Constants.h"  // CCNuPionIncConsts, CCNuPionIncShifts
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/MinervaUniverse.h"

class CVUniverse : public PlotUtils::MinervaUniverse {
 private:
 public:
#include "PlotUtils/LowRecoilPionFunctions.h"
#include "PlotUtils/MuonFunctions.h"
#include "PlotUtils/RecoilEnergyFunctions.h"
#include "PlotUtils/TruthFunctions.h"
#include "PlotUtils/WeightFunctions.h"
  // CTOR
  CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma = 0);

  // DTOR
  virtual ~CVUniverse(){};

  // Print arachne link
  void PrintArachneLink() const;

  // Dummy access for variable constructors
  virtual double GetDummyVar() const;
  virtual double GetDummyHadVar(const int x) const;

  // No stale cache!
  virtual void OnNewEntry() override {
    // RDF: there was some pion business here, I figured I'd preserve
    // the structure in case it ends up being useful
  }

  virtual bool IsVerticalOnly() const override { return true; }

  //==============================================================================
  // Analysis Variables
  //==============================================================================
  // muon
  virtual double GetPTmu() const;
  virtual double GetPXmu() const;
  virtual double GetPYmu() const;
  virtual double GetPZmu() const;
  virtual double GetThetamuDeg() const;

  // event-wide
  virtual double GetEhad() const;
  virtual double GetEnu() const;
  virtual double GetQ2() const;
  virtual double GetWexp() const;
  virtual double Getq0() const;
  virtual double Getq3() const;

  //==============================================================================
  // Truth
  //==============================================================================
  virtual double GetAllTrackEnergyTrue() const;
  virtual double GetEmuTrue() const;
  virtual double GetIntVtxXTrue() const;
  virtual double GetIntVtxYTrue() const;
  virtual double GetIntVtxZTrue() const;
  virtual double GetPTmuTrue() const;
  virtual double GetPXmuTrue() const;
  virtual double GetPYmuTrue() const;
  virtual double GetPZmuTrue() const;
  virtual double GetPmuTrue() const;
  virtual double GetThetamuTrue() const;
  virtual double GetThetamuTrueDeg() const;
  virtual double GetWexpTrue() const;

  //==============================
  // Ehad (GetErecoil) Variables
  //==============================
  // ehad and related variables
  virtual double GetCalRecoilEnergy() const;
  virtual double GetCalRecoilEnergy_CCPiSpline() const;
  virtual double GetCalRecoilEnergy_DefaultSpline() const;
  virtual double GetNonCalRecoilEnergy() const;

  // ehad old variables
  virtual double GetCalRecoilEnergy_CCIncSpline() const;

  // ehad truth variables
  virtual double GetEhadTrue() const;

  //==============================================================================
  // Cuts, Systematics, Studies
  //==============================================================================
  virtual bool IsInHexagon(double x, double y, double apothem) const;
  virtual bool IsInPlastic() const;
  virtual double GetFitVtxX() const;
  virtual double GetFitVtxY() const;
  virtual double GetFitVtxZ() const;
  virtual int GetNhadrons() const;

  //==============================================================================
  // Weights
  //==============================================================================
  virtual double GetAnisoDeltaDecayWarpWeight() const;
  virtual double GetDiffractiveWeight() const;
  virtual double GetGenieWarpWeight() const;
  virtual double GetLowQ2PiWeight(double q2, std::string channel) const;
  virtual double GetWeight() const;

  //=============================================================================3
  // Physics Calculations
  //==============================================================================
  TVector3 AdlerAngle(int RefSystemDef, double dmumom, double dpimom,
                      TVector3 NeuDir, TVector3 MuDir, TVector3 PiDir,
                      double Enu) const;
  double CalcQ2(const double Enu, const double Emu, const double Thetamu) const;
  double CalcWexp(const double Q2, const double Ehad) const;
  double Calcq0(const double Enu, const double Emu) const;
  double Calcq3(const double Q2, const double Enu, const double Emu) const;
  double Calct(const double epi, const double emu, const double pzpi,
               const double pzmu, const double pxpi, const double pxmu,
               const double pypi, const double pymu) const;

  //==============================================================================
  // Untracked pions functions
  //==============================================================================
  // Mehreen's full fit, circa 2022-02:
  // Tpi = p + q*range + r*sqrt(range) with
  // p = -2.93015 +- 4.44962    // Yikes, BTW
  // r = 0.132851 +- 0.0199247
  // q = 3.95884  +- 0.657313
  virtual double thetaWRTBeam(double x, double y, double z) const {
    double pyp = -1.0 * sin(MinervaUnits::numi_beam_angle_rad) * z +
                 cos(MinervaUnits::numi_beam_angle_rad) * y;
    double pzp = cos(MinervaUnits::numi_beam_angle_rad) * z +
                 sin(MinervaUnits::numi_beam_angle_rad) * y;
    double denom2 = pow(x, 2) + pow(pyp, 2) + pow(pzp, 2);
    if (0. == denom2)
      return -9999.;
    else
      return acos(pzp / sqrt(denom2));
  }
};

#endif  // CVUniverse_H
