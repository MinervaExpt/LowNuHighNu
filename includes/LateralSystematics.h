#ifndef LateralSystematics_H
#define LateralSystematics_H

#include <iostream>

#include "CVUniverse.h"
#include "Constants.h"         // CCNuPionIncShifts
#include "TRandom.h"           // gRandom
#include "common_functions.h"  // FixAngle

class TrackAngleShiftCVUniverse : public CVUniverse {
 public:
  TrackAngleShiftCVUniverse(PlotUtils::ChainWrapper* chw, double nsigma)
      : CVUniverse(chw, nsigma) {}

  virtual double GetThetamu() const {
    double shift_val =
        (CCNuPionIncShifts::muon_angle_res) * (gRandom->Gaus(0.0, 1.0));
    return FixAngle(shift_val + CVUniverse::GetThetamu());
  }

  virtual std::string ShortName() const { return "TrackAngle"; }
  virtual std::string LatexName() const { return "Track Angle"; }
};

class BeamAngleShiftCVUniverse : public CVUniverse {
 public:
  BeamAngleShiftCVUniverse(PlotUtils::ChainWrapper* chw, double nsigma)
      : CVUniverse(chw, nsigma) {}

  virtual double GetThetamu() const {
    double shift_val = m_nsigma > 0
                           ? GetDouble("MasterAnaDev_muon_theta_biasUp")
                           : GetDouble("MasterAnaDev_muon_theta_biasDown");
    return FixAngle(shift_val + CVUniverse::GetThetamu());
  }

  virtual std::string ShortName() const { return "BeamAngle"; }
  virtual std::string LatexName() const { return "Beam Angle"; }
};

#endif  // LateralSystematics_H
