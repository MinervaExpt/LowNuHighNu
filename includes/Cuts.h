//==============================================================================
// This file contains all event selection cuts definitions.
//
// As well as generic PassesCut and PassesCuts functions.
//==============================================================================
#ifndef Cuts_H
#define Cuts_H

#include <tuple>
#include <vector>

#include "CVUniverse.h"
#include "Constants.h"        // enum ECuts, CCNuPionIncConsts, PassesCutsInfo
#include "CutUtils.h"         // kCutsVector
#include "SignalDefinition.h"

//==============================================================================
// Generic Pass Cut(s) Functions
//      * PassesCutsInfo(passes) = PassesCuts()
//      * tuple(passes) = PassesCut(cut)
//      * PassedCuts <-- just an event counter
//==============================================================================
PassesCutsInfo PassesCuts(CVUniverse&, const bool is_mc, const SignalDefinition,
                          const std::vector<ECuts> cuts = kCutsVector);

// Event Counter
EventCount PassedCuts(const CVUniverse&, std::vector<int>& pion_candidate_idxs,
                      bool is_mc, SignalDefinition,
                      std::vector<ECuts> cuts = kCutsVector);

// Passes Single, Given Cut
bool PassesCut(const CVUniverse& univ, const ECuts cut, const bool is_mc,
          const SignalDefinition);

//==============================================================================
// Cuts Definitions
//==============================================================================
// Gaudi tool cuts -- read from truth tuple.
// (won't work if we pass a reco mc universe)
bool GoodObjectsCut(const CVUniverse&);
bool GoodVertexCut(const CVUniverse&);
bool FiducialVolumeCut(const CVUniverse&);
bool MinosActivityCut(const CVUniverse&);

// Cut Definitions -- eventwide
bool MinosMatchCut(const CVUniverse&);
bool MinosChargeCut(const CVUniverse&);
bool vtxCut(const CVUniverse& univ);
bool zVertexCut(const CVUniverse& univ, const double upZ, const double downZ);
bool XYVertexCut(const CVUniverse& univ, const double a);
bool PmuCut(const CVUniverse& univ);


#endif
