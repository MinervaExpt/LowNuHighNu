#ifndef utilities_h
#define utilities_h

#include "Constants.h"  // Constants::PI
#include "TString.h"

TString uniq() {
  static int i = 0;
  return TString::Format("uniq%d", i++);
}

//======================================================================
// Erase element from container if element satisfies condition
// "if item in container satisfies property X then remove it"
// Example:
// ContainerEraser::erase_if(variables, [](Variable* v) {
//                       return v->Name() == "TpiMBR";
//                     });
namespace ContainerEraser {
template <class T, class predicate>
void erase_if(T& items, const predicate& p) {
  for (auto it = items.begin(); it != items.end();) {
    if (p(*it))
      it = items.erase(it);
    else
      ++it;
  }
};
}  // namespace ContainerEraser

// Angle and Geometry Functions
// Restrict to [0,pi]
double FixAngle(double angle) {
  double ret = angle;
  if (ret < 0.0) ret = -1.0 * ret;
  if (ret > Constants::PI) ret = 2.0 * Constants::PI - ret;
  return ret;
}

double ConvertRadToDeg(double rad) { return rad * 180. / Constants::PI; }

#endif  // utilities_h
