#ifndef ABSORPTION_H
#define ABSORPTION_H

#include "NuPath.h"

namespace OscProb {

  class Absorption {
    public:
      Absorption();
      virtual ~Absorption();

      virtual double Trans(double xsec);

      virtual void SetPath(std::vector<NuPath> paths); ///< Set a path sequence

    protected:
      std::vector<NuPath> fNuPaths; ///< Vector of neutrino paths

      static const double kU; ///< Atomic mass unit
  };

} // namespace OscProb

#endif
