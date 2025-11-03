///////////////////////////////////////////////////////////////////////////////
/// \struct OscProb::SKMatricesPoint
///
/// \brief Struct to organise S and K matrices for caching
///
/// A completer 
///
/// \author Joao Coelho - jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#ifndef SKMatricesPoint_H
#define SKMatricesPoint_H

#include "Definitions.h"
#include "NuPath.h"

namespace OscProb {

  struct SKMatricesPoint {
      /// Constructor
      SKMatricesPoint(int numNus = 3, double e = 0, std::vector<NuPath> p = {} ,
                 bool n = false);

      /// Set eigensystem parameters
      void SetVars(double e = 0, std::vector<NuPath> p = {}, bool n = false);

      double fEnergy; ///< Neutrino energy
      //NuPath fPath;   ///< Neutrino pathstd::vector<NuPath> fNuPaths; ///< Vector of neutrino paths
      std::vector<NuPath> fNuPaths; ///< Vector of neutrino paths
      bool   fNubar;  ///< Nu-Nubar flag
      double fNE;     ///< Energy-density

      //void SetNE(); ///< Set energy-density

      //bool operator<(const SKMatricesPoint& rhs) const;  ///< Comparison operator
      //bool operator==(const SKMatricesPoint& rhs) const; ///< Identity operator

      matrixC fS; ///< S matrix to be cached
      matrixC fK; ///< K matrix to be cached
  };

} // namespace OscProb

/*namespace std {

  template <> struct hash<OscProb::SKMatricesPoint> {
      auto operator()(const OscProb::SKMatricesPoint& ep) const -> size_t
      {
        return hash<double>()(ep.fNE) ^ hash<double>()(ep.fPath.zoa);
      }
  };

} // namespace std */

#endif
