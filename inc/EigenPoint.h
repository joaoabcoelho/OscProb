///////////////////////////////////////////////////////////////////////////////
/// \struct OscProb::EigenPoint
///
/// \brief Struct to organise eigensystems for caching
///
/// This struct allows for comparisons of eigensystems based on
/// the neutrino energy, nu-nubar status, path density and Z/A.
///
/// \author Joao Coelho - jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#ifndef EIGENPOINT_H
#define EIGENPOINT_H

#include "Definitions.h"
#include "NuPath.h"

namespace OscProb {

  struct EigenPoint {
      /// Constructor
      EigenPoint(int numNus = 3, double e = 0, NuPath p = NuPath(0, 0),
                 bool n = false);

      /// Set eigensystem parameters
      void SetVars(double e = 0, NuPath p = NuPath(0, 0), bool n = false);

      double fEnergy; ///< Neutrino energy
      NuPath fPath;   ///< Neutrino path
      bool   fNubar;  ///< Nu-Nubar flag
      double fNE;     ///< Energy-density

      void SetNE(); ///< Set energy-density

      bool operator<(const EigenPoint& rhs) const;  ///< Comparison operator
      bool operator==(const EigenPoint& rhs) const; ///< Identity operator

      vectorD fEval; ///< Eigenvalues to be cached
      matrixC fEvec; ///< Eigenvectors to be cached
  };

} // namespace OscProb

namespace std {

  template <> struct hash<OscProb::EigenPoint> {
      auto operator()(const OscProb::EigenPoint& ep) const -> size_t
      {
        return hash<double>()(ep.fNE) ^ hash<double>()(ep.fPath.zoa);
      }
  };

} // namespace std

#endif
