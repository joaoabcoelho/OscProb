///////////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_Deco
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework with decoherence.
///
/// This class expands the PMNS_Fast class including effects from
/// decoherence in an increasing entropy and energy conserving model.
///
/// The model assumes a power law energy dependence of the decoherence
/// parameters and that decoherence occurs in the effective mass basis.
///
/// Reference: https://doi.org/10.1140/epjc/s10052-013-2434-6
///
/// \sa PMNS_Fast
///
/// \author jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_Deco_H
#define PMNS_Deco_H

#include "PMNS_DensityMatrix.h"

namespace OscProb {

  class PMNS_Deco : public PMNS_DensityMatrix {
    public:
      PMNS_Deco();          ///< Constructor
      virtual ~PMNS_Deco(); ///< Destructor

      /// Set any given decoherence parameter
      virtual void SetGamma(int j, double val);

      /// Set the \f$\Gamma_{32}\f$ parameter
      virtual void SetGamma32(double val);

      /// Set the decoherence angle
      virtual void SetDecoAngle(double th);

      /// Set the power index
      virtual void SetPower(double n);

      /// Get any given decoherence parameter
      virtual double GetGamma(int i, int j);

      /// Get the decoherence angle
      virtual double GetDecoAngle();

      /// Get the power index
      virtual double GetPower();

    protected:
      /// Propagate neutrino through a single path
      virtual void PropagatePath(NuPath p);

      double fGamma[3]; ///< Stores each decoherence parameter
      double fPower;    ///< Stores the power index parameter
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
