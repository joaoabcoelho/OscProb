///////////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_GQD
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework with gravitationally induced decoherence.
///
/// This class expands the PMNS_Base class including a effects from
/// decoherence in an increasing entropy and energy conserving model.
///
/// The quantum mechanical model assumes a neutrino coupling, eta,
/// with the thermal gravitational waves environment
/// with a temperature T.
/// The decoherence occurs in the effective mass basis.
///
/// This developement is part of the QGRANT project with
/// ID: 101068013,
/// founded by the HORIZON-MSCA-2021-PF-01-01 programme.
///
/// Reference: https://doi.org/10.48550/arXiv.2403.03106
///
/// \sa PMNS_Base
///
/// \author Joao Coelho - jcoelho\@apc.in2p3.fr
/// \author Alba Domi - alba.domi\@fau.de
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_GQD_H
#define PMNS_GQD_H

#include "PMNS_Deco.h"

namespace OscProb {

  class PMNS_GQD : public PMNS_Deco {
    public:
      PMNS_GQD();          ///< Constructor
      virtual ~PMNS_GQD(); ///< Destructor

      /// Set neutrino coupling with the environment in [s]
      virtual void SetEta(double val);

      /// Set temperature T of the thermal gravitational waves
      /// environment in [K]
      virtual void SetTemperature(double val);

      /// Get neutrino coupling with the environment in [s]
      virtual double GetEta();

      /// Get temperature T in [K]
      virtual double GetTemperature();

    protected:
      /// Propagate neutrino through a single path
      virtual void PropagatePath(NuPath p);

      double fEta; ///< The neutrino coupling with the environment in [s]
      double fT;   ///< The temperature of the thermal gravitational waves
                   ///< environment in [K]

      // Unit conversion constants
      static const double kkB;   ///< Boltzmann constant [eV/K]
      static const double khbar; ///< Planck constant [eV.s]
      static const double kc;    ///< Speed of light [km/s]
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
