///////////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_Decay
///
/// \brief Implementation of neutrino decay in a three-neutrino framework.
///
/// This class expands the PMNS_Fast class including the decay of the
/// second and third mass state of the neutrino through a decay constant
/// \f$\alpha_i=m_i/\tau_i\ (eV^2)\f$, where \f$m_i\f$ is the mass in the
/// restframe and \f$tau_i\f$ is the lifetime in the restframe.
///
/// Reference: https://doi.org/10.1007/JHEP04(2023)090
///
/// Propagation is computed by exponentiating the Hamiltonian directly
/// instead of solving the eigensystem. This is done with the
/// [Eigen](https://eigen.tuxfamily.org/) library.
///
/// \author Victor Carretero - vcarretero\@km3net.de
/// \author Joao Coelho - jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_Decay_H
#define PMNS_Decay_H

#include <Eigen/Core>

#include "PMNS_Base.h"

namespace OscProb {

  class PMNS_Decay : public PMNS_Base {
    public:
      PMNS_Decay();          ///< Constructor
      virtual ~PMNS_Decay(); ///< Destructor

      /// Set the all mixing parameters at once
      virtual void SetMix(double th12, double th23, double th13,
                          double deltacp);

      /// Set both mass-splittings at once
      virtual void SetDeltaMsqrs(double dm21, double dm32);

      virtual void   SetAlpha3(double alpha3);
      virtual void   SetAlpha2(double alpha2);
      virtual double GetAlpha3();
      virtual double GetAlpha2();
      virtual void   SetIsNuBar(bool isNuBar);

    protected:
      /// Build the Hms Hamiltonian
      virtual void BuildHms();

      /// Build the full Hamiltonian
      virtual void UpdateHam();

      /// Solve the full Hamiltonian for eigenvalues
      virtual void SolveHam();

      /// Propagation with Decay
      virtual void PropagatePath(NuPath p);

      matrixC          fHd;  ///< Decay hamiltonian
      Eigen::Matrix3cd fHam; ///< Final hamiltonian

      vectorD fAlpha; ///< alpha parameters
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
