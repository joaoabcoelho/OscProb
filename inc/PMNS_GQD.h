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

#include "PMNS_Base.h"

#include <Eigen/Core>

namespace OscProb {

  class PMNS_GQD : public PMNS_Base {
  public:
    PMNS_GQD();          ///< Constructor
    virtual ~PMNS_GQD(); ///< Destructor
    
    /// Set neutrino coupling with the environment in [s]
    virtual void SetEta(double val);
    
    /// Set temperature T of the thermal gravitational waves environment in [K]
    virtual void SetTemperature(double val);

    /// Set cutoff frequency in [Hz]
    virtual void SetOmega(double val);
    
    /// Compute the probability matrix
    using PMNS_Base::ProbMatrix;
    virtual matrixD ProbMatrix(int nflvi, int nflvf);
    
  protected:
    
    // Resetting and propagating
    virtual void ResetToFlavour(int flv); ///< Reset neutrino state to pure flavour flv
    virtual void SetPureState(vectorC nu_in); ///< Set the density matrix from a pure state
    virtual void SetVacuumEigensystem();
    virtual void BuildHms();
    virtual void UpdateHam();
    virtual void RotateHam(bool to_mass, matrixC& Ham);
    virtual void SolveHam();
    
    template <typename T> void SolveEigenSystem();
    
    virtual void PropagatePath(NuPath p); ///< Propagate neutrino through a single path

    virtual double P(int flv); ///< Return the probability of final state in flavour flv

    virtual void RotateState(bool to_mass); ///< Rotate rho to/from mass basis

    double fEta; ///< The neutrino coupling with the environment in [s]
    double fT; /// The temperature of the thermal gravitational waves environment in [K]
    double fOmega; ///< The cutoff frequency in [Hz]

    complexD fHam[3][3]; ///< The full hamiltonian
    
    Eigen::Matrix3cd fHmatter;
    
    matrixC fRho; ///< The neutrino density matrix state

    matrixC fMBuffer; ///< Some memory buffer for matrix operations
  };
  
} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
