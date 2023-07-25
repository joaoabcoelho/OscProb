////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_Fast
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework.
///
/// Two optimizations are relevant:\n
/// - The construction of the Hamiltonian avoids computing null terms\n
/// - The eigensystem determination is based on the following reference:\n
///
/// 
/// <pre>
/// 
///......................................................................
///
/// Int. J. Mod. Phys. C       VOLUME 19, NUMBER 03            MARCH 2008
///
///     Efficient numerical diagonalization of hermitian 3x3 matrices
///
///                            Joachim Kopp
///                  Max-Planck-Institut für Kernphysik 
///             Postfach 10 39 80, 69029 Heidelberg, Germany
///                    (Received 19 October 2007)
///
///                                523
///......................................................................
/// </pre>
///
/// \author Joao Coelho - jcoelho\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

#ifndef PMNS_Fast_H
#define PMNS_Fast_H
#include "PMNS_Base.h"

namespace OscProb {

  class PMNS_Fast : public PMNS_Base {

  public:

    PMNS_Fast();          ///< Constructor
    virtual ~PMNS_Fast(); ///< Destructor
    
    /// Set the all mixing parameters at once
    virtual void SetMix(double th12, double th23, double th13, double deltacp);
    
    /// Set both mass-splittings at once
    virtual void SetDeltaMsqrs(double dm21, double dm32);
    
  protected:

    /// Build the full Hamiltonian
    virtual void UpdateHam();

    /// Solve the full Hamiltonian for eigenvectors and eigenvalues
    virtual void SolveHam();
    
    /// Set the eigensystem to the analytic solution of the vacuum Hamiltonian
    virtual void SetVacuumEigensystem();
    
    complexD fHam[3][3]; ///< The full hamiltonian
    
  };

}
#endif
////////////////////////////////////////////////////////////////////////
