////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_Decay
///
/// \brief Implementation of neutrino decay of the third mass state 
///        \f$\nu_3\f$ with oscillations of neutrinos in matter in a
///        three-neutrino framework.
///
///                               
///......................................................................
/// 
///
/// \author Victor Carretero - vcarretero@km3net.de
////////////////////////////////////////////////////////////////////////

#ifndef PMNS_Decay_H
#define PMNS_Decay_H
#include "PMNS_Base.h"

namespace OscProb {

  class PMNS_Decay : public PMNS_Base {

  public:

    PMNS_Decay();          ///< Constructor
    virtual ~PMNS_Decay(); ///< Destructor
    
    /// Set the all mixing parameters at once
    virtual void SetMix(double th12, double th23, double th13, double deltacp);
    
    /// Set both mass-splittings at once
    virtual void SetDeltaMsqrs(double dm21, double dm32);
    //Set Alpha 3
    vectorD falpha;
    matrixC fEvecinv;
    vectorD fEvalI;
    
    virtual void SetAlpha3(double alpha3);
    virtual void SetAlpha2(double alpha2);
    virtual double GetAlpha3();
    virtual double GetAlpha2();
    virtual void SetIsNuBar(bool isNuBar);
    
  protected:
    virtual void BuildHam();
   
    /// Build the full Hamiltonian
    virtual void UpdateHam();
   

    /// Solve the full Hamiltonian for eigenvectors and eigenvalues
    virtual void SolveHam();
    
    ///Propagation with D
    virtual void PropagatePath(NuPath p);    ///< Propagate neutrino through a single path
    
    
    matrixC fHd;  //Decay hamiltonian
    matrixC fHam; //Final hamiltonian
    matrixC fHt;  //Standard+Decay hamiltonian
   
    vectorC fEvalEigen;
    matrixC fEvecEigen;
    matrixC fEvecEigeninv;
  };

}
#endif
////////////////////////////////////////////////////////////////////////
