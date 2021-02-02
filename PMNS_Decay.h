////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_Decay
///
/// \brief Implementation of neutrino decay of the third mass state \nu_3 with oscillations of neutrinos in matter in a
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
    std::vector<double>               	 falpha;
    std::vector< std::vector<complexD> > fEvecinv; 
    std::vector<double>                  fEvalI;
    
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
    virtual void PropagatePath(OscProb::NuPath p);    ///< Propagate neutrino through a single path
    
    
    std::vector< std::vector<complexD> > fHd; //Decay hamiltonian
    complexD fHam[3][3]; //Final hamiltonian 
    complexD fHt[3][3]; //Standard+Decay hamiltonian
    
  };

}
#endif
////////////////////////////////////////////////////////////////////////
