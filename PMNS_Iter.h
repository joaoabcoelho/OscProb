////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_Iter
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework.
///
/// Based on an approximate iterative model from:\n
/// 
/// <pre>
/// 
///......................................................................
///
///                      Universe 6 (2020) no.1, 16
///
///    Effects of Atomic-Scale Electron Density Profile and a Fast and 
/// Efficient Iteration Algorithm for Matter Effect of Neutrino Oscillation
///
///            Mihai Horoi, Adam Zettel (Central Michigan U.)
///
///......................................................................
/// </pre>
///
/// \author Joao Coelho - jcoelho\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

#ifndef PMNS_Iter_H
#define PMNS_Iter_H
#include "PMNS_Fast.h"

namespace OscProb {

  class PMNS_Iter : public PMNS_Fast {

  public:

    PMNS_Iter();          ///< Constructor
    virtual ~PMNS_Iter(); ///< Destructor
    
    virtual void SetPrec(double prec);
    
  protected:

    /// Just use the vacuum to start
    virtual void SolveHam();
    
    /// Propagate through matter part
    virtual void PropMatter();

    virtual void SetExpVL(NuPath p);
    
    /// Reimplement full propagation
    virtual void Propagate();
    
    virtual void SplitPropagate(NuPath p);
    
    virtual void SetVacuumEigensystem();
    
    double fPrec;
    
    double fVL;
    complexD fExpVL;
    
  };

}
#endif
////////////////////////////////////////////////////////////////////////
