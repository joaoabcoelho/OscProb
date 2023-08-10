///////////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_Iter
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework.
///
/// Based on an approximate iterative model from:
///
/// https://doi.org/10.3390/universe6010016
///
/// <pre>
///
///............................................................................
///
///                      Universe 6 (2020) no.1, 16
///
///    Effects of Atomic-Scale Electron Density Profile and a Fast and
/// Efficient Iteration Algorithm for Matter Effect of Neutrino Oscillation
///
///            Mihai Horoi, Adam Zettel (Central Michigan U.)
///
///............................................................................
/// </pre>
///
/// \author Joao Coelho - jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_Iter_H
#define PMNS_Iter_H

#include "PMNS_Fast.h"

namespace OscProb {

  class PMNS_Iter : public PMNS_Fast {

    public:

      PMNS_Iter();          ///< Constructor
      virtual ~PMNS_Iter(); ///< Destructor

      /// Set the iterative precision
      virtual void SetPrec(double prec);

    protected:

      /// Just use the vacuum to start
      virtual void SolveHam();

      /// Propagate through matter part
      virtual void PropMatter();

      /// Set the matter propagation term
      virtual void SetExpVL(NuPath p);

      /// Reimplement propagation
      virtual void PropagatePath(NuPath p);

      /// Iterative precision
      double fPrec;

      double fVL; ///< Matter potential
      complexD fExpVL; ///< Matter phase shift

      double fPrevEnergy;

  };

}

#endif

///////////////////////////////////////////////////////////////////////////////
