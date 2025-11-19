///////////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_DensityMatrix
///
/// \brief Base class for methods based on density matrices
///
/// This class expands the PMNS_Fast class to allow the use of density matrices
///
/// \sa PMNS_Fast
///
/// \author jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_DensityMatrix_H
#define PMNS_DensityMatrix_H

#include "PMNS_Fast.h"

namespace OscProb {

  class PMNS_DensityMatrix : public PMNS_Fast {
    public:
      PMNS_DensityMatrix();          ///< Constructor
      virtual ~PMNS_DensityMatrix(); ///< Destructor

      /// Compute the probability matrix
      using PMNS_Base::ProbMatrix;
      virtual matrixD ProbMatrix(int nflvi, int nflvf);

      virtual void SetIsOscProbAvg(bool isOscProbAvg)
      {
        fIsOscProbAvg = true;
      } ///< Deactivate Maltoni

    protected:
      // Resetting and propagating
      /// Reset neutrino state to pure flavour flv
      virtual void ResetToFlavour(int flv);
      /// Set the density matrix from a pure state
      virtual void SetPureState(vectorC nu_in);
      /// Set the density matrix from an arbitrary state
      virtual void SetInitialRho(matrixC rho_in);

      /// Propagate neutrino through a single path
      virtual void PropagatePath(NuPath p) = 0;

      /// Return the probability of final state in flavour flv
      virtual double P(int flv);

      /// Rotate rho to/from a given basis
      virtual void RotateState(bool to_basis, matrixC U);

      matrixC fRho; ///< The neutrino density matrix state

      matrixC fMBuffer; ///< Some memory buffer for matrix operations
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
