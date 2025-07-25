#ifndef PMNS_OQS_H
#define PMNS_OQS_H

#include <Eigen/Core>

#include "PMNS_DensityMatrix.h"

namespace OscProb {

  class PMNS_OQS : public PMNS_DensityMatrix {
    public:
      PMNS_OQS();          ///< Constructor
      virtual ~PMNS_OQS(); ///< Destructor

      virtual void SetParameterisation(int par);
      virtual void SetPhi(int i, double val);
      virtual void SetDissipatorElement(int i, int j, double val,
                                        bool print = false);

    protected:
      virtual void InitializeVectors();
      virtual void SetHeff(NuPath p);
      virtual void SetHGM();
      virtual void SetM();
      virtual void RotateState(bool to_mass); ///< Rotate rho to/from mass basis
      virtual void ChangeBaseToGM();
      virtual void ChangeBaseToSU3();

      /// Specialized solver for NxN matrices
      template <typename T> void SolveEigenSystem();

      virtual void Diagonalise();

      /// Propagate neutrino through a single path
      virtual void PropagatePath(NuPath p);

      int      fParameterisation;
      double   fPhi[2]; ///< Majorana phases
      complexD fR[9];
      complexD fRt[9];

      matrixC fHeff;
      matrixC fHGM;

      std::vector<matrixC> fGM; ///< 3x3 Gell-Mann matrices: they are 9

      matrixC fD; ///< Off-diagonal, 9x9 dissipator
      matrixC fM; ///< M

      Eigen::MatrixXcd fMd;
      Eigen::MatrixXcd fMd8;
      Eigen::MatrixXcd fMEvec;

      vectorC fEvalC; ///< Complex eigenvalues
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
