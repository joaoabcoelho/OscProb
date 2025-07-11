#ifndef PMNS_OQS_H
#define PMNS_OQS_H

#include <Eigen/Core>

#include "PMNS_Fast.h"

namespace OscProb {

  class PMNS_OQS : public PMNS_Fast {
    public:
      PMNS_OQS();          ///< Constructor
      virtual ~PMNS_OQS(); ///< Destructor

      virtual void SetParameterisation(int par);
      virtual void SetPhi(int i, double val);
      virtual void SetDissipatorElement(int i, int j, double val,
                                        bool print = false);

      /// Compute the probability matrix
      using PMNS_Base::ProbMatrix;
      virtual matrixD ProbMatrix(int nflvi, int nflvf);
      virtual double  Prob(int flvi, int flvf);
      virtual double  Prob(int flvi, int flvf, double energy);

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

      // Resetting and propagating
      virtual void ResetToFlavour(
          int flv); ///< Reset neutrino state to pure flavour flv
      virtual void SetPureState(
          vectorC nu_in); ///< Set the density matrix from a pure state
      virtual void PropagatePath(
          NuPath p); ///< Propagate neutrino through a single path
      virtual double P(
          int flv); ///< Return the probability of final state in flavour flv

      int      fParameterisation;
      double   fPhi[2]; ///< Majorana phases
      complexD fR[9];
      complexD fRt[9];

      matrixC fRho; ///< The neutrino density matrix state
      matrixC fHeff;
      matrixC fHGM;

      std::vector<matrixC> fGM; ///< 3x3 Gell-Mann matrices: they are 9

      matrixC fD; ///< Off-diagonal, 9x9 dissipator
      matrixC fM; ///< M

      Eigen::MatrixXcd fMd;
      Eigen::MatrixXcd fMd8;
      Eigen::MatrixXcd fMEvec;

      vectorC fEval; ///< Eigenvalues of the Hamiltonian
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
