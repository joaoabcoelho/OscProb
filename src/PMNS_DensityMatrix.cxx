///////////////////////////////////////////////////////////////////////////////
//
// Base class for methods using density matrices
//
// This  class inherits from the PMNS_Base class
//
// jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

//#include <cassert>

#include "MatrixDecomp/zheevh3.h"

#include "PMNS_DensityMatrix.h"
#include "exceptions.h"

using namespace OscProb;

using namespace std;

//.............................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_DensityMatrix::PMNS_DensityMatrix()
    : PMNS_Fast(), fRho(3, vectorC(3, 0)), fMBuffer(3, vectorC(3, 0))
{
}

//.............................................................................
///
/// Nothing to clean.
///
PMNS_DensityMatrix::~PMNS_DensityMatrix() {}

//.............................................................................
///
/// Rotate the density matrix to or from a basis given by U
///
/// @param to_basis - true if to the given basis
/// @param U - the rotation matrix to the given basis
///
void PMNS_DensityMatrix::RotateState(bool to_basis, matrixC U)
{
  // buffer = rho . U
  for (int i = 0; i < fNumNus; i++) {
    for (int j = 0; j < fNumNus; j++) {
      fMBuffer[i][j] = 0;
      for (int k = 0; k < fNumNus; k++) {
        if (to_basis)
          fMBuffer[i][j] += fRho[i][k] * U[k][j];
        else
          fMBuffer[i][j] += fRho[i][k] * conj(U[j][k]);
      }
    }
  }

  // rho = U^\dagger . buffer = U^\dagger . rho . U
  // Final matrix is Hermitian, so copy upper to lower triangle
  for (int i = 0; i < fNumNus; i++) {
    for (int j = i; j < fNumNus; j++) {
      fRho[i][j] = 0;
      for (int k = 0; k < fNumNus; k++) {
        if (to_basis)
          fRho[i][j] += conj(U[k][i]) * fMBuffer[k][j];
        else
          fRho[i][j] += U[i][k] * fMBuffer[k][j];
      }
      if (j > i) fRho[j][i] = conj(fRho[i][j]);
    }
  }
}

//.............................................................................
///
/// Reset the neutrino state back to a pure flavour where it starts
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
/// @param flv - The neutrino starting flavour.
///
void PMNS_DensityMatrix::ResetToFlavour(int flv)
{
  THROW_ON_LOGIC_ERR(0 <= flv && flv < fNumNus, "flavour not valid,", flv);

  PMNS_Base::ResetToFlavour(flv);

  for (int i = 0; i < fNumNus; ++i) {
    for (int j = 0; j < fNumNus; ++j) {
      if (i == flv && i == j)
        fRho[i][j] = one;
      else
        fRho[i][j] = zero;
    }
  }
}

//.............................................................................
///
/// Compute oscillation probability of flavour flv from current state
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
/// @param flv - The neutrino final flavour.
///
/// @return Neutrino oscillation probability
///
double PMNS_DensityMatrix::P(int flv)
{
  THROW_ON_LOGIC_ERR(0 <= flv && flv < fNumNus, "flavour not valid.", flv);
  return abs(fRho[flv][flv]);
}

//.............................................................................
///
/// Set the density matrix from a pure state
///
/// @param nu_in - The neutrino initial state in flavour basis.
///
void PMNS_DensityMatrix::SetPureState(vectorC nu_in)
{
  THROW_ON_LOGIC_ERR(nu_in.size() == fNumNus, "Invalid state vector size.",
                     nu_in.size());

  for (int i = 0; i < fNumNus; i++) {
    for (int j = 0; j < fNumNus; j++) {
      fRho[i][j] = conj(nu_in[i]) * nu_in[j];
    }
  }
}

//.............................................................................
///
/// Set the density matrix from an arbitrary initial matrix
///
/// The initial matrix needs to satisfy the basic conditions:
/// * Hermitian matrix
/// * Trace = 1
/// * Positive semidefinite
///
/// @param rho_in - The initial density matrix in flavour basis.
///
void PMNS_DensityMatrix::SetInitialRho(matrixC rho_in)
{
  THROW_ON_LOGIC_ERR(rho_in.size() == fNumNus && rho_in[0].size() == fNumNus,
                     "Invalid density matrix size", rho_in.size(),
                     rho_in[0].size());

  double eps = 1e-12;

  double   trace = 0;
  complexD rho_array[3][3];
  for (int i = 0; i < fNumNus; i++) {
    trace += abs(rho_in[i][i]);
    for (int j = i; j < fNumNus; j++) {
      // Ensure rho_in is hermitian
      THROW_ON_LOGIC_ERR(abs(rho_in[i][j] - conj(rho_in[j][i])) < eps,
                         "Matrix not hermitian", rho_in[i][j], rho_in[j][i]);
      rho_array[i][j] = rho_in[i][j];
    }
  }

  // Ensure rho_in has trace 1
  THROW_ON_LOGIC_ERR(abs(trace - 1) < eps, "Total probability not 1", trace);

  double   fEvalGLoBES[3];
  complexD fEvecGLoBES[3][3];

  // Find the eigenvalues of rho_in
  zheevh3(rho_array, fEvecGLoBES, fEvalGLoBES);

  // Ensure rho_in is positive definite
  for (int i = 0; i < fNumNus; i++)
    THROW_ON_LOGIC_ERR(fEvalGLoBES[i] >= -eps, "Matrix not positive definite",
                       fEvalGLoBES[i], i);

  fRho = rho_in;
}

//.............................................................................
///
/// Compute the probability matrix for the first nflvi and nflvf states.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
/// @param nflvi - The number of initial flavours in the matrix.
/// @param nflvf - The number of final flavours in the matrix.
///
/// @return Neutrino oscillation probabilities
///
matrixD PMNS_DensityMatrix::ProbMatrix(int nflvi, int nflvf)
{
  THROW_ON_INVALID_ARG(nflvi >= 0, "", nflvi);
  THROW_ON_INVALID_ARG(nflvi <= fNumNus, "", nflvi, fNumNus);
  THROW_ON_INVALID_ARG(nflvf >= 0, "", nflvf);
  THROW_ON_INVALID_ARG(nflvf <= fNumNus, "", nflvf, fNumNus);

  // Output probabilities
  matrixD probs(nflvi, vectorD(nflvf));

  // List of states
  vector<matrixC> allstates(nflvi, matrixC(fNumNus, vectorC(fNumNus)));

  // Reset all initial states
  for (int i = 0; i < nflvi; i++) {
    ResetToFlavour(i);
    allstates[i] = fRho;
  }

  // Propagate all states in parallel
  for (int i = 0; i < int(fNuPaths.size()); i++) {
    for (int flvi = 0; flvi < nflvi; flvi++) {
      fRho = allstates[flvi];
      PropagatePath(fNuPaths[i]);
      allstates[flvi] = fRho;
    }
  }

  // Get all probabilities
  for (int flvi = 0; flvi < nflvi; flvi++) {
    for (int flvj = 0; flvj < nflvf; flvj++) {
      probs[flvi][flvj] = abs(allstates[flvi][flvj][flvj]);
    }
  }

  return probs;
}

////////////////////////////////////////////////////////////////////////
