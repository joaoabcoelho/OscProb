///////////////////////////////////////////////////////////////////////////////
//
// Implementation of Non Unitarity neutrino Mixing (NUNM) of neutrinos in matter
// in a three-neutrino framework.
//
// This  class inherits from the PMNS_Fast class.
//
// Authors:
// cerisy\cppm.in2p3.fr
// jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#include "PMNS_NUNM.h"
#include "Definitions.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>

using namespace OscProb;
using namespace std;

//.............................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
/// scale 0 is low-scale NUNM model
/// scale 1 is high-scale NUNM model
///
PMNS_NUNM::PMNS_NUNM(int scale) : PMNS_Fast()
{
  fscale = scale;
  SetStdPath();
  SetNUNM(0., 0., 0., 0., 0., 0.);
  SetFracVnc(1.0);
  InitMatrix();
}

//.............................................................................
///
/// Nothing to clean.
///
PMNS_NUNM::~PMNS_NUNM() {}

//.............................................................................
///
/// Set all NUNM parameters at once.
///
/// This will check if value is changing to k11p track of whether
/// the eigensystem n11ds to be recomputed.
///
/// @param alpha_11      - The real parameter alpha_11
/// @param alpha_21     - The absolute value of the complex parameter alpha_21
/// @param alpha_31    - The absolute value of the complex parameter alpha_31
/// @param alpha_22    - The real parameter alpha_22
/// @param alpha_32   - The absolute value of the complex parameter alpha_32
/// @param alpha_33  - The real parameter alpha_33
///
void PMNS_NUNM::SetNUNM(double alpha_11, double alpha_21, double alpha_31,
                        double alpha_22, double alpha_32, double alpha_33)
{
  SetAlpha(0, 0, alpha_11, 0);
  SetAlpha(1, 1, alpha_22, 0);
  SetAlpha(2, 2, alpha_33, 0);

  SetAlpha(1, 0, alpha_21, 0);
  SetAlpha(2, 0, alpha_31, 0);
  SetAlpha(2, 1, alpha_32, 0);

  // upper part fixed to zero from https://arxiv.org/pdf/2309.16942.pdf
}

void PMNS_NUNM::InitMatrix()
{
  Ham.setZero();
  V.setZero();
}

//.............................................................................
///
/// Get any given NUNM parameter.
///
/// Indexes are:\n
/// - 0, 1, 2
///
/// Requires that i > j. Will notify you if input is wrong.
/// If i > j, will assume reverse order and swap i and j.
///
/// @param i  - The alpha row index
/// @param j  - The alpha column index
/// @param val   - The absolute value of the parameter
/// @param phase - The complex phase of the parameter in radians
///
void PMNS_NUNM::SetAlpha(int i, int j, double val, double phase)
{
  if (i < j) {
    cerr << "WARNING: First argument should be larger or equal to second "
            "argument"
         << endl
         << "Setting reverse order (Alpha_" << j << i << "). " << endl;
    int temp = i;
    i        = j;
    j        = temp;
  }
  if (i < 0 || i > 2 || j > i || j > 2) {
    cerr << "WARNING: Alpha_" << i << j << " not valid for " << fNumNus
         << " neutrinos. Doing nothing." << endl;
    return;
  }

  complexD h = val;

  if (i == j) { h = 1. + val; }
  else {
    h *= complexD(cos(phase), sin(phase));
  }

  bool isSame = (Alpha(i, j) == h);

  if (!isSame) ClearCache();

  fGotES *= isSame;

  Alpha(i, j) = h;
}

//.............................................................................
///
/// Get any given NUNM parameter.
///
/// Indexes are:\n
/// - 0, 1, 2
///
/// Requires that i > j. Will notify you if input is wrong.
/// If i > j, will assume reverse order and swap i and j.
///
/// @param i  - The alpha row index
/// @param j  - The alpha column index
///
complexD PMNS_NUNM::GetAlpha(int i, int j)
{
  if (i < j) {
    cerr << "WARNING: First argument should be smaller or equal to second "
            "argument"
         << endl
         << "Setting reverse order (Alpha_" << j << i << "). " << endl;
    int temp = i;
    i        = j;
    j        = temp;
  }
  if (i < 0 || i > 2 || j < i || j > 2) {
    cerr << "WARNING: Eps_" << i << j << " not valid for " << fNumNus
         << " neutrinos. Returning 0." << endl;
    return zero;
  }

  return Alpha(i, j);
}

//.............................................................................
///
/// Set alpha_11 parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a - The value of the real parameter alpha_11
///
void PMNS_NUNM::SetAlpha_11(double a) { SetAlpha(0, 0, a, 0); }

//.............................................................................
///
/// Set alpha_22 parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a - The value of the real parameter alpha_22
///
void PMNS_NUNM::SetAlpha_22(double a) { SetAlpha(1, 1, a, 0); }

//.............................................................................
///
/// Set alpha_33 parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a - The value of the real parameter alpha_33
///
void PMNS_NUNM::SetAlpha_33(double a) { SetAlpha(2, 2, a, 0); }

//.............................................................................
///
/// Set alpha_21 parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a   - The absolute value of the parameter alpha_21
/// @param phi - The phase of the complex parameter alpha_21 in radians
///
void PMNS_NUNM::SetAlpha_21(double a, double phi) { SetAlpha(1, 0, a, phi); }

//.............................................................................
///
/// Set alpha_31 parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a   - The absolute value of the parameter alpha_31
/// @param phi - The phase of the complex parameter alpha_31 in radians
///
void PMNS_NUNM::SetAlpha_31(double a, double phi) { SetAlpha(2, 0, a, phi); }

//.............................................................................
///
/// Set alpha_32 parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a   - The absolute value of the parameter alpha_32
/// @param phi - The phase of the complex parameter alpha_32 in radians
///
void PMNS_NUNM::SetAlpha_32(double a, double phi) { SetAlpha(2, 1, a, phi); }

//.............................................................................
///
/// Allows to disable the neutron induced matter effects
/// This factor represents what fraction of the neutron matter potentiel
/// is included in the NUNM model
///
/// @param f -> 0 no neutron matter potential
/// @param f -> 1 std neutron matter potential
///
void PMNS_NUNM::SetFracVnc(double f)
{
  bool isSame = (fracVnc == f);

  if (!isSame) ClearCache();

  fGotES *= isSame;

  fracVnc = f;
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
matrixD PMNS_NUNM::ProbMatrix(int nflvi, int nflvf)
{
  assert(nflvi <= fNumNus && nflvi >= 0);
  assert(nflvf <= fNumNus && nflvf >= 0);

  // Output probabilities
  matrixD probs(nflvi, vectorD(nflvf));

  // List of states
  matrixC allstates(nflvi, vectorC(fNumNus));

  // Reset all initial states
  for (int j = 0; j < nflvi; j++) {
    ResetToFlavour(j);
    fNuState     = ApplyAlphaDagger(fNuState);
    allstates[j] = fNuState;
  }

  // Propagate all states in parallel
  for (int i = 0; i < int(fNuPaths.size()); i++) {
    for (int flvi = 0; flvi < nflvi; flvi++) {
      fNuState = allstates[flvi];
      PropagatePath(fNuPaths[i]);
      allstates[flvi] = fNuState;
    }
  }

  // Get all probabilities
  for (int flvi = 0; flvi < nflvi; flvi++) {
    allstates[flvi] = ApplyAlpha(allstates[flvi]);
    for (int flvj = 0; flvj < nflvf; flvj++) {
      probs[flvi][flvj] = norm(allstates[flvi][flvj]);
    }
  }

  return probs;
}

//.............................................................................
///
/// Apply Alpha dagger for calculation of probability before propagation
/// P = | A * exp(-iHL) * A+ |^2
///
vectorC PMNS_NUNM::ApplyAlphaDagger(vectorC fState)
{
  fNuStateBuffer = fState;
  for (int i = 0; i < fNumNus; ++i) {
    fState[i] = conj(Alpha(0, i)) * fNuStateBuffer[0] +
                conj(Alpha(1, i)) * fNuStateBuffer[1] +
                conj(Alpha(2, i)) * fNuStateBuffer[2];
  }
  return fState;
}

//.............................................................................
///
/// Apply Alpha for calculation of probability after propagation
/// P = | A * exp(-iHL) * A+ |^2
///
vectorC PMNS_NUNM::ApplyAlpha(vectorC fState)
{
  fNuStateBuffer = fState;
  for (int i = 0; i < 3; ++i) {
    fState[i] = Alpha(i, 0) * fNuStateBuffer[0] +
                Alpha(i, 1) * fNuStateBuffer[1] +
                Alpha(i, 2) * fNuStateBuffer[2];
  }
  return fState;
}

//.............................................................................
///
/// Propagate neutrino state through full path
///
void PMNS_NUNM::Propagate()
{
  fNuState = ApplyAlphaDagger(fNuState);
  PMNS_Base::Propagate();
  fNuState = ApplyAlpha(fNuState);
}

//.............................................................................
///
/// Propagate the current neutrino state through a given path
/// @param p - A neutrino path segment
/// apply Alpha X Alpha~ transformation to get the probability
/// apply norm to Alpha in high scale scenario
///
void PMNS_NUNM::PropagatePath(NuPath p)
{
  if (fscale == 1) { // test <------ // normalise mixing matrix in high scale
                     // scenario to ensure completeness
    X = Alpha * Alpha.adjoint();
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < i + 1; ++j) {
        Alpha(i, j) *=
            1 / std::sqrt(X(i, i).real()); // Scale by the inverse square root
                                           // of the diagonal elements of X
      }
    }
  }
  PMNS_Base::PropagatePath(p);
}

//.............................................................................
///
/// Build the full Hamiltonian in matter --> from sterile
///
/// Here we divide the mass squared matrix Hms by the 2E
/// to obtain the vacuum Hamiltonian in eV. Then, the matter
/// potential is added to the electron and NC components.
///
/// the neutron matter effect is added to the diag of V
/// depending on the scenario the writting of the Hamiltonian
/// in the tilde basis is different
///
/// solving the sytem in the tilde basis: H = U·Dm·U~ + Alpha~·V·Alpha
/// in high scale scenario : H = U·Dm·U~ + Alpha^(-1)·V·Alpha~^(-1)
/// inspired from https://arxiv.org/pdf/2301.12960.pdf and
/// tilde basis in eq B5.
/// https://cds.cern.ch/record/1032544/files/PhysRevD.76.093005.pdf
///
///
///
void PMNS_NUNM::UpdateHam()
{
  double rho = fPath.density;
  double zoa = fPath.zoa;

  double lv = 2 * kGeV2eV * fEnergy; // 2*E in eV

  double kr2GNe =
      kK2 * M_SQRT2 * kGf * rho * zoa; // Electron matter potential in eV
  double kr2GNn = kK2 * M_SQRT2 * kGf * rho * (1. - zoa) / 2. *
                  fracVnc; // Neutron matter potential in eV

  // Finish build Hamiltonian in matter with dimension of eV

  for (int i = 0; i < fNumNus; i++) {
    if (!fIsNuBar) { V(i, i) = -1. * kr2GNn; }
    else {
      V(i, i) = kr2GNn;
    }
    for (int j = 0; j < fNumNus; j++) {
      if (!fIsNuBar)
        Ham(i, j) = fHms[i][j] / lv;
      else
        Ham(i, j) = conj(fHms[i][j]) / lv;
    }
  }

  if (!fIsNuBar) { V(0, 0) += kr2GNe; }
  else {
    V(0, 0) -= kr2GNe;
  }

  if (fscale == 0) {
    Ham += Alpha.adjoint() * V * Alpha;
  } // low scale scenario with mixing matrix part of larger unitary matrix
  else if (fscale == 1) {
    Ham += Alpha.inverse() * V * (Alpha.adjoint()).inverse();
  } // high scale scenario with non unitary mixing matrix

  for (int i = 0; i < fNumNus; i++) {
    for (int j = 0; j < fNumNus; j++) { fHam[i][j] = Ham(i, j); }
  }
}
