///////////////////////////////////////////////////////////////////////////////
/// Implementation of three-flavor neutrino oscillations in vacuum and matter,
/// incorporating quantum dissipation effects.
///
/// @author Alba Domi - alba.domi@fau.de
/// @author Joao Coelho - jcoelho@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#include <unsupported/Eigen/MatrixFunctions>

#include "PMNS_OQS.h"

#include <iostream>

using namespace std;
using namespace OscProb;

constexpr int SU3_DIM = PMNS_OQS::SU3_DIM;

//.............................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_OQS::PMNS_OQS()
    : PMNS_DensityMatrix(), fUM(3, vectorC(3, 0)), fHVMB(3, vectorC(3, 0)),
      fHGM(SU3_DIM, vectorD(SU3_DIM, 0)), fD(SU3_DIM, vectorD(SU3_DIM, 0)),
      fcos(SU3_DIM, vectorD(SU3_DIM, 1)), fa(SU3_DIM, 0), fR(SU3_DIM),
      fRt(SU3_DIM), fM(8, 8)
{
  SetPower(0);
  fBuiltDissipator = false;
}

//.............................................................................
///
/// Nothing to clean.
///
PMNS_OQS::~PMNS_OQS() {}

//.............................................................................
///
/// Set the power-law index for the energy dependence of decoherence parameters.
///
/// @param n - Power-law index.
///
void PMNS_OQS::SetPower(int n) { fPower = n; }

//.............................................................................
///
/// Set the value of decoherence parameter |a_i| in Gell-Mann representation.
///
/// @param i   - Index of the parameter in range [1, 8].
/// @param val - Value to assign.
///
void PMNS_OQS::SetDecoElement(int i, double val)
{
  if (1 > i || i >= SU3_DIM) {
    cerr << "WARNING: a_" << i << " is not valid. Doing nothing." << endl;
    return;
  }

  fBuiltDissipator *= (fa[i] == abs(val));
  fa[i] = abs(val);
}

//.............................................................................
///
/// Set angle between two decoherence vectors a_i and a_j.
///
/// @param i  - First index in GM representation in range [1, 8].
/// @param j  - Second index in GM representation in range [2, 8].
/// @param th - Angle in radians.
///
void PMNS_OQS::SetDecoAngle(int i, int j, double th)
{
  if (1 > i || i >= SU3_DIM || 1 > j || j >= SU3_DIM || i == j) {
    cerr << "WARNING: deco angle " << i << j << " not valid. Doing nothing."
         << endl;
    return;
  }

  double val = cos(th);
  fBuiltDissipator *= (fcos[i][j] == val);
  fcos[i][j] = val;
  fcos[j][i] = val;
}

//.............................................................................
///
/// Get the currently set power-law index for decoherence parameters.
///
/// @return Power-law index.
///
int PMNS_OQS::GetPower() { return fPower; }

//.............................................................................
///
/// Get the value of a_i decoherence parameter in GM representation.
///
/// @param i - Index of the parameter in range [1, 8].
///
/// @return Decoherence parameter a_i.
///
double PMNS_OQS::GetDecoElement(int i)
{
  assert(0 < i && i < SU3_DIM);
  return fa[i];
}

//.............................................................................
///
/// Get mixing angle between two decoherence vectors a_i and a_j.
///
/// @param i - First index in GM representation in range [1, 8].
/// @param j - Second index in GM representation in range [2, 8].
///
/// @return Decoherence angle between a_i and a_j.
///
double PMNS_OQS::GetDecoAngle(int i, int j)
{
  assert(0 < i && i < SU3_DIM && 0 < j && j < SU3_DIM && i != j);
  return acos(fcos[i][j]);
}

//.............................................................................
///
/// Get element of the Hamiltonian in GM representation.
///
/// @param i - Row index in range [0, 8].
/// @param j - Column index in range [0, 8].
///
/// @return Hamiltonian element i,j in GM representation.
///
double PMNS_OQS::GetHGM(int i, int j)
{
  assert(0 <= i && i < SU3_DIM && 0 <= j && j < SU3_DIM);
  return fHGM[i][j];
}

//.............................................................................
///
/// Get specific element of the dissipator matrix in the GM representation.
///
/// @param i Row index in range [0, 8].
/// @param j Column index in range [0, 8].
///
/// @return Dissipator element i,j.
///
double PMNS_OQS::GetDissipatorElement(int i, int j)
{
  assert(0 <= i && i < SU3_DIM && 0 <= j && j < SU3_DIM);
  return fD[i][j];
}

//.............................................................................
///
/// Reimplement from PMNS_Base flagging to rebuild Hamiltonian.
///
/// @param isNuBar - Flag for anti-neutrinos.
///
void PMNS_OQS::SetIsNuBar(bool isNuBar)
{
  fBuiltHms *= (fIsNuBar == isNuBar);
  PMNS_Base::SetIsNuBar(isNuBar);
}

//.............................................................................
///
/// Reimplement from PMNS_Base and save the PMNS matrix for later use.
///
void PMNS_OQS::BuildHms()
{
  if (fBuiltHms) return;
  SetVacuumEigensystem();
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) { fUM[i][j] = conj(fEvec[i][j]); }
  }
  PMNS_Base::BuildHms();
}

//.............................................................................
///
/// Build the effective Hamiltonian in VMB for a given propagation path.
///
/// @param p - Neutrino path segment.
///
void PMNS_OQS::BuildHVMB(NuPath p)
{
  // Set the neutrino path
  SetCurPath(p);

  double kr2GNe = kK2 * M_SQRT2 * kGf;
  kr2GNe *= fPath.density * fPath.zoa; // Matter potential in eV

  double Ve = 0;

  if (!fIsNuBar)
    Ve = kr2GNe;
  else
    Ve = -kr2GNe;

  BuildHms();

  for (int i = 0; i < 3; i++) {
    for (int j = i; j < 3; j++) {
      fHVMB[i][j] = conj(fUM[0][i]) * Ve * fUM[0][j];
      if (i > j) fHVMB[j][i] = conj(fHVMB[i][j]);
    }
  }

  double lv = 2 * kGeV2eV * fEnergy; // 2E in eV

  // add mass terms
  fHVMB[1][1] += fDm[1] / lv;
  fHVMB[2][2] += fDm[2] / lv;
}

//.............................................................................
///
/// Auxiliary function to convert Hermitian 3x3 matrices
/// to GM representation vectors.
///
/// @param A - The input Hermitian 3x3 matrix.
/// @param v - The output real 9-vector.
///
void get_GM(const matrixC& A, vectorD& v)
{
  v[0] = real(A[0][0] + A[1][1] + A[2][2]) / sqrt(6);
  v[3] = real(A[0][0] - A[1][1]) / 2;
  v[8] = real(A[0][0] + A[1][1] - 2. * A[2][2]) / sqrt(12);

  v[1] = real(A[0][1]);
  v[4] = real(A[0][2]);
  v[6] = real(A[1][2]);

  v[2] = -imag(A[0][1]);
  v[5] = -imag(A[0][2]);
  v[7] = -imag(A[1][2]);
}

//.............................................................................
///
/// Auxiliary function to convert GM representation vectors
/// to Hermitian 3x3 matrices.
///
/// @param v - The input real 9-vector.
/// @param A - The output Hermitian 3x3 matrix.
///
void get_SU3(const vectorD& v, matrixC& A)
{
  A[0][0] = v[0] * sqrt(2 / 3.) + v[3] + v[8] / sqrt(3);
  A[1][1] = v[0] * sqrt(2 / 3.) - v[3] + v[8] / sqrt(3);
  A[2][2] = v[0] * sqrt(2 / 3.) - 2 * v[8] / sqrt(3);

  A[1][0] = complexD(v[1], v[2]);
  A[2][0] = complexD(v[4], v[5]);
  A[2][1] = complexD(v[6], v[7]);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < i; j++) { A[j][i] = conj(A[i][j]); }
  }
}

//.............................................................................
///
/// Auxiliary function to convert Hermitian 3x3 matrices to
/// GM representation 9x9 operators.
///
/// @param A - The input Hermitian 3x3 matrix.
/// @param B - The output real 9x9 matrix.
///
void get_GMOP(const matrixC& A, matrixD& B)
{
  B[1][2] = real(A[0][0] - A[1][1]);
  B[1][3] = 2 * imag(A[0][1]);
  B[1][4] = -imag(A[1][2]);
  B[1][5] = -real(A[1][2]);
  B[1][6] = -imag(A[0][2]);
  B[1][7] = -real(A[0][2]);

  B[2][3] = 2 * real(A[0][1]);
  B[2][4] = real(A[1][2]);
  B[2][5] = -imag(A[1][2]);
  B[2][6] = -real(A[0][2]);
  B[2][7] = imag(A[0][2]);

  B[3][4] = -imag(A[0][2]);
  B[3][5] = -real(A[0][2]);
  B[3][6] = imag(A[1][2]);
  B[3][7] = real(A[1][2]);

  B[4][5] = real(A[0][0] - A[2][2]);
  B[4][6] = -imag(A[0][1]);
  B[4][7] = real(A[0][1]);
  B[4][8] = sqrt(3) * imag(A[0][2]);

  B[5][6] = -real(A[0][1]);
  B[5][7] = -imag(A[0][1]);
  B[5][8] = sqrt(3) * real(A[0][2]);

  B[6][7] = real(A[1][1] - A[2][2]);
  B[6][8] = sqrt(3) * imag(A[1][2]);

  B[7][8] = sqrt(3) * real(A[1][2]);

  for (int i = 1; i < SU3_DIM; ++i) {
    for (int j = i + 1; j < SU3_DIM; ++j) { B[j][i] = -B[i][j]; }
  }
}

//.............................................................................
///
/// Build the effective Hamiltonian in Gell-Mann representation.
///
/// @param p - Neutrino path segment.
///
void PMNS_OQS::BuildHGM(NuPath p)
{
  BuildHVMB(p);
  get_GMOP(fHVMB, fHGM);
}

//.............................................................................
///
/// Build dissipator in Gell-Mann representation.
///
void PMNS_OQS::BuildDissipator()
{
  if (fBuiltDissipator) return;

  double aa[SU3_DIM][SU3_DIM];
  for (int i = 1; i < SU3_DIM; i++) {
    for (int j = i; j < SU3_DIM; j++) {
      aa[i][j] = fa[i] * fa[j];
      if (i == 8) aa[i][j] *= sqrt(3);
      if (j == 8) aa[i][j] *= sqrt(3);
      if (i < j) aa[i][j] *= fcos[i][j];
    }
  }
  double sum12   = aa[1][1] + aa[2][2];
  double sum45   = aa[4][4] + aa[5][5];
  double sum67   = aa[6][6] + aa[7][7];
  double gamma21 = aa[3][3];
  double gamma31 = (aa[3][3] + aa[8][8] + 2 * aa[3][8]) / 4;
  double gamma32 = (aa[3][3] + aa[8][8] - 2 * aa[3][8]) / 4;

  fD[1][1] = gamma21 + aa[2][2] + (sum45 + sum67) / 4;
  fD[2][2] = gamma21 + aa[1][1] + (sum45 + sum67) / 4;
  fD[3][3] = sum12 + (sum45 + sum67) / 4;

  fD[4][4] = gamma31 + aa[5][5] + (sum12 + sum67) / 4;
  fD[5][5] = gamma31 + aa[4][4] + (sum12 + sum67) / 4;
  fD[6][6] = gamma32 + aa[7][7] + (sum12 + sum45) / 4;
  fD[7][7] = gamma32 + aa[6][6] + (sum12 + sum45) / 4;

  fD[8][8] = (sum45 + sum67) * 3 / 4;
  fD[3][8] = (sum45 - sum67) * sqrt(3) / 4;

  fD[1][2] = -aa[1][2];
  fD[1][3] = -aa[1][3];
  fD[2][3] = -aa[2][3];
  fD[4][5] = -aa[4][5];
  fD[6][7] = -aa[6][7];

  fD[1][8] = (aa[4][6] + aa[5][7]) * sqrt(3) / 2;
  fD[2][8] = (aa[5][6] - aa[4][7]) * sqrt(3) / 2;

  fD[4][6] = -(aa[4][6] - 2 * aa[1][8] - 3 * aa[5][7]) / 4;
  fD[4][7] = -(aa[4][7] + 2 * aa[2][8] + 3 * aa[5][6]) / 4;
  fD[5][6] = -(aa[5][6] - 2 * aa[2][8] + 3 * aa[4][7]) / 4;
  fD[5][7] = -(aa[5][7] - 2 * aa[1][8] - 3 * aa[4][6]) / 4;

  fD[1][4] = -(aa[1][4] - 3 * (aa[2][5] - aa[3][6]) + aa[6][8]) / 4;
  fD[1][5] = -(aa[1][5] + 3 * (aa[2][4] + aa[3][7]) + aa[7][8]) / 4;
  fD[1][6] = -(aa[1][6] + 3 * (aa[2][7] - aa[3][4]) + aa[4][8]) / 4;
  fD[1][7] = -(aa[1][7] - 3 * (aa[2][6] + aa[3][5]) + aa[5][8]) / 4;
  fD[2][4] = -(aa[2][4] + 3 * (aa[1][5] - aa[3][7]) - aa[7][8]) / 4;
  fD[2][5] = -(aa[2][5] - 3 * (aa[1][4] - aa[3][6]) + aa[6][8]) / 4;
  fD[2][6] = -(aa[2][6] - 3 * (aa[1][7] + aa[3][5]) + aa[5][8]) / 4;
  fD[2][7] = -(aa[2][7] + 3 * (aa[1][6] + aa[3][4]) - aa[4][8]) / 4;
  fD[3][4] = -(aa[3][4] - 3 * (aa[1][6] - aa[2][7]) + aa[4][8]) / 4;
  fD[3][5] = -(aa[3][5] - 3 * (aa[1][7] + aa[2][6]) + aa[5][8]) / 4;
  fD[3][6] = -(aa[3][6] + 3 * (aa[1][4] + aa[2][5]) - aa[6][8]) / 4;
  fD[3][7] = -(aa[3][7] + 3 * (aa[1][5] - aa[2][4]) - aa[7][8]) / 4;

  fD[4][8] = -(aa[1][6] - aa[2][7] + aa[3][4] + aa[4][8]) * sqrt(3) / 4;
  fD[5][8] = -(aa[1][7] + aa[2][6] + aa[3][5] + aa[5][8]) * sqrt(3) / 4;
  fD[6][8] = -(aa[1][4] + aa[2][5] - aa[3][6] + aa[6][8]) * sqrt(3) / 4;
  fD[7][8] = -(aa[1][5] - aa[2][4] - aa[3][7] + aa[7][8]) * sqrt(3) / 4;

  for (int j = 1; j < SU3_DIM; j++) {
    for (int k = j; k < SU3_DIM; k++) {
      fD[j][k] = -fD[j][k] * kGeV2eV;
      fD[k][j] = fD[j][k];
    }
  }

  fBuiltDissipator = true;
}

//.............................................................................
///
/// Build the matrix M used in evolution equations.
///
/// @param p - Neutrino path segment.
///
void PMNS_OQS::BuildM(NuPath p)
{
  BuildHGM(p);
  BuildDissipator();
  double energyCorr = pow(fEnergy, fPower);
  for (int k = 1; k < SU3_DIM; ++k) {
    for (int j = 1; j < SU3_DIM; ++j) {
      fM(k - 1, j - 1) = fHGM[k][j] + fD[k][j] * energyCorr;
    }
  }
}

//.............................................................................
///
/// Build Gell-Mann representation of density matrix.
///
void PMNS_OQS::BuildR() { get_GM(fRho, fR); }

//.............................................................................
///
/// Update density matrix from Gell-Mann representation.
///
void PMNS_OQS::UpdateRho() { get_SU3(fR, fRho); }

//.............................................................................
///
/// Rotate the density matrix to or from the mass basis.
///
/// @param to_mass - true if to mass basis.
///
void PMNS_OQS::RotateState(bool to_mass)
{
  BuildHms();
  if (!to_mass) UpdateRho();
  PMNS_DensityMatrix::RotateState(to_mass, fUM);
  if (to_mass) BuildR();
}

//.............................................................................
///
/// Reimplemented from PMNS_Base to rotate to mass basis before propagation.
///
void PMNS_OQS::Propagate()
{
  RotateState(true);
  PMNS_Base::Propagate();
  RotateState(false);
}

//.............................................................................
///
/// Propagate the current neutrino state through a given path.
///
/// @param p - A neutrino path segment.
///
void PMNS_OQS::PropagatePath(NuPath p)
{
  BuildM(p);

  // Solve for eigensystem
  SolveHam();

  // Convert path length to eV
  double lengthIneV = kKm2eV * p.length;

  fM *= lengthIneV;
  fM = fM.exp();

  fRt[0] = fR[0];

  for (int i = 1; i < SU3_DIM; ++i) {
    fRt[i] = 0;
    for (int j = 1; j < SU3_DIM; ++j) { fRt[i] += fM(i - 1, j - 1) * fR[j]; }
  }
  fR = fRt;
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
/// @return Neutrino oscillation probabilities.
///
matrixD PMNS_OQS::ProbMatrix(int nflvi, int nflvf)
{
  assert(nflvi <= fNumNus && nflvi >= 0);
  assert(nflvf <= fNumNus && nflvf >= 0);

  // Output probabilities
  matrixD probs(nflvi, vectorD(nflvf));

  // List of states
  vector<vectorD> allstates(nflvi, vectorD(9));

  // Reset all initial states
  for (int i = 0; i < nflvi; i++) {
    ResetToFlavour(i);
    RotateState(true);
    allstates[i] = fR;
  }

  // Propagate all states in parallel
  for (int i = 0; i < int(fNuPaths.size()); i++) {
    for (int flvi = 0; flvi < nflvi; flvi++) {
      fR = allstates[flvi];
      PropagatePath(fNuPaths[i]);
      allstates[flvi] = fR;
    }
  }

  // Get all probabilities
  for (int flvi = 0; flvi < nflvi; flvi++) {
    fR = allstates[flvi];
    RotateState(false);
    for (int flvj = 0; flvj < nflvf; flvj++) { probs[flvi][flvj] = P(flvj); }
  }

  return probs;
}

////////////////////////////////////////////////////////////////////////
