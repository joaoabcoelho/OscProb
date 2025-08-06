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
    : PMNS_DensityMatrix(), fPhi{}, fR(SU3_DIM), fRt(SU3_DIM), fa(SU3_DIM, 0), fM(8, 8),
      fD(SU3_DIM, vectorD(SU3_DIM, 0)), fcos(SU3_DIM, vectorD(SU3_DIM, 1)), fHGM(SU3_DIM, vectorD(SU3_DIM, 0)),
      fHeff(3, vectorC(3, 0)), fUM(3, vectorC(3, 0))
{
  SetParameterisation(1);
  SetPower(0);
  SetPhi(1, 0);
  SetPhi(2, 0);
  fBuiltDissipator = false;
}

//.............................................................................
///
/// Nothing to clean.
///
PMNS_OQS::~PMNS_OQS() {}

void PMNS_OQS::SetIsNuBar(bool isNuBar)
{
  fBuiltHms *= (fIsNuBar == isNuBar);
  PMNS_Base::SetIsNuBar(isNuBar);
}

void PMNS_OQS::SetParameterisation(int param) { fParameterisation = param; }

// set Heff in vacuum-mass basis
void PMNS_OQS::SetHeff(NuPath p)
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
      fHeff[i][j] = conj(fUM[0][i]) * Ve * fUM[0][j];
      if (i > j) fHeff[j][i] = conj(fHeff[i][j]);
    }
  }

  double lv = 2 * kGeV2eV * fEnergy; // 2E in eV

  // add mass terms
  fHeff[1][1] += fDm[1] / lv;
  fHeff[2][2] += fDm[2] / lv;
}

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

// Set Heff in Gell-Mann basis: set only right-triangle as the left part is
// -right
void PMNS_OQS::SetHGM() { get_GMOP(fHeff, fHGM); }

void PMNS_OQS::SetDissipator()
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

  for (int j = 0; j < SU3_DIM; j++) {
    for (int k = j; k < SU3_DIM; k++) {
      if (j == 0 || k == 0) {
        fD[j][k] = 0;
        continue;
      }
      fD[j][k] = -fD[j][k] * kGeV2eV;
      fD[k][j] = fD[j][k];
    }
  }
  
  fBuiltDissipator = true;
}


double PMNS_OQS::GetDissipatorElement(int i, int j){

  return fD[i][j];
  
}


void PMNS_OQS::SetDecoElement(int i, double val)
{
  fBuiltDissipator *= (fa[i] == abs(val));
  fa[i] = abs(val);
}

void PMNS_OQS::SetPower(int n) { fPower = n; }

void PMNS_OQS::SetDecoAngle(int i, int j, double th)
{
  double val = cos(th);
  fBuiltDissipator *= (fcos[i][j] == val);

  if (i == j){
    fcos[i][j] = 1;
  } else {
    fcos[i][j] = val;
    fcos[j][i] = val;
  }
}

void PMNS_OQS::SetM()
{
  SetDissipator();
  double energyCorr = pow(fEnergy, fPower);
  for (int k = 1; k < SU3_DIM; ++k) {
    for (int j = 1; j < SU3_DIM; ++j) {
      fM(k - 1, j - 1) = fHGM[k][j] + fD[k][j] * energyCorr;
    }
  }
}

void PMNS_OQS::SetPhi(int i, double val) { fPhi[i - 1] = val; }

void PMNS_OQS::BuildHms()
{
  if (fBuiltHms) return;
  BuildUM();
  PMNS_Base::BuildHms();
}

void PMNS_OQS::BuildUM()
{
  SetVacuumEigensystem();
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) { fUM[i][j] = conj(fEvec[i][j]); }
  }

  complexD iphi1(0.0, fPhi[0]);
  complexD iphi2(0.0, fPhi[1]);

  fUM[0][1] *= exp(iphi1);
  fUM[0][2] *= exp(iphi2);

  if (fParameterisation == 1) {
    fUM[1][0] *= exp(-iphi1);
    fUM[1][2] *= exp(iphi2 - iphi1);

    fUM[2][0] *= exp(-iphi2);
    fUM[2][1] *= exp(iphi1 - iphi2);
  }
  else {
    fUM[1][1] *= exp(iphi1);
    fUM[1][2] *= exp(iphi2);

    fUM[2][1] *= exp(iphi1);
    fUM[2][2] *= exp(iphi2);
  }
}

//.............................................................................
///
/// Rotate the density matrix to or from the mass basis
///
/// @param to_mass - true if to mass basis
///
void PMNS_OQS::RotateState(bool to_mass)
{
  BuildHms();
  if (!to_mass) ChangeBaseToSU3();
  PMNS_DensityMatrix::RotateState(to_mass, fUM);
  if (to_mass) ChangeBaseToGM();
}

void PMNS_OQS::ChangeBaseToGM() { get_GM(fRho, fR); }

void PMNS_OQS::ChangeBaseToSU3() { get_SU3(fR, fRho); }

void PMNS_OQS::Propagate()
{
  RotateState(true);
  PMNS_Base::Propagate();
  RotateState(false);
}

//.............................................................................
///
/// Propagate the current neutrino state through a given path
///
/// @param p - A neutrino path segment
///
void PMNS_OQS::PropagatePath(NuPath p)
{
  SetHeff(p);

  SetHGM();

  SetM();

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
/// @return Neutrino oscillation probabilities
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
