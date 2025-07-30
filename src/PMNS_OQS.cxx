#include <cassert>
#include <complex>
#include <iostream>
#include <math.h>

#include <unsupported/Eigen/MatrixFunctions>

#include "PMNS_OQS.h"

using namespace std;
using namespace OscProb;

//.............................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_OQS::PMNS_OQS()
    : PMNS_DensityMatrix(), fPhi(), fR(), fRt(), fa(9, 0), fMe(8, 8),
      fD(9, vectorD(9, 0)), fM(9, vectorD(9, 0)), fcos(9, vectorD(9, 1)),
      fHGM(9, vectorD(9, 0)), fHeff(3, vectorC(3, 0)), fUM(3, vectorC(3, 0))
{
  InitializeVectors();
  SetParameterisation(1);
}

//.............................................................................
///
/// Nothing to clean.
///
PMNS_OQS::~PMNS_OQS() {}

void PMNS_OQS::InitializeVectors()
{
  SetPhi(1, 0);
  SetPhi(2, 0);
}

void PMNS_OQS::SetIsNuBar(bool isNuBar)
{
  fBuiltHms *= (fIsNuBar == isNuBar);
  PMNS_Base::SetIsNuBar(isNuBar);
}

void PMNS_OQS::SetParameterisation(int param = 1) { fParameterisation = param; }

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

  double lv = 2. * kGeV2eV * fEnergy; // 2E in eV

  // add mass terms
  fHeff[1][1] += fDm[1] / lv;
  fHeff[2][2] += fDm[2] / lv;
}

// Set Heff in Gell-Mann basis: set only right-triangle as the left part is
// -right
void PMNS_OQS::SetHGM()
{
  fHGM[1][2] = real(fHeff[0][0] - fHeff[1][1]);
  fHGM[1][3] = 2. * imag(fHeff[0][1]);
  fHGM[1][4] = -imag(fHeff[1][2]);
  fHGM[1][5] = -real(fHeff[1][2]);
  fHGM[1][6] = -imag(fHeff[0][2]);
  fHGM[1][7] = -real(fHeff[0][2]);

  fHGM[2][3] = 2. * real(fHeff[0][1]);
  fHGM[2][4] = real(fHeff[1][2]);
  fHGM[2][5] = -imag(fHeff[1][2]);
  fHGM[2][6] = -real(fHeff[0][2]);
  fHGM[2][7] = imag(fHeff[0][2]);

  fHGM[3][4] = -imag(fHeff[0][2]);
  fHGM[3][5] = -real(fHeff[0][2]);
  fHGM[3][6] = imag(fHeff[1][2]);
  fHGM[3][7] = real(fHeff[1][2]);

  fHGM[4][5] = real(fHeff[0][0] - fHeff[2][2]);
  fHGM[4][6] = -imag(fHeff[0][1]);
  fHGM[4][7] = real(fHeff[0][1]);
  fHGM[4][8] = sqrt(3.) * imag(fHeff[0][2]);

  fHGM[5][6] = -real(fHeff[0][1]);
  fHGM[5][7] = -imag(fHeff[0][1]);
  fHGM[5][8] = sqrt(3.) * real(fHeff[0][2]);

  fHGM[6][7] = real(fHeff[1][1] - fHeff[2][2]);
  fHGM[6][8] = sqrt(3.) * imag(fHeff[1][2]);

  fHGM[7][8] = sqrt(3.) * real(fHeff[1][2]);

  for (int i = 1; i < 9; ++i) {
    for (int j = i + 1; j < 9; ++j) { fHGM[j][i] = -fHGM[i][j]; }
  }
}

bool is_same_array(array<int, 3> input, array<int, 3> test)
{
  sort(input.begin(), input.end());
  sort(test.begin(), test.end());
  return input == test;
}

int get_permutation(array<int, 3> input, array<int, 3> test)
{
  if (!is_same_array(input, test)) return 0;

  if (input[0] == test[0] && input[1] == test[1]) return 1;
  if (input[0] == test[1] && input[1] == test[2]) return 1;
  if (input[0] == test[2] && input[1] == test[0]) return 1;

  return -1;
}

double get_f(array<int, 3> indices)
{
  int perm = 0;

  if ((perm = get_permutation(indices, {1, 2, 3}))) return perm;

  double value = 0.5;

  if ((perm = get_permutation(indices, {1, 4, 7}))) return perm * value;
  if ((perm = get_permutation(indices, {1, 6, 5}))) return perm * value;
  if ((perm = get_permutation(indices, {2, 4, 6}))) return perm * value;
  if ((perm = get_permutation(indices, {2, 5, 7}))) return perm * value;
  if ((perm = get_permutation(indices, {3, 4, 5}))) return perm * value;
  if ((perm = get_permutation(indices, {3, 7, 6}))) return perm * value;

  value = sqrt(3) / 2;

  if ((perm = get_permutation(indices, {4, 5, 8}))) return perm * value;
  if ((perm = get_permutation(indices, {6, 7, 8}))) return perm * value;

  return 0;
}

void PMNS_OQS::SetDissipatorElement(int j, int k)
{
  double val = 0;

  for (int l = 0; l < 9; l++) {
    for (int m = l; m < 9; m++) {
      double f = 0;
      for (int n = 0; n < 9; n++) {
        f += get_f({l, k, n}) * get_f({n, m, j});
        if (m > l) f += get_f({m, k, n}) * get_f({n, l, j});
      }
      val += fa[l] * fa[m] * fcos[l][m] * f;
    }
  }

  val *= kGeV2eV; // convert to eV

  fD[j][k] = -val;
  fD[k][j] = -val;
}

void PMNS_OQS::SetDissipator()
{
  if (fBuiltDissipator) return;

  for (int j = 0; j < 9; j++) {
    for (int k = j; k < 9; k++) { SetDissipatorElement(j, k); }
  }

  fBuiltDissipator = true;
}

void PMNS_OQS::Seta(int i, double val)
{
  fBuiltDissipator *= (fa[i] == val);
  fa[i] = abs(val);
}

void PMNS_OQS::Setcos(int i, int j, double val)
{
  fBuiltDissipator *= (fcos[i][j] == val);
  if (i == j) fcos[i][j] = 1;
  if (val > 1) val = 1;
  if (val < -1) val = -1;
  fcos[i][j] = val;
  fcos[j][i] = val;
}

void PMNS_OQS::SetM()
{
  SetDissipator();
  for (int k = 0; k < 9; ++k) {
    for (int j = 0; j < 9; ++j) { fM[k][j] = fHGM[k][j] + fD[k][j]; }
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
  PMNS_DensityMatrix::RotateState(to_mass, fUM);
}

void PMNS_OQS::ChangeBaseToGM()
{
  fR[0] = real(fRho[0][0] + fRho[1][1] + fRho[2][2]) / sqrt(6.);
  fR[1] = real(fRho[0][1]);
  fR[2] = -imag(fRho[0][1]);
  fR[3] = real(fRho[0][0] - fRho[1][1]) / 2.;
  fR[4] = real(fRho[0][2]);
  fR[5] = -imag(fRho[0][2]);
  // fR[5] =  0;
  fR[6] = real(fRho[1][2]);
  fR[7] = -imag(fRho[1][2]);
  //  fR[7] =  0;
  fR[8] = real(fRho[0][0] + fRho[1][1] - 2. * fRho[2][2]) / (2. * sqrt(3.));
}

void PMNS_OQS::ChangeBaseToSU3()
{
  fRho[0][0] = (sqrt(2.) * fRt[0] + sqrt(3.) * fRt[3] + fRt[8]) / sqrt(3.);
  fRho[0][1] = (fRt[1] - complexD(0.0, 1.0) * fRt[2]);
  fRho[0][2] = (fRt[4] - complexD(0.0, 1.0) * fRt[5]);
  fRho[1][0] = (fRt[1] + complexD(0.0, 1.0) * fRt[2]);
  fRho[1][1] = (sqrt(2. / 3.) * fRt[0] - fRt[3] + fRt[8] / sqrt(3.));
  fRho[1][2] = (fRt[6] - complexD(0.0, 1.0) * fRt[7]);
  fRho[2][0] = (fRt[4] + complexD(0.0, 1.0) * fRt[5]);
  fRho[2][1] = (fRt[6] + complexD(0.0, 1.0) * fRt[7]);
  fRho[2][2] = (1. / sqrt(3.) * (sqrt(2.) * fRt[0] - 2. * fRt[8]));
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
  for (int i = 1; i < 9; i++) {
    for (int j = 1; j < 9; j++) { fMe(i - 1, j - 1) = fM[i][j]; }
  }

  // Convert path length to eV
  double lengthIneV = kKm2eV * p.length;

  fMe *= lengthIneV;
  fMe = fMe.exp();

  RotateState(true);

  ChangeBaseToGM();

  fRt[0] = fR[0];

  for (int i = 1; i < 9; ++i) {
    fRt[i] = 0;
    for (int j = 1; j < 9; ++j) { fRt[i] += fMe(i - 1, j - 1) * fR[j]; }
  }

  ChangeBaseToSU3();

  RotateState(false); // go back to flavour basis
}

////////////////////////////////////////////////////////////////////////
