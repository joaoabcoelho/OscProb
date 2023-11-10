///////////////////////////////////////////////////////////////////////////////
//
// Implementation of oscillations of neutrinos with scalar NSI
// three-neutrino framework.
//
// jcoelho@apc.in2p3.fr and urahaman@km3net.de
///////////////////////////////////////////////////////////////////////////////

#include "PMNS_SNSI.h"

using namespace OscProb;

//.............................................................................
///
/// Constructor. \sa PMNS_NSI::PMNS_NSI
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_SNSI::PMNS_SNSI() : PMNS_NSI() { SetLowestMass(0); }

//.............................................................................
///
/// Nothing to clean.
///
PMNS_SNSI::~PMNS_SNSI() {}

//.............................................................................
///
/// Set lightest neutrino mass
///
/// @param m - The lightest mass in eV
///
void PMNS_SNSI::SetLowestMass(double m)
{
  // Check if value is actually changing
  fBuiltHms *= (fM == m);

  fM = m;
}

//.............................................................................
///
/// Get lightest neutrino mass
///
/// @return The lightest mass in eV
///
double PMNS_SNSI::GetLowestMass() { return fM; }

//.............................................................................
///
/// Reimplement vacuum Hms update
///
/// \sa PMNS_Base::BuildHms
///
void PMNS_SNSI::BuildHms()
{
  // Check if anything changed
  if (fBuiltHms) return;

  // Tag to recompute eigensystem
  fGotES = false;

  // Start with m1 = fM
  double m1 = fM;
  // Make sure all masses are larger than fM
  for (int i = 1; i < fNumNus; i++) {
    if (fDm[i] + m1 * m1 < fM * fM) m1 = sqrt(fabs(fM * fM - fDm[i]));
  }

  // Build M - m1
  fHms[0][0] = 0;
  for (int j = 1; j < fNumNus; j++) {
    // Set mass difference m_j - m_1
    fHms[j][j] = sqrt(fabs(fDm[j] + m1 * m1)) - m1;
    // Reset off-diagonal elements
    for (int i = 0; i < j; i++) { fHms[i][j] = 0; }
    // Rotate j neutrinos
    for (int i = 0; i < j; i++) { RotateH(i, j, fHms); }
  }

  // Add back m1
  for (int i = 0; i < fNumNus; i++) { fHms[i][i] += m1; }

  // Tag as built
  fBuiltHms = true;
}

//.............................................................................
///
/// Square a Hermitian matrix
///
/// @param A - input matrix A
///
void HermitianSquare(complexD (&A)[3][3])
{
  complexD tmp[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = i; j < 3; j++) {
      tmp[i][j] = A[i][j];
      if (i < j) tmp[j][i] = conj(A[i][j]);
      A[i][j] = 0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = i; j < 3; j++) {
      for (int k = 0; k < 3; k++) { A[i][j] += tmp[i][k] * tmp[k][j]; }
    }
  }
}

//.............................................................................
///
/// Reimplement Hamiltonian update for scalar NSI
///
/// Here we divide the mass squared matrix Hms by the 2E
/// to obtain the vacuum Hamiltonian in eV. Then, the matter
/// potential is added to the electron component.
///
void PMNS_SNSI::UpdateHam()
{
  double sqrtlv = sqrt(2 * kGeV2eV * fEnergy); // sqrt(2*E) in eV^1/2

  // Effective density of fermions in eV * MeV^2
  double nsiCoup = 1e6 * kK2 * fPath.density * GetZoACoup();

  // Add the dM matrix
  for (int i = 0; i < fNumNus; i++) {
    for (int j = i; j < fNumNus; j++) {
      fHam[i][j] = (fHms[i][j] + nsiCoup * fEps[i][j]) / sqrtlv;
      if (fIsNuBar) fHam[i][j] = conj(fHam[i][j]);
    }
  }

  // Square fHam
  HermitianSquare(fHam);

  // Add matter potential in eV
  double kr2GNe = kK2 * M_SQRT2 * kGf * fPath.density * fPath.zoa;

  if (!fIsNuBar)
    fHam[0][0] += kr2GNe;
  else
    fHam[0][0] -= kr2GNe;
}

///////////////////////////////////////////////////////////////////////////////
