////////////////////////////////////////////////////////////////////////
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework. 
//
// coelho@lal.in2p3.fr
////////////////////////////////////////////////////////////////////////

#include "PMNS_ScalarNSI.h"

using namespace OscProb;

//......................................................................
///
/// Constructor. \sa PMNS_NSI::PMNS_NSI
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_ScalarNSI::PMNS_ScalarNSI() : PMNS_NSI()
{
  SetLowestMass(0);
}

//......................................................................
///
/// Nothing to clean.
///
PMNS_ScalarNSI::~PMNS_ScalarNSI(){}

//......................................................................
///
/// Set lightest neutrino mass
///
/// @param m - The lightest mass in eV
///
void PMNS_ScalarNSI::SetLowestMass(double m)
{
  fM = m;
}

//......................................................................
///
/// Get lightest neutrino mass
///
/// @return The lightest mass in eV
///
double PMNS_ScalarNSI::GetLowestMass()
{
  return fM;
}

//......................................................................
///
/// Reimplement vacuum Hms update
///
/// \sa PMNS_Base::BuildHms
///
void PMNS_ScalarNSI::BuildHms()
{

  // Check if anything changed
  if(fBuiltHms) return;
  
  // Tag to recompute eigensystem
  fGotES = false;
 
  // Start with m1 = 0
  double m1 = 0;
  // Make sure all masses are positive
  for(int i=1; i<fNumNus; i++){
    if(fDm[i] + m1*m1 < 0) m1 = sqrt(fabs(fDm[i]));
  }
  // Add minimum mass
  m1 += fM;

  // Build M - m1
  fHms[0][0] = 0;
  for(int j=1; j<fNumNus; j++){
    // Set mass difference m_j - m_1
    fHms[j][j] = sqrt(fDm[j] + m1*m1) - m1;
    // Reset off-diagonal elements
    for(int i=0; i<j; i++){
      fHms[i][j] = 0;
    }
    // Rotate j neutrinos
    for(int i=0; i<j; i++){
      RotateH(i,j,fHms);
    }
  }
 
  // Add back m1
  for(int i=0; i<fNumNus; i++){
    fHms[i][i] += m1;
  }

  ClearCache();
 
  // Tag as built
  fBuiltHms = true;

}

//.....................................................................
///
/// Square a Hermitian matrix
///
/// @param A - input matrix A
///
void HermitianSquare(complexD (&A)[3][3])
{

  complexD tmp[3][3];
  for(int i=0; i<3; i++){
  for(int j=i; j<3; j++){
    tmp[i][j] = A[i][j];
    if(i<j) tmp[j][i] = conj(A[i][j]);
    A[i][j] = 0;
  }}

  for(int i=0; i<3; i++){
  for(int j=i; j<3; j++){
  for(int k=0; k<3; k++){

    A[i][j] += tmp[i][k] * tmp[k][j];

  }}}

}

//......................................................................
///
/// Reimplement Hamiltonian update for scalar NSI
///
/// Here we divide the mass squared matrix Hms by the 2E
/// to obtain the vacuum Hamiltonian in eV. Then, the matter
/// potential is added to the electron component.
///
void PMNS_ScalarNSI::UpdateHam()
{

  double lv = 2 * kGeV2eV*fEnergy; // 2*E in eV

  // Effective density of fermions
  double nsiCoup = fPath.density * GetZoACoup();

  // Add the dM matrix
  for(int i=0; i<fNumNus; i++){
  for(int j=i; j<fNumNus; j++){
    fHam[i][j] = (fHms[i][j] + nsiCoup * fEps[i][j]) / lv;
    if(i<j) fHam[j][i] = conj(fHam[i][j]);
  }}

  // Square fHam
  HermitianSquare(fHam); 

  // Add matter potential
  double kr2GNe   = kK2*M_SQRT2*kGf * fPath.density;

  if(!fIsNuBar) fHam[0][0] += kr2GNe;
  else          fHam[0][0] -= kr2GNe;

}


////////////////////////////////////////////////////////////////////////
