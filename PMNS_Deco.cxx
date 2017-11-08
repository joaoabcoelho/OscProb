////////////////////////////////////////////////////////////////////////
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework with decoherence.
//
// This  class inherits from the PMNS_Fast class
//
// jcoelho@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cassert>
#include <stdlib.h>

#include "PMNS_Deco.h"

using namespace OscProb;

using namespace std;


//......................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_Deco::PMNS_Deco() : PMNS_Fast(), fGamma(),
fRho(3, row(3,0))
{
  SetStdPath();
  SetGamma(2,0);
  SetGamma(3,0);
  SetGH(true);
}

//......................................................................
///
/// Nothing to clean.
///
PMNS_Deco::~PMNS_Deco(){}

//......................................................................
///
/// Set any given decoherence parameter.
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// Requires that j > 0. Will notify you if input is wrong.
///
/// @param j   - The second mass index
/// @param val - The absolute value of the parameter
///
void PMNS_Deco::SetGamma(int j, double val){

  if(j < 2 || j > 3){
    cout << "Gamma_" << j << 1 << " not valid for " << fNumNus;
    cout << " neutrinos. Doing nothing." << endl;
    return;
  }

  fGotES *= (fGamma[j-1] == val);
  
  fGamma[j-1] = val;
  
}

//......................................................................
///
/// Set the decoherence hierarchy. This will define the relationship:
///
///   \f$\Gamma_{32} = (\sqrt{\Gamma_{31}} \pm \sqrt{\Gamma_{21}})^{2}\f$
///
/// The + (-) sign will be refered to as IH (NH)
///
/// @param isNH - Is the hierarchy normal?
///
void PMNS_Deco::SetGH(bool isNH)
{

  fGamma[0] = isNH ? 1 : -1;
  
}

//......................................................................
///
/// Get any given decoherence parameter.
///
/// Requires that i > j. Will notify you if input is wrong.
/// If i < j, will assume reverse order and swap i and j.
///
/// @param i  - The first mass index
/// @param j  - The second mass index
///
double PMNS_Deco::GetGamma(int i, int j)
{

  if(i < j){
    cout << "First argument should be larger than second argument" << endl;
    cout << "Setting reverse order (Gamma_" << j << i << "). " << endl;
    int temp = i;
    i = j;
    j = temp;
  }
  if(i<1 || i>3 || i <= j || j < 1){
    cout << "Gamma_" << i << j << " not valid for " << fNumNus;
    cout << " neutrinos. Returning 0." << endl;
    return 0;
  }

  if(j == 1){ 
    return fGamma[i-1];
  }
  else {
    return pow( sqrt(fGamma[2]) - fGamma[0]*sqrt(fGamma[1]), 2);
  }

}

//.....................................................................
/// Dot product of two matrices.
///
/// @param A - input matrix A
/// @param B - input matrix B
///
/// @return - matrix A.B
///
PMNS_Deco::matrix PMNS_Deco::Dot(matrix A, matrix B)
{

  matrix out(3, row(3,0));
  
  for(int i=0; i<3; i++){
  for(int j=0; j<3; j++){
  for(int k=0; k<3; k++){

    out[i][j] += A[i][k] * B[k][j];

  }}}
  
  return out;

}

//.....................................................................
/// Product of elements of two matrices.
///
/// @param A - input matrix A
/// @param B - input matrix B
///
/// @return - matrix A * B
///
PMNS_Deco::matrix PMNS_Deco::Mult(matrix A, matrix B)
{

  matrix out(3, row(3,0));
  
  for(int i=0; i<3; i++){
  for(int j=0; j<3; j++){
  
    out[i][j] += A[i][j] * B[i][j];

  }}
  
  return out;

}

//.....................................................................
/// Conjugate transpose matrix.
///
/// @param A - input matrix A
///
/// @return - \f$A^{\dagger}\f$
///
PMNS_Deco::matrix PMNS_Deco::CTransp(matrix A)
{

  matrix out(3, row(3,0));
  
  for(int i=0; i<3; i++){
  for(int j=0; j<3; j++){
  
    out[i][j] += conj(A[j][i]);

  }}
  
  return out;

}

//.....................................................................
///
/// Propagate the current neutrino state through a given path
///
/// @param p - A neutrino path segment
///
void PMNS_Deco::PropagatePath(NuPath p)
{

  // Set the neutrino path
  SetCurPath(p);

  // Solve for eigensystem
  SolveHam();
  
  // The code below is not pretty or optimised at all. Lots of matrix
  // multiplications are going on which could be otpimised due to 
  // symmetry. The index sorting is also terrible. Needs improving.

  // Store rotation matrices
  matrix U = fEvec;
  matrix Ut = CTransp(fEvec);

  // Rotate to effective mass basis
  fRho = Dot(Ut, Dot(fRho, U));
  
  // Compute evolution matrix
  matrix evolve(3, row(3,1));
  
  // Some ugly way of matching gamma and dmsqr indices
  int maxdm = 0; 
  if(fDm[1] > fDm[maxdm]) maxdm = 1;
  if(fDm[2] > fDm[maxdm]) maxdm = 2;

  int mindm = 0; 
  if(fDm[1] < fDm[mindm]) mindm = 1;
  if(fDm[2] < fDm[mindm]) mindm = 2;

  int middm = 0;
  if(middm == maxdm || middm == mindm) middm = 1; 
  if(middm == maxdm || middm == mindm) middm = 2; 

  int maxev = 0; 
  if(fEval[1] > fEval[maxev]) maxev = 1;
  if(fEval[2] > fEval[maxev]) maxev = 2;

  int minev = 0; 
  if(fEval[1] < fEval[minev]) minev = 1;
  if(fEval[2] < fEval[minev]) minev = 2;

  int midev = 0;
  if(midev == maxev || midev == minev) midev = 1; 
  if(midev == maxev || midev == minev) midev = 2; 
  
  int idx[3];
  idx[minev] = mindm;
  idx[midev] = middm;
  idx[maxev] = maxdm;
  
  for(int j=0; j<3; j++){ 
  for(int i=0; i<j; i++){
    complex arg = complex(GetGamma(max(idx[i],idx[j])+1,min(idx[i],idx[j])+1) * kGeV2eV, fEval[i] - fEval[j]) * kKm2eV * p.length;
    evolve[i][j] = exp(-arg);
    evolve[j][i] = exp(-conj(arg));
  }}
  
  // Evolve density matrix
  fRho = Mult(fRho, evolve);
  
  // Rotate back to flavour basis
  fRho = Dot(U, Dot(fRho, Ut));

}

//......................................................................
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
void PMNS_Deco::ResetToFlavour(int flv)
{
  assert(flv>=0 && flv<fNumNus);
  for (int i=0; i<fNumNus; ++i){
  for (int j=0; j<fNumNus; ++j){
    if (i==flv && i==j) fRho[i][j] = one;
    else                fRho[i][j] = zero;
  }}
}

//......................................................................
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
double PMNS_Deco::P(int flv)
{
  assert(flv>=0 && flv<fNumNus);
  return sqrt(norm(fRho[flv][flv]));
}

////////////////////////////////////////////////////////////////////////
