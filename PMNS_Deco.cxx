////////////////////////////////////////////////////////////////////////
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework with decoherence.
//
// This  class inherits from the PMNS_Fast class
//
// coelho@lal.in2p3.fr
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
  SetDecoAngle(0);
  SetPower(0);
}

//......................................................................
///
/// Nothing to clean.
///
PMNS_Deco::~PMNS_Deco(){}

//......................................................................
///
/// Set the decoherence parameter \f$\Gamma_{j1}\f$.
///
/// Requires that j = 2 or 3. Will notify you if input is wrong.
///
/// @param j   - The first mass index
/// @param val - The absolute value of the parameter
///
void PMNS_Deco::SetGamma(int j, double val){

  if(j < 2 || j > 3){
    cout << "Gamma_" << j << 1 << " not valid for " << fNumNus;
    cout << " neutrinos. Doing nothing." << endl;
    return;
  }

  fGamma[j-1] = val;
  
}

//......................................................................
///
/// Set the decoherence parameter \f$\Gamma_{32}\f$.
///
/// In practice, this sets \f$\Gamma_{31}\f$ internally via the formula:
///
///   \f$
///   \Gamma_{31} = \Gamma_{32} + \Gamma_{21} \cos^2\theta +
///   \cos\theta \sqrt{\Gamma_{21} (4\Gamma_{32} - \Gamma_{21} (1 - \cos^2\theta))}
///   \f$
///
/// IMPORTANT: Note this needs to be used AFTER defining \f$\Gamma_{21}\f$ and \f$\theta\f$.
///
/// @param val - The absolute value of the parameter
///
void PMNS_Deco::SetGamma32(double val){

  double min32 = 0.25*fGamma[1];

  if(fGamma[0] >= 0) min32 *= 1 - pow(fGamma[0],2);
  else               min32 *= 1 + 3*pow(fGamma[0],2);
  
  if(val < min32){
    
    if(fGamma[0]>=0 || 4*val/fGamma[1] < 1){
      fGamma[0] = sqrt(1 - 4*val/fGamma[1]);
      fGamma[2] = 0.25*fGamma[1]*(1 + 3*pow(fGamma[0],2));
    }
    else{
      fGamma[0] = -sqrt((4*val/fGamma[1] - 1)/3);
      fGamma[2] = 0.25*fGamma[1]*(1 - pow(fGamma[0],2));
    }
    
    cout << "WARNING: Impossible to have Gamma32 = "        << val
         << " with current Gamma21 and theta parameters."   << endl
         << "         Changing the value of cos(theta) to " << fGamma[0] 
         << endl;
         
    return;
         
  }

  double arg = fGamma[1] * ( 4*val - fGamma[1]*(1 - pow(fGamma[0],2)) );

  if(arg < 0){
    arg = 0;
    fGamma[0] = sqrt(1 - 4*val/fGamma[1]);
    cout << "WARNING: Imaginary Gamma31 found. Changing the value of cos(theta) to " 
         << fGamma[0] << endl;
  }
  
  double gamma31 = val + fGamma[1]*pow(fGamma[0],2) + fGamma[0] * sqrt(arg);
  
  fGamma[2] = gamma31;
  
  if(fabs(val - GetGamma(3,2)) > 1e-6){
    cout << "ERROR: Failed sanity check: GetGamma(3,2) = "
         << GetGamma(3,2) << " != " << val << endl;
  }
  
}

//......................................................................
///
/// Set the decoherence angle. This will define the relationship:
///
///   \f$ 
///   \Gamma_{32} = \Gamma_{31} + \Gamma_{21} \cos^2\theta - 
///   \cos\theta \sqrt{\Gamma_{21} (4\Gamma_{31} - \Gamma_{21} (1 - \cos^2\theta))} 
///   \f$
///
/// @param th - decoherence angle
///
void PMNS_Deco::SetDecoAngle(double th)
{

  fGamma[0] = cos(th);
  
}

//......................................................................
///
/// Set the power index of the decoherence energy dependence. 
/// This will multiply the Gammas such that the final parameters are:
///
///   \f$ 
///   \Gamma_{ij} \rightarrow \Gamma_{ij} \times \left(\frac{E}{[1 \mbox{GeV}]}\right)^n
///   \f$
///
/// @param n - power index
///
void PMNS_Deco::SetPower(double n)
{

  fPower = n;
  
}

//......................................................................
///
/// Get the decoherence angle value.
///
double PMNS_Deco::GetDecoAngle()
{

  return acos(fGamma[0]);
  
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
    // Combine Gamma31, Gamma21 and an angle (fGamma[0] = cos(th)) to make Gamma32
    double arg = fGamma[1] * ( 4*fGamma[2] - fGamma[1]*(1 - pow(fGamma[0],2)) );
    
    if(arg < 0) return fGamma[1] - 3*fGamma[2];

    return fGamma[2] + fGamma[1]*pow(fGamma[0],2) - fGamma[0] * sqrt(arg);

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
  
  double energyCorr = pow(fEnergy, fPower);

  for(int j=0; j<3; j++){ 
  for(int i=0; i<j; i++){
    double gamma_ij = GetGamma(max(idx[i],idx[j])+1, min(idx[i],idx[j])+1) * energyCorr;
    complex arg = complex(gamma_ij * kGeV2eV, fEval[i] - fEval[j]) * kKm2eV * p.length;
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
  return std::abs(fRho[flv][flv]);
}

////////////////////////////////////////////////////////////////////////
