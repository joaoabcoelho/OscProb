///////////////////////////////////////////////////////////////////////////////
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework with gravitationally induced decoherence.
//
// This  class inherits from the PMNS_Fast class
//
// This developement is part of the QGRANT project with                                   
// ID: 101068013,                                                                           
// founded by the HORIZON-MSCA-2021-PF-01-01 programme. 
//
// \author Joao Coelho - jcoelho\@apc.in2p3.fr                                           
// \author Alba Domi - alba.domi\@fau.de      
///////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <cmath>
#include <complex>
#include <limits>
#include <iomanip>

#include <Eigen/Eigenvalues>

#include "PMNS_GQD.h"

using namespace OscProb;

using namespace std;

double hbar = 6.582119569e-16; // eV*s

//.............................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_GQD::PMNS_GQD()
  : PMNS_Fast(), fRho(3, vectorC(3, 0)), fMBuffer(3, vectorC(3, 0))
{
  SetStdPath();
  SetEta(0);
  SetTemperature(0);
  SetOmega(0);
 }

//.............................................................................
///
/// Nothing to clean.
///
PMNS_GQD::~PMNS_GQD() {}

//.............................................................................
///
/// Set the coupling parameter \f$\Eta\f$ in [s].
///
/// @param val - The absolute value of the parameter
///
void PMNS_GQD::SetEta(double val)
{
  if (val < 0) {
    cerr << "WARNING: Eta must be positive."
         << "Setting it to absolute value of input: " << -val << endl;
    val = -val;
  }

  fEta = val;
  
  return;
}

//.............................................................................
///
/// Set the parameter \f$\Temperature\f$ in [K].
///
/// @param val - The absolute value of the parameter
///
void PMNS_GQD::SetTemperature(double val)
{
  if (val < 0) {
    cerr << "WARNING: Temperature must be positive. "
         << "Setting it to absolute value of input: " << -val << endl;
    val = abs(val);
  }

  fT = val;
  
  return;
}

//.............................................................................
///
/// Set the cutoff frequency \f$\Omega\f$ in [Hz].
///
/// @param val - The absolute value of the parameter
///
void PMNS_GQD::SetOmega(double val)
{
  if (val < 0) {
    cerr << "WARNING: Omega must be positive. "
         << "Setting it to absolute value of input: " << -val << endl;
    val = -val;
  }

  fOmega = val;
  
  return;
}


void PMNS_GQD::SetVacuumEigensystem()
{  
  double   s12, s23, s13, c12, c23, c13;
  complexD idelta(0.0, fDelta[0][2]);
  if (fIsNuBar) idelta = conj(idelta);

  s12 = sin(fTheta[0][1]);
  s23 = sin(fTheta[1][2]);
  s13 = sin(fTheta[0][2]);
  c12 = cos(fTheta[0][1]);
  c23 = cos(fTheta[1][2]);
  c13 = cos(fTheta[0][2]);

  fEvec[0][0] =  c12 * c13;
  fEvec[0][1] =  s12 * c13;
  fEvec[0][2] =  s13 * exp(-idelta);

  fEvec[1][0] = -s12 * c23 - c12 * s23 * s13 * exp(idelta);
  fEvec[1][1] =  c12 * c23 - s12 * s23 * s13 * exp(idelta);
  fEvec[1][2] =  s23 * c13;

  fEvec[2][0] =  s12 * s23 - c12 * c23 * s13 * exp(idelta);
  fEvec[2][1] = -c12 * s23 - s12 * c23 * s13 * exp(idelta);
  fEvec[2][2] =  c23 * c13;
  
  double lv = 6 * kGeV2eV * fEnergy; // 6E in eV     
  
  fEval[0] = -(fDm[1] + fDm[2])/lv;
  fEval[1] =  (fDm[1] - (fDm[2] - fDm[1]))/lv;
  fEval[2] =  (fDm[2] + (fDm[2] - fDm[1]))/lv;  
}

//.............................................................................           
///                                                                                          
/// Build the neutrino Hamiltonian in vacuum:
/// here, the diagonal elements are the neutrino masses
/// instead of the typically used neutrino mass squared differences.
///                                                                                         
void PMNS_GQD::BuildHms()
{
  // Check if anything changed                                                           
  if (fBuiltHms) return;

  // Tag to recompute eigensystem                                                        
  fGotES = false;

  // Set neutrino masses from mass splittings                     
  fHms[0][0] = -(fDm[1] + fDm[2])/3 + fDm[1];
  fHms[1][1] =  (fDm[1] - (fDm[2] - fDm[1]))/3 + fDm[1];
  fHms[2][2] =  (fDm[2] + (fDm[2] - fDm[1]))/3 + fDm[1];

  for (int j = 0; j < fNumNus; j++) {
    // Reset off-diagonal elements                                                       
    for (int i = 0; i < j; i++) { fHms[i][j] = 0; }
  }

  RotateHam(false, fHms);
    
  ClearCache();

  // Tag as built                                                                        
  fBuiltHms = true;
}

//.............................................................................           
///                                                                                          
/// Rotate the Hamiltonian to or from the mass basis                                    
///                                                                                         
/// @param to_mass - true if to mass basis                                                     
/// @param Ham - Hamiltonian to be rotated
///   
void PMNS_GQD::RotateHam(bool to_mass, matrixC& Ham)
{
  
  matrixC U(3, vectorC(3, 0));
  
  double s12 = sin(fTheta[0][1]);
  double s13 = sin(fTheta[0][2]);
  double s23 = sin(fTheta[1][2]);

  double c12 = cos(fTheta[0][1]);
  double c13 = cos(fTheta[0][2]);
  double c23 = cos(fTheta[1][2]);

  complexD idelta(0.0, fDelta[0][2]);
  if (fIsNuBar){
    idelta = conj(idelta);
  }

  U[0][0] =  c12 * c13;
  U[0][1] =  s12 * c13;
  U[0][2] =  s13 * exp(-idelta);

  U[1][0] = -s12 * c23 - c12 * s23 * s13 * exp(idelta);
  U[1][1] =  c12 * c23 - s12 * s23 * s13 * exp(idelta);
  U[1][2] =  s23 * c13;

  U[2][0] =  s12 * s23 - c12 * c23 * s13 * exp(idelta);
  U[2][1] = -c12 * s23 - s12 * c23 * s13 * exp(idelta);
  U[2][2] =  c23 * c13;

  matrixC mult(3, vectorC(3, 0));

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      mult[i][j] = 0;
      for (int k = 0; k < 3; k++) {
        if (to_mass)
          mult[i][j] += Ham[i][k] * U[k][j];
        else
          mult[i][j] += Ham[i][k] * conj(U[j][k]);
      }
    }
   }

  for (int i = 0; i < 3; i++) {
    for (int j = i; j < 3; j++) {
      Ham[i][j] = 0;
      for (int k = 0; k < 3; k++) {
        if (to_mass)
          Ham[i][j] += conj(U[k][i]) * mult[k][j];
        else
          Ham[i][j] += U[i][k] * mult[k][j];
      }
      if (j > i) Ham[j][i] = conj(Ham[i][j]);
    }
  }
}

//.............................................................................          
///                                                                                        
/// Solve the full Hamiltonian for eigenvectors and eigenvalues.                             
///                                                                                           
/// This is using a method from the GLoBES software available at                               
/// http://www.mpi-hd.mpg.de/personalhomes/globes/3x3/ \n                                    
/// We should cite them accordingly                                                         
///   
void PMNS_GQD::SolveHam()
{
  // Do vacuum oscillation in low density                                                
  if (fPath.density < 1.0e-6) {
    SetVacuumEigensystem();
    return;
  }

  // Build Hamiltonian                                             
  BuildHms();
  
  // Check if anything changed                                                       
  if (fGotES) return;

  // Try caching if activated                                                        
  if (TryCache()) return;

  UpdateHam();
  
  double   fEvalGLoBES[3];
  complexD fEvecGLoBES[3][3];

  // Solve Hamiltonian for eigensystem using the GLoBES method
  zheevh3(fHam, fEvecGLoBES, fEvalGLoBES);

  // Fill fEval and fEvec vectors from GLoBES arrays
  for (int i = 0; i < fNumNus; i++) {
    fEval[i] = fEvalGLoBES[i];
    for (int j = 0; j < fNumNus; j++) {
      fEvec[i][j] = fEvecGLoBES[i][j];
    }
  }
  
  // Mark eigensystem as solved      
  fGotES = true;

  // Fill cache if activated                                              
  FillCache();
}

//.............................................................................
///
/// Rotate the density matrix to or from the mass basis
///
/// @param to_mass - true if to mass basis
///
void PMNS_GQD::RotateState(bool to_mass)
{

  // buffer = rho . U
  for (int i = 0; i < fNumNus; i++) {
    for (int j = 0; j < fNumNus; j++) {
      fMBuffer[i][j] = 0;
      for (int k = 0; k < fNumNus; k++) {
        if (to_mass){
	  fMBuffer[i][j] += fRho[i][k] * fEvec[k][j];
	}else
          fMBuffer[i][j] += fRho[i][k] * conj(fEvec[j][k]);
      }
    }
  }

  // rho = U^\dagger . buffer = U^\dagger . rho . U
  // Final matrix is Hermitian, so copy upper to lower triangle
  for (int i = 0; i < fNumNus; i++) {
    for (int j = i; j < fNumNus; j++) {
      fRho[i][j] = 0;
      for (int k = 0; k < fNumNus; k++) {
        if (to_mass)
          fRho[i][j] += conj(fEvec[k][i]) * fMBuffer[k][j];
        else
          fRho[i][j] += fEvec[i][k] * fMBuffer[k][j];
      }
      if (j > i) fRho[j][i] = conj(fRho[i][j]);
    }
  }
}

//.............................................................................
///
/// Propagate the current neutrino state through a given path
///
/// @param p - A neutrino path segment
///
void PMNS_GQD::PropagatePath(NuPath p)
{
  // Set the neutrino path
  SetCurPath(p);

  // Solve for eigensystem
  SolveHam();

  // Rotate to effective mass basis
  RotateState(true);

  // Convert path length to eV
  double lengthIneV = kKm2eV * p.length;
  double kB = 8.6173*pow(10, -5); // eV*K^-1
  double beta = 1. / (kB*fT);
  
  // Apply phase rotation to off-diagonal elements
  for (int j = 0; j < fNumNus; j++) {
    for (int i = 0; i < j; i++) {
      
      // This is the standard oscillation phase
      double arg = (fEval[i] - fEval[j]) * lengthIneV;

      double arg2 = (2./hbar) * fEta*fEta * fOmega * (fEval[i]*fEval[i]
						      - fEval[j]*fEval[j]
						      + 2*fEval[i]*kGeV2eV*fEnergy
						      - 2*fEval[j]*kGeV2eV*fEnergy) * lengthIneV; 
      
      fRho[i][j] *= exp(-4. * (1./(hbar*hbar)) * (fEta*fEta/beta) * (fEval[i]-fEval[j])*(fEval[i]-fEval[j]) * lengthIneV) * complexD(cos(arg), -sin(arg)) * complexD(cos(arg2), sin(arg2));
      
      fRho[j][i] = conj(fRho[i][j]);
    }
  }

  // Rotate back to flavour basis
  RotateState(false);
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
void PMNS_GQD::ResetToFlavour(int flv)
{
  PMNS_Base::ResetToFlavour(flv);

  assert(flv >= 0 && flv < fNumNus);
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
double PMNS_GQD::P(int flv)
{
  assert(flv >= 0 && flv < fNumNus);
  return abs(fRho[flv][flv]);
}

//.............................................................................
///
/// Set the density matrix from a pure state
///
/// @param nu_in - The neutrino initial state in flavour basis.
///
void PMNS_GQD::SetPureState(vectorC nu_in)
{
  assert(nu_in.size() == fNumNus);
  
  for (int i = 0; i < fNumNus; i++) {
    for (int j = 0; j < fNumNus; j++) {
      fRho[i][j] = conj(nu_in[i]) * nu_in[j];
    }
  }
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
matrixD PMNS_GQD::ProbMatrix(int nflvi, int nflvf)
{
  assert(nflvi <= fNumNus && nflvi >= 0);
  assert(nflvf <= fNumNus && nflvf >= 0);

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


//.............................................................................          
///                                                                                      
/// Solve the Homiltonian using the appropriate Eigen::MatrixType.\n
/// This will implement faster methods for N<5
///                                                                                      
template <typename T> void PMNS_GQD::SolveEigenSystem()
{
  
  for(int j = 1; j < 3; ++j){
    for(int i = 0; i < j; ++i){
      fHmatter(j, i) = conj(fHmatter(i, j));
    }
  }
  
  Eigen::Ref<T> H(fHmatter);

  Eigen::SelfAdjointEigenSolver<T> eigensolver(H);
  
  // Fill fEval and fEvec vectors from GSL objects                            
  for (int i = 0; i < 3; i++) {
    fEval[i] = eigensolver.eigenvalues()(i);
    cout << "EVALS \n";
    cout << fixed
    	 << setprecision( 
              numeric_limits<double>::max_digits10) 
	 << fEval[i] << " "; 
    for (int j = 0; j < 3; j++) {
      fEvec[i][j] = eigensolver.eigenvectors()(i, j);
      cout <<"evec: i " << i << " j " << j << " " << fEvec[i][j] << endl;
    }
  }
  cout << endl;

}
