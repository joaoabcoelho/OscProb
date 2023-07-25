////////////////////////////////////////////////////////////////////////
//
// Implementation of oscillations of neutrinos in matter in a
// n-neutrino framework. 
//
//......................................................................
//
// Throughout I have taken:
//   - L to be the neutrino flight distance in km
//   - E to be the neutrino energy in GeV
//   - Dmij to be the differences between the mass-squares in eV^2
//   - Rho to be the matter density in g/cm^3
//   - theta_ij and delta_ij to be in radians
//
// jcoelho\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

#include "PMNS_Sterile.h"

#include <iostream>
#include <cassert>
#include <stdlib.h>

#include <gsl/gsl_complex_math.h>

using namespace std;

using namespace OscProb;

//......................................................................
///
/// GSL Eigensystem Constructor.
///
/// @param numNus - the number of neutrino flavours
///
GSL_EinSys::GSL_EinSys(int numNus) : 
fNumNus(numNus), 
fEvalGSL(0), fEvecGSL(0), H_GSL(0), W_GSL(0)
{

  // Allocate memory for the GSL objects
  fEvalGSL = gsl_vector_alloc(fNumNus);
  fEvecGSL = gsl_matrix_complex_alloc(fNumNus, fNumNus);
  H_GSL = gsl_matrix_complex_alloc(fNumNus, fNumNus);
  W_GSL = gsl_eigen_hermv_alloc(fNumNus);

}

//......................................................................
///
/// GSL Copy-constructor.
///
/// Allocate new memory for GSL objects when copying.
///
/// @param other - GSL object to copy
///
GSL_EinSys::GSL_EinSys(const GSL_EinSys &other) :
fNumNus(other.fNumNus), 
fEvalGSL(0), fEvecGSL(0), H_GSL(0), W_GSL(0)
{

  // Allocate memory for the GSL objects
  fEvalGSL = gsl_vector_alloc(fNumNus);
  fEvecGSL = gsl_matrix_complex_alloc(fNumNus, fNumNus);
  H_GSL = gsl_matrix_complex_alloc(fNumNus, fNumNus);
  W_GSL = gsl_eigen_hermv_alloc(fNumNus);

  gsl_vector_memcpy(fEvalGSL, other.fEvalGSL);
  gsl_matrix_complex_memcpy(fEvecGSL, other.fEvecGSL);
  gsl_matrix_complex_memcpy(H_GSL, other.H_GSL);

}

//......................................................................
///
/// Copy assignment operator.
///
/// Creates a new deep copy of object..
///
/// @param other - GSL object to assign
///
GSL_EinSys& GSL_EinSys::operator=(const GSL_EinSys &other)
{

  if(fNumNus != other.fNumNus){
 
    // Free memory from GSL objects
    gsl_matrix_complex_free(H_GSL); H_GSL = 0;
    gsl_eigen_hermv_free(W_GSL); W_GSL = 0;
    gsl_matrix_complex_free(fEvecGSL); fEvecGSL = 0;
    gsl_vector_free(fEvalGSL); fEvalGSL = 0;

    fNumNus = other.fNumNus;

    // Allocate memory for the GSL objects
    fEvalGSL = gsl_vector_alloc(fNumNus);
    fEvecGSL = gsl_matrix_complex_alloc(fNumNus, fNumNus);
    H_GSL = gsl_matrix_complex_alloc(fNumNus, fNumNus);
    W_GSL = gsl_eigen_hermv_alloc(fNumNus);
 
  }

  gsl_vector_memcpy(fEvalGSL, other.fEvalGSL);
  gsl_matrix_complex_memcpy(fEvecGSL, other.fEvecGSL);
  gsl_matrix_complex_memcpy(H_GSL, other.H_GSL);
  
  return *this;

}

//......................................................................
///
/// Destructor.
///
/// Clear GSL objects.
///
GSL_EinSys::~GSL_EinSys(){

  // Free memory from GSL objects
  gsl_matrix_complex_free(H_GSL); H_GSL = 0;
  gsl_eigen_hermv_free(W_GSL); W_GSL = 0;
  gsl_matrix_complex_free(fEvecGSL); fEvecGSL = 0;
  gsl_vector_free(fEvalGSL); fEvalGSL = 0;

}

//......................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
/// This class can implement any number of neutrino flavours.
///
/// @param numNus - the number of neutrino flavours
///
PMNS_Sterile::PMNS_Sterile(int numNus) :
PMNS_Base(numNus), fGSL(numNus) {}

//......................................................................
///
/// Build and solve the full Hamiltonian for eigenvectors and eigenvalues.
///
/// This uses GSL to solve the eigensystem.
///
/// Here we divide the mass squared matrix Hms by the 2E
/// to obtain the vacuum Hamiltonian in eV. Then, the matter
/// potential is added to the electron and sterile components.
///
/// The matter potential is implemented as in arXiv:hep-ph/0606054v3 (eq. 3.15).
///
void PMNS_Sterile::SolveHam()
{

  // Build Hamiltonian
  BuildHms();

  // Check if anything changed
  if(fGotES) return;
  
  // Try caching if activated  
  if(TryCache()) return;

  double rho = fPath.density;
  double zoa = fPath.zoa;

  double lv = 2 * kGeV2eV*fEnergy;                 // 2E in eV 
  double kr2GNe = kK2*M_SQRT2*kGf * rho*zoa;       // Electron matter potential in eV
  double kr2GNn = kK2*M_SQRT2*kGf * rho*(1-zoa)/2; // Neutron matter potential in eV

  // Finish building Hamiltonian in matter with dimension of eV
  for(size_t i=0;i<size_t(fNumNus);i++){
    complexD buf = fHms[i][i]/lv;
    *gsl_matrix_complex_ptr(fGSL.H_GSL, i, i) = gsl_complex_rect(real(buf),0);
    for(size_t j=i+1;j<size_t(fNumNus);j++){
      buf = fHms[i][j]/lv;
      if(!fIsNuBar) *gsl_matrix_complex_ptr(fGSL.H_GSL, j, i) = gsl_complex_rect(real(buf),-imag(buf));
      else          *gsl_matrix_complex_ptr(fGSL.H_GSL, j, i) = gsl_complex_rect(real(buf), imag(buf));
    }
    if(i>2){
      // Subtract NC coherent forward scattering from sterile neutrinos. See arXiv:hep-ph/0606054v3, eq. 3.15, for example. 
      if(!fIsNuBar) *gsl_matrix_complex_ptr(fGSL.H_GSL, i, i) = gsl_complex_add_real(gsl_matrix_complex_get(fGSL.H_GSL,i,i) ,  kr2GNn);
      else          *gsl_matrix_complex_ptr(fGSL.H_GSL, i, i) = gsl_complex_add_real(gsl_matrix_complex_get(fGSL.H_GSL,i,i) , -kr2GNn);;
    }
  }
  // Add nue CC coherent forward scattering from sterile neutrinos. 
  if(!fIsNuBar) *gsl_matrix_complex_ptr(fGSL.H_GSL, 0, 0) = gsl_complex_add_real(gsl_matrix_complex_get(fGSL.H_GSL,0,0) ,  kr2GNe);
  else          *gsl_matrix_complex_ptr(fGSL.H_GSL, 0, 0) = gsl_complex_add_real(gsl_matrix_complex_get(fGSL.H_GSL,0,0) , -kr2GNe);

  // Solve Hamiltonian for eigensystem
  gsl_eigen_hermv(fGSL.H_GSL, fGSL.fEvalGSL, fGSL.fEvecGSL, fGSL.W_GSL);

  // Fill fEval and fEvec vectors from GSL objects
  for(int i=0;i<fNumNus;i++){
    fEval[i] = gsl_vector_get(fGSL.fEvalGSL,i);
    for(int j=0;j<fNumNus;j++){
      gsl_complex buf = gsl_matrix_complex_get(fGSL.fEvecGSL,i,j);
      fEvec[i][j] = complexD( GSL_REAL(buf), GSL_IMAG(buf) );
    }
  }
  
  fGotES = true;
  
  // Fill cache if activated  
  FillCache();
  
}

////////////////////////////////////////////////////////////////////////
