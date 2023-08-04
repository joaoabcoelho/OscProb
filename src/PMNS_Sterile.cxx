///////////////////////////////////////////////////////////////////////////////
///
/// Implementation of oscillations of neutrinos in matter in a
/// n-neutrino framework.
///
///............................................................................
///
/// Throughout I have taken:
///   - L to be the neutrino flight distance in km
///   - E to be the neutrino energy in GeV
///   - Dmij to be the differences between the mass-squares in eV^2
///   - Rho to be the matter density in g/cm^3
///   - theta_ij and delta_ij to be in radians
///
/// jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#include <Eigen/Eigenvalues>

#include "PMNS_Sterile.h"

using namespace std;

using namespace OscProb;

//.............................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
/// This class can implement any number of neutrino flavours.
///
/// @param numNus - the number of neutrino flavours
///
PMNS_Sterile::PMNS_Sterile(int numNus) :
PMNS_Base(numNus), fHam(numNus, numNus) {}


//.............................................................................
///
/// Build the full Hamiltonian in matter.
///
/// Here we divide the mass squared matrix Hms by the 2E
/// to obtain the vacuum Hamiltonian in eV. Then, the matter
/// potential is added to the electron and NC components.
///
void PMNS_Sterile::UpdateHam(){

  double rho = fPath.density;
  double zoa = fPath.zoa;

  double lv = 2 * kGeV2eV*fEnergy;                 // 2E in eV
  double kr2GNe = kK2*M_SQRT2*kGf * rho*zoa;       // Electron matter potential in eV
  double kr2GNn = kK2*M_SQRT2*kGf * rho*(1-zoa)/2; // Neutron matter potential in eV

  // Finish building Hamiltonian in matter with dimension of eV
  for(int i=0; i<fNumNus; i++){

    // Set diagonal elements which are real
    fHam(i,i) = fHms[i][i]/lv;

    // Set off-diagonal elements and flip deltaCP if needed
    for(int j=i+1; j<fNumNus; j++){
      // Obs: fHam is lower tirangular while fHms is upper triangular
      if(!fIsNuBar) fHam(j,i) = fHms[i][j]/lv;
      else          fHam(j,i) = conj(fHms[i][j])/lv;
    }

    // Subtract NC coherent forward scattering from sterile neutrinos.
    // See arXiv:hep-ph/0606054v3, eq. 3.15, for example.
    if(i>2){
      if(!fIsNuBar) fHam(i,i) += kr2GNn;
      else          fHam(i,i) -= kr2GNn;
    }

  }

  // Add nue CC coherent forward scattering.
  if(!fIsNuBar) fHam(0,0) += kr2GNe;
  else          fHam(0,0) -= kr2GNe;

}

//.............................................................................
///
/// Solve the Homiltonian using the appropriate Eigen::MatrixType.\n
/// This will implement faster methods for N<5
///
template <typename T>
void PMNS_Sterile::SolveEigenSystem()
{

  Eigen::Ref<T> H(fHam);

  Eigen::SelfAdjointEigenSolver<T> eigensolver(H);

  // Fill fEval and fEvec vectors from GSL objects
  for(int i=0; i<fNumNus; i++){
    fEval[i] = eigensolver.eigenvalues()(i);
    for(int j=0; j<fNumNus; j++){
      fEvec[i][j] = eigensolver.eigenvectors()(i,j);
    }
  }

}

//.............................................................................
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

  UpdateHam();

  // Choose the appropriate MatrixType for the solver
  if     (fNumNus==4) SolveEigenSystem<Eigen::Matrix4cd>();
  else if(fNumNus==3) SolveEigenSystem<Eigen::Matrix3cd>();
  else if(fNumNus==2) SolveEigenSystem<Eigen::Matrix2cd>();
  else                SolveEigenSystem<Eigen::MatrixXcd>();

  // Mark eigensystem as solved
  fGotES = true;

  // Fill cache if activated
  FillCache();

}

///////////////////////////////////////////////////////////////////////////////
