///////////////////////////////////////////////////////////////////////////////
///
/// Implementation of oscillations of neutrinos in matter in a
/// three-neutrino framework with Neutrino Decay in the 3rd state.
///
/// This class expands the PMNS_Fast class including the decay of the 
/// second and third mass state of the neutrino through a decay constant
/// alpha_i=m_i/tau_i (eV^2), where m_i is the mass in the restframe and
/// tau_i is the lifetime in the restframe.
///
/// Diagonalization of the matrix is done using the Eigen library.
///
/// \author Victor Carretero - vcarretero\@km3net.de
/// \author Joao Coelho - jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <unsupported/Eigen/MatrixFunctions>

#include "complexsolver.h"
#include "PMNS_Decay.h"

using namespace OscProb;
using namespace std;

//.............................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_Decay::PMNS_Decay() : PMNS_Base() {

  fAlpha         = vectorD(fNumNus, 0);
  fHam           = Eigen::MatrixXcd(fNumNus, fNumNus);
  fHd            = matrixC(fNumNus, vectorC(fNumNus,zero));

}

//.............................................................................
///
/// Nothing to clean.
///
PMNS_Decay::~PMNS_Decay(){}

//.............................................................................
///
/// Set all mixing parameters at once.
///
/// @param th12    - The value of the mixing angle theta_12
/// @param th23    - The value of the mixing angle theta_23
/// @param th13    - The value of the mixing angle theta_13
/// @param deltacp - The value of the CP phase delta_13
///
void PMNS_Decay::SetMix(double th12, double th23, double th13, double deltacp)
{

  SetAngle(1,2, th12);
  SetAngle(1,3, th13);
  SetAngle(2,3, th23);
  SetDelta(1,3, deltacp);

}

//.............................................................................
///
/// Setting Alpha3 parameter, it must be possitive, and units are eV^2.
///
/// @param alpha3 = m3/tau3, mass and lifetime of the 3rd state in the restframe
///
void PMNS_Decay::SetAlpha3(double alpha3)
{

  if(alpha3<0){
    cerr << "WARNING: Alpha3 must be positive. Doing nothing." << endl;
    return;
  }

  fBuiltHms *= (fAlpha[2] == alpha3);
  fAlpha[2] = alpha3;

}

//.............................................................................
///
/// Setting Alpha2 parameter, it must be possitive, and units are eV^2.
///
/// @param alpha2 = m2/tau2, mass and lifetime of the 2nd state in the restframe
///
void PMNS_Decay::SetAlpha2(double alpha2)
{

  if(alpha2<0){
    cerr << "WARNING: Alpha2 must be positive. Doing nothing." << endl;
    return;
  }

  fBuiltHms *= (fAlpha[1] == alpha2);
  fAlpha[1] = alpha2;

}

//.............................................................................
///
/// Reimplement SetIsNuBar to also rebuild hamiltonian.
///
/// @param isNuBar - true if antineutrinos
///
void PMNS_Decay::SetIsNuBar(bool isNuBar)
{

  // Check if value is actually changing
  fGotES *= (fIsNuBar == isNuBar);
  fBuiltHms *= (fIsNuBar == isNuBar);
  fIsNuBar = isNuBar;

}

//.............................................................................
///
/// Set both mass-splittings at once.
///
/// These are Dm_21 and Dm_32 in eV^2.\n
/// The corresponding Dm_31 is set in the class attributes.
///
/// @param dm21 - The solar mass-splitting Dm_21 
/// @param dm32 - The atmospheric mass-splitting Dm_32 
///
void PMNS_Decay::SetDeltaMsqrs(double dm21, double dm32)
{

  SetDm(2, dm21);
  SetDm(3, dm32 + dm21);

}

//.............................................................................
///
/// Get alpha3 parameter
///
/// @return alpha3
///
double PMNS_Decay::GetAlpha3()
{
  return fAlpha[2];
}

//.............................................................................
///
/// Get alpha2 parameter
///
/// @return alpha2
///
double PMNS_Decay::GetAlpha2()
{
  return fAlpha[1];
}


//.............................................................................
///
///  Implement building decay hamiltonian
///
void PMNS_Decay::BuildHms()
{

  // Check if anything changed
  if(fBuiltHms) return;

  // Tag to recompute eigensystem
  fGotES = false;

  for(int j=0; j<fNumNus; j++){
    // Set mass splitting
    fHms[j][j] = fDm[j];
    fHd[j][j] = fAlpha[j];
    // Reset off-diagonal elements
    for(int i=0; i<j; i++){
      fHms[i][j] = 0;
      fHd[i][j] = 0;
    }
    // Rotate j neutrinos
    for(int i=0; i<j; i++){
      RotateH(i,j,fHms);
      RotateH(i,j,fHd);
    }
  }

  // Taking care of antineutrinos delta-> -delta and filling the upper triangle
  // because final hamiltonian will be non-hermitian.
  for(int i=0; i<fNumNus; i++){
    for(int j=i+1; j<fNumNus; j++){
      if(fIsNuBar){
        fHms[i][j] = conj(fHms[i][j]);
        fHd[i][j] = conj(fHd[i][j]);
      }
      fHms[j][i] = conj(fHms[i][j]);
      fHd[j][i] = conj(fHd[i][j]);
    }
  }

  const complexD numi(0.0,1.0);
  ///Construct the total Hms+Hd
  for(int j=0; j<fNumNus; j++){
    for(int i=0; i<fNumNus; i++){
      fHms[i][j] = fHms[i][j] - numi*fHd[i][j];
    }
  }

  // Clear eigensystem cache
  //ClearCache();

  // Tag as built
  fBuiltHms = true;

}


//.............................................................................
///
/// Build the full Hamiltonian in matter.
///
/// Here we divide the mass squared matrix Hms by the 2E
/// to obtain the vacuum Hamiltonian in eV. Then, the matter
/// potential is added to the electron component.
///
void PMNS_Decay::UpdateHam()
{

  double lv = 2 * kGeV2eV*fEnergy;     // 2E in eV 

  double kr2GNe = kK2*M_SQRT2*kGf;
  kr2GNe *= fPath.density * fPath.zoa; // Matter potential in eV

  // Finish building Hamiltonian in matter with dimension of eV
  for(int i=0; i<fNumNus; i++){
    for(int j=0; j<fNumNus; j++){
      fHam(i,j) = fHms[i][j]/lv;
    }
  }

  if(!fIsNuBar) fHam(0,0) += kr2GNe;
  else          fHam(0,0) -= kr2GNe;

}

//.............................................................................
///
/// Solve the full Hamiltonian for eigenvalues.
///
/// This is needed only when estimating the oscillation frequency
/// for computing the AvgProb method.
///
/// This is using a method from the Eigen library
///
void PMNS_Decay::SolveHam()
{

  // Build Hamiltonian
  BuildHms();

  // Check if anything changed  
  if(fGotES) return;

  // Try caching if activated
  //if(TryCache()) return;

  UpdateHam();

  // Solve Hamiltonian for eigenvalues using the Eigen library method 
  complexsolver(fHam, fEval);

  fGotES = true;

  // Fill cache if activated
  //FillCache();

}

//.............................................................................
///
/// Propagate the neutrino through a constant density path
///
/// @param p - The neutrino path
///
void PMNS_Decay::PropagatePath(NuPath p)
{

  // Set the neutrino path
  SetCurPath(p);

  // Build Hamiltonian
  BuildHms();
  UpdateHam();

  double LengthIneV = kKm2eV * p.length;

  // Compute evolution operator exp(-I*H*L)
  Eigen::MatrixXcd H = fHam;
  H *= complexD(0, -LengthIneV);
  H = H.exp();

  // Propagate using evolution operator
  for(int i=0;i<fNumNus;i++){
    fBuffer[i] = 0;
    for(int j=0;j<fNumNus;j++){
      fBuffer[i] += H(i,j) * fNuState[j];
    }
  }

  // Update neutrino state
  for(int i=0;i<fNumNus;i++){
    fNuState[i] = fBuffer[i];
  }

}


////////////////////////////////////////////////////////////////////////
