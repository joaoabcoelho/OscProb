///////////////////////////////////////////////////////////////////////////////
//
// Implementation of Non Unitarity neutrino Mixing (NUNM) of neutrinos in matter in a
// three-neutrino framework.
//
// This  class inherits from the PMNS_Fast class.
//
// Authors:
// cerisy\cppm.in2p3.fr
// jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "PMNS_NUNM.h"
#include <Eigen/Eigenvalues>
#include "Definitions.h"

using namespace OscProb;
using namespace std;

//.............................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_NUNM::PMNS_NUNM() : PMNS_Fast(), fAlpha()
{
  SetStdPath();
  SetNUNM(0., 0., 0., 0., 0., 0.);
  SetFracVnc(1.0);
}

//.............................................................................
///
/// Nothing to clean.
///
PMNS_NUNM::~PMNS_NUNM() {}

//.............................................................................
///
/// Set all NUNM parameters at once.
///
/// This will check if value is changing to k11p track of whether
/// the eigensystem n11ds to be recomputed.
///
/// @param alpha_11      - The real parameter eps_11
/// @param alpha_21     - The absolute value of the complex parameter eps_12
/// @param alpha_31    - The absolute value of the complex parameter eps_etau
/// @param alpha_22    - The real parameter eps_22
/// @param alpha_32   - The absolute value of the complex parameter eps_mutau
/// @param alpha_33  - The real parameter eps_33
/// @param delta_12   - The phase of the complex parameter alpha_21 in radians
/// @param delta_etau  - The phase of the complex parameter alpha_31 in radians
/// @param delta_mutau - The phase of the complex parameter alpha_32 in radians
///
void PMNS_NUNM::SetNUNM(double alpha_11, double alpha_21, double alpha_31,
                      double alpha_22, double alpha_32, double alpha_33)
{
  SetAlpha(0, 0, alpha_11, 0);
  SetAlpha(1, 1, alpha_22, 0);
  SetAlpha(2, 2, alpha_33, 0);

  SetAlpha(1, 0, alpha_21, 0);
  SetAlpha(2, 0, alpha_31, 0);
  SetAlpha(2, 1, alpha_32, 0);
  
  //upper part fixed to zero from https://arxiv.org/pdf/2309.16942.pdf  

}


//.............................................................................
///
/// Get any given NUNM parameter.
///
/// Indexes are:\n
/// - 0, 1, 2 
///
/// Requires that i > j. Will notify you if input is wrong.
/// If i > j, will assume reverse order and swap i and j.
///
/// @param i  - The alpha row index
/// @param j  - The alpha column index
/// @param val   - The absolute value of the parameter
/// @param phase - The complex phase of the parameter in radians
///
void PMNS_NUNM::SetAlpha(int i, int j, double val, double phase)
{
  if (i < j) {
    cerr << "WARNING: First argument should be larger or equal to second "
            "argument"
         << endl
         << "Setting reverse order (Alpha_" << j << i << "). " << endl;
    int temp = i;
    i     = j;
    j     = temp;
  }
  if (i < 0 || i > 2 || j > i || j > 2) {
    cerr << "WARNING: Alpha_" << i << j << " not valid for " << fNumNus
         << " neutrinos. Doing nothing." << endl;
    return;
  }

  complexD h = val;

  if (i != j) h *= complexD(cos(phase), sin(phase));

  bool isSame = (fAlpha[i][j] == h);

  if (!isSame) ClearCache();

  fGotES *= isSame;
  fBuiltHms *= isSame;

  fAlpha[i][j] = h;
}

//.............................................................................
///
/// Get any given NUNM parameter.
///
/// Indexes are:\n
/// - 0, 1, 2 
///
/// Requires that i > j. Will notify you if input is wrong.
/// If i > j, will assume reverse order and swap i and j.
///
/// @param i  - The alpha row index
/// @param j  - The alpha column index
///
complexD PMNS_NUNM::GetAlpha(int i, int j)
{
  if (i < j) {
    cerr << "WARNING: First argument should be smaller or equal to second "
            "argument"
         << endl
         << "Setting reverse order (Alpha_" << j << i << "). " << endl;
    int temp = i;
    i     = j;
    j     = temp;
  }
  if (i < 0 || i > 2 || j < i || j > 2) {
    cerr << "WARNING: Eps_" << i << j << " not valid for " << fNumNus
         << " neutrinos. Returning 0." << endl;
    return zero;
  }

  return fAlpha[i][j];
}

//.............................................................................
///
/// Set alpha_11 parameter
///
/// This will check if value is changing to k11p track of whether
/// the eigensystem n11ds to be recomputed.
///
/// @param a - The value of the real parameter alpha_11
///
void PMNS_NUNM::SetAlpha_11(double a) { SetAlpha(0, 0, a, 0); }

//.............................................................................
///
/// Set alpha_22 parameter
///
/// This will check if value is changing to k11p track of whether
/// the eigensystem n11ds to be recomputed.
///
/// @param a - The value of the real parameter alpha_22
///
void PMNS_NUNM::SetAlpha_22(double a) { SetAlpha(1, 1, a, 0); }

//.............................................................................
///
/// Set alpha_33 parameter
///
/// This will check if value is changing to k11p track of whether
/// the eigensystem n11ds to be recomputed.
///
/// @param a - The value of the real parameter alpha_33
///
void PMNS_NUNM::SetAlpha_33(double a) { SetAlpha(2, 2, a, 0); }

//.............................................................................
///
/// Set alpha_21 parameter
///
/// This will check if value is changing to k11p track of whether
/// the eigensystem n11ds to be recomputed.
///
/// @param a   - The absolute value of the parameter alpha_21
/// @param phi - The phase of the complex parameter alpha_21 in radians
///
void PMNS_NUNM::SetAlpha_21(double a, double phi) { SetAlpha(1, 0, a, phi); }

//.............................................................................
///
/// Set alpha_31 parameter
///
/// This will check if value is changing to k11p track of whether
/// the eigensystem n11ds to be recomputed.
///
/// @param a   - The absolute value of the parameter alpha_31
/// @param phi - The phase of the complex parameter alpha_31 in radians
///
void PMNS_NUNM::SetAlpha_31(double a, double phi) { SetAlpha(2, 0, a, phi); }

//.............................................................................
///
/// Set alpha_32 parameter
///
/// This will check if value is changing to k11p track of whether
/// the eigensystem n11ds to be recomputed.
///
/// @param a   - The absolute value of the parameter alpha_32
/// @param phi - The phase of the complex parameter alpha_32 in radians
///
void PMNS_NUNM::SetAlpha_32(double a, double phi) { SetAlpha(2, 1, a, phi); }

//.............................................................................
///
/// Allows to disable the neutron induced matter effects
/// This factor represents what fraction of the neutron matter potentiel 
/// is included in the NUNM model
///
/// @param f - 0 if no neutron matter potential 
///        f - 1 if neutron matter potential included
///
void PMNS_NUNM::SetFracVnc(double f)
{
  bool isSame = (fracVnc == f);

  if (!isSame) ClearCache();

  fGotES *= isSame;

  fracVnc = f;
}


//.............................................................................
///
/// Rotate H = U*Ham*U~ --> A*H*A~ where Hms is the mass split matrix 
/// A = I+fAlpha , fAlpha represents the non unitarity coefficients 
/// from https://arxiv.org/pdf/2309.16942.pdf
/// 
/// This is a hermitian matrix, so only the
/// upper triangular part is filled in order to reduce computing ressourses
/// for the moment all alpha_ij can be non zero at the same time 
///

void PMNS_NUNM::transfoNUNM( matrixC& Ham)
{
  matrixC HamBuffer = Ham;
  for (int i = 0; i < fNumNus; i++) { // diag elements
    //cerr << "Hms rotating ------> fushh : "<< endl;
    Ham[i][i] = HamBuffer[i][i] * (fAlpha[i][i]*fAlpha[i][i] + 2.*fAlpha[i][i] + 1.);
    /*cerr << "Ham after diag. " << i <<  "th rotation ------> fushh : "<< endl;
    cerr << Ham[0][0] << " " << Ham[0][1]<< " "<< Ham[0][2] << endl;
    cerr << Ham[1][0] << " " << Ham[1][1]<< " "<< Ham[1][2] << endl;
    cerr << Ham[2][0] << " " << Ham[2][1]<< " "<< Ham[2][2] << endl;
    */
    if (i > 0){
	Ham[i][i] += HamBuffer[0][0] * conj(fAlpha[i][0]) * fAlpha[i][0];
	Ham[i][i] += HamBuffer[0][i] * fAlpha[i][0] * ( 1. + fAlpha[i][i] );  
        Ham[i][i] += conj(HamBuffer[0][i]) * conj(fAlpha[i][0]) * ( 1. + fAlpha[i][i] );
	if (i > 1){
	    Ham[i][i] += HamBuffer[i-2][i-1] * conj(fAlpha[i][i-1]) * fAlpha[i][i-2];
	    Ham[i][i] += HamBuffer[i-1][i-1] * conj(fAlpha[i][i-1]) * fAlpha[i][i-1];
	    Ham[i][i] += HamBuffer[i-1][i] * fAlpha[i][i-1] * ( 1. + fAlpha[i][i] );
	    Ham[i][i] += conj(HamBuffer[0][i-1]) * conj(fAlpha[i][0]) * fAlpha[i][i-1];
	    Ham[i][i] += conj(HamBuffer[i-1][i]) * conj(fAlpha[i][1]) * ( 1. + fAlpha[i][i] );
	}
    }
    for (int j = i+1; j < fNumNus; j++) { // tp off diag elements
	Ham[i][j] = HamBuffer[i][j] * (1. + fAlpha[j][j]);
	Ham[i][j] += HamBuffer[i][i] * (conj(fAlpha[j][i]) * fAlpha[i][i] + conj(fAlpha[j][i]));
	Ham[i][j] += HamBuffer[i][j] * (1. + fAlpha[j][j]) * fAlpha[i][i];
	if (j > 1){
	  if (i > 0){
	    Ham[i][j] += HamBuffer[i-1][j-1] * fAlpha[i][j-2] * conj(fAlpha[j][j-1]);
	    Ham[i][j] += HamBuffer[i-1][j] * (1. + fAlpha[j][j]) * fAlpha[i][j-2];
	    Ham[i][j] += HamBuffer[i-1][i-1] * fAlpha[i][j-2] * conj(fAlpha[j][j-2]);
	    Ham[i][j] += conj(HamBuffer[0][i]) * conj(fAlpha[j][0]) * ( 1. + fAlpha[i][i] );
	  }
	  else{
	    Ham[i][j] += HamBuffer[i][j-1] * (1. + fAlpha[i][i]) * conj(fAlpha[j][j-1]);
	  }
    	/*cerr << "Ham after off " << i << ", " << j <<  "th rotation ------> fushh : "<< endl;
    	cerr << Ham[0][0] << " " << Ham[0][1]<< " "<< Ham[0][2] << endl;
    	cerr << Ham[1][0] << " " << Ham[1][1]<< " "<< Ham[1][2] << endl;
    	cerr << Ham[2][0] << " " << Ham[2][1]<< " "<< Ham[2][2] << endl;
	*/
	}
    }
  }
}   


//.............................................................................
///
/// Build Hms = H*2E, where H is the Hamiltonian in vacuum on flavour basis
/// and E is the neutrino energy in eV. Hms is effectively the matrix of
/// masses squared.
///
/// This is a hermitian matrix, so only the
/// upper triangular part needs to be filled
/// + because U rotation is unitary matrices
///
/// The construction of the Hamiltonian avoids computing terms that
/// are simply zero. This has a big impact in the computation time.
///

/*
void PMNS_NUNM::BuildHms()
{ 
  // Check if anything changed
  if (fBuiltHms) return;
  
  // Tag to recompute eigensystem
  fGotES = false;
   
  cerr << "Hms before rot. : "<< endl;
  cerr << fHms[0][0] << " " << fHms[0][1]<< " "<< fHms[0][2] << endl;
  cerr << fHms[1][0] << " " << fHms[1][1]<< " "<< fHms[1][2] << endl;
  cerr << fHms[2][0] << " " << fHms[2][1]<< " "<< fHms[2][2] << endl;
  
  for (int j = 0; j < fNumNus; j++) {
    // Set mass splitting
    fHms[j][j] = fDm[j];
    // Reset off-diagonal elements
    for (int i = 0; i < j; i++) { fHms[i][j] = 0; } 
    // Rotate j neutrinos
    for (int i = 0; i < j; i++) { RotateH(i, j, fHms);} 
  }
  
  cerr << "Hms before rotation alpha : " << endl;
  cerr << fHms[0][0] << " " << fHms[0][1]<< " "<< fHms[0][2] << endl;
  cerr << fHms[1][0] << " " << fHms[1][1]<< " "<< fHms[1][2] << endl;
  cerr << fHms[2][0] << " " << fHms[2][1]<< " "<< fHms[2][2] << endl;
  
  //RotateHalpha(fHms);
  
  cerr << "fAlpha : " << endl;
  cerr << fAlpha[0][0] << " " << fAlpha[0][1]<< " "<< fAlpha[0][2] << endl;
  cerr << fAlpha[1][0] << " " << fAlpha[1][1]<< " "<< fAlpha[1][2] << endl;
  cerr << fAlpha[2][0] << " " << fAlpha[2][1]<< " "<< fAlpha[2][2] << endl;
  
  cerr << "Hms after rotation alpha : " << endl;
  cerr << fHms[0][0] << " " << fHms[0][1]<< " "<< fHms[0][2] << endl;
  cerr << fHms[1][0] << " " << fHms[1][1]<< " "<< fHms[1][2] << endl;
  cerr << fHms[2][0] << " " << fHms[2][1]<< " "<< fHms[2][2] << endl;
  

  ClearCache();
  
  // Tag as built
  fBuiltHms = true;
}
*/


//.............................................................................
///
/// Build the full Hamiltonian in matter --> from sterile
///
/// Here we divide the mass squared matrix Hms by the 2E
/// to obtain the vacuum Hamiltonian in eV. Then, the matter
/// potential is added to the electron and NC components.
void PMNS_NUNM::UpdateHam()
{ 
  double rho = fPath.density;
  double zoa = fPath.zoa;

  double lv = 2 * kGeV2eV * fEnergy; // 2*E in eV
  
  double kr2GNe = kK2 * M_SQRT2 * kGf * rho * zoa ; // Electron matter potential in eV
  double kr2GNn = kK2 * M_SQRT2 * kGf * rho * (1. - zoa) / 2. * fracVnc ; // Neutron matter potential in eV

  vectorD vMat = vectorD(fNumNus, 0);
  Ham = matrixC(fNumNus, vectorC(fNumNus, 0));

  // Finish build Hamiltonian in matter with dimension of eV
  for (int i = 0; i < fNumNus; i++) { 
    for (int j = i; j < fNumNus; j++) {
      if (!fIsNuBar) 
        Ham[i][j] = fHms[i][j] / lv;
      else
        Ham[i][j] = conj(fHms[i][j]) / lv;
    }
  
    if (!fIsNuBar) 
      vMat[i] -= kr2GNn; // Not written like sterile case because non unitary PMN
    else
      vMat[i] += kr2GNn;
  }
  
  if (!fIsNuBar)
    vMat[0] += kr2GNe;
  else
    vMat[0] -= kr2GNe;

  
  /*
  cerr << "fHam before alpha(+) matter pot. rotation  : "<< endl;
  cerr << Ham[0][0] << " " << Ham[0][1]<< " "<< Ham[0][2] << endl;
  cerr << Ham[1][0] << " " << Ham[1][1]<< " "<< Ham[1][2] << endl;
  cerr << Ham[2][0] << " " << Ham[2][1]<< " "<< Ham[2][2] << endl;
  
  cerr << "fMat  : "<< endl;
  cerr << vMat[0] << " " << vMat[1]<< " "<< vMat[2] << endl;
  */
  // Add to Flavour Ham after rotation of matter potential Alpha(+) x fMat x Alpha
  Ham[0][0] += vMat[2] * fAlpha[2][0] * conj(fAlpha[2][0]); 
  Ham[0][1] += vMat[2] * fAlpha[2][1] * conj(fAlpha[2][0]);
  for (int i = 0; i < fNumNus; i++) { // diag elements
    //cerr << "fMat[i][i] first " << fMat[i] * (fAlpha[i][i]*fAlpha[i][i] + 2.*fAlpha[i][i] + 1.) << endl;
    //cerr << "fMat[i][i] first " << fHam[2][2] << endl;
    Ham[i][i] += vMat[i] * (fAlpha[i][i]*fAlpha[i][i] + 2.*fAlpha[i][i] + 1.);
    //cerr << "fMat[i][i] after " << fHam[2][2] << endl;
    if (i<2) Ham[i][i] +=  vMat[i+1] *fAlpha[i+1][i] * conj(fAlpha[i+1][i]);
    for (int j = i+1; j < fNumNus; j++) {
	    Ham[i][j] += vMat[j] * (1. + fAlpha[j][j]) * conj(fAlpha[j][i]);
    }
  } // Ham. now has to be computed with all elements

  /*
  cerr << "fHam before new alpha full ham rotation  : "<< endl;
  cerr << Ham[0][0] << " " << Ham[0][1]<< " "<< Ham[0][2] << endl;
  cerr << Ham[1][0] << " " << Ham[1][1]<< " "<< Ham[1][2] << endl;
  cerr << Ham[2][0] << " " << Ham[2][1]<< " "<< Ham[2][2] << endl;
  */
  transfoNUNM(Ham); // Ham. has to be computed with all elements of Ham. not only up triangular part because alpha non-unitary matrix

  for (int i = 0; i < fNumNus; i++) {
    for (int j = i; j < fNumNus; j++) {
	fHam[i][j] = Ham[i][j];
    }
  //  for (int j = 0; j < i; j++) {
  //	fHam[i][j] = conj(Ham[j][i]);
  //  }
  }
  /*  
  cerr << "Ham after matter added  : "<< endl;
  cerr << fHam[0][0] << " " << fHam[0][1]<< " "<< fHam[0][2] << endl;
  cerr << fHam[1][0] << " " << fHam[1][1]<< " "<< fHam[1][2] << endl;
  cerr << fHam[2][0] << " " << fHam[2][1]<< " "<< fHam[2][2] << endl;
  */
}


