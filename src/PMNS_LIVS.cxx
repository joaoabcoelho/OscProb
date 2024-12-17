///////////////////////////////////////////////////////////////////////////////
/// class OscProb::PMNS_LIVS
///
/// Implementation of neutrino oscillations in matter in a
/// three-neutrino framework with LIVS as modelled by the SME.
/// The SME coefficients are included up to the 8th order,
/// following the approach described in
/// https://doi.org/10.1103/PhysRevD.85.096005.
///
/// This developement is part of the QGRANT project with
/// ID: 101068013,
/// founded by the HORIZON-MSCA-2021-PF-01-01 programme.
///
/// \author Alba Domi - alba.domi\@fau.de
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "MatrixDecomp/zheevh3.h"
#include "PMNS_LIVS.h"

using namespace OscProb;

using namespace std;

//.............................................................................
/// Constructor.
/// This class is restricted to 3 neutrino flavours.
/// By default, all Sidereal LIV coefficients are set to zero.
///
PMNS_LIVS::PMNS_LIVS() : PMNS_Fast()
{
  SetStdPath();

  N[0] = 0;
  N[1] = 0;
  N[2] = 0;
  
  zenith = 0;
  azimuth = 0;
  chi = 0;
  
  omega = 2*M_PI/23.933;
  T = 0; 
  
  
  for (int flvi = 0; flvi < 3; flvi++) {
    for (int flvj = flvi; flvj < 3; flvj++) {
      for (int coord1 = 0; coord1 < 3; coord1++) {
	SetA(flvi, flvj, coord1, 0);
	for (int coord2 = 0; coord2 < 3; coord2++) {
	  SetC(flvi, flvj, coord1, coord2, 0);
	}
      }
    }
  }
}

//.............................................................................
///
/// Destructor.
///
PMNS_LIVS::~PMNS_LIVS() {}

//.............................................................................
/// Set any given aT parameter of a chosen dimension.
///
/// Flavours are:\n
/// - 0 = nue, 1 = numu, 2 = nutau
///
/// Requires that flvi < flvj. Will notify you if input is wrong.
/// If flvi > flvj, will assume reverse order and swap flvi and flvj.
///
/// @param flvi  - First flavour index
/// @param flvj  - Second flavour index
/// @param coord - Space coordinate: 0=X, 1=Y, 2=Z
/// @param val   - Absolute value of the parameter
///
void PMNS_LIVS::SetA(int flvi, int flvj, int coord, double val)
{
  if (flvi > flvj) {
    cerr << "WARNING: First argument should be smaller or equal to second "
            "argument"
         << endl
         << "Setting reverse order (aT_" << flvj << flvi << "). " << endl;
    int temp = flvi;
    flvi     = flvj;
    flvj     = temp;
  }
  if (flvi < 0 || flvi > 2 || flvj < flvi || flvj > 2) {
    cerr << "WARNING: aT_" << flvi << flvj << " not valid for " << fNumNus
         << " neutrinos. Doing nothing." << endl;
    return;
  }

  fa[flvi][flvj][coord] = val;
}

//.............................................................................
/// Set any given cT parameter of a chosen dimension.
///
/// Flavours are:\n
/// - 0 = nue, 1 = numu, 2 = nutau
///
/// Requires that flvi < flvj. Will notify you if input is wrong.
/// If flvi > flvj, will assume reverse order and swap flvi and flvj.
///
/// @param flvi   - First flavour index
/// @param flvj   - Second flavour index
/// @param coord1 - Space coordinate: 0=X, 1=Y, 2=Z
/// @param coord2 - Space coordinate: 0=X, 1=Y, 2=Z
/// @param val    - Absolute value of the parameter
///
void PMNS_LIVS::SetC(int flvi, int flvj, int coord1, int coord2, double val)
{
  if (flvi > flvj) {
    cerr << "WARNING: First argument should be smaller or equal to second "
            "argument"
         << endl
         << "Setting reverse order (cT_" << flvj << flvi << "). " << endl;
    int temp = flvi;
    flvi     = flvj;
    flvj     = temp;
  }
  if (flvi < 0 || flvi > 2 || flvj < flvi || flvj > 2) {
    cerr << "WARNING: cT_" << flvi << flvj << " not valid for " << fNumNus
         << " neutrinos. Doing nothing." << endl;
    return;
  }

  fc[flvi][flvj][coord1][coord2] = val;
}

//.............................................................................
///
/// Get any given aT parameter.
///
/// Flavours are:\n
/// - 0 = nue, 1 = numu, 2 = nutau
///
/// Requires that flvi < flvj. Will notify you if input is wrong.
/// If flvi > flvj, will assume reverse order and swap flvi and flvj.
///
/// @param flvi  - The first flavour index
/// @param flvj  - The second flavour index
/// @param coord - Space coordinate: 0=X, 1=Y, 2=Z
///
/// @return - The aT parameter value
///
double PMNS_LIVS::GetA(int flvi, int flvj, int coord)
{
  if (flvi > flvj) {
    cerr << "WARNING: First argument should be smaller or equal to second "
            "argument"
         << endl
         << "Setting reverse order (aT_" << flvj << flvi << "). " << endl;
    int temp = flvi;
    flvi     = flvj;
    flvj     = temp;
  }
  if (flvi < 0 || flvi > 2 || flvj < flvi || flvj > 2) {
    cerr << "WARNING: aT_" << flvi << flvj << " not valid for " << fNumNus
         << " neutrinos. Returning 0." << endl;
    return 0;
  }

  return fa[flvi][flvj][coord];
}

//.............................................................................
///
/// Get any given cT parameter.
///
/// Flavours are:\n
/// - 0 = nue, 1 = numu, 2 = nutau
///
/// Requires that flvi < flvj. Will notify you if input is wrong.
/// If flvi > flvj, will assume reverse order and swap flvi and flvj.
///
/// @param flvi   - The first flavour index
/// @param flvj   - The second flavour index
/// @param coord1 - Space coordinate: 0=X, 1=Y, 2=Z
/// @param coord2 - Space coordinate: 0=X, 1=Y, 2=Z
///
/// @return - The cT parameter value
///
double PMNS_LIVS::GetC(int flvi, int flvj, int coord1, int coord2)
{
  if (flvi > flvj) {
    cerr << "WARNING: First argument should be smaller or equal to second "
            "argument"
         << endl
         << "Setting reverse order (cT_" << flvj << flvi << "). " << endl;
    int temp = flvi;
    flvi     = flvj;
    flvj     = temp;
  }
  if (flvi < 0 || flvi > 2 || flvj < flvi || flvj > 2) {
    cerr << "WARNING: cT_" << flvi << flvj << " not valid for " << fNumNus
         << " neutrinos. Returning 0." << endl;
    return 0;
  }

  return fc[flvi][flvj][coord1][coord2];
}


//.............................................................................
///
/// Set the neutrino direction.
///
/// @param theta - zenith angle in (rad) 
/// @param phi   - azimuth angle in (rad) 
/// @param ch    - colatitude of the detector in (rad) 
///
void PMNS_LIVS::SetNeutrinoDirection(double theta, double phi, double ch){

  zenith = theta;
  azimuth = phi;
  chi = ch;
  
  N[0] = cos(chi)*sin(zenith)*cos(azimuth) + sin(chi)*cos(zenith);
  N[1] = sin(zenith)*sin(azimuth);
  N[2] = -sin(chi)*sin(zenith)*cos(azimuth) + cos(chi)*cos(zenith);
}


void PMNS_LIVS::SetTime(double hours){
  T = hours;
}


void PMNS_LIVS::SolveHam()
{
  // Build Hamiltonian                                                        
  BuildHms();

  UpdateHam();

  double   fEvalGLoBES[3];
  complexD fEvecGLoBES[3][3];

  // Solve Hamiltonian for eigensystem using the GLoBES method               
  zheevh3(fHam, fEvecGLoBES, fEvalGLoBES);

  // Fill fEval and fEvec vectors from GLoBES arrays                        
  for (int i = 0; i < fNumNus; i++) {
    fEval[i] = fEvalGLoBES[i];
    for (int j = 0; j < fNumNus; j++) { fEvec[i][j] = fEvecGLoBES[i][j]; }
  }

}


//.............................................................................
///
/// Build the full LIVS Hamiltonian in matter
///
void PMNS_LIVS::UpdateHam()
{
  double lv = 2 * kGeV2eV * fEnergy; // 2*E in eV

  // Set the vacuum Hamiltonian
  for (int i = 0; i < fNumNus; i++) {
    for (int j = i; j < fNumNus; j++) { fHam[i][j] = fHms[i][j] / lv; }
  }

  // Add matter potential
  double kr2GNe = kK2 * M_SQRT2 * kGf * fPath.density * fPath.zoa;
  if (!fIsNuBar)
    fHam[0][0] += kr2GNe;
  else
    fHam[0][0] -= kr2GNe;

  double phi_orientation = atan(N[1]/N[0]);

  double alpha = azimuth - M_PI;

  double C0 = 0;
  double As0 = 0, Ac0 = 0;
  double As1 = 0, Ac1 = 0;
  double Bs1 = 0, Bs = 0;
  double Bc1 = 0, Bc = 0;
  
  for (int i = 0; i < fNumNus; i++) {
    for (int j = i; j < fNumNus; j++) {

      double sign = 0;
      
      if (fIsNuBar)
	sign = -kGeV2eV;
      else
	sign = kGeV2eV;      

      C0  = -sign * fa[i][j][2]*N[2];
      
      As0 =  sign * (fa[i][j][0]*N[1] - fa[i][j][1]*N[0]);
      Ac0 = -sign * (fa[i][j][0]*N[0] - fa[i][j][1]*N[1]);

      As1 =  2*N[1]*N[2]*fc[i][j][0][2] - 2*N[0]*N[2]*fc[i][j][1][2];
      Ac1 = -2*N[0]*N[2]*fc[i][j][0][2] - 2*N[1]*N[2]*fc[i][j][1][2];

      Bs1 = N[0]*N[1]*(fc[i][j][0][0] - fc[i][j][1][1]) - (N[0]*N[0] - N[1]*N[1])*fc[i][j][0][1];
      Bs = fEnergy * Bs1; 

      Bc1 = -1./2.*(N[0]*N[0] - N[1]*N[1])*(fc[i][j][0][0] - fc[i][j][1][1]) - 2.*N[0]*N[1]*fc[i][j][0][1];

      Bc = fEnergy * Bc1;
      
      double As = As0 + fEnergy*As1;
      double Ac = Ac0 + fEnergy*Ac1;

      double liv_term = C0 + (As*sin(omega*T) +
			      Ac*cos(omega*T) +
			      Bs*2*sin(omega*T)*cos(omega*T) +
			      Bc*(cos(omega*T)*cos(omega*T) - sin(omega*T)*sin(omega*T)));

      fHam[i][j] += liv_term;
    }
  }
  
  // Conjugate Hamiltonian for antineutrinos
  if (fIsNuBar) {
    for (int i = 0; i < fNumNus; i++) {
      for (int j = i + 1; j < fNumNus; j++) { fHam[i][j] = conj(fHam[i][j]); }
    }
  }
  
}

///////////////////////////////////////////////////////////////////////////////
