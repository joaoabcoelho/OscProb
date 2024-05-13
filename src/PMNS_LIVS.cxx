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
/// By default, all LIVS coefficients are set to zero.
///
PMNS_LIVS::PMNS_LIVS() : PMNS_Fast()
{
  SetStdPath();

  N[0] = 0;
  N[1] = 0;
  N[2] = 0;

  omega = 2*M_PI/23.933;
  T = 0; 
  
  
  for (int flvi = 0; flvi < 3; flvi++) {
    for (int flvj = flvi; flvj < 3; flvj++) {
      SetaT(flvi, flvj, 0);
      SetcT(flvi, flvj, 0);
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
/// @param dim   - Dimension of the coefficient: can be 3,5,7.
/// @param val   - Absolute value of the parameter
/// @param phase - Complex phase of the parameter in radians
///
void PMNS_LIVS::SetaT(int flvi, int flvj, double val)
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

  faT[flvi][flvj] = val;
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
/// @param flvi  - First flavour index
/// @param flvj  - Second flavour index
/// @param dim   - Dimension of the coefficient: can be 4,6,8.
/// @param val   - Absolute value of the parameter
/// @param phase - Complex phase of the parameter in radians
///
void PMNS_LIVS::SetcT(int flvi, int flvj, double val)
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

  fcT[flvi][flvj] = val;
}

//.............................................................................
///
/// Get any given aT parameter of a chosen dimension.
///
/// Flavours are:\n
/// - 0 = nue, 1 = numu, 2 = nutau
///
/// Requires that flvi < flvj. Will notify you if input is wrong.
/// If flvi > flvj, will assume reverse order and swap flvi and flvj.
///
/// @param flvi  - The first flavour index
/// @param flvj  - The second flavour index
/// @param dim   - Dimension of the coefficient: can be 3,5,7.
///
/// @return - The aT parameter value
///
double PMNS_LIVS::GetaT(int flvi, int flvj)
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

  return faT[flvi][flvj];
}

//.............................................................................
///
/// Get any given cT parameter of a chosen dimension.
///
/// Flavours are:\n
/// - 0 = nue, 1 = numu, 2 = nutau
///
/// Requires that flvi < flvj. Will notify you if input is wrong.
/// If flvi > flvj, will assume reverse order and swap flvi and flvj.
///
/// @param flvi  - The first flavour index
/// @param flvj  - The second flavour index
/// @param dim   - Dimension of the coefficient: can be 4,6,8.
///
/// @return - The cT parameter value
///
double PMNS_LIVS::GetcT(int flvi, int flvj)
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

  return fcT[flvi][flvj];
}



void PMNS_LIVS::SetNeutrinoDirection(double zenith, double azimuth, double chi){
  N[0] = cos(chi)*sin(zenith)*cos(azimuth) + sin(chi)*cos(zenith);
  N[1] = sin(zenith)*sin(azimuth);
  N[3] = -sin(chi)*sin(zenith)*cos(azimuth) + cos(chi)*cos(zenith);
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


  double R = sqrt(N[0]*N[0] + N[1]*N[1]);
  double phi_orientation = atan(N[1]/N[0]);

  omega = 2*M_PI/23.933;
  
  for (int i = 0; i < fNumNus; i++) {
    for (int j = i; j < fNumNus; j++) {

      double liv_term = (N[1]*sin(omega*T-phi_orientation) - N[0]*cos(omega*T-phi_orientation));

      if (fIsNuBar)
	liv_term *= -faT[i][j]*kGeV2eV;
      else
	liv_term *= faT[i][j]*kGeV2eV;

      cout << "liv_term " << liv_term << ", ham " << fHam[i][j] << endl;
      
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
