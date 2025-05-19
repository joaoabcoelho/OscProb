///////////////////////////////////////////////////////////////////////////////
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework with gravitationally induced decoherence.
//
// This  class inherits from the PMNS_Base class
//
// This developement is part of the QGRANT project with
// ID: 101068013,
// founded by the HORIZON-MSCA-2021-PF-01-01 programme.
//
// \author Joao Coelho - jcoelho\@apc.in2p3.fr
// \author Alba Domi - alba.domi\@fau.de
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "PMNS_GQD.h"

using namespace OscProb;

using namespace std;

// Define some constants from PDG 2015
const double PMNS_GQD::kkB   = 8.6173303e-5;      // Boltzmann constant [eV/K]
const double PMNS_GQD::kc    = 299792.458;        // Speed of light [km/s]
const double PMNS_GQD::khbar = 1 / (kKm2eV * kc); // Planck constant [eV.s]

//.............................................................................
///
/// Constructor. \sa PMNS_Deco::PMNS_Deco
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_GQD::PMNS_GQD() : PMNS_Deco()
{
  SetStdPath();
  SetEta(0);
  SetTemperature(0);
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
}

//.............................................................................
///
/// Get the coupling parameter \f$\Eta\f$ in [s].
///
/// @return Eta parameter
///
double PMNS_GQD::GetEta() { return fEta; }

//.............................................................................
///
/// Get the parameter \f$\Temperature\f$ in [K].
///
/// @return Environment temperature
///
double PMNS_GQD::GetTemperature() { return fT; }

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

  // Coupling factor in eV^-2
  double coup = 4 * pow(fEta / khbar, 2) * (kkB * fT) * lengthIneV;

  // Apply phase rotation to off-diagonal elements
  for (int j = 0; j < fNumNus; j++) {
    for (int i = 0; i < j; i++) {
      // Get the appropriate gamma parameters based on sorted indices
      double gamma_ij = coup * pow(fEval[i] - fEval[j], 2);
      // This is the standard oscillation phase
      double arg = (fEval[i] - fEval[j]) * lengthIneV;
      // apply phase to rho
      fRho[i][j] *= exp(-gamma_ij) * complexD(cos(arg), -sin(arg));
      fRho[j][i] = conj(fRho[i][j]);
    }
  }

  // Rotate back to flavour basis
  RotateState(false);
}
