///////////////////////////////////////////////////////////////////////////////
/// class OscProb::PMNS_SiderealLIV
///
/// Implementation of neutrino oscillations in matter in a
/// three-neutrino framework with Sidereal LIV as modelled by the SME.
/// The SME coefficients are included following the approach described in
/// https://doi.org/10.1103/PhysRevD.70.076002.
///
/// This developement is part of the QGRANT project with
/// ID: 101068013,
/// founded by the HORIZON-MSCA-2021-PF-01-01 programme.
///
/// \author Alba Domi - alba.domi\@fau.de
///////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cmath>
#include <iostream>

#include "MatrixDecomp/zheevh3.h"
#include "PMNS_SiderealLIV.h"

using namespace OscProb;

using namespace std;

//.............................................................................
/// Constructor.
/// This class is restricted to 3 neutrino flavours.
/// By default, all Sidereal LIV coefficients are set to zero.
///
PMNS_SiderealLIV::PMNS_SiderealLIV()
    : PMNS_Fast(),
      fa{},
      fc{},
      fN{0.0, 0.0, 0.0},
      fTime(0.0),
      fChi(0.0)
{
  SetStdPath();
}

//.............................................................................
///
/// Destructor.
///
PMNS_SiderealLIV::~PMNS_SiderealLIV() {}

//.............................................................................
/// Validate and reorder a flavour pair so that flvi <= flvj.
///
bool PMNS_SiderealLIV::ValidateAndOrder(int& flvi, int& flvj, const char* name)
{
  if (flvi > flvj) {
    cerr << "WARNING: First argument should be smaller or equal to second "
            "argument"
         << endl
         << "Setting reverse order (" << name << "_" << flvj << flvi << "). "
         << endl;
    std::swap(flvi, flvj);
  }
  if (flvi < 0 || flvi > 2 || flvj < flvi || flvj > 2) {
    cerr << "WARNING: " << name << "_" << flvi << flvj
         << " not valid for 3 neutrinos. Doing nothing." << endl;
    return false;
  }
  return true;
}

//.............................................................................
/// Validate that a Lorentz spatial index is in [0, 2].
///
bool PMNS_SiderealLIV::ValidateCoord(int coord, const char* coordName,
                              const char* paramName)
{
  if (coord < 0 || coord > 2) {
    cerr << "WARNING: " << coordName << " = " << coord
         << " out of range [0,2] for " << paramName
         << ". Doing nothing." << endl;
    return false;
  }
  return true;
}

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
void PMNS_SiderealLIV::SetA(int flvi, int flvj, int coord, double val)
{
  if (!ValidateAndOrder(flvi, flvj, "aT")) return;
  if (!ValidateCoord(coord, "coord", "aT")) return;

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
void PMNS_SiderealLIV::SetC(int flvi, int flvj, int coord1, int coord2, double val)
{
  if (!ValidateAndOrder(flvi, flvj, "cT")) return;
  if (!ValidateCoord(coord1, "coord1", "cT")) return;
  if (!ValidateCoord(coord2, "coord2", "cT")) return;

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
double PMNS_SiderealLIV::GetA(int flvi, int flvj, int coord)
{
  if (!ValidateAndOrder(flvi, flvj, "aT")) return 0;
  if (!ValidateCoord(coord, "coord", "aT")) return 0;

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
double PMNS_SiderealLIV::GetC(int flvi, int flvj, int coord1, int coord2)
{
  if (!ValidateAndOrder(flvi, flvj, "cT")) return 0;
  if (!ValidateCoord(coord1, "coord1", "cT")) return 0;
  if (!ValidateCoord(coord2, "coord2", "cT")) return 0;

  return fc[flvi][flvj][coord1][coord2];
}


//.............................................................................
///
/// Set the neutrino direction using the stored colatitude (fChi).
///
/// @param zenith  - zenith angle in radians
/// @param azimuth - azimuth angle in radians
///
void PMNS_SiderealLIV::SetNeutrinoDirection(double zenith, double azimuth)
{
  double chi = fChi * M_PI / 180.0;
  double N0  =  cos(chi) * sin(zenith) * cos(azimuth) + sin(chi) * cos(zenith);
  double N1  =  sin(zenith) * sin(azimuth);
  double N2  = -sin(chi) * sin(zenith) * cos(azimuth) + cos(chi) * cos(zenith);

  fGotES *= (fN[0] == N0 && fN[1] == N1 && fN[2] == N2);

  fN[0] = N0;
  fN[1] = N1;
  fN[2] = N2;
}


//.............................................................................
///
/// Set the local sidereal time at which the Hamiltonian is evaluated.
///
/// @param hours - Local sidereal time, in hours.
///
void PMNS_SiderealLIV::SetTimeHours(double hours)
{
  fGotES *= (fTime == hours);
  fTime = hours;
}


//.............................................................................
///
/// Set the detector colatitude (chi = 90 - latitude), in degrees.
///
/// @param chi - Colatitude of the detector, in degrees.
///
void PMNS_SiderealLIV::SetColatitude(double chi)
{
  fChi = chi;
}

//.............................................................................
///
/// Set the detector colatitude from geographic coordinates.
/// All three arguments carry the sign: negative for South, positive for North.
/// chi = 90 - latitude.
///
/// @param deg - Degrees of latitude
/// @param min - Arcminutes of latitude
/// @param sec - Arcseconds of latitude (default 0)
///
void PMNS_SiderealLIV::SetColatitude(double deg, double min, double sec)
{
  double latitude = deg + min / 60.0 + sec / 3600.0;
  fChi            = 90.0 - latitude;
}

//.............................................................................
///
/// Get the detector colatitude (chi), in degrees.
///
/// @return - The colatitude of the detector, in degrees.
///
double PMNS_SiderealLIV::GetColatitude() const
{
  return fChi;
}


//.............................................................................
///
/// Build the full LIVS Hamiltonian in matter
///
void PMNS_SiderealLIV::UpdateHam()
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

  for (int i = 0; i < fNumNus; i++) {
    for (int j = i; j < fNumNus; j++) {

      double sign = fIsNuBar ? -kGeV2eV : kGeV2eV;

      double C0  = -sign * fa[i][j][2] * fN[2];

      double As0 =  sign * (fa[i][j][0]*fN[1] - fa[i][j][1]*fN[0]);
      double Ac0 = -sign * (fa[i][j][0]*fN[0] + fa[i][j][1]*fN[1]);

      double As1 =  2*fN[1]*fN[2]*fc[i][j][0][2]
                  - 2*fN[0]*fN[2]*fc[i][j][1][2];
      double Ac1 = -2*fN[0]*fN[2]*fc[i][j][0][2]
                  - 2*fN[1]*fN[2]*fc[i][j][1][2];

      double Bs1 =  fN[0]*fN[1]*(fc[i][j][0][0] - fc[i][j][1][1])
                  - (fN[0]*fN[0] - fN[1]*fN[1])*fc[i][j][0][1];
      double Bs  = fEnergy * Bs1;

      double Bc1 = -0.5*(fN[0]*fN[0] - fN[1]*fN[1])
                       *(fc[i][j][0][0] - fc[i][j][1][1])
                  - 2.0*fN[0]*fN[1]*fc[i][j][0][1];
      double Bc  = fEnergy * Bc1;

      double As = As0 + fEnergy*As1;
      double Ac = Ac0 + fEnergy*Ac1;

      double liv_term = C0
                      + As * sin(kOmegaSidereal*fTime)
                      + Ac * cos(kOmegaSidereal*fTime)
                      + Bs * 2*sin(kOmegaSidereal*fTime)*cos(kOmegaSidereal*fTime)
                      + Bc * (cos(kOmegaSidereal*fTime)*cos(kOmegaSidereal*fTime)
                            - sin(kOmegaSidereal*fTime)*sin(kOmegaSidereal*fTime));

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
