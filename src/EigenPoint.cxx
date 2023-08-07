////////////////////////////////////////////////////////////////////////
//
// Struct to organise eigensystems for caching
//
//.................................................
//
// jcoelho\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

#include "EigenPoint.h"

using namespace OscProb;

//......................................................................
///
/// Constructor.
///
/// Uses number of neutrinos to fix eigensystem size.
///
/// @param numNus - the number of neutrino flavours
/// @param e      - the neutrino energy
/// @param p      - the neutrino path
/// @param n      - nu-nubar flag
///
EigenPoint::EigenPoint(int numNus, double e, NuPath p, bool n) :
fEval(numNus,0), fEvec(numNus, vectorC(numNus,0))
{
  SetVars(e,p,n);
}

//......................................................................
///
/// Set the eigensystem properties to new values
///
/// @param e      - the neutrino energy
/// @param p      - the neutrino path
/// @param n      - nu-nubar flag
///
void EigenPoint::SetVars(double e, NuPath p, bool n)
{
  fEnergy = e;
  fPath = p;
  fNubar = n;
  SetNE();
}

//......................................................................
///
/// Compute the combined energy-density parameter
///
void EigenPoint::SetNE()
{
  fNE = fEnergy * fPath.density * fPath.zoa;
  if(fNE < 1e-12) fNE = 1e-12;
  if(fNubar) fNE = -fNE;
}

//......................................................................
///
/// Comparison operator used for sorting into set
///
bool EigenPoint::operator < (const EigenPoint &rhs) const {
  if(fNE == rhs.fNE) return fPath.zoa < rhs.fPath.zoa;
  return fNE < rhs.fNE;
}

//......................................................................
///
/// Identity operator used for finding existing eigensystems
///
bool EigenPoint::operator == (const EigenPoint &rhs) const {
  return (fNE == rhs.fNE) && (fPath.zoa == rhs.fPath.zoa);
}
