///////////////////////////////////////////////////////////////////////////////
//
// Struct to organise S and K matrices for caching
//
//.................................................
//
// jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#include "SKMatricesPoint.h"

using namespace OscProb;

//.............................................................................
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
SKMatricesPoint::SKMatricesPoint(int numNus, double e, std::vector<NuPath>p, bool n)
    : fS(numNus, vectorC(numNus, 0)), fK(numNus, vectorC(numNus, 0))
{
  SetVars(e, p, n);
}

//.............................................................................
///
/// Set the eigensystem properties to new values
///
/// @param e      - the neutrino energy
/// @param p      - the neutrino path
/// @param n      - nu-nubar flag
///
void SKMatricesPoint::SetVars(double e, std::vector<NuPath>p, bool n)
{
  fEnergy = e;
  fNuPaths = p;
  fNubar  = n;
  //SetNE();
}

//.............................................................................
///
/// Compute the combined energy-density parameter
///
/*void SKMatricesPoint::SetNE()
{
  fNE = fEnergy * fPath.density * fPath.zoa;
  if (fNE < 1e-12) fNE = 1e-12;
  if (fNubar) fNE = -fNE;
}

//.............................................................................
///
/// Comparison operator used for sorting into set
///
bool SKMatricesPoint::operator<(const SKMatricesPoint& rhs) const
{
  if (fNE == rhs.fNE) return fPath.zoa < rhs.fPath.zoa;
  return fNE < rhs.fNE;
}

//.............................................................................
///
/// Identity operator used for finding existing eigensystems
///
bool SKMatricesPoint::operator==(const SKMatricesPoint& rhs) const
{
  return (fNE == rhs.fNE) && (fPath.zoa == rhs.fPath.zoa);
}*/
