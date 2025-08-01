///////////////////////////////////////////////////////////////////////////////
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework with decoherence.
//
// This  class inherits from the PMNS_DensityMatrix class
//
// jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "PMNS_Deco.h"

using namespace OscProb;

using namespace std;

//.............................................................................
///
/// Constructor. \sa PMNS_DensityMatrix::PMNS_DensityMatrix
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_Deco::PMNS_Deco() : PMNS_DensityMatrix(), fGamma()
{
  SetGamma(2, 0);
  SetGamma(3, 0);
  SetDecoAngle(0);
  SetPower(0);
}

//.............................................................................
///
/// Nothing to clean.
///
PMNS_Deco::~PMNS_Deco() {}

//.............................................................................
///
/// Set the decoherence parameter \f$\Gamma_{j1}\f$.
///
/// Requires that j = 2 or 3. Will notify you if input is wrong.
///
/// @param j   - The first mass index
/// @param val - The absolute value of the parameter
///
void PMNS_Deco::SetGamma(int j, double val)
{
  if (j < 2 || j > 3) {
    cerr << "WARNING: Gamma_" << j << 1 << " not valid for " << fNumNus
         << " neutrinos. Doing nothing." << endl;
    return;
  }

  if (val < 0) {
    cerr << "WARNING: Gamma_" << j << 1 << " must be positive. "
         << "Setting it to absolute value of input: " << -val << endl;
    val = -val;
  }

  fGamma[j - 1] = val;
}

//.............................................................................
///
/// Set the decoherence parameter \f$\Gamma_{32}\f$.
///
/// In practice, this sets \f$\Gamma_{31}\f$ internally via the formula:
///
///   \f$
///   \Gamma_{31} = \Gamma_{32} + \Gamma_{21} \cos^2\theta +
///   \cos\theta \sqrt{\Gamma_{21} (4\Gamma_{32} - \Gamma_{21} (1 -
///   \cos^2\theta))} \f$
///
/// IMPORTANT: Note this needs to be used AFTER defining \f$\Gamma_{21}\f$ and
/// \f$\theta\f$.
///
/// @param val - The absolute value of the parameter
///
void PMNS_Deco::SetGamma32(double val)
{
  if (val < 0) {
    cerr << "WARNING: Gamma_32 must be positive. "
         << "Setting it to absolute value of input: " << -val << endl;
    val = -val;
  }

  double min32 = 0.25 * fGamma[1];

  if (fGamma[0] >= 0)
    min32 *= 1 - pow(fGamma[0], 2);
  else
    min32 *= 1 + 3 * pow(fGamma[0], 2);

  if (val < min32) {
    if (fGamma[0] >= 0 || 4 * val / fGamma[1] < 1) {
      fGamma[0] = sqrt(1 - 4 * val / fGamma[1]);
      fGamma[2] = 0.25 * fGamma[1] * (1 + 3 * pow(fGamma[0], 2));
    }
    else {
      fGamma[0] = -sqrt((4 * val / fGamma[1] - 1) / 3);
      fGamma[2] = 0.25 * fGamma[1] * (1 - pow(fGamma[0], 2));
    }

    cerr << "WARNING: Impossible to have Gamma32 = " << val
         << " with current Gamma21 and theta parameters." << endl
         << "         Changing the value of cos(theta) to " << fGamma[0]
         << endl;

    return;
  }

  double arg = fGamma[1] * (4 * val - fGamma[1] * (1 - pow(fGamma[0], 2)));

  if (arg < 0) {
    arg       = 0;
    fGamma[0] = sqrt(1 - 4 * val / fGamma[1]);
    cerr << "WARNING: Imaginary Gamma31 found. Changing the value of "
            "cos(theta) to "
         << fGamma[0] << endl;
  }

  double gamma31 = val + fGamma[1] * pow(fGamma[0], 2) + fGamma[0] * sqrt(arg);

  fGamma[2] = gamma31;

  if (fabs(val - GetGamma(3, 2)) > 1e-6) {
    cerr << "ERROR: Failed sanity check: GetGamma(3,2) = " << GetGamma(3, 2)
         << " != " << val << endl;
  }
}

//.............................................................................
///
/// Set the decoherence angle. This will define the relationship:
///
///   \f$
///   \Gamma_{32} = \Gamma_{31} + \Gamma_{21} \cos^2\theta -
///   \cos\theta \sqrt{\Gamma_{21} (4\Gamma_{31} - \Gamma_{21} (1 -
///   \cos^2\theta))} \f$
///
/// @param th - decoherence angle
///
void PMNS_Deco::SetDecoAngle(double th) { fGamma[0] = cos(th); }

//.............................................................................
///
/// Set the power index of the decoherence energy dependence.
/// This will multiply the Gammas such that the final parameters are:
///
///   \f$
///   \Gamma_{ij} \rightarrow \Gamma_{ij} \times \left(\frac{E}{[1
///   \mbox{GeV}]}\right)^n \f$
///
/// @param n - power index
///
void PMNS_Deco::SetPower(double n) { fPower = n; }

//.............................................................................
///
/// Get the decoherence angle value.
///
double PMNS_Deco::GetDecoAngle() { return acos(fGamma[0]); }

//.............................................................................
///
/// Get the power index for energy dependence.
///
double PMNS_Deco::GetPower() { return fPower; }

//.............................................................................
///
/// Get any given decoherence parameter.
///
/// Requires that i > j. Will notify you if input is wrong.
/// If i < j, will assume reverse order and swap i and j.
///
/// @param i  - The first mass index
/// @param j  - The second mass index
///
double PMNS_Deco::GetGamma(int i, int j)
{
  if (i < j) {
    cerr << "WARNING: First argument should be larger than second argument"
         << endl
         << "Setting reverse order (Gamma_" << j << i << "). " << endl;
    int temp = i;
    i        = j;
    j        = temp;
  }
  if (i < 1 || i > 3 || i <= j || j < 1) {
    cerr << "WARNING: Gamma_" << i << j << " not valid for " << fNumNus
         << " neutrinos. Returning 0." << endl;
    return 0;
  }

  if (j == 1) { return fGamma[i - 1]; }
  else {
    // Combine Gamma31, Gamma21 and an angle (fGamma[0] = cos(th)) to make
    // Gamma32
    double arg =
        fGamma[1] * (4 * fGamma[2] - fGamma[1] * (1 - pow(fGamma[0], 2)));

    if (arg < 0) return fGamma[1] - 3 * fGamma[2];

    return fGamma[2] + fGamma[1] * pow(fGamma[0], 2) - fGamma[0] * sqrt(arg);
  }
}

//.............................................................................
///
/// Simple index sorting of 3-vector
///
vectorI sort3(const vectorD& x)
{
  vectorI out(3, 0);

  // 1st element is smallest
  if (x[0] < x[1] && x[0] < x[2]) {
    // 3rd element is largest
    if (x[1] < x[2]) {
      out[1] = 1;
      out[2] = 2;
    }
    // 2nd element is largest
    else {
      out[1] = 2;
      out[2] = 1;
    }
  }
  // 2nd element is smallest
  else if (x[1] < x[2]) {
    out[0] = 1;
    // 3rd element is largest
    if (x[0] < x[2]) out[2] = 2;
    // 1st element is largest
    else
      out[1] = 2;
  }
  // 3rd element is smallest
  else {
    out[0] = 2;
    // 2nd element is largest
    if (x[0] < x[1]) out[2] = 1;
    // 1st element is largest
    else
      out[1] = 1;
  }

  return out;
}

//.............................................................................
///
/// Propagate the current neutrino state through a given path
///
/// @param p - A neutrino path segment
///
void PMNS_Deco::PropagatePath(NuPath p)
{
  // Set the neutrino path
  SetCurPath(p);

  // Solve for eigensystem
  SolveHam();

  // Rotate to effective mass basis
  RotateState(true, fEvec);

  // Some ugly way of matching gamma and dmsqr indices
  vectorI dm_idx = sort3(fDm);
  vectorI ev_idx = sort3(fEval);

  int idx[3];
  for (int i = 0; i < fNumNus; i++) idx[ev_idx[i]] = dm_idx[i];

  // Power law dependency of gamma parameters
  double energyCorr = pow(fEnergy, fPower);

  // Convert path length to eV
  double lengthIneV = kKm2eV * p.length;

  // Apply phase rotation to off-diagonal elements
  for (int j = 0; j < fNumNus; j++) {
    for (int i = 0; i < j; i++) {
      // Get the appropriate gamma parameters based on sorted indices
      double gamma_ij =
          GetGamma(max(idx[i], idx[j]) + 1, min(idx[i], idx[j]) + 1) *
          energyCorr;
      // Multiply by path length
      gamma_ij *= kGeV2eV * lengthIneV;
      // This is the standard oscillation phase
      double arg = (fEval[i] - fEval[j]) * lengthIneV;
      // apply phase to rho
      fRho[i][j] *= exp(-gamma_ij) * complexD(cos(arg), -sin(arg));
      fRho[j][i] = conj(fRho[i][j]);
    }
  }

  // Rotate back to flavour basis
  RotateState(false, fEvec);
}

////////////////////////////////////////////////////////////////////////
