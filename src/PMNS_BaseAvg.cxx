///////////////////////////////////////////////////////////////////////////////
//
// Base class for faster oscillations of neutrinos in matter in a
// n-neutrino framework using a Taylor Expansion.
//
//.............................................................................
//
// jcoelho\@apc.in2p3.fr
//  
///////////////////////////////////////////////////////////////////////////////

#include "PMNS_BaseAvg.h"
#include <Eigen/Eigenvalues>
#include "PremModel.h"

using namespace std;

using namespace OscProb;


//.............................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
PMNS_BaseAvg::PMNS_BaseAvg(int numNus) : PMNS_Base(numNus)
{
  fPrem.LoadModel("");
  InitializeTaylorsVectors();

}

//.............................................................................
/// Functions for averaging over oscillations
/// Introduces speed up to code
///
///
/// Set vector sizes and initialize elements to zero.
/// Initialize diagonal elements of S to one
///
void PMNS_BaseAvg::InitializeTaylorsVectors()
{
  fdensityMatrix = matrixC(fNumNus, vectorC(fNumNus, 0));

  flambdaInvE = vectorD(fNumNus, 0);
  fVInvE      = matrixC(fNumNus, vectorC(fNumNus, 0));
  fKInvE      = Eigen::MatrixXcd::Zero(fNumNus, fNumNus);

  flambdaCosT = vectorD(fNumNus, 0);
  fVcosT      = matrixC(fNumNus, vectorC(fNumNus, 0));
  fKcosT      = Eigen::MatrixXcd::Zero(fNumNus, fNumNus);

  fevolutionMatrixS = matrixC(fNumNus, vectorC(fNumNus, 0));

  fSflavor = matrixC(fNumNus, vectorC(fNumNus, 0));
  fKmass   = matrixC(fNumNus, vectorC(fNumNus, 0));
  fKflavor = matrixC(fNumNus, vectorC(fNumNus, 0));

  for (int i = 0; i < fNumNus; i++) { fevolutionMatrixS[i][i] = 1; }

  fLayer = fPrem.GetPremLayers().size() - 1;

  fdl = -1;

  fDetRadius = fPrem.GetDetRadius();
}

////////////////////////////////////////////////////////////////////////
///
/// Copy the earth model used
///
/// This is done to get access to the PremLayer list to be used in the
/// LnDerivative() function.
///
/// @param prem - The earth model used
///
void PMNS_BaseAvg::SetPremModel(OscProb::PremModel& prem) { fPrem = prem; }

//.............................................................................
//.............................................................................
///
/// Set neutrino angle.
///
/// @param cosT - The cosine of the neutrino angle
///
void PMNS_BaseAvg::SetCosT(double cosT) { fcosT = cosT; }

///
/// Set bin's widths.
///
/// @param dE - The width of the bin for energy in GeV
/// @param dcosT - The width of the bin for angle
///
void PMNS_BaseAvg::SetwidthBin(double dE, double dcosT)
{
  fdInvE = dE;
  fdcosT = dcosT;
}


//.............................................................................
///
/// Build K matrix for the inverse of energy in mass basis.
///
/// The variable for which a Taylor expansion is done here is not directly the
/// energy but the inverse of it. This change of variable allow to express the
/// hamiltonien as linear with respect to this new variable.
///
/// @param L - The length of the layer in km
///
void PMNS_BaseAvg::BuildKE(double L)
{
  double lenghtEV = L * kKm2eV;     // L in eV-1
  double bufK     = lenghtEV * 0.5; // L/2 in eV-1

  complexD buffer[fNumNus];

  for (int i = 0; i < fNumNus; i++) {
    complexD Hms_kl;

    for (int l = 0; l < fNumNus; l++) {
      buffer[l] = 0;

      for (int k = 0; k < fNumNus; k++) {
        if (k <= l)
          Hms_kl = fHms[k][l];
        else
          Hms_kl = conj(fHms[l][k]);

        if (fIsNuBar && k != l) Hms_kl = conj(Hms_kl);

        buffer[l] += conj(fEvec[k][i]) * Hms_kl;
      }
    }

    for (int j = 0; j <= i; j++) {
      fKmass[i][j] = 0;

      for (int l = 0; l < fNumNus; l++) {
        fKmass[i][j] += buffer[l] * fEvec[l][j];
      }

      complexD C;

      if (i == j) { C = complexD(1, 0); }
      else {
        double arg = (fEval[i] - fEval[j]) * lenghtEV;

        C = (complexD(cos(arg), sin(arg)) - complexD(1, 0.0)) /
            complexD(0.0, arg);
      }

      fKmass[i][j] *= bufK * C;

      if (i != j) fKmass[j][i] = conj(fKmass[i][j]);
    }
  }
}


//.............................................................................
///
/// Compute the derivation of one layer's length depending on the angle
///
double PMNS_BaseAvg::LnDerivative()
{
  double dL = 0;

  double L1 = pow(fPrem.GetPremLayers()[fLayer].radius, 2) - fminRsq;

  double L2 = -fminRsq;
  if (fLayer > 0) L2 += pow(fPrem.GetPremLayers()[fLayer - 1].radius, 2);

  bool ismin = (L2 <= 0 && fcosT < 0);

  if (ismin)
    dL = 2 * pow(fDetRadius, 2) * fcosT * pow(L1, -0.5);
  else
    dL = pow(fDetRadius, 2) * fcosT * (pow(L1, -0.5) - pow(L2, -0.5));

  if (ismin) fdl = 1;

  fLayer += fdl;

  return dL;
}

