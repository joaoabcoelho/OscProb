///////////////////////////////////////////////////////////////////////////////
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework with a first order Taylor expansion.
//
// This  class inherits from the PMNS_Fast class
//
///////////////////////////////////////////////////////////////////////////////

#include "PMNS_TaylorExp.h"
#include "MatrixDecomp/zheevh3.h"
#include "PremModel.h"
#include <algorithm>
#include <iostream>

using namespace OscProb;

using namespace std;

//.............................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_TaylorExp::PMNS_TaylorExp() : PMNS_Fast()
{
  fPrem.LoadModel("");

  InitializeTaylorsVectors();

  SetwidthBin(0, 0);

  SetAvgProbPrec(1e-4);
}

//.............................................................................
///
/// Nothing to clean.
///
PMNS_TaylorExp::~PMNS_TaylorExp() {}

//.............................................................................
///
/// Set vector sizes and initialize elements to zero.
/// Initialize diagonal elements of S to one
///
void PMNS_TaylorExp::InitializeTaylorsVectors()
{
  fdensityMatrix = matrixC(fNumNus, vectorC(fNumNus, 0));

  flambdaInvE = vectorD(fNumNus, 0);
  fVInvE      = matrixC(fNumNus, vectorC(fNumNus, 0));

  flambdaCosT = vectorD(fNumNus, 0);
  fVcosT      = matrixC(fNumNus, vectorC(fNumNus, 0));

  fevolutionMatrixS = matrixC(fNumNus, vectorC(fNumNus, 0));

  fSflavor = matrixC(fNumNus, vectorC(fNumNus, 0));
  fKmass   = matrixC(fNumNus, vectorC(fNumNus, 0));
  fKflavor = matrixC(fNumNus, vectorC(fNumNus, 0));

  for (int i = 0; i < fNumNus; i++) {
    fevolutionMatrixS[i][i] = 1;

    for (int j = 0; j < fNumNus; j++) {
      fKInvE[i][j] = 0;
      fKcosT[i][j] = 0;
    }
  }

  fLayer = fPrem.GetPremLayers().size() - 1;

  fdl = -1;

  fDetRadius = fPrem.GetDetRadius();
}

//.............................................................................
///
/// Copy the earth model used
///
/// This is done to get access to the PremLayer list to be used in the
/// LnDerivative() function.
///
/// @param prem - The earth model used
///
void PMNS_TaylorExp::SetPremModel(OscProb::PremModel& prem) { fPrem = prem; }

//.............................................................................
///
/// Set neutrino angle.
///
/// @param cosT - The cosine of the neutrino angle
///
void PMNS_TaylorExp::SetCosT(double cosT) { fcosT = cosT; }

//.............................................................................
///
/// Set bin's widths.
///
/// @param dE - The width of the bin for energy in GeV
/// @param dcosT - The width of the bin for angle
///
void PMNS_TaylorExp::SetwidthBin(double dE, double dcosT)
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
void PMNS_TaylorExp::BuildKE(double L)
{
  double lenghtEV = L * kKm2eV;     // L in eV-1
  double bufK     = lenghtEV * 0.5; // L/2 in eV-1

  complexD buffer[3];

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
/// Build K matrix for angle in flavor basis
///
/// The variable for which a Taylor expansion is done here is not directly the
/// angle but the cosine of the angle
///
/// @param L - The length of the layer in km
///
void PMNS_TaylorExp::BuildKcosT(double L)
{
  UpdateHam();

  double dL = LnDerivative() * kKm2eV;

  for (int j = 0; j < fNumNus; j++) {
    for (int i = 0; i <= j; i++) {
      fKflavor[i][j] = dL * fHam[i][j];

      if (i != j) { fKflavor[j][i] = conj(fKflavor[i][j]); }
    }
  }
}

//.............................................................................
///
/// Compute the derivation of one layer's length depending on the angle
///
double PMNS_TaylorExp::LnDerivative()
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

//.............................................................................
///
/// Rotate the S matrix for the current layer from mass to flavor basis
///
void PMNS_TaylorExp::rotateS()
{
  complexD buffer[3];

  for (int j = 0; j < fNumNus; j++) {
    for (int k = 0; k < fNumNus; k++) {
      buffer[k] = fPhases[k] * conj(fEvec[j][k]);
    }

    for (int i = 0; i < fNumNus; i++) {
      fSflavor[i][j] = 0;
      for (int k = 0; k < fNumNus; k++) {
        fSflavor[i][j] += fEvec[i][k] * buffer[k];
      }
    }
  }
}

//.............................................................................
///
/// Rotate the K matrix from mass to flavor basis
///
void PMNS_TaylorExp::rotateK()
{
  complexD buffer[3];

  for (int j = 0; j < fNumNus; j++) {
    for (int k = 0; k < fNumNus; k++) {
      buffer[k] = 0;

      for (int l = 0; l < fNumNus; l++) {
        buffer[k] += fKmass[k][l] * conj(fEvec[j][l]);
      }
    }

    for (int i = 0; i <= j; i++) {
      fKflavor[i][j] = 0;

      for (int k = 0; k < fNumNus; k++) {
        fKflavor[i][j] += fEvec[i][k] * buffer[k];
      }

      if (i != j) { fKflavor[j][i] = conj(fKflavor[i][j]); }
    }
  }
}

//.............................................................................
///
/// Product between two S matrices.
///
/// This is used to calculate the matrix S corresponding to the propagation
/// between the beginning of the path and the end of the current layer.
///
/// The matrix fevolutionMatrixS represent the propagation between the beginning
/// of the path and the beginning of the current layer. This matrix is updated
/// after every layer with this function. The matrix fSflavor represent the
/// propagation between the beginning and the end of the layer.
///
void PMNS_TaylorExp::MultiplicationRuleS()
{
  complexD save[3];

  for (int j = 0; j < fNumNus; j++) {
    for (int n = 0; n < fNumNus; n++) { save[n] = fevolutionMatrixS[n][j]; }

    for (int i = 0; i < fNumNus; i++) {
      fevolutionMatrixS[i][j] = 0;

      for (int k = 0; k < fNumNus; k++) {
        fevolutionMatrixS[i][j] += fSflavor[i][k] * save[k];
      }
    }
  }
}

//.............................................................................
///
/// Product between two K matrices.
///
/// This is used to calculate the matrix K corresponding to the propagation
/// between the beginning of the path and the end of the current layer.
///
/// The matrix fKflavor correspond to the propagation between the beginning and
/// the end of the layer.
///
/// @param K - The K matrix corresponding to the propagation between the
/// beginning
///            of the path and the beginning of the current layer
///
void PMNS_TaylorExp::MultiplicationRuleK(complexD K[3][3])
{
  for (int i = 0; i < fNumNus; i++) {
    complexD buffer[3];

    for (int l = 0; l < fNumNus; l++) {
      for (int k = 0; k < fNumNus; k++) {
        buffer[l] += conj(fevolutionMatrixS[k][i]) * fKflavor[k][l];
      }
    }

    for (int j = 0; j <= i; j++) {
      for (int l = 0; l < fNumNus; l++) {
        K[i][j] += buffer[l] * fevolutionMatrixS[l][j];
      }

      if (i != j) { K[j][i] = conj(K[i][j]); }
    }
  }
}

//.............................................................................
///
/// Propagate neutrino state through full path
///
void PMNS_TaylorExp::PropagateTaylor()
{
  for (int i = 0; i < int(fNuPaths.size()); i++) {
    PropagatePathTaylor(fNuPaths[i]);
  }
}

//.............................................................................
///
/// Propagate the current neutrino state through a given path
/// @param p - A neutrino path segment
///
void PMNS_TaylorExp::PropagatePathTaylor(NuPath p)
{
  // Set the neutrino path
  SetCurPath(p);

  // Solve for eigensystem
  SolveHam();

  // Get the evolution matrix in mass basis
  double LengthIneV = kKm2eV * p.length;
  for (int i = 0; i < fNumNus; i++) {
    double arg = fEval[i] * LengthIneV;
    fPhases[i] = complexD(cos(arg), -sin(arg));
  }

  // Rotate S in flavor basis
  rotateS();

  // if avg on E
  if (fdInvE != 0) {
    // Build KE in mass basis
    BuildKE(p.length);

    // Rotate KE in flavor basis
    rotateK();

    // Multiply this layer K's with the previous path K's
    MultiplicationRuleK(fKInvE);
  }

  // if avg on cosT
  if (fdcosT != 0) {
    // Build KcosT in mass basis
    BuildKcosT(p.length);

    // Multiply this layer K's with the previous path K's
    MultiplicationRuleK(fKcosT);
  }

  // Multiply this layer S's with the previous path S's
  MultiplicationRuleS();
}

//.............................................................................
///
/// Solve one K matrix for eigenvectors and eigenvalues.
///
/// This is using a method from the GLoBES software available at
/// http://www.mpi-hd.mpg.de/personalhomes/globes/3x3/ \n
/// We should cite them accordingly
///
/// @param K - The K matrix
/// @param lambda - The eigenvalues of K
/// @param V - The eigenvectors of K
///
void PMNS_TaylorExp::SolveK(complexD K[3][3], vectorD& lambda, matrixC& V)
{
  double   fEvalGLoBES[3];
  complexD fEvecGLoBES[3][3];

  // Solve K for eigensystem using the GLoBES method
  zheevh3(K, fEvecGLoBES, fEvalGLoBES);

  // Fill flambdaInvE and fVInvE vectors from GLoBES arrays
  for (int i = 0; i < fNumNus; i++) {
    lambda[i] = fEvalGLoBES[i];
    for (int j = 0; j < fNumNus; j++) { V[i][j] = fEvecGLoBES[i][j]; }
  }
}

//.............................................................................
///
/// Formula for the average probability of flvi going to flvf over
/// a bin of width dbin with a Taylor expansion.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param dbin - The width of the bin
/// @param lambda - The eigenvalues of K
/// @param V - The eigenvectors of K
///
/// @return Average neutrino oscillation probability
///
double PMNS_TaylorExp::AvgFormula(int flvi, int flvf, double dbin,
                                  vectorD lambda, matrixC V)
{
  vectorC SV = vectorC(fNumNus, 0);

  for (int i = 0; i < fNumNus; i++) {
    for (int k = 0; k < fNumNus; k++) {
      SV[i] += fevolutionMatrixS[flvf][k] * V[k][i];
    }
  }

  complexD buffer[3];

  for (int n = 0; n < fNumNus; n++) { buffer[n] = SV[n] * conj(V[flvi][n]); }

  complexD sinc[fNumNus][fNumNus];

  for (int j = 0; j < fNumNus; j++) {
    sinc[j][j] = 1;

    for (int i = 0; i < j; i++) {
      double arg = (lambda[i] - lambda[j]) * dbin / 2;
      sinc[i][j] = sin(arg) / arg;
      sinc[j][i] = sinc[i][j];
    }
  }

  complexD P = 0;

  for (int j = 0; j < fNumNus; j++) {
    for (int i = 0; i < fNumNus; i++) {
      P += buffer[i] * conj(buffer[j]) * sinc[j][i];
    }
  }

  return real(P);
}

//.............................................................................
///
/// Compute the average probability of flvi going to flvf over
/// a bin of energy E with width dE with a Taylor expansion.
///
/// This gets transformed into L/E, since the oscillation terms
/// have arguments linear in 1/E and not E.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param E - The neutrino energy in GeV
/// @param dE - The energy bin width in GeV
///
/// @return Average neutrino oscillation probability
///
double PMNS_TaylorExp::AvgProb(int flvi, int flvf, double E, double dE)
{
  if (E <= 0) return 0;

  if (fNuPaths.empty()) return 0;

  // Don't average zero width
  if (dE <= 0) return Prob(flvi, flvf, E);

  vectorD Ebin = ConvertEtoLoE(E, dE);

  // return fct avr proba
  return AvgProbLoE(flvi, flvf, Ebin[0], Ebin[1]);
}

//.............................................................................
///
/// Compute the average probability of flvi going to flvf over
/// a bin of energy L/E with width dLoE with a Taylor expansion.
///
/// The probabilities are weighted by (L/E)^-2 so that event
/// density is flat in energy. This avoids giving too much
/// weight to low energies. Better approximations would be
/// achieved if we used an interpolated event density.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param LoE - The neutrino  L/E value in the bin center in km/GeV
/// @param dLoE - The L/E bin width in km/GeV
///
/// @return Average neutrino oscillation probability
///
double PMNS_TaylorExp::AvgProbLoE(int flvi, int flvf, double LoE, double dLoE)
{
  if (LoE <= 0) return 0;

  if (fNuPaths.empty()) return 0;

  // Don't average zero width
  if (dLoE <= 0) return Prob(flvi, flvf, fPath.length / LoE);

  // Get sample points for this bin
  vectorD samples = GetSamplePoints(LoE, dLoE);

  double avgprob = 0;
  double L       = fPath.length;
  double sumw    = 0;

  // Loop over all sample points
  for (int j = 1; j < int(samples.size()); j++) {
    double w = 1. / pow(samples[j], 2);

    avgprob += w * AvgAlgo(flvi, flvf, samples[j], samples[0], L);

    sumw += w;
  }

  // Return average of probabilities
  return avgprob / sumw;
}

//.............................................................................
///
/// Algorithm for the compute of the average probability over a bin of LoE
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param LoE - The neutrino  L/E value in the bin center in km/GeV
/// @param dLoE - The L/E bin width in km/GeV
/// @param L - The length of the path in km
///
/// @return Average neutrino oscillation probability
///
double PMNS_TaylorExp::AvgAlgo(int flvi, int flvf, double LoE, double dLoE,
                               double L)
{
  // Set width bin as 1/E
  double d1oE = dLoE / L;

  // reset K et S et Ve et lambdaE
  InitializeTaylorsVectors();

  SetEnergy(L / LoE);
  SetwidthBin(d1oE, 0);

  // Propagate -> get S and K matrix (on the whole path)
  PropagateTaylor();

  // DiagolK -> get VE and lambdaE
  SolveK(fKInvE, flambdaInvE, fVInvE);

  // return fct avr proba
  return AvgFormula(flvi, flvf, d1oE / kGeV2eV, flambdaInvE, fVInvE);
}

//.............................................................................
///
/// Compute the average probability of flvi going to flvf over
/// a bin of angle cost with width dcosT with a Taylor expansion.
///
/// IMPORTANT: The PremModel object used must be copied by this
/// class using the function SetPremModel.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param E - The neutrino energy in GeV
/// @param cosT - The cosine of the neutrino angle
/// @param dcosT - The cosT bin width
///
/// @return Average neutrino oscillation probability
///
double PMNS_TaylorExp::AvgProb(int flvi, int flvf, double E, double cosT,
                               double dcosT)
{
  if (cosT > 0) return 0;

  // if (fNuPaths.empty()) return 0; //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  // Don't average zero width
  if (dcosT == 0) return Prob(flvi, flvf, E);

  vectorD samples = GetSamplePoints(E, cosT, dcosT);

  double avgprob = 0;

  // Loop over all sample points
  for (int j = 1; j < int(samples.size()); j++) {
    avgprob += AvgAlgoCosT(flvi, flvf, E, samples[j], samples[0]);
  }

  // Return average of probabilities
  return avgprob / (samples.size() - 1);
}

//.............................................................................
///
/// Algorithm for the compute of the average probability over a bin of cosT
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param E - The neutrino energy in GeV
/// @param cosT - The cosine of the neutrino angle
/// @param dcosT - The cosT bin width
///
/// @return Average neutrino oscillation probability
///
double PMNS_TaylorExp::AvgAlgoCosT(int flvi, int flvf, double E, double cosT,
                                   double dcosT)
{
  // reset K et S et Ve et lambdaE
  InitializeTaylorsVectors();

  SetEnergy(E);
  SetCosT(cosT);
  SetwidthBin(0, dcosT);

  fPrem.FillPath(cosT);
  SetPath(fPrem.GetNuPath());

  fminRsq = pow(fDetRadius * sqrt(1 - cosT * cosT), 2);

  // Propagate -> get S and K matrix (on the whole path)
  PropagateTaylor();

  // DiagolK -> get VE and lambdaE
  SolveK(fKcosT, flambdaCosT, fVcosT);

  // return fct avr proba
  return AvgFormula(flvi, flvf, fdcosT, flambdaCosT, fVcosT);
}

//.............................................................................
///
/// Compute the average probability of flvi going to flvf over
/// a bin of energy E and angle cosT with width dE and dcosT
/// with a Taylor expansion.
///
/// This gets transformed into L/E, since the oscillation terms
/// have arguments linear in 1/E and not E.
///
/// IMPORTANT: The function SetPremLayers must be used in the
/// macro file to make this function work. The argument for
/// SetPremLayers must be premModel.GetPremLayers().
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param E - The neutrino energy in GeV
/// @param dE - The energy bin width in GeV
/// @param cosT - The cosine of the neutrino angle
/// @param dcosT - The cosT bin width
///
/// @return Average neutrino oscillation probability
///
double PMNS_TaylorExp::AvgProb(int flvi, int flvf, double E, double dE,
                               double cosT, double dcosT)
{
  if (E <= 0) return 0;
  if (cosT > 0) return 0;

  if (fNuPaths.empty()) return 0;

  // Don't average zero width
  if (dE <= 0 && dcosT == 0) return Prob(flvi, flvf, E);
  if (dE <= 0) return AvgProb(flvi, flvf, E, cosT, dcosT);
  if (dcosT == 0) return AvgProb(flvi, flvf, E, dE);

  vectorD Ebin = ConvertEtoLoE(E, dE);

  return AvgProbLoE(flvi, flvf, Ebin[0], Ebin[1], cosT, dcosT);
}

//.............................................................................
///
/// Compute the average probability of flvi going to flvf over
/// a bin of energy L/E and cosT with width dLoE and dcosT
/// with a Taylor expansion.
///
/// The probabilities are weighted by (L/E)^-2 so that event
/// density is flat in energy. This avoids giving too much
/// weight to low energies. Better approximations would be
/// achieved if we used an interpolated event density.
///
/// IMPORTANT: The PremModel object used must be copied by this
/// class using the function SetPremModel.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param LoE - The neutrino  L/E value in the bin center in km/GeV
/// @param dLoE - The L/E bin width in km/GeV
/// @param cosT - The cosine of the neutrino angle
/// @param dcosT - The cosT bin width
///
/// @return Average neutrino oscillation probability
///
double PMNS_TaylorExp::AvgProbLoE(int flvi, int flvf, double LoE, double dLoE,
                                  double cosT, double dcosT)
{
  if (LoE <= 0) return 0;
  if (cosT > 0) return 0;

  // if (fNuPaths.empty()) return 0;

  // Don't average zero width
  if (dLoE <= 0 && dcosT == 0) return Prob(flvi, flvf, fPath.length / LoE);
  if (dLoE <= 0)
    return AvgProb(flvi, flvf, fPath.length / LoE, cosT, dcosT); /// chg ici
  if (dcosT == 0) return AvgProbLoE(flvi, flvf, LoE, dLoE);

  // Make sample with 1oE and not LoE
  matrixC samples =
      GetSamplePoints(LoE / fPath.length, dLoE / fPath.length, cosT, dcosT);

  int rows = samples.size();
  int cols = samples[0].size();

  double avgprob = 0;
  double sumw    = 0;

  // Loop over all sample points
  for (int k = 1; k < int(rows); k++) {
    for (int l = 1; l < int(cols); l++) {
      double w = 1. / pow(real(samples[k][l]), 2);

      avgprob +=
          w * AvgAlgo(flvi, flvf, real(samples[k][l]), real(samples[0][0]),
                      imag(samples[k][l]), imag(samples[0][0]));

      sumw += w;
    }
  }

  // Return average of probabilities
  return avgprob / ((sumw));
}

//.............................................................................
///
/// Algorithm for the compute of the average probability
/// over a bin of 1oE and cosT
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param InvE - The neutrino  1/E value in the bin center in GeV-1
/// @param dInvE - The 1/E bin width in GeV-1
/// @param cosT - The cosine of the neutrino angle
/// @param dcosT - The cosT bin width
///
/// @return Average neutrino oscillation probability
///
double PMNS_TaylorExp::AvgAlgo(int flvi, int flvf, double InvE, double dInvE,
                               double cosT, double dcosT)
{
  fPrem.FillPath(cosT);
  SetPath(fPrem.GetNuPath());

  // reset K et S et Ve et lambdaE
  InitializeTaylorsVectors();

  SetEnergy(1 / InvE);
  SetCosT(cosT);
  SetwidthBin(dInvE, dcosT);

  fminRsq = pow(fDetRadius * sqrt(1 - cosT * cosT), 2);

  // Propagate -> get S and K matrix (on the whole path)
  PropagateTaylor();

  // DiagolK -> get VE and lambdaE
  SolveK(fKInvE, flambdaInvE, fVInvE);
  SolveK(fKcosT, flambdaCosT, fVcosT);

  return AlgorithmDensityMatrix(flvi, flvf);
}

//.............................................................................
///
/// Algorithm for the transformations on the density matrix
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
///
/// @return Average neutrino oscillation probability
///
double PMNS_TaylorExp::AlgorithmDensityMatrix(int flvi, int flvf)
{
  fdensityMatrix[flvi][flvi] = 1;

  RotateDensityM(true, fVcosT);
  HadamardProduct(flambdaCosT, fdcosT);
  RotateDensityM(false, fVcosT);

  RotateDensityM(true, fVInvE);
  HadamardProduct(flambdaInvE, fdInvE / kGeV2eV);
  RotateDensityM(false, fVInvE);

  RotateDensityM(false, fevolutionMatrixS);

  return real(fdensityMatrix[flvf][flvf]);
}

//.............................................................................
///
/// Apply rotation to the density matrix from or to the basis
/// where V is diagonal
///
/// @param to_basis - Rotation from (false) or to (true)
/// @param V - The matrix used for the denisty matrix rotation.
///
void PMNS_TaylorExp::RotateDensityM(bool to_basis, matrixC V)
{
  matrixC Buffer = matrixC(fNumNus, vectorC(fNumNus, 0));

  for (int i = 0; i < fNumNus; i++) {
    for (int j = 0; j < fNumNus; j++) {
      for (int k = 0; k < fNumNus; k++) {
        if (to_basis)
          Buffer[i][j] += fdensityMatrix[i][k] * V[k][j];
        else
          Buffer[i][j] += fdensityMatrix[i][k] * conj(V[j][k]);
      }
    }
  }

  for (int i = 0; i < fNumNus; i++) {
    for (int j = i; j < fNumNus; j++) {
      fdensityMatrix[i][j] = 0;
      for (int k = 0; k < fNumNus; k++) {
        if (to_basis)
          fdensityMatrix[i][j] += conj(V[k][i]) * Buffer[k][j];
        else
          fdensityMatrix[i][j] += V[i][k] * Buffer[k][j];
      }
      if (j > i) fdensityMatrix[j][i] = conj(fdensityMatrix[i][j]);
    }
  }
}

//.............................................................................
///
/// Apply an Hadamard Product to the density matrix
///
/// @param lambda - Eigenvalues of the K matrix
/// @param dbin - Width of the bin
///
void PMNS_TaylorExp::HadamardProduct(vectorD lambda, double dbin)
{
  matrixC sinc = matrixC(fNumNus, vectorC(fNumNus, 0));
  for (int j = 0; j < fNumNus; j++) {
    for (int i = 0; i < j; i++) {
      double arg = (lambda[i] - lambda[j]) * dbin;
      sinc[i][j] = sin(arg) / arg;

      sinc[j][i] = sinc[i][j];
    }
  }

  for (int i = 0; i < fNumNus; i++) { sinc[i][i] = 1; }

  for (int j = 0; j < fNumNus; j++) {
    for (int i = 0; i < fNumNus; i++) {
      fdensityMatrix[i][j] = fdensityMatrix[i][j] * sinc[i][j];
    }
  }
}

//.............................................................................
///
/// Compute the sample points for a bin of L/E with width dLoE
///
/// This is used to increase the average probability over a bin of L/E,
/// calculated with a Taylor expansion
///
/// @param LoE  - The neutrino  L/E value in the bin center in km/GeV
/// @param dLoE   - The L/E bin width in km/GeV
///
vectorD PMNS_TaylorExp::GetSamplePoints(double LoE, double dLoE)
{
  // Set a number of sub-divisions to achieve "good" accuracy
  // This needs to be studied better
  int n_div = ceil(3 * pow(dLoE / LoE, 0.8) * pow(LoE, 0.3) /
                   sqrt(fAvgProbPrec / 1e-4));

  // A vector to store sample points
  vectorD Samples;

  // Define sub-division width
  Samples.push_back(dLoE / n_div);

  // Loop over sub-divisions
  for (int k = 0; k < n_div; k++) {
    // Define sub-division center
    double bctr = LoE - dLoE / 2 + (k + 0.5) * dLoE / n_div;

    Samples.push_back(bctr);

  } // End of loop over sub-divisions

  // Return sample points
  return Samples;
}

//.............................................................................
///
/// Compute the sample points for a bin of cosTheta with width dcosTheta
///
/// This is used to increase the average probability over a bin of cosT,
/// calculated with a Taylor expansion
///
/// @param E  - The neutrino  Energy value GeV
/// @param cosT   - The neutrino  cosT value in the bin center
/// @param dcosT   - The cosT bin width
///
vectorD PMNS_TaylorExp::GetSamplePoints(double E, double cosT, double dcosT)
{
  // Set a number of sub-divisions to achieve "good" accuracy
  // This needs to be studied better
  int n_div = ceil(35 * pow(E, -0.4) * pow(abs(dcosT / cosT), 0.8) /
                   sqrt(fAvgProbPrec / 1e-4));

  // A vector to store sample points
  vectorD Samples;

  // Define sub-division width
  Samples.push_back(dcosT / n_div);

  // Loop over sub-divisions
  for (int k = 0; k < n_div; k++) {
    // Define sub-division center
    double bctr = cosT - dcosT / 2 + (k + 0.5) * dcosT / n_div;

    Samples.push_back(bctr);

  } // End of loop over sub-divisions

  // Return sample points
  return Samples;
}

//.............................................................................
///
/// Compute the sample points for a bin of 1oE and cosTheta
/// with width d1oE and dcosTheta
///
/// This is used to increase the average probability over a bin of L/E
/// and cosT, calculated with a Taylor expansion
///
/// @param InvE  - The neutrino  1/E value in the bin center in GeV-1
/// @param dInvE   - The 1/E bin width in GeV-1
/// @param cosT   - The neutrino  cosT value in the bin center
/// @param dcosT   - The cosT bin width
///
matrixC PMNS_TaylorExp::GetSamplePoints(double InvE, double dInvE, double cosT,
                                        double dcosT)
{
  // Set a number of sub-divisions to achieve "good" accuracy
  // This needs to be studied better
  int n_divCosT = ceil(380 * pow(InvE, 0.5) * pow(abs(dcosT / cosT), 0.8) /
                       sqrt(fAvgProbPrec / 1e-4));
  int n_divE    = ceil(260 * pow(dInvE / InvE, 0.8) * pow(InvE, 0.6) /
                       sqrt(fAvgProbPrec / 1e-4));

  // A matrix to store sample points
  matrixC Samples = matrixC(n_divE + 1, vectorC(n_divCosT + 1, 0));

  // Define sub-division width
  Samples[0][0] = complexD(dInvE / n_divE, dcosT / n_divCosT);

  // Loop over sub-divisions
  for (int k = 0; k < n_divE; k++) {
    // Define sub-division center for energy
    double bctr_InvE = InvE - dInvE / 2 + (k + 0.5) * dInvE / n_divE;

    for (int l = 0; l < n_divCosT; l++) {
      // Define sub-division center for angle
      double bctr_CosT = cosT - dcosT / 2 + (l + 0.5) * dcosT / n_divCosT;

      Samples[k + 1][l + 1] = complexD(bctr_InvE, bctr_CosT);
    }

  } // End of loop over sub-divisions

  // Return sample points
  return Samples;
}

//.............................................................................
///
/// Fomula for the propability for flvi going to flvf for an energy E+dE
/// using a first order Taylor expansion from a reference energy E.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param dE - The energy variation in GeV
/// @param lambda - The eigenvalues of K
/// @param V - The eigenvectors of K
///
/// @return Neutrino oscillation probability
///
double PMNS_TaylorExp::AvgFormulaExtrapolation(int flvi, int flvf, double dbin,
                                               vectorD lambda, matrixC V)
{
  vectorC SV = vectorC(fNumNus, 0);

  for (int i = 0; i < fNumNus; i++) {
    for (int k = 0; k < fNumNus; k++) {
      SV[i] += fevolutionMatrixS[flvf][k] * V[k][i];
    }
  }

  complexD buffer[3];

  for (int n = 0; n < fNumNus; n++) { buffer[n] = SV[n] * conj(V[flvi][n]); }

  complexD expo[fNumNus][fNumNus];

  for (int j = 0; j < fNumNus; j++) {
    for (int i = 0; i < fNumNus; i++) {
      double arg = (lambda[j] - lambda[i]) * dbin;
      expo[j][i] = exp(complexD(0.0, arg));
    }
  }

  complexD P = 0;

  for (int j = 0; j < fNumNus; j++) {
    for (int i = 0; i < fNumNus; i++) {
      P += buffer[i] * conj(buffer[j]) * expo[j][i];
    }
  }

  return real(P);
}

//.............................................................................
///
/// Compute the propability for flvi going to flvf for an energy E+dE
/// using a first order Taylor expansion from a reference energy E.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param E - The reference energy in GeV
/// @param dE - The energy variation in GeV
///
/// @return Neutrino oscillation probability
///
double PMNS_TaylorExp::ExtrapolationProb(int flvi, int flvf, double E,
                                         double dE)
{
  // reset K et S et Ve et lambdaE
  InitializeTaylorsVectors();

  SetEnergy(E);
  SetwidthBin(1 / (E + dE) - 1 / E, 0);

  // Propagate -> get S and K matrix (on the whole path)
  PropagateTaylor();

  // DiagolK -> get VE and lambdaE
  SolveK(fKInvE, flambdaInvE, fVInvE);

  return AvgFormulaExtrapolation(flvi, flvf, fdInvE / kGeV2eV, flambdaInvE,
                                 fVInvE);
}

//.............................................................................
///
/// Compute the propability for flvi going to flvf for an energy LoE+dLoE
/// using a first order Taylor expansion from a reference value LoE.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param LoE - The reference energy in GeV
/// @param dLoE - The energy variation in GeV
///
/// @return Neutrino oscillation probability
///
double PMNS_TaylorExp::ExtrapolationProbLoE(int flvi, int flvf, double LoE,
                                            double dLoE)
{
  // reset K et S et Ve et lambdaE
  InitializeTaylorsVectors();

  SetCurPath(AvgPath(fNuPaths));
  double L = fPath.length;

  SetEnergy(L / LoE);
  SetwidthBin(dLoE / L, 0);

  // Propagate -> get S and K matrix (on the whole path)
  PropagateTaylor();

  // DiagolK -> get VE and lambdaE
  SolveK(fKInvE, flambdaInvE, fVInvE);

  return AvgFormulaExtrapolation(flvi, flvf, dLoE * kGeV2eV / L, flambdaInvE,
                                 fVInvE);
}

//.............................................................................
///
/// Compute the propability for flvi going to flvf for an angle cosT+dcosT
/// using a first order Taylor expansion from a reference angle cosT.
///
/// IMPORTANT: The PremModel object used must be copied by this
/// class using the function SetPremModel.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
///
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param cosT - The reference angle
/// @param dcosT - The angle variation
///
/// @return Neutrino oscillation probability
///
double PMNS_TaylorExp::ExtrapolationProbCosT(int flvi, int flvf, double cosT,
                                             double dcosT)
{
  fPrem.FillPath(cosT);
  SetPath(fPrem.GetNuPath());

  // reset K et S et Ve et lambdaE
  InitializeTaylorsVectors();

  // SetEnergy(E);
  SetCosT(cosT);
  SetwidthBin(0, dcosT);

  fminRsq = pow(fDetRadius * sqrt(1 - cosT * cosT), 2);

  // Propagate -> get S and K matrix (on the whole path)
  PropagateTaylor();

  // DiagolK -> get VE and lambdaE
  SolveK(fKcosT, flambdaCosT, fVcosT);

  return AvgFormulaExtrapolation(flvi, flvf, fdcosT, flambdaCosT, fVcosT);
}