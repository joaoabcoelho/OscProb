///////////////////////////////////////////////////////////////////////////////
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework with a first order Taylor expansion.
//
// This  class inherits from the PMNS_Fast class
//
///////////////////////////////////////////////////////////////////////////////

#include "PMNS_Avg.h"
#include <algorithm>
#include <iostream>
#include <Eigen/Eigenvalues>

using namespace OscProb;

using namespace std;

//.............................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
///
///
/// This bit would copy over to PMNS_MyClass, as the myOsc = ROOT.OscProb.PMNS_MyClass(4)
PMNS_Avg::PMNS_Avg(int numNus) : PMNS_Sterile(numNus)
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
PMNS_Avg::~PMNS_Avg() {}

//.............................................................................
///
/// Build K matrix for angle in flavor basis
///
/// The variable for which a Taylor expansion is done here is not directly the
/// angle but the cosine of the angle
///
/// @param L - The length of the layer in km
///
void PMNS_Avg::BuildKcosT(double L)
{
  UpdateHam();

  double dL = LnDerivative() * kKm2eV;

  for (int j = 0; j < fNumNus; j++) {
    for (int i = 0; i <= j; i++) {
      fKflavor[i][j] = dL * fHam(i, j);

      if (i != j) { fKflavor[j][i] = conj(fKflavor[i][j]); }
    }
  }
}

//.............................................................................
///
/// Compute the derivation of one layer's length depending on the angle
///
double PMNS_Avg::LnDerivative()
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
void PMNS_Avg::rotateS()
{
  complexD buffer[fNumNus];

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
void PMNS_Avg::rotateK()
{
  complexD buffer[fNumNus];

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
void PMNS_Avg::MultiplicationRuleS()
{
  complexD save[fNumNus];

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
void PMNS_Avg::MultiplicationRuleK(Eigen::MatrixXcd& K)
{
  for (int i = 0; i < fNumNus; i++) {
    complexD buffer[fNumNus];

    for (int l = 0; l < fNumNus; l++) {
      for (int k = 0; k < fNumNus; k++) {
        buffer[l] += conj(fevolutionMatrixS[k][i]) * fKflavor[k][l];
      }
    }

    for (int j = 0; j <= i; j++) {
      for (int l = 0; l < fNumNus; l++) {
        K(i, j) += buffer[l] * fevolutionMatrixS[l][j];
      }

      if (i != j) { K(j, i) = conj(K(i, j)); }
    }
  }
}

//.............................................................................
///
/// Propagate neutrino state through full path
///
void PMNS_Avg::PropagateTaylor()
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
void PMNS_Avg::PropagatePathTaylor(NuPath p)
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
void PMNS_Avg::SolveK(Eigen::MatrixXcd& K, vectorD& lambda, matrixC& V)
{
  if (fNumNus == 4) TemplateSolver<Eigen::Matrix4cd>(K, lambda, V);
  if (fNumNus == 3) TemplateSolver<Eigen::Matrix3cd>(K, lambda, V);
  if (fNumNus == 2)
    TemplateSolver<Eigen::Matrix2cd>(K, lambda, V);
  else
    TemplateSolver<Eigen::MatrixXcd>(K, lambda, V);
}

template <typename T>
void PMNS_Avg::TemplateSolver(Eigen::MatrixXcd& K, vectorD& lambda, matrixC& V)
{
  Eigen::Ref<T> R(K);

  Eigen::SelfAdjointEigenSolver<T> eigensolver(R);

  // Fill flambdaInvE and fVInvE vectors from GLoBES arrays
  for (int i = 0; i < fNumNus; i++) {
    lambda[i] = eigensolver.eigenvalues()(i);
    for (int j = 0; j < fNumNus; j++) {
      V[i][j] = eigensolver.eigenvectors()(i, j);
    }
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
double PMNS_Avg::AvgFormula(int flvi, int flvf, double dbin, vectorD lambda,
                            matrixC V)
{
  vectorC SV = vectorC(fNumNus, 0);

  for (int i = 0; i < fNumNus; i++) {
    for (int k = 0; k < fNumNus; k++) {
      SV[i] += fevolutionMatrixS[flvf][k] * V[k][i];
    }
  }

  complexD buffer[fNumNus];

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
double PMNS_Avg::AvgProb(int flvi, int flvf, double E, double dE)
{
  // Do nothing if energy is not positive
  if (E <= 0) return 0;

  if (fNuPaths.empty()) return 0;

  // Don't average zero width
  if (dE <= 0) return Prob(flvi, flvf, E);

  vectorD LoEbin = ConvertEtoLoE(E, dE);

  // Compute average in LoE
  return AvgProbLoE(flvi, flvf, LoEbin[0], LoEbin[1]);
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
double PMNS_Avg::AvgProbLoE(int flvi, int flvf, double LoE, double dLoE)
{
  if (LoE <= 0) return 0;

  if (fNuPaths.empty()) return 0;

  // Don't average zero width
  if (dLoE <= 0) return Prob(flvi, flvf, fPath.length / LoE);

  // Get sample points for this bin
  vectorD samples = GetSamplePointsAvgClass(LoE, dLoE);

  double avgprob = 0;
  double L       = fPath.length;
  double sumw    = 0;

  // Loop over all sample points
  for (int j = 1; j < int(samples.size()); j++) {
    // Set (L/E)^-2 weights
    double w = 1. / pow(samples[j], 2);

    avgprob += w * AvgAlgo(flvi, flvf, samples[j], samples[0], L);

    // Increment sum of weights
    sumw += w;
  }

  // Return weighted average of probabilities
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
double PMNS_Avg::AvgAlgo(int flvi, int flvf, double LoE, double dLoE, double L)
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
///
///
vectorD PMNS_Avg::AvgProbVector(int flvi, double E, double dE)
{
  vectorD probs(fNumNus, 0);

  // Do nothing if energy is not positive
  if (E <= 0) return probs;

  if (fNuPaths.empty()) return probs;

  // Don't average zero width
  if (dE <= 0) return ProbVector(flvi, E);

  vectorD LoEbin = ConvertEtoLoE(E, dE);

  // Compute average in LoE
  return AvgProbVectorLoE(flvi, LoEbin[0], LoEbin[1]);
}

//.............................................................................
///
///
///
vectorD PMNS_Avg::AvgProbVectorLoE(int flvi, double LoE, double dLoE)
{
  vectorD probs(fNumNus, 0);

  if (LoE <= 0) return probs;

  if (fNuPaths.empty()) return probs;

  double L = fPath.length;

  // Don't average zero width
  if (dLoE <= 0) return ProbVector(flvi, L / LoE);

  // Get sample points for this bin
  vectorD samples = GetSamplePointsAvgClass(LoE, dLoE);

  double sumw = 0;

  // Loop over all sample points
  for (int j = 1; j < int(samples.size()); j++) {
    // Set (L/E)^-2 weights
    double w = 1. / pow(samples[j], 2);

    for (int flvf = 0; flvf < fNumNus; flvf++) {
      if (flvf == 0)
        probs[flvf] += w * AvgAlgo(flvi, flvf, samples[j], samples[0], L);
      else
        probs[flvf] +=
            w * AvgFormula(flvi, flvf, fdInvE / kGeV2eV, flambdaInvE, fVInvE);
    }

    // Increment sum of weights
    sumw += w;
  }

  for (int flvf = 0; flvf < fNumNus; flvf++) {
    // Divide by total sampling weight
    probs[flvf] /= sumw;
  }

  // Return weighted average of probabilities
  return probs;
}

//.............................................................................
///
///
///
matrixD PMNS_Avg::AvgProbMatrix(int nflvi, int nflvf, double E, double dE)
{
  matrixD probs(nflvi, vectorD(nflvf, 0));

  // Do nothing if energy is not positive
  if (E <= 0) return probs;

  if (fNuPaths.empty()) return probs;

  // Don't average zero width
  if (dE <= 0) return ProbMatrix(nflvi, nflvf, E);

  vectorD LoEbin = ConvertEtoLoE(E, dE);

  // Compute average in LoE
  return AvgProbMatrixLoE(nflvi, nflvf, LoEbin[0], LoEbin[1]);
}

//.............................................................................
///
///
///
matrixD PMNS_Avg::AvgProbMatrixLoE(int nflvi, int nflvf, double LoE,
                                   double dLoE)
{
  matrixD probs(nflvi, vectorD(nflvf, 0));

  if (LoE <= 0) return probs;

  if (fNuPaths.empty()) return probs;

  double L = fPath.length;

  // Don't average zero width
  if (dLoE <= 0) return ProbMatrix(nflvi, nflvf, L / LoE);

  // Get sample points for this bin
  vectorD samples = GetSamplePointsAvgClass(LoE, dLoE);

  double sumw = 0;

  // Loop over all sample points
  for (int j = 1; j < int(samples.size()); j++) {
    // Set (L/E)^-2 weights
    double w = 1. / pow(samples[j], 2);

    for (int flvi = 0; flvi < nflvi; flvi++) {
      for (int flvf = 0; flvf < nflvf; flvf++) {
        // if (!TryCacheK())
        if (flvi == 0 && flvf == 0)
          probs[flvi][flvf] +=
              w * AvgAlgo(flvi, flvf, samples[j], samples[0], L);
        else
          probs[flvi][flvf] +=
              w * AvgFormula(flvi, flvf, fdInvE / kGeV2eV, flambdaInvE, fVInvE);
      }
    }

    // Increment sum of weights
    sumw += w;
  }

  for (int flvi = 0; flvi < nflvi; flvi++) {
    for (int flvf = 0; flvf < nflvf; flvf++) {
      // Divide by total sampling weight
      probs[flvi][flvf] /= sumw;
    }
  }

  // Return weighted average of probabilities
  return probs;
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
double PMNS_Avg::AvgProb(int flvi, int flvf, double E, double cosT,
                         double dcosT)
{
  if (cosT > 0) return 0;

  // if (fNuPaths.empty()) return 0; //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  // Don't average zero width
  if (dcosT == 0) return Prob(flvi, flvf, E);

  vectorD samples = GetSamplePointsAvgClass(E, cosT, dcosT);

  double avgprob = 0;

  // Loop over all sample points
  for (int j = 1; j < int(samples.size()); j++) {
    avgprob += AvgAlgoCosT(flvi, flvf, E, samples[j], samples[0]);
  }

  // Compute average
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
double PMNS_Avg::AvgAlgoCosT(int flvi, int flvf, double E, double cosT,
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

  // Compute average
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
double PMNS_Avg::AvgProb(int flvi, int flvf, double E, double dE, double cosT,
                         double dcosT)
{
  if (E <= 0) return 0;
  if (cosT > 0) return 0;

  if (fNuPaths.empty()) return 0;

  // Don't average zero width
  if (dE <= 0 && dcosT == 0) return Prob(flvi, flvf, E);
  if (dE <= 0) return AvgProb(flvi, flvf, E, cosT, dcosT);
  if (dcosT == 0) return AvgProb(flvi, flvf, E, dE);

  vectorD Ebin = ConvertEtoLoE(E, dE);

  // Compute average in LoE
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
double PMNS_Avg::AvgProbLoE(int flvi, int flvf, double LoE, double dLoE,
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
  matrixC samples = GetSamplePointsAvgClass(LoE / fPath.length,
                                            dLoE / fPath.length, cosT, dcosT);

  int rows = samples.size();
  int cols = samples[0].size();

  double avgprob = 0;
  double sumw    = 0;

  // Loop over all sample points
  for (int k = 1; k < int(rows); k++) {
    for (int l = 1; l < int(cols); l++) {
      // Set (L/E)^-2 weights
      double w = 1. / pow(real(samples[k][l]), 2);

      avgprob +=
          w * AvgAlgo(flvi, flvf, real(samples[k][l]), real(samples[0][0]),
                      imag(samples[k][l]), imag(samples[0][0]));

      // Increment sum of weights
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
double PMNS_Avg::AvgAlgo(int flvi, int flvf, double InvE, double dInvE,
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
double PMNS_Avg::AlgorithmDensityMatrix(int flvi, int flvf)
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
void PMNS_Avg::RotateDensityM(bool to_basis, matrixC V)
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
void PMNS_Avg::HadamardProduct(vectorD lambda, double dbin)
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
vectorD PMNS_Avg::GetSamplePointsAvgClass(double LoE, double dLoE)
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
vectorD PMNS_Avg::GetSamplePointsAvgClass(double E, double cosT, double dcosT)
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
matrixC PMNS_Avg::GetSamplePointsAvgClass(double InvE, double dInvE,
                                          double cosT, double dcosT)
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
double PMNS_Avg::AvgFormulaExtrapolation(int flvi, int flvf, double dbin,
                                         vectorD lambda, matrixC V)
{
  vectorC SV = vectorC(fNumNus, 0);

  for (int i = 0; i < fNumNus; i++) {
    for (int k = 0; k < fNumNus; k++) {
      SV[i] += fevolutionMatrixS[flvf][k] * V[k][i];
    }
  }

  complexD buffer[fNumNus];

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
double PMNS_Avg::ExtrapolationProb(int flvi, int flvf, double E, double dE)
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
double PMNS_Avg::ExtrapolationProbLoE(int flvi, int flvf, double LoE,
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
double PMNS_Avg::ExtrapolationProbCosT(int flvi, int flvf, double cosT,
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
