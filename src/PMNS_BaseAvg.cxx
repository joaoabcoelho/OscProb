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
/// Build K matrix for angle in flavor basis
///
/// The variable for which a Taylor expansion is done here is not directly the
/// angle but the cosine of the angle
///
/// @param L - The length of the layer in km
///
void PMNS_BaseAvg::BuildKcosT(double L)
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


//.............................................................................
///
/// Compute the derivation of one layer's length depending on the angle

//.............................................................................
///
/// Rotate the S matrix for the current layer from mass to flavor basis
///
void PMNS_BaseAvg::rotateS()
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
void PMNS_BaseAvg::rotateK()
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
void PMNS_BaseAvg::MultiplicationRuleS()
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
void PMNS_BaseAvg::MultiplicationRuleK(Eigen::MatrixXcd& K)
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
void PMNS_BaseAvg::PropagateTaylor()
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
void PMNS_BaseAvg::PropagatePathTaylor(NuPath p)
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
void PMNS_BaseAvg::SolveK(Eigen::MatrixXcd& K, vectorD& lambda, matrixC& V)
{
  if (fNumNus == 4) TemplateSolver<Eigen::Matrix4cd>(K, lambda, V);
  if (fNumNus == 3) TemplateSolver<Eigen::Matrix3cd>(K, lambda, V);
  if (fNumNus == 2)
    TemplateSolver<Eigen::Matrix2cd>(K, lambda, V);
  else
    TemplateSolver<Eigen::MatrixXcd>(K, lambda, V);
}

template <typename T>
void PMNS_BaseAvg::TemplateSolver(Eigen::MatrixXcd& K, vectorD& lambda, matrixC& V)
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
double PMNS_BaseAvg::AvgFormula(int flvi, int flvf, double dbin, vectorD lambda,
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
double PMNS_BaseAvg::AvgProb(int flvi, int flvf, double E, double dE)
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
double PMNS_BaseAvg::AvgProbLoE(int flvi, int flvf, double LoE, double dLoE)
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
double PMNS_BaseAvg::AvgAlgo(int flvi, int flvf, double LoE, double dLoE, double L)
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
/// Compute the sample points for a bin of L/E with width dLoE
///
/// This is used to increase the average probability over a bin of L/E,
/// calculated with a Taylor expansion
///
/// @param LoE  - The neutrino  L/E value in the bin center in km/GeV
/// @param dLoE   - The L/E bin width in km/GeV
///
vectorD PMNS_BaseAvg::GetSamplePointsAvgClass(double LoE, double dLoE)
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
///
///
vectorD PMNS_BaseAvg::AvgProbVector(int flvi, double E, double dE)
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
vectorD PMNS_BaseAvg::AvgProbVectorLoE(int flvi, double LoE, double dLoE)
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
matrixD PMNS_BaseAvg::AvgProbMatrix(int nflvi, int nflvf, double E, double dE)
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
matrixD PMNS_BaseAvg::AvgProbMatrixLoE(int nflvi, int nflvf, double LoE,
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


