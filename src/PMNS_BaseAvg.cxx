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

