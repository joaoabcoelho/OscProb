///////////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_Avg
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework with a first order Taylor expansion.
///
/// This class expands the PMNS_Fast class including the use of a first order
/// Taylor expansion to calculate the average on bins faster.
///
/// The model assumes a first order expansion over neutrino energy and angle
/// for both dynamical variables at the same time or for only one.
///
/// This is the first version of this class. A second version will be release
/// with a better implementation with the other classes.
///
/// Reference: https://doi.org/10.48550/arXiv.2308.00037
///
/// \sa PMNS_Fast
///
/// \author jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_Avg_H
#define PMNS_Avg_H

#include "PMNS_Fast.h"
#include "PremModel.h"

namespace OscProb {

  class PMNS_Avg : public PMNS_Fast {
    public:
      PMNS_Avg();          ///< Constructor
      virtual ~PMNS_Avg(); ///< Destructor

      virtual void SetPremModel(OscProb::PremModel& prem);

      // Get probability averaged over a bin
      virtual double AvgProb(
          int flvi, int flvf, double E,
          double dE); ///< Compute the average probability over
                      ///< a bin of energy with a Taylor expansion

      virtual double AvgProbLoE(
          int flvi, int flvf, double LoE,
          double dLoE); ///< Compute the average probability over
                        ///< a bin of LoE with a Taylor expansion

      virtual double AvgProb(
          int flvi, int flvf, double E, double cosT,
          double dcosT); ///< Compute the average probability over a
                         ///< bin of cosTheta with a Taylor expansion

      virtual double AvgProb(int flvi, int flvf, double E, double dE,
                             double cosT,
                             double dcosT); ///< Compute the average probability
                                            ///< over a bin of cosTheta and
                                            ///< energy with a Taylor expansion

      virtual double AvgProbLoE(
          int flvi, int flvf, double LoE, double dLoE, double cosT,
          double dcosT); ///< Compute the average probability over a
                         ///< bin of cosTheta and LoE with a Taylor
                         ///< expansion

            // Get probability vector averaged over a bin
      virtual vectorD AvgProbVector(
          int flvi, double E,
          double dE = 0); ///< Compute the average probability vector over a bin
                          ///< of energy using a Taylor expansion
      virtual vectorD AvgProbVectorLoE(
          int flvi, double LoE,
          double dLoE =
              0); ///< Compute the average probability vector over a bin of L/E using a Taylor expansion

      virtual matrixD AvgProbMatrix(int nflvi, int nflvf, double E, double dE);

      virtual matrixD AvgProbMatrixLoE(int nflvi, int nflvf, double LoE, double dLoE);

      // Get probability
      virtual double ExtrapolationProb(
          int flvi, int flvf, double E,
          double dE); ///< Compute the probability of flvi going to
                      ///< flvf for an energy E+dE

      virtual double ExtrapolationProbLoE(
          int flvi, int flvf, double LoE,
          double dLoE); ///< Compute the probability of flvi going to
                        ///< flvf at LoE+dLoE

      virtual double ExtrapolationProbCosT(
          int flvi, int flvf, double cosT,
          double dcosT); ///< Compute the probability of flvi going to
                         ///< flvf for an angle cosT+dcosT

    protected:
      virtual void InitializeTaylorsVectors(); ///< Initialize all member
                                               ///< vectors with zeros

      virtual void SetwidthBin(double dE,
                               double dcosT); ///< Set bin's widths for
                                              ///< energy and angle

      virtual void SetCosT(double cosT); ///< Set neutrino angle.

      virtual vectorD GetSamplePoints(
          double LoE,
          double dLoE); ///< Compute the sample points for
                        ///< a bin of L/E with width dLoE

      virtual vectorD GetSamplePoints(
          double E, double cosT,
          double dcosT); ///< Compute the sample points for
                         ///< a bin of cosT with width dcosT

      virtual matrixC GetSamplePoints(double InvE, double dInvE, double cosT,
                                      double dcosT); ///< Compute the sample
                                                     ///< points for a bin
                                                     ///< of E and cosT with
                                                     ///< width dE and dcosT

      // Construction of the K matrices
      virtual void BuildKE(
          double L); ///< build K matrix for the inverse of energy in mass basis
      virtual void BuildKcosT(
          double L); ///< build K matrix for angle in flavor basis

      virtual double LnDerivative(); ///< Compute the derivation of one layer's
                                     ///< length

      // Rotation from mass to flavor basis
      virtual void rotateS(); ///< Rotate the S matrix
      virtual void rotateK(); ///< Rotate one K matrix

      // Multiplication rule between two pairs of (S,K)
      virtual void MultiplicationRuleS(); ///< Product between two S matrices
      virtual void MultiplicationRuleK(
          complexD K[3][3]); ///< Product between two K matrices

      /// Solve the K matrix
      virtual void SolveK(complexD K[3][3], vectorD& lambda,
                          matrixC& V); ///< Solve the K matrix for
                                       ///< eigenvectors and eigenvalues

      // Propagating
      virtual void PropagatePathTaylor(
          NuPath p); ///< Propagate neutrino through a single path
      virtual void PropagateTaylor(); ///< Propagate neutrino through full path

      virtual double AvgFormula(
          int flvi, int flvf, double dbin, vectorD flambda,
          matrixC fV); ///< Formula for the average probability over a bin
                       ///< of width dbin

      virtual double AvgAlgo(
          int flvi, int flvf, double LoE, double dLoE,
          double L); ///< Algorithm for the compute of the average
                     ///< probability over a bin of LoE

      virtual double AvgAlgoCosT(
          int flvi, int flvf, double E, double cosT,
          double dcosT); ///< Algorithm for the compute of the average
                         ///< probability over a bin of cosT

      virtual double AvgAlgo(
          int flvi, int flvf, double InvE, double dInvE, double cosT,
          double dcosT); ///< Algorithm for the compute of the average
                         ///< probability over a bin of 1oE and cosT

      virtual double AlgorithmDensityMatrix(
          int flvi, int flvf); ///< Algorithm for the transformations on the
                               ///< density matrix

      virtual void RotateDensityM(
          bool to_basis, matrixC V); ///< Apply rotation to the density matrix

      virtual void HadamardProduct(vectorD lambda,
                                   double  dbin); ///< Apply an Hadamard Product
                                                 ///< to the density matrix

      virtual double AvgFormulaExtrapolation(
          int flvi, int flvf, double dbin, vectorD flambda,
          matrixC fV); ///< Formula for the extrapolation of probability

      // Attributes

      matrixC fevolutionMatrixS; ///< Evolution matrix S for reference energy
                                 ///< and angle for the entire path

      matrixC fSflavor; ///< S matrix for one layer
      matrixC fKmass;   ///< K matrix in mass basis for one layer
      matrixC fKflavor; ///< K matrix in flavor basis for one layer

      double fdInvE; ///< Bin's width for the inverse of energy in GeV-1

      complexD fKInvE[3][3]; ///< K matrix for the inverse of energy in GeV for
                             ///< the entire path
      vectorD flambdaInvE;   ///< Eigenvectors of K_invE
      matrixC fVInvE;        ///< Eigenvalues of K_invE

      double fcosT;  ///<  Cosine of neutrino angle
      double fdcosT; ///< Bin's width for angle

      complexD fKcosT[3]
                     [3];  ///< K matrix for neutrino angle for the entire path
      vectorD flambdaCosT; ///< Eigenvectors of K_cosTheta
      matrixC fVcosT;      ///< Eigenvalues of K_cosTheta

      matrixC fdensityMatrix; ///< The neutrino density matrix state

      // Variables for the compute of the derivation of one layer's length
      int    fLayer;
      int    fdl;
      double fDetRadius;
      double fminRsq;

      // Copy of the earth model used
      OscProb::PremModel fPrem;
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
