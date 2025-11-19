///////////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_Maltoni
///
/// \brief Implementation of oscillations of neutrinos in matter
///        in a framework with a Taylor expansion.
///
/// This class expands the PMNS_Base class including the use of a
/// Taylor expansion to calculate the average on bins faster.
///
/// The model assumes a first order expansion over energy or direction
/// for both dynamical variables at the same time or for only one.
///
/// This version is limited to classes where the vacuum Hamiltonian
/// is inversely proportional to energy and the matter potential is
/// independent of energy.
///
/// Reference: https://doi.org/10.48550/arXiv.2308.00037
///
/// \sa PMNS_Base
///
/// \author mloup@apc.in2p3.fr rogan.clark@uclouvain.be
///////////////////////////////////////////////////////////////////////////////

#include "PMNS_Base.h"
#include <Eigen/Core>

#ifndef PMNS_BASEAVG_H
#define PMNS_BASEAVG_H

#include "PremModel.h"

namespace OscProb {

  class PMNS_Maltoni : public PMNS_Base {
    public:
      PMNS_Maltoni(int numNus); ///< Constructor
      virtual void SetPremModel(OscProb::PremModel& prem);

      using PMNS_Base::AvgProb;
      using PMNS_Base::AvgProbLoE;

      /// Compute the average probability over a
      /// bin of energy with a Taylor expansion
      virtual double AvgProb(int flvi, int flvf, double E, double dE);

      /// Compute the average probability over a
      /// bin of cosTheta with a Taylor expansion
      virtual double AvgProb(int flvi, int flvf, double E, double cosT,
                             double dcosT);

      /// Compute the average probabilit over a
      /// bin of cosTheta and energy with a Taylor expansion
      virtual double AvgProb(int flvi, int flvf, double E, double dE,
                             double cosT, double dcosT);

      /// Compute the average probability over a
      /// bin of LoE with a Taylor expansion
      virtual double AvgProbLoE(int flvi, int flvf, double LoE, double dLoE);

      /// Compute the average probability over a
      /// bin of cosTheta and LoE with a Taylor expansion
      virtual double AvgProbLoE(int flvi, int flvf, double LoE, double dLoE,
                                double cosT, double dcosT);

      /// Compute the average probability vector over a
      /// bin of energy using a Taylor expansion
      virtual vectorD AvgProbVector(int flvi, double E, double dE);

      /// Compute the average probability vector over a
      /// bin of L/E using a Taylor expansion
      virtual vectorD AvgProbVectorLoE(int flvi, double LoE, double dLoE);

      /// Compute the average probability matrix over a
      /// bin of energy using a Taylor expansion
      virtual matrixD AvgProbMatrix(int nflvi, int nflvf, double E, double dE);

      /// Compute the average probability matrix over a
      /// bin of L/E using a Taylor expansion
      virtual matrixD AvgProbMatrixLoE(int nflvi, int nflvf, double LoE,
                                       double dLoE);

      /// Compute the probability of flvi going to flvf for an energy E+dE
      virtual double ExtrapolationProb(int flvi, int flvf, double E, double dE);

      /// Compute the probability of flvi going to flvf at LoE+dLoE
      virtual double ExtrapolationProbLoE(int flvi, int flvf, double LoE,
                                          double dLoE);

      /// Compute the probability of flvi going to flvf for an angle cosT+dcosT
      virtual double ExtrapolationProbCosT(int flvi, int flvf, double cosT,
                                           double dcosT);

      /// Set flag for which averaging to use
      virtual void SetIsOscProbAvg(bool isOscProbAvg);

    protected:
      /// Compute the sample points fo a bin of L/E with width dLoE
      virtual vectorD GetSamplePointsAvgClass(double LoE, double dLoE);

      /// Compute the sample points for a bin of cosT with width dcosT
      virtual vectorD GetSamplePointsAvgClass(double E, double cosT,
                                              double dcosT);

      /// Compute the sample points for a bin of
      /// E and cosT with width dE and dcosT
      virtual matrixC GetSamplePointsAvgClass(double InvE, double dInvE,
                                              double cosT, double dcosT);

      virtual void UpdateHam() = 0; ///< Build the full Hamiltonian

      /// Initialize all member vectors with zeros
      virtual void InitializeTaylorsVectors();

      /// Set bin's widths for energy and angle
      virtual void SetwidthBin(double dE, double dcosT);

      virtual void SetCosT(double cosT); ///< Set neutrino angle.

      /// build K matrix for the inverse of energy in mass basis
      virtual void BuildKE(double L);

      /// build K matrix for angle in flavor basis
      virtual void BuildKcosT();

      /// Compute the derivation of one layer's length
      virtual double LnDerivative();

      /// Propagate neutrino through a single path
      virtual void PropagatePathTaylor(NuPath p);

      /// Propagate neutrino through full path
      virtual void PropagateTaylor();

      // Rotation from mass to flavor basis
      virtual void rotateS(); ///< Rotate the S matrix
      virtual void rotateK(); ///< Rotate one K matrix

      /// Product between two S matrices
      virtual void MultiplicationRuleS();

      /// Product between two K matrices
      virtual void MultiplicationRuleK(Eigen::MatrixXcd& K);

      /// Solve the K matrix for eigenvectors and eigenvalues
      void SolveK(Eigen::MatrixXcd& K, vectorD& lambda, matrixC& V);

      /// Auxiliary function to choose eigensystem method
      template <typename T>
      void TemplateSolver(Eigen::MatrixXcd& K, vectorD& lambda, matrixC& V);

      // Propagating
      /// Formula for the average probability over a bin of width dbin
      virtual double AvgFormula(int flvi, int flvf, double dbin,
                                vectorD flambda, matrixC fV);

      /// Formula for the extrapolation of probability
      virtual double AvgFormulaExtrapolation(int flvi, int flvf, double dE,
                                             vectorD flambda, matrixC fV);

      /// Algorithm for the compute of the average
      /// probability over a bin of LoE
      virtual double AvgAlgo(int flvi, int flvf, double LoE, double dLoE,
                             double L);

      /// Algorithm for the compute of the average
      /// probability over a bin of 1oE and cosT
      virtual double AvgAlgo(int flvi, int flvf, double InvE, double dInvE,
                             double cosT, double dcosT);

      /// Algorithm for the compute of the average
      /// probability over a bin of cosT
      virtual double AvgAlgoCosT(int flvi, int flvf, double E, double cosT,
                                 double dcosT);

      /// Algorithm for the transformations on the density matrix
      virtual double AlgorithmDensityMatrix(int flvi, int flvf);

      /// Apply rotation to the density matrix
      virtual void RotateDensityM(bool to_basis, matrixC V);

      /// Apply an Hadamard Product to the density matrix
      virtual void HadamardProduct(vectorD lambda, double dbin);

      double fcosT; ///< Cosine of zenith angle

      double fdInvE; ///< Bin's width for the inverse of energy
      double fdcosT; ///< Bin's width for zenith angle

      vectorD flambdaInvE; ///< Eigenvectors of K_invE
      vectorD flambdaCosT; ///< Eigenvectors of K_cosTheta
      matrixC fVInvE;      ///< Eigenvalues of K_invE

      /// K matrix for the inverse of energy for the entire path
      Eigen::MatrixXcd fKInvE;
      /// K matrix for neutrino angle for the entire path
      Eigen::MatrixXcd fKcosT;

      /// Evolution matrix S for reference energy
      /// and angle for the entire path
      matrixC fevolutionMatrixS;
      matrixC fVcosT;   ///< Eigenvalues of K_cosTheta
      matrixC fSflavor; ///< S matrix for one layer
      matrixC fKmass;   ///< K matrix in mass basis for one layer
      matrixC fKflavor; ///< K matrix in flavor basis for one layer

      matrixC fdensityMatrix; ///< The neutrino density matrix state

      // Variables for the compute of the derivation of one layer's length
      int    fLayer;     ///< Layer index
      int    fdl;        ///< Length derivative
      double fDetRadius; ///< Detector radius
      double fminRsq;    ///< Minimum square radius

      /// Earth model used
      OscProb::PremModel fPrem;

      Eigen::MatrixXcd fHam; ///< The full Hamiltonian

      bool fIsOscProbAvg; ///< Flag to call OscProb default or Maltoni average
  };

} // namespace OscProb

#endif
