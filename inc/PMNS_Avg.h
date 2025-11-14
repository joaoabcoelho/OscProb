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

#include "PMNS_Sterile.h"

namespace OscProb {

  class PMNS_Avg : public PMNS_Sterile {
    public:
      PMNS_Avg(int numNus); ///< Constructor
      virtual ~PMNS_Avg();  ///< Destructor
      // Get probability averaged over a bin
      using PMNS_Base::AvgProb;
      using PMNS_Base::AvgProbLoE;
      using PMNS_Base::AvgProbVector;
      using PMNS_Base::AvgProbVectorLoE;

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
          double dE); ///< Compute the average probability vector over a bin
                      ///< of energy using a Taylor expansion
      virtual vectorD AvgProbVectorLoE(
          int flvi, double LoE,
          double dLoE); ///< Compute the average probability vector over a
                        ///< bin of L/E using a Taylor expansion

      virtual matrixD AvgProbMatrix(int nflvi, int nflvf, double E, double dE);

      virtual matrixD AvgProbMatrixLoE(int nflvi, int nflvf, double LoE,
                                       double dLoE);

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
      virtual vectorD GetSamplePointsAvgClass(
          double LoE,
          double dLoE); ///< Compute the sample points for
                        ///< a bin of L/E with width dLoE

      virtual vectorD GetSamplePointsAvgClass(
          double E, double cosT,
          double dcosT); ///< Compute the sample points for
                         ///< a bin of cosT with width dcosT

      virtual matrixC GetSamplePointsAvgClass(
          double InvE, double dInvE, double cosT,
          double dcosT); ///< Compute the sample
                         ///< points for a bin
                         ///< of E and cosT with
                         ///< width dE and dcosT

      // Construction of the K matrices
      virtual void BuildKcosT(
          double L); ///< build K matrix for angle in flavor basis


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

                                 ///< and angle for the entire path



  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
