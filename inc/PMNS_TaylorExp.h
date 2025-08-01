///////////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_TaylorExp
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
/// Reference: https://doi.org/10.48550/arXiv.2308.00037
///
/// \sa PMNS_Fast
///
/// \author jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_TaylorExp_H
#define PMNS_TaylorExp_H

#include "PMNS_Fast.h"
#include "MatrixDecomp/zheevh3.h"

namespace OscProb {

  class PMNS_TaylorExp : public PMNS_Fast {
    public:
      PMNS_TaylorExp();          ///< Constructor
      virtual ~PMNS_TaylorExp(); ///< Destructor

      //@@@@@@
      virtual void printMatrix1(matrixC M);
      virtual void printMatrix2(complexD M[3][3]);
      virtual void LenghtLayer();
      //@@@@@@

      virtual void SetwidthBin(double dE, 
        double dcosT); ///< Set bin's widths for 
                       ///< energy and angle

      virtual void SetCosT(double cosT); ///< Set neutrino angle.

      virtual vectorD ConvertEto1oE(
        double E, double dE); ///< Convert a bin of energy into a bin of 1/E

      // Get probability averaged over a bin
      virtual double avgProbTaylor(int flvi, int flvf, 
        double E, double dE); ///< Compute the average probability over 
                              ///< a bin of energy with a Taylor expansion

      virtual double avgProbTaylorLoE(int flvi, int flvf, 
        double LoE, double dLoE);  ///< Compute the average probability over 
                                   ///< a bin of LoE with a Taylor expansion

      virtual double avgProbTaylor1oE(int flvi, int flvf, 
        double ONEoE, double d1oE);  ///< Compute the average probability over 
                                     ///< a bin of 1oE with a Taylor expansion

      virtual double avgProbTaylorAngle(int flvi, int flvf, double E, 
        double cosT, double dcosT); ///< Compute the average probability over a 
                                    ///< bin of cosTheta with a Taylor expansion

      virtual double avgProbTaylor(int flvi, int flvf, double E, 
        double dE, 
        double cosT, double dcosT); ///< Compute the average probability over a 
                                    ///< bin of cosTheta and energy with a Taylor
                                    ///< expansion

      // blablabla
      virtual double interpolationEnergy(int flvi, int flvf, 
        double E , double dE);      /// < 

      virtual double interpolationEnergyLoE(int flvi, int flvf, 
        double LoE , double dLoE);

      virtual double interpolationCosT(int flvi, int flvf, 
        double cosT , double dcosT);


    protected:

      virtual void InitializeTaylorsVectors(); ///< Initialize all member vectors 
                                               ///< with zeros

      // Construction of the K matrices                                   
      virtual void BuildKE(
        double L, matrixC& K); ///< build K matrix for the inverse of energy in mass basis
      virtual void BuildKcosT(
        double L, matrixC& K); ///< build K matrix for angle in flavor basis

      // Rotation from mass to flavor basis
      virtual void rotateS(vectorC fPhases,matrixC& S); ///< Rotate the S matrix 
      virtual void rotateK(matrixC Kmass,matrixC& Kflavor); ///< Rotate one K matrix 

      // Multiplication rule between two pairs of (S,K) 
      virtual void MultiplicationRuleS(matrixC SLayer); /// < Product between two S matrices
      virtual void MultiplicationRuleK(matrixC KLayer, complexD K[3][3]); ///< Product between 
                                                                          ///< two K matrices

      /// Solve one K matrix
      virtual void SolveK(complexD K[3][3], 
        vectorD& lambda, matrixC& V); ///< Solve one K matrix for
                                      ///< eigenvectors and eigenvalues

      // Propagating
      virtual void PropagatePathTaylor(
        NuPath p );                   ///< Propagate neutrino through a single path
      virtual void PropagateTaylor(); ///< Propagate neutrino through full path

      // Avg for only one dynamical variable
      virtual double avgFormula(int flvi, int flvf, double dbin, 
        vectorD flambda, matrixC fV); ///< Formula for the average probability over a bin

      // Avg on energy and cosT at the same time 
      virtual void RotateDensityM(bool to_mass, matrixC V, matrixC& densityMatrix);
      virtual void HadamardProduct(vectorD lambda, matrixC& densityMatrix, double dbin);
      virtual double avgAlgorithm(int flvi, int flvf);

      // Interpolation for only one dynamical variable
      virtual double avgFormulaExtrapolation(int flvi, int flvf, double dbin, vectorD flambda, matrixC fV);


      // Attributes

      matrixC fevolutionMatrixS;  ///< Evolution matrix S for reference energy and angle 

      double fdInvE; ///< Bin's width for the inverse of energy in GeV-1 

      complexD fKInvE[3][3];  ///< K matrix for the inverse of energy in GeV
      vectorD flambdaInvE;    ///< Eigenvectors of K_invE
      matrixC fVInvE;         ///< Eigenvalues of K_invE

      double fcosT;   ///<  Cosine of neutrino angle 
      double fdcosT;  ///< Bin's width for angle

      complexD fKcosT[3][3];  ///< K matrix for neutrino angle 
      vectorD flambdaCosT;    ///< Eigenvectors of K_cosTheta
      matrixC fVcosT;         ///< Eigenvalues of K_cosTheta

      matrixC densityMatrix; ///< The neutrino density matrix state

      std::vector<NuPath> fNuPathsVariation ;

  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
