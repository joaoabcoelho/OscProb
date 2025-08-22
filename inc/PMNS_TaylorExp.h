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

      virtual void GetPremLayers(std::vector<PremLayer> PremLayers);

      // Get probability averaged over a bin
      using PMNS_Base::AvgProb;
      virtual double AvgProb(int flvi, int flvf, 
        double E, double dE); ///< Compute the average probability over 
                              ///< a bin of energy with a Taylor expansion

      using PMNS_Base::AvgProbLoE;
      virtual double AvgProbLoE(int flvi, int flvf, 
        double LoE, double dLoE);  ///< Compute the average probability over 
                                   ///< a bin of LoE with a Taylor expansion

      virtual double AvgProb(int flvi, int flvf, double E, 
        double cosT, double dcosT); ///< Compute the average probability over a 
                                    ///< bin of cosTheta with a Taylor expansion

      virtual double AvgProb(int flvi, int flvf, double E, 
        double dE, 
        double cosT, double dcosT); ///< Compute the average probability over a 
                                    ///< bin of cosTheta and energy with a Taylor
                                    ///< expansion

      virtual double AvgProbLoE(int flvi, int flvf, double LoE, double dLoE, double cosT , double dcosT);

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

      virtual vectorD GetSamplePoints(double LoE, double dLoE);
      virtual vectorD GetSamplePoints(double E, double cosT, double dcosT);
      virtual vector<vectorD> GetSamplePoints(double InvE, double dInvE, double cosT, double dcosT);


      // Construction of the K matrices                                   
      virtual void BuildKE(
        double L); ///< build K matrix for the inverse of energy in mass basis
      virtual void BuildKcosT(
        double L); ///< build K matrix for angle in flavor basis

      virtual double LnDerivative();


      // Rotation from mass to flavor basis
      virtual void rotateS(); ///< Rotate the S matrix 
      virtual void rotateK(); ///< Rotate one K matrix 

      // Multiplication rule between two pairs of (S,K) 
      virtual void MultiplicationRuleS(); /// < Product between two S matrices
      virtual void MultiplicationRuleK( complexD K[3][3]); ///< Product between 
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
      virtual double AvgFormula(int flvi, int flvf, double dbin, 
        vectorD flambda, matrixC fV); ///< Formula for the average probability over a bin

      virtual double AvgAlgo(int flvi, int flvf, double LoE , double dLoE, double L);
      virtual double AvgAlgoCosT(int flvi, int flvf, double E , double cosT, double dcosT);
      virtual double AvgAlgo(int flvi, int flvf, double InvE , double dInvE, double cosT , double dcosT);

      // Avg on energy and cosT at the same time 
      virtual void RotateDensityM(bool to_mass, matrixC V);
      virtual void HadamardProduct(vectorD lambda, double dbin);
      virtual double AlgorithmDensityMatrix(int flvi, int flvf);

      // Interpolation for only one dynamical variable
      virtual double AvgFormulaExtrapolation(int flvi, int flvf, double dbin, vectorD flambda, matrixC fV);


      // Attributes

      matrixC fevolutionMatrixS;  ///< Evolution matrix S for reference energy and angle for the entire path 

      matrixC fSflavor; ///< S matrix for one layer 
      matrixC fKmass;   ///< K matrix in mass basis for one layer 
      matrixC fKflavor; ///< K matrix in flavor basis for one layer 

      double fdInvE; ///< Bin's width for the inverse of energy in GeV-1 

      complexD fKInvE[3][3];  ///< K matrix for the inverse of energy in GeV for the entire path 
      vectorD flambdaInvE;    ///< Eigenvectors of K_invE
      matrixC fVInvE;         ///< Eigenvalues of K_invE

      double fcosT;   ///<  Cosine of neutrino angle 
      double fdcosT;  ///< Bin's width for angle

      complexD fKcosT[3][3];  ///< K matrix for neutrino angle for the entire path 
      vectorD flambdaCosT;    ///< Eigenvectors of K_cosTheta
      matrixC fVcosT;         ///< Eigenvalues of K_cosTheta

      matrixC fdensityMatrix; ///< The neutrino density matrix state

      int flayer;
      int fdl;
      double fDetRadius;
      double fminRsq ;
      std::vector<PremLayer> fPremLayers; ///< The layers in the earth model

  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
