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

      virtual void SetwidthBin(double dE , double dcosT);

      virtual void SetCosT(double cosT);

      virtual double avgProbTaylor(int flvi, int flvf, double E , double dE);

      virtual double avgProbTaylorLoE(int flvi, int flvf, double LoE , double dLoE);

      virtual double avgProbTaylor1oE(int flvi, int flvf, double ONEoE , double d1oE);

      virtual double avgProbTaylor(int flvi, int flvf, double E , double dE, double cosT , double dcosT);

      virtual double avgProbTaylorAngle(int flvi, int flvf, double E, double cosT , double dcosT);

      virtual vectorD ConvertEto1oE(double E, double dE);

      virtual double interpolationEnergy(int flvi, int flvf, double E , double dE);

      virtual double interpolationEnergyLoE(int flvi, int flvf, double LoE , double dLoE);

      virtual double interpolationCosT(int flvi, int flvf, double cosT , double dcosT);

    protected:
      virtual void InitializeTaylorsVectors(); ///< Initialize all member vectors with
      ///< zeros

      virtual void BuildKE(double L , matrixC& K);

      virtual void BuildKcosT(double L, matrixC& K);

      virtual void SolveK(complexD K[3][3], vectorD& lambda, matrixC& V);

      virtual void PropagatePathTaylor(
        NuPath p );            ///< Propagate neutrino through a single path
      virtual void PropagateTaylor(); ///< Propagate neutrino through full path

      virtual void rotateS(vectorC fPhases,matrixC& S);

      virtual void rotateK(matrixC Kmass,matrixC& Kflavor);

      virtual void RotateDensityM(bool to_mass, matrixC V, matrixC& densityMatrix);

      virtual void HadamardProduct(vectorD lambda, matrixC& densityMatrix, double dbin);

      virtual void MultiplicationRuleK(matrixC KLayer, complexD K[3][3]);

      virtual void MultiplicationRuleS(matrixC SLayer);

      virtual double avgFormula(int flvi, int flvf, double dbin, vectorD flambda, matrixC fV); 

      virtual double avgFormulaExtrapolation(int flvi, int flvf, double dbin, vectorD flambda, matrixC fV);

      virtual double avgAlgorithm(int flvi, int flvf);

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
