///////////////////////////////////////////////////////////////////////////////
/// mettre info et biblio
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

      virtual void printMatrix1(matrixC M);
      virtual void printMatrix2(complexD M[3][3]);

      virtual void SetwidthBin(double dE , double dcosT);

      virtual void SetCosT(double cosT);

      virtual double avgProbTaylor(int flvi, int flvf, double E , double dE);

      virtual double avgProbTaylorLoE(int flvi, int flvf, double LoE , double dLoE);

      virtual double avgProbTaylor1oE(int flvi, int flvf, double ONEoE , double d1oE);

      virtual double avgProbTaylor(int flvi, int flvf, double E , double dE, double cosT , double dcosT);

      virtual double avgProbTaylorAngle(int flvi, int flvf, double E, double cosT , double dcosT);

      virtual vectorD ConvertLoEtoE(double LoE, double dLoE);

      virtual vectorD ConvertEto1oE(double E, double dE);

      virtual void LenghtLayer();

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

      matrixC fevolutionMatrixS;  

      complexD fKE[3][3];
      vectorD flambdaE;
      matrixC fVE;
      double fdE;

      complexD fKcosT[3][3];
      vectorD flambdaCosT;
      matrixC fVcosT;
      double fdcosT;

      double fcosT;

      std::vector<NuPath> fNuPathsVariation ;
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
