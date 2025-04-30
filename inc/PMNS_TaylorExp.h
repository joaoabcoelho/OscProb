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

      virtual void SetwidthBin(double dE , double dcosT);

      virtual void SetCosT(double cosT);

      virtual double avgProbTaylor(int flvi, int flvf, double E , double dE);

      virtual double avgProbTaylorLoE(int flvi, int flvf, double LoE, double dLoE);

      virtual double avgProbTaylorAngle(int flvi, int flvf, double cosT , double dcosT);

      virtual vectorD ConvertLoEtoE(double LoE, double dLoE);

    protected:
      virtual void InitializeTaylorsVectors(); ///< Initialize all member vectors with
      ///< zeros

      virtual void BuildKE(double L , matrixC& K);

      virtual void SolveK(complexD K[3][3], vectorD& lambda, matrixC& V);

      virtual void PropagatePathTaylor(
        NuPath p , complexD K[3][3]);            ///< Propagate neutrino through a single path
      virtual void PropagateTaylor(complexD K[3][3]); ///< Propagate neutrino through full path

      virtual void rotateS(vectorC fPhases,matrixC& S);

      virtual void rotateK(matrixC Kmass,matrixC& Kflavor);

      virtual void MultiplicationRule(matrixC SLayer,matrixC KLayer, complexD K[3][3]);

      virtual double avgFormula(int flvi, int flvf, double dbin, vectorD flambda, matrixC fV); 

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
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
