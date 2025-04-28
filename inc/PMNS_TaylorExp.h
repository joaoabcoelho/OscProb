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

      virtual void SetwidthBin(double dE , double dcosTheta);

      virtual double avgProbTaylor(int flvi, int flvf, double E , double dE);

      virtual double avgProbTaylorLoE(int flvi, int flvf, double LoE, double dLoE);

      virtual double avgProbTaylorAngle(int flvi, int flvf, double cosTheta , double dcosTheta);

      virtual vectorD ConvertLoEtoE(double LoE, double dLoE);

    protected:
      virtual void InitializeTaylorsVectors(); ///< Initialize all member vectors with
      ///< zeros

      virtual void BuildKE(double L , matrixC& K);

      virtual void SolveK();

      virtual void PropagatePathTaylor(
        NuPath p);            ///< Propagate neutrino through a single path
      virtual void PropagateTaylor(); ///< Propagate neutrino through full path

      virtual void rotateS(vectorC fPhases,matrixC& S);

      virtual void rotateK(matrixC Kmass,matrixC& Kflavor);

      virtual void MultiplicationRule(matrixC SLayer,matrixC KLayer);

      virtual double avgFormula(int flvi, int flvf); 

      // Attributes

      matrixC fevolutionMatrixS;  

      complexD fKE[3][3];
      vectorD flambdaE;
      matrixC fVE;

      complexD fKcosTheta[3][3];
      vectorD flambdacosTheta;
      matrixC fVcosTheta;

      double fdE;
      double fdcosTheta;
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
