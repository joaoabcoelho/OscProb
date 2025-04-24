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

      virtual void SetwidthBin(double widthBin);

      virtual double avgProbTaylor(int flvi, int flvf, double E , double widthBin);

      virtual double avgProbTaylorLoE(int flvi, int flvf, double ELoE, double widthBin);

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

      complexD fKE[3][3]; 
      matrixC fevolutionMatrixS;  

      vectorD flambdaE;
      matrixC fVE;

      double fwidthBin;
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
