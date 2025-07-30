#ifndef PMNS_OQS_H
#define PMNS_OQS_H

#include <Eigen/Core>

#include "PMNS_DensityMatrix.h"

namespace OscProb {

  class PMNS_OQS : public PMNS_DensityMatrix {
    public:
      PMNS_OQS();          ///< Constructor
      virtual ~PMNS_OQS(); ///< Destructor

      virtual void SetParameterisation(int par);
      virtual void SetPhi(int i, double val);
      virtual void Seta(int i, double val);
      virtual void Setcos(int i, int j, double val);

    protected:
      virtual void SetDissipatorElement(int i, int j);
      virtual void SetDissipator();
      virtual void InitializeVectors();
      virtual void SetHeff(NuPath p);
      virtual void SetHGM();
      virtual void SetM();
      virtual void RotateState(bool to_mass); ///< Rotate rho to/from mass basis
      virtual void ChangeBaseToGM();
      virtual void ChangeBaseToSU3();

      virtual void FillCache() {} ///< Deactivate cache

      /// Propagate neutrino through a single path
      virtual void PropagatePath(NuPath p);

      int      fParameterisation;
      double   fPhi[2]; ///< Majorana phases
      complexD fR[9];
      complexD fRt[9];

      matrixC fHeff;
      matrixC fHGM;

      matrixD fD;   ///< Off-diagonal, 9x9 dissipator
      matrixC fM;   ///< M
      vectorD fa;   ///< a vector
      matrixD fcos; ///< cosines ai . aj

      Eigen::MatrixXcd fMe; ///< Buffer matrix for exponential

      bool fBuiltDissipator; ///< Flag to rebuilt D
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
