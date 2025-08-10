#ifndef PMNS_OQS_H
#define PMNS_OQS_H

#include <Eigen/Core>

#include "PMNS_DensityMatrix.h"

namespace OscProb {

  class PMNS_OQS : public PMNS_DensityMatrix {
    public:
      PMNS_OQS();          ///< Constructor
      virtual ~PMNS_OQS(); ///< Destructor

      static constexpr int SU3_DIM = 9;

      virtual void SetPower(int n);
      virtual void SetDecoElement(int i, double val);
      virtual void SetDecoAngle(int i, int j, double th);

      virtual int    GetPower();
      virtual double GetDecoElement(int i);
      virtual double GetDecoAngle(int i, int j);

      virtual double GetHGM(int i, int j);
      virtual double GetDissipatorElement(int i, int j);

      virtual void    SetIsNuBar(bool isNuBar);
      virtual matrixD ProbMatrix(int nflvi, int nflvf);

    protected:
      virtual void BuildUM();
      virtual void BuildHms();
      virtual void SetHeff(NuPath p);
      virtual void SetHGM();
      virtual void SetDissipator();
      virtual void SetM();
      virtual void RotateState(bool to_mass); ///< Rotate rho to/from mass basis
      virtual void ChangeBaseToGM();
      virtual void ChangeBaseToSU3();

      virtual void Propagate();
      /// Propagate neutrino through a single path
      virtual void PropagatePath(NuPath p);

      int     fPower;
      vectorD fa;   ///< a vector
      matrixD fcos; ///< cosines ai . aj

      vectorD fR;
      vectorD fRt;
      matrixC fUM; ///< PMNS Matrix
      matrixC fHeff;
      matrixD fHGM;
      matrixD fD; ///< Off-diagonal, 9x9 dissipator

      Eigen::MatrixXd fM; ///< Buffer matrix for exponential

      bool fBuiltDissipator; ///< Flag to rebuilt D

  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
