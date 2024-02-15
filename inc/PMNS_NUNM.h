///////////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_NUNM
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework with Non unitary Neutrino Mixing (NUNM).
///
/// This class expands the PMNS_Fast class including a general NU mixing matrix
/// 
///
/// The non unitarity effect is parametrized by dimensionless quantities
/// alpha which quantify the deviation from unitarity with respect to the
/// standard mixing
///
/// Reference: https://arxiv.org/pdf/2111.00329.pdf
///
/// \sa PMNS_Fast
/// 
/// \author cerisy@cppm.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_NUNM_H
#define PMNS_NUNM_H

#include <Eigen/Dense>
#include "PMNS_Fast.h"
#include "Definitions.h"

namespace OscProb {

  class PMNS_NUNM : public PMNS_Fast {
    public:
      PMNS_NUNM(int scale = 0);          ///< Constructor
      virtual ~PMNS_NUNM(); ///< Destructor

      /// Set any given NUNM parameter
      virtual void SetAlpha(int i, int j, double val, double phase); // i, j between 0 and 2

      /// Get any given NUNM parameter
      virtual complexD GetAlpha(int i, int j);

      /// Set the NUNM parameters all at once
      void SetNUNM(double alpha_11, double alpha_21, double alpha_31,
                  double alpha_22, double alpha_32, double alpha_33);

      // Set the diagonal real NUNM parameters
      virtual void SetAlpha_11(double a);     ///< Set alpha_11 parameter
      virtual void SetAlpha_22(double a);   ///< Set alpha_22 parameter
      virtual void SetAlpha_33(double a); ///< Set alpha_33 parameter

      // Set the off-diagonal complex NUNM parameters
      virtual void SetAlpha_21(double a, double phi); ///< Set alpha_21 parameter
      virtual void SetAlpha_31(double a,
                               double phi); ///< Set alpha_31 parameter
      virtual void SetAlpha_32(double a,
                                double phi); ///< Set alpha_32 parameter
      virtual void SetFracVnc(double f);

      //virtual matrixD ProbMatrix(int nflvi, int nflvf);
    
    protected:
      int fscale;
      virtual void UpdateHam();
      //virtual matrixD ProbMatrix(int nflvi, int nflvf);
      //virtual void Propagate();
      virtual void PropagatePath(NuPath p);
      //virtual void getBar(Eigen::Matrix3cd& M); // apply normalisation to a matrix
      double fracVnc; // set fraction of matter potential affecting NC
      void InitMatrix();
      Eigen::Matrix<std::complex<double>, 3, 3> X;
      Eigen::Matrix<std::complex<double>, 3, 3> Alpha;
      Eigen::Matrix<std::complex<double>, 3, 3> V;
      Eigen::Matrix<std::complex<double>, 3, 3> Ham;
      Eigen::Matrix<std::complex<double>, 3, 3> Evec0;
      Eigen::Matrix<std::complex<double>, 3, 3> Evec;
      Eigen::Matrix<std::complex<double>, 3, 3> EvecA;
  };

} // namespace OscProb

#endif

////////////////////////////////////////////////////////////////////////
