///////////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_NUNM
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework with NUNM.
///
/// This class expands the PMNS_Fast class including a general
/// matter potential matrix describing Non-Standard Interactions (NUNM).
///
/// The non unitarity effect is parametrized by dimensionless quantities
/// alpha which quantify the intensity of the NUNM with respect to the
/// standard mixing
///
/// Reference: https://arxiv.org/pdf/2309.16942.pdf
///
/// \sa PMNS_Fast
/// 
/// \author cerisy@cppm.in2p3.fr
/// \author jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_NUNM_H
#define PMNS_NUNM_H

#include "PMNS_Fast.h"

namespace OscProb {

  class PMNS_NUNM : public PMNS_Fast {
    public:
      PMNS_NUNM();          ///< Constructor
      virtual ~PMNS_NUNM(); ///< Destructor

      /// Set any given NUNM parameter
      virtual void SetAlpha(int i, int j, double val, double phase);

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

    protected:
      /// Build the full Hamiltonian
      virtual void UpdateHam();
      //virtual void BuildHms();
      virtual void transfoNUNM(matrixC& Ham);
      complexD fAlpha[3][3]; ///< Stores each NUNM parameter
      matrixC Ham;
      double fracVnc;
  };

} // namespace OscProb

#endif

////////////////////////////////////////////////////////////////////////
