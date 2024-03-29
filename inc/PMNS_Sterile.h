///////////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_Sterile
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        N-neutrino framework.
///
/// This class implements neutrino oscillations for any number of
/// neutrino falvours. The extra flavours above 3 are treated as
/// sterile neutrinos, i.e. it implements a 3+N model. One can also
/// run a 2-neutrino model in principle, but an analytical solution
/// would be feasible and more efficient in that case.
///
/// Reference: https://doi.org/10.1007/JHEP12(2013)014
///
/// The extended mixing matrix follows the recommendations from:
/// https://doi.org/10.1016/0370-2693(86)91268-2
///
/// The eigensystem solution is performed by
/// [Eigen](https://eigen.tuxfamily.org/), which makes it slower
/// than the PMNS_Fast class, which only works for 3 neutrinos.
///
/// \author jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_STERILE_H
#define PMNS_STERILE_H

#include <Eigen/Core>

#include "PMNS_Base.h"

namespace OscProb {

  class PMNS_Sterile : public PMNS_Base {
    public:
      PMNS_Sterile(int numNus); ///< Constructor

    protected:
      /// Build the full Hamiltonian
      virtual void UpdateHam();

      /// Solve the full Hamiltonian for eigenvectors and eigenvalues
      virtual void SolveHam();

      /// Specialized solver for NxN matrices
      template <typename T> void SolveEigenSystem();

      Eigen::MatrixXcd fHam; ///< The full Hamiltonian
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
