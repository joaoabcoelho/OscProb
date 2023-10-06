
#include <Eigen/Eigenvalues>

#include "complexsolver.h"

//.............................................................................
///
/// Wrapper to solve non-hermitian matrix eigenvalues.
///
/// @param A    - Input matrix
/// @param w    - Output eigenvalues
///
void complexsolver(const Eigen::Matrix3cd& A, OscProb::vectorD& w)
{
  Eigen::ComplexEigenSolver<Eigen::Matrix3cd> eigensolver;
  eigensolver.compute(A);

  for (int t = 0; t < w.size(); t++) {
    w[t] = eigensolver.eigenvalues()(t).real();
  }
}
