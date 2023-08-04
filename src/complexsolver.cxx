
#include <iostream>

#include <Eigen/Eigenvalues>

#include "complexsolver.h"

//.............................................................................
///
/// Wrapper to solve non-hermitian matrix eigenvalues.
///
/// @param A    - Input matrix
/// @param w    - Output eigenvalues
///
void complexsolver(const Eigen::MatrixXcd& A,
                   OscProb::vectorD& w)
{

  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigensolver;
  eigensolver.compute(A);

  if(eigensolver.info() != Eigen::Success){
    std::cerr << "ERROR: The diagonalization is failing" << std::endl;
    abort();
  }

  for(int t=0; t<w.size(); t++){
    w[t] = eigensolver.eigenvalues()(t).real();
  }

}
