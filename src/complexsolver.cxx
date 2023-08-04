
#include <iostream>

#include "Eigen/Eigenvalues"

#include "complexsolver.h"

//.............................................................................
///
/// Wrapper to solve non-hermitian matrix eigensystem.
///
/// @param A    - Input matrix
/// @param Q    - Output eigenvectors
/// @param Qinv - Inverse eigenvectors
/// @param w    - Output eigenvalues
///
void complexsolver(OscProb::matrixC  A,
                   OscProb::vectorD& w)
{

  int size=w.size();
  Eigen::MatrixXcd l(size,size);
  Eigen::VectorXcd evl(size);
  
  for(int n=0; n<size; n++){
    for(int m=0; m<size; m++){
      l(n,m)=A[n][m];
    } 
  }
                  
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigensolver;
  eigensolver.compute(l);

  evl = eigensolver.eigenvalues();

  for(int t=0; t<size; t++){
    w[t] = evl(t).real();
  }
 
}
