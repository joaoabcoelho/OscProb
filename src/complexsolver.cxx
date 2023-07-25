
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
                   OscProb::matrixC& Q,
                   OscProb::matrixC& Qinv,
                   OscProb::vectorC& w)
{

  int size=w.size();
  Eigen::MatrixXcd l(size,size);
  Eigen::MatrixXcd evc(size,size), evcinv(size,size);
  Eigen::VectorXcd evl(size);
  
  for(int n=0; n<size; n++){
    for(int m=0; m<size; m++){
      l(n,m)=A[n][m];
    } 
  }  
                  
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigensolver;
  eigensolver.compute(l);
 
  if(eigensolver.info() != Eigen::Success){
    std::cout << "The diagonalization is failing" << std::endl;
    abort();
  }
  
  evl = eigensolver.eigenvalues();
  evc= eigensolver.eigenvectors();
  evcinv=evc.inverse();
  for(int n=0; n<size; n++){
    for(int m=0; m<size; m++){
      Q[n][m] = evc(n,m);
      Qinv[n][m]= evcinv(n,m);
    }
  }  
  for(int t=0; t<size; t++){
    w[t] = evl(t);
  }
              
}
