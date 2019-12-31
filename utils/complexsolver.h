#include <iostream>
#include "Eigenvalues"

using namespace Eigen;
using namespace std;

void complexsolver(std::complex<double> A[3][3], std::complex<double> Q[3][3], std::complex<double> Qinv[3][3], std::complex<double> w[3])
{

  Matrix3cd l;
  Matrix3cd evc, evcinv;
  Vector3cd evl;
  
  for(int n=0; n<3; n++){
    for(int m=0; m<3; m++){
      l(n,m)=A[n][m];
    } 
  }  
                  
  ComplexEigenSolver<Matrix3cd> eigensolver;
  eigensolver.compute(l);
 
  if(eigensolver.info() != Success){ 
    cout << "The diagonalization is failing" << endl;
    abort();
  }
  
  evl = eigensolver.eigenvalues();
  evc= eigensolver.eigenvectors();
  evcinv=evc.inverse();
  for(int n=0; n<3; n++){
    for(int m=0; m<3; m++){
      Q[n][m] = evc(n,m);
      Qinv[n][m]= evcinv(n,m);
    } 
  }  
  for(int t=0; t<3; t++){
    w[t] = evl(t);
  }
              
}
