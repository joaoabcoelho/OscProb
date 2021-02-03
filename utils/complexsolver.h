#include <iostream>
#include "Eigenvalues"
#include <complex>
#include <vector>

using namespace Eigen;
using namespace std;


void complexsolver(std::vector< std::vector<complexD> > A, std::vector< std::vector<complexD> >& Q, std::vector< std::vector<complexD> >& Qinv, std::vector< complexD >& w)  
{
  int size=w.size();
  MatrixXcd l(size,size);
  MatrixXcd evc(size,size), evcinv(size,size);
  VectorXcd evl(size);
  
  for(int n=0; n<size; n++){
    for(int m=0; m<size; m++){
      l(n,m)=A[n][m];
    } 
  }  
                  
  ComplexEigenSolver<MatrixXcd> eigensolver;
  eigensolver.compute(l);
 
  if(eigensolver.info() != Success){ 
    cout << "The diagonalization is failing" << endl;
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
