////////////////////////////////////////////////////////////////////////
///
/// \brief Some useful general definitions
///
/// \author Joao Coelho - jcoelho\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <complex>
#include <vector>

namespace OscProb {

  typedef std::vector<int> vectorI;

  typedef std::vector<double>  vectorD;
  typedef std::vector<vectorD> matrixD;

  typedef std::complex<double> complexD;
  typedef std::vector<complexD> vectorC;
  typedef std::vector<vectorC>  matrixC;
  
}

#endif
