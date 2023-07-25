
#include <cmath>
#include <numeric>

#include "Absorption.h"

using namespace std;

using namespace OscProb;

const double Absorption::kU  = 1.660539066e-24; // atomic mass unit [g]

Absorption::Absorption(){}
Absorption::~Absorption(){}

//.............................................................................
double Absorption::Trans(double xsec)
{

  vector<double> p_trans_vec;

  for(int i=0; i<int(fNuPaths.size()); i++){
    double n = fNuPaths[i].density / kU;
    double l = fNuPaths[i].length * 1e5; //l in km -> cm
    double p_trans = exp( -l * n * xsec );
    p_trans_vec.push_back(p_trans); //probability of NO absorption!
  }

  //returns probability of transmission
  return accumulate(p_trans_vec.begin(),
                    p_trans_vec.end(),
                    1.,
                    multiplies<double>());

}

//.............................................................................
void Absorption::SetPath(vector<NuPath> paths)
{
  fNuPaths = paths;
}
