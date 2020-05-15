#include "Absorption.h"
#include "NuPath.h"
#include "PremModel.h"
//#include "../EarthModel/EarhtModel.h"

#include <vector>
#include <numeric>
#include <math.h>

using namespace OscProb;

//Absorption::Absorption();

const double Absorption::kNA     = 6.022140857e23;             // Avogadro constant (N_A)

double Absorption::Trans(double xsec, double E){
  return 1;
  std::vector<double> p_trans_vec;
  for(int i=0; i<int(fNuPaths.size()); i++){
    double u = 55.845; // atomic weight of Iron [g/mol] 
    double n = fNuPaths[i].density/(u*kNA);
    double l = fNuPaths[i].length;
    double p_trans = exp( - l* (n*xsec) ); 
    p_trans_vec.push_back(1 - p_trans); //probability of NO absorption!
  }
  return 1 - std::accumulate(p_trans_vec.begin(),p_trans_vec.end(),1,std::multiplies<double>()); //returns probability of absorption
}

void Absorption::SetPath(std::vector<OscProb::NuPath> paths){
  fNuPaths = paths;
}
