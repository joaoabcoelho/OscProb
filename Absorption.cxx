#include "Absorption.h"
#include "NuPath.h"
#include "PremModel.h"
//#include "../EarthModel/EarhtModel.h"

#include <vector>
#include <numeric>
#include <math.h>
#include <iostream>
//extern ostream cout;

using std::cout;
using std::endl;
using std::accumulate;
using std::multiplies;
//using namespace std;
using namespace OscProb;

Absorption::Absorption(){}

Absorption::~Absorption(){}

const double Absorption::kNA     = 6.022140857e+23;             // Avogadro constant (N_A)
const double Absorption::kU      = 1.660539066e-24;             // atomic mass unit [g] 

double Absorption::Trans(double xsec){
  //cout << "test" << std::endl;
  //return 0.1;
  std::vector<double> p_trans_vec;
  //cout << "test" << std::endl;
  for(int i=0; i<int(fNuPaths.size()); i++){
    double n = fNuPaths[i].density/(kU);
    double l = fNuPaths[i].length * 100000; //l in km -> cm
    double p_trans = exp( - l* (n*xsec) ); 
    //cout << 'l' << '\t' << 'n' << '\t' << "xsec" << '\t' << "p_trans" << endl;
    //cout << l << '\t' << n << '\t' << xsec << '\t' << p_trans << endl;
    p_trans_vec.push_back(p_trans); //probability of NO absorption!
  }
  return 1.-accumulate(p_trans_vec.begin(),p_trans_vec.end(),1.,multiplies<double>()); //returns probability of absorption
}

void Absorption::SetPath(std::vector<OscProb::NuPath> paths){
  fNuPaths = paths;
}
