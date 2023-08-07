
#include <iostream>

#include "Utils.h"

#include "PremModel.h"
#include "../tutorial/SetNiceStyle.C"

using namespace std;

//.............................................................................
///
/// Run benchmark test of AvgProb
///
void StressTest(OscProb::PMNS_Base* p, string model){

  // Use a PremModel to make paths
  // through the earth
  OscProb::PremModel prem;

  // Chose an angle for the neutrino
  // and fill the paths with cosTheta
  // e.g. cosTheta = -1 (vertical up-going)
  prem.FillPath(-1);

  // Give the path to the PMNS object
  // and get the probability
  p->SetPath(prem.GetNuPath());

  // Define energy sample points
  int nbins = 100;
  vector<double> xbins = GetLogAxis(nbins, 1, 100);

  cout << "PMNS_" << model << ": Starting..." << endl;

  for(int k=0; k<1; k++){
    p->ClearCache(); // Clear cache at each iteration
    // Compute AvgProb for the energy range
    for(int i=0; i<nbins; i++){
      double energy = 0.5*(xbins[i] + xbins[i+1]);
      double dE = xbins[i+1] - xbins[i];
      p->AvgProb(1,0, energy, dE);
    }
  }

  cout << "PMNS_" << model << ": Done!" << endl;

}

//.............................................................................
///
/// Run the time test over all models
///
int main(int argc, char* argv[]){

  string model = "Fast";

  if(argc>1) model = argv[1];

  OscProb::PMNS_Base* p = GetModel(model);

  StressTest(p, model);

  return 0;

}
