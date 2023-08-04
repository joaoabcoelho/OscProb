
#include "Utils.h"

using hrc = std::chrono::high_resolution_clock;
using std::chrono::nanoseconds;
using std::chrono::duration_cast;

//.............................................................................
///
/// Timing performance utility
///
struct TimeIt {

  // Constructor
  TimeIt() { reset(); }

  // Reset the clock
  void reset(){ count = 0; start = hrc::now(); }

  // Return elapsed time in seconds
  double time(){
    return 1e-9 * duration_cast<nanoseconds>(hrc::now() - start).count();
  }

  // Print performance metric
  void Print(){
    double tpi = time() / count;
    string scale = "s";
    if(tpi<1){ tpi *= 1e3; scale = "ms"; }
    if(tpi<1){ tpi *= 1e3; scale = "µs"; }
    if(tpi<1){ tpi *= 1e3; scale = "ns"; }
    cout << "Performance = " 
         << TString::Format(tpi>10 ? "%.0f" : "%.1f",tpi) 
         << " " << scale << "/iteration" << endl;
  }

  // Member attributes
  hrc::time_point start; // start time
  int count; // number of iterations

};


//.............................................................................
///
/// Test the speed performance of a given PMNS object
///
void TimeTest(OscProb::PMNS_Base* p, string model, int max_length){

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

  // Time the AvgProb performance for at most 1s
  TimeIt time;
  for(int k=0; k<100 && time.time()<1; k++){
    p->ClearCache(); // Clear cache at each iteration
    // Compute AvgProb for the energy range
    for(int i=0; i<=nbins; i++){
      double energy = xbins[i];
      p->AvgProb(1,0, energy, 0.2*energy);
      time.count++;
    }
  }
  
  // Print performance estimate
  cout << "PMNS_" << model << ": " 
       << string(max_length - model.size(), ' ');
  time.Print();

}

//.............................................................................
///
/// Run the time test over all models
///
int StressTest(){

  vector<string> models = {"Fast", "Iter", "Sterile", "NSI",
                           "Deco", "Decay", "LIV", "SNSI"};

  // Keep track of the longest name
  int max_length = 0;
  for(int i=0; i<models.size(); i++){
    max_length = max(int(models[i].size()), max_length);
  }

  for(int i=0; i<models.size(); i++){

    OscProb::PMNS_Base* p = GetModel(models[i]);

    TimeTest(p, models[i], max_length);

  }

  return 0;
  
}
