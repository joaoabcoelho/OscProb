
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

  double get_tpi(){
    if(!count) count = 1;
    return time() / count;
  }

  // Print performance metric
  void Print(double ref_tpi=0){
    double tpi = get_tpi();
    string scale = "s";
    if(tpi<1){ tpi *= 1e3; scale = "ms"; }
    if(tpi<1){ tpi *= 1e3; scale = "µs"; }
    if(tpi<1){ tpi *= 1e3; scale = "ns"; }
    TString stpi = TString::Format(tpi>9.95 ? "%.0f" : "%.1f",tpi);
    if(stpi.Length()<3) stpi = " " + stpi;
    cout << "Performance = " << stpi << " " << scale << "/iteration";
    if(ref_tpi) printf(" (%.1fx)", get_tpi()/ref_tpi);
    cout << endl;
  }

  // Member attributes
  hrc::time_point start; // start time
  int count; // number of iterations

};


//.............................................................................
///
/// Test the speed performance of a given PMNS object
///
double TimeTest(OscProb::PMNS_Base* p, string model, int max_length, double ref_tpi=0){

  // Use a PremModel to make paths
  // through the earth
  OscProb::PremModel prem;

  // Chose an angle for the neutrino
  // and fill the paths with cosTheta
  // e.g. cosTheta = -1 (vertical up-going)

  // Give the path to the PMNS object
  // and get the probability

  // Define energy sample points
  int nbins = 100;
  vector<double> xbins = GetLogAxis(nbins, 1, 100);

  // Time the AvgProb performance for at most 1s
  TimeIt time;
  for(int k=0; k<100 && time.time()<1; k++){
    p->ClearCache(); // Clear cache at each iteration
    // Compute oscillogram
    for(double cosZ=-0.95; cosZ<1; cosZ+=0.1){
      prem.FillPath(cosZ);
      p->SetPath(prem.GetNuPath());
      for(int i=0; i<nbins; i++){
        double energy = 0.5*(xbins[i] + xbins[i+1]);
        double dE = xbins[i+1] - xbins[i];
        p->AvgProbMatrix(2,3, energy, dE);
        time.count++;
      }
    }
  }

  // Print performance estimate
  if(!ref_tpi){
    cout << string(20, '=') << " Reference " << string(20, '=') << endl;
  }
  cout << "PMNS_" << model << ": "
       << string(max_length - model.size(), ' ');
  time.Print(ref_tpi);
  if(!ref_tpi) cout << string(51, '=') << endl;

  return time.get_tpi();

}

//.............................................................................
///
/// Run the time test over all models
///
int StressTest(vector<string> models = {}){

  if(!models.size()) models = GetListOfModels();

  // Keep track of the longest name
  int max_length = 4;
  for(int i=0; i<models.size(); i++){
    max_length = max(int(models[i].size()), max_length);
  }

  OscProb::PMNS_Base* p0 = GetModel("Fast");
  double ref_tpi = TimeTest(p0, "Fast", max_length);

  for(int i=0; i<models.size(); i++){

    if(models[i]=="Fast") continue;

    OscProb::PMNS_Base* p = GetModel(models[i]);

    TimeTest(p, models[i], max_length, ref_tpi);

  }

  return 0;

}
