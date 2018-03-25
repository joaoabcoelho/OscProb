
// Macro to load OscProb library
#include "LoadOscProb.C"

void simpleExamples(){
 
  // Load the library.
  // You can also include this in your rootrc file
  // The library must be in your path
  LoadOscProb();



  // Initialize your objects
  OscProb::PMNS_Fast myPMNS;
  
  // Other examples
  OscProb::PMNS_NSI  myNSI; // NSI extension
  OscProb::PMNS_Sterile mySterile(5); // 3+2 Model
  


  // Set oscillation paramaters
  // By default, they are set to
  // PDG values with SetStdPars()
  // in the constructor
  myPMNS.SetDm(3, -2.4e-3);   // Set Dm31 in eV^2
  myPMNS.SetAngle(1,3, 0.15); // Set Theta13 in radians
  myPMNS.SetDelta(1,3, 1.57); // Set Delta13 in radians

  // Set other quantities
  myPMNS.SetDensity(2.5); // Set the matter density in g/cm^3
  myPMNS.SetLength(1300); // Set baseline in km
  myPMNS.SetEnergy(3.0);  // Set neutrino energy in GeV
  
  // Say whether you want antineutrinos
  // Default is false, i.e. neutrinos
  myPMNS.SetIsNuBar(true); 
  


  // Calculate probability of 
  // numubar to nuebar transition
  double Pme = myPMNS.Prob(1,0);
  
  // Calculate again at 10 GeV
  Pme = myPMNS.Prob(1,0, 10.0);
  
  // Calculate average probability
  // between 9.5 and 10.5 GeV
  // Assumes constant event rate
  // across energy range
  Pme = myPMNS.AvgProb(1,0, 10.0, 1.0);



  // Set 3 paths in order
  myPMNS.SetPath(1000, 2.50); // Set initial path with L = 1000 km and rho = 2.50 g/cm^3
  myPMNS.AddPath(2000, 4.25); // Add a second path with L = 2000 km and rho = 4.25 g/cm^3
  myPMNS.AddPath(1500, 6.00); // Add a third path with L = 1500 km and rho = 6.00 g/cm^3

  // Get numubar to nuebar probability
  Pme = myPMNS.Prob(1,0);
  
  
  
  // Use a PremModel to make paths 
  // through the earth
  OscProb::PremModel prem;
  
  // Chose an angle for the neutrino
  // and fill the paths with cosTheta
  // e.g. cosTheta = -1 (vertical up-going)
  prem.FillPath(-1);
  
  // Give the path to the PMNS object
  // and get the probability
  myPMNS.SetPath(prem.GetNuPath());
  Pme = myPMNS.Prob(1,0);


}
