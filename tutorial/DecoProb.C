#include <string>

#include "SetNiceStyle.C"
#include "LoadOscProb.C"
#include "../prem_default.hpp"

#ifndef __CINT__
#include "PremModel.h"
#include "PMNS_Deco.h"
#include "TMath.h"
bool isCINT = false;
#else
bool isCINT = true;
#endif


void DecoProb(double cosZ = -1, bool isNuBar = false, int mh = 1, int flvi = 1, int flvf = 2, double power = 0, double gamma = 1e-23){

  // Name files
  std::string MH = "NH";
  std::string FLVI = "Mu";
  std::string FLVF = "Mu";
  std::string POWER = "0";

  if(isNuBar == true){
    switch(flvi){
      case 0: FLVI = "Ebar";
              break;
      case 1: FLVI = "Mubar";
              break;
      case 2: FLVI = "Taubar";
              break;
    }
  }
  else{
      switch(flvi){
      case 0: FLVI = "E";
              break;
      case 1: FLVI = "Mu";
              break;
      case 2: FLVI = "Tau";
              break;
      }
  }
  
  if(isNuBar == true){
    switch(flvf){
      case 0: FLVF = "Ebar";
              break;
      case 1: FLVF = "Mubar";
              break;
      case 2: FLVF = "Taubar";
              break;
    }
  }
  else{
      switch(flvf){
      case 0: FLVF = "E";
              break;
      case 1: FLVF = "Mu";
              break;
      case 2: FLVF = "Tau";
              break;
      }
  }
  
  if(power == -2) POWER = "Min2";
  else if(power == -1) POWER = "Min1";
  else if(power == 0) POWER = "0";
  else if(power == 1) POWER = "1";
  else POWER = "2";

  if(mh == -1) MH = "IH";
  
  std::string filename = "P_" + FLVI + FLVF + "_" + MH + "_Power" + POWER + "_Theta23_45deg.txt";
  
  std::ofstream output_file(filename,std::ios::out);


  // Load OscProb library
  if(isCINT) LoadOscProb();

  // Probability Calculator with decoherence 
  OscProb::PMNS_Deco p;

  // Set neutrino or antineutrino
  p.SetIsNuBar(isNuBar);
  
  // PREM Model
  OscProb::PremModel prem;
  
  // Fill path for given angle
  prem.FillPath(cosZ);
  
  // Give path to calculator
  p.SetPath(prem.GetNuPath());
  
  // Make some histograms
  int nbins = 1000;
  double xmin = 1;
  double xmax = 1000;
  vector<double> xbins = GetLogAxis(nbins, xmin, xmax);
  
  TH1D* hMuMu_Std = new TH1D("","",nbins,&xbins[0]);

  // Set oscillation parameters
  double dm31{};
  
  if(mh == -1){
    dm31 = -2.498e-3;
  }
  else{
    dm31 = 2.515e-3;
  }
  p.SetDm(3, dm31);
  p.SetDm(2, 7.42e-5);
  p.SetAngle(1, 3, 8.57  * TMath::DegToRad());
  p.SetAngle(1, 2, 33.44 * TMath::DegToRad());
  p.SetAngle(2, 3, 49.2  * TMath::DegToRad());

  // Fill histograms
  for(int i=1; i<=nbins; i++){
  
  // Create probabilities 
    double prob_Std{};
    double prob_Atm{};
    double prob_Sol1{};
    double prob_Sol2{};
    
    // Set energy from bin center
    double energy = hMuMu_Std->GetBinCenter(i);
    
    // Set standard oscillations
    p.SetGamma(3, 0);
    p.SetGamma(2, 0);
    p.SetPower(power);
    
    // Fill standard oscillation probability
    prob_Std = p.Prob(flvi,flvf, energy);
    
    // Set Atmospheric limit
    p.SetGamma(2, 0);
    p.SetGamma(3, gamma);
    p.SetPower(power);
    
    prob_Atm = p.Prob(flvi,flvf, energy);
 
    // Set Solar limit 1
    p.SetGamma(3, gamma);
    p.SetGamma(2, gamma);
    p.SetPower(power);
    
    prob_Sol1 = p.Prob(flvi,flvf, energy);

  // Set Solar limit 2
    p.SetGamma(3, 0);
    p.SetGamma(2, gamma);
    p.SetPower(power);
    
    prob_Sol2 = p.Prob(flvi,flvf, energy);
    
    output_file << energy << "\t" << prob_Std << "\t" << prob_Atm << "\t" << prob_Sol1 << "\t" << prob_Sol2 << endl;
  }
}