#include <iostream>

#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TSystem.h"

#ifndef __CINT__
#include "PremModel.h"
#include "PMNS_Fast.h"
#include "PMNS_Deco.h"
bool isCINT = false;
#else
bool isCINT = true;
#endif

// Macro to load OscProb library
#include "LoadOscProb.C"

// Some functions to make nice plots
#include "SetNiceStyle.C"

void DecoProb2D(bool isNuBar = true, int mh = -1, int flvi = 0, int flvf = 0, double power = 0, double gamma = 1e-23){
  
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
  
  std::string filename = "P2D_" + FLVI + FLVF + "_" + MH + "_Power" + POWER + ".txt";
  
  std::ofstream output_file(filename,std::ios::out);
  
  // Load OscProb library
  if(isCINT) LoadOscProb();
  
  // Probability calculator with decoherence
  OscProb::PMNS_Deco p;
  
  // Set neutrino or antineutrino
  p.SetIsNuBar(isNuBar);
  
  // PREM Model
  OscProb::PremModel prem;
  
  // Set binning 
  int nbinsx = 120; // 200 previuosly
  int nbinsy = 120; // 100 previously
  double xmin = 5;
  double xmax = 400;
  vector<double> xbins = GetLogAxis(nbinsx, xmin, xmax);
  
  TH2D* h2 = new TH2D("","",nbinsx,&xbins[0],nbinsy,-1,0);

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
  // Loop over cos(theta_z) and E
  for(int ct=0; ct<=nbinsy; ct++){
  
    // Initialize probabilities 
    double prob_Std{};
    double prob_Atm{};
    double prob_Sol1{};
    double prob_Sol2{};

    // Set cos(theta_z) from bin center
    double cosZ = h2->GetYaxis()->GetBinCenter(ct);
    
    // Skip if cosZ is unphysical  
    if(cosZ < -1 || cosZ > 1) continue;
    
    // Fill paths from PREM model
    prem.FillPath(cosZ);
      
    // Set paths in OscProb  
    p.SetPath(prem.GetNuPath());

    // Loop of E  
    for(int e=1; e<=nbinsx; e++){
  
      // Set E from bin center
      double energy  = h2->GetXaxis()->GetBinCenter(e);
      
      // Set standard oscillations
      p.SetGamma(3, 0);
      p.SetGamma(2, 0);
      p.SetPower(power); 
         
      prob_Std = p.Prob(flvi, flvf, energy);
      
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
    
      output_file << energy << "\t" << cosZ << "\t" << prob_Std << "\t" << prob_Atm << "\t" << prob_Sol1 << "\t" << prob_Sol2 << endl;
          
    }// energy loop
  }// cosZ loop
}



