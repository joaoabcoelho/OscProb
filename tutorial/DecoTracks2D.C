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

void DecoTracks2D(int mh = 1, double power = 0, double gamma = 1e-23){
  
    // Name files
  std::string MH = "NH";
  std::string POWER = "0";

  if(mh == -1) MH = "IH";
  
  if(power == -2) POWER = "Min2";
  else if(power == -1) POWER = "Min1";
  else if(power == 0) POWER = "0";
  else if(power == 1) POWER = "1";
  else POWER = "2";
  
  std::string filename = "Tracks_" + MH + "_Power" + POWER + "_gamma2.txt";
  
  std::ofstream output_file(filename,std::ios::out);
  
  // Load OscProb library
  if(isCINT) LoadOscProb();
  
  // Probability calculator with decoherence
  OscProb::PMNS_Deco p;
  
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
  

    // Set cos(theta_z) from bin center
    double cosZ = h2->GetYaxis()->GetBinCenter(ct);
    
    // Skip if cosZ is unphysical  
    if(cosZ < -1 || cosZ > 1) continue;
    
    // Fill paths from PREM model
    prem.FillPath(cosZ);
      
    // Set paths in OscProb  
    p.SetPath(prem.GetNuPath());

    // Set muon nuetrinos (tracks) as final states
    int flvf = 1;
    
    // Loop of E  
    for(int e=1; e<=nbinsx; e++){
  
      // Set E from bin center
      double energy  = h2->GetXaxis()->GetBinCenter(e);
      
      // Initialize probabilities 
      double prob_Std{0};
      double prob_Atm{0};
      double prob_Sol1{0};
      double prob_Sol2{0};
      
      // Loop over initial flavour and nu or nubar
      for(int flvi = 0; flvi<2; flvi++){
        for(int nunubar = -1; nunubar<2; nunubar+=2){
          
          // Define some basic weights for nue/numu and nu/nubar
          double weight = (0.75 + 0.25*nunubar) * (0.5 + 0.5*flvi);

          // Add probabilities from OscProb
          p.SetIsNuBar(nunubar <= 0);
          
          // Set standard oscillations
          p.SetGamma(3, 0);
          p.SetGamma(2, 0);
          p.SetPower(power); 
          prob_Std += weight*p.Prob(flvi, flvf, energy);
          //cout << "Prob: "  << p.Prob(flvi, flvf, energy) << endl;
          //cout << "Prob total: " <<  prob_Std << endl;
          
          // Set Atmospheric limit
          p.SetGamma(2, 0);
          p.SetGamma(3, gamma);
          p.SetPower(power);
    
          prob_Atm += weight*p.Prob(flvi, flvf, energy);
      
          // Set Solar limit 1
          p.SetGamma(3, gamma);
          p.SetGamma(2, gamma);
          p.SetPower(power);
    
          prob_Sol1 += weight*p.Prob(flvi, flvf, energy);

          // Set Solar limit 2
          p.SetGamma(3, 0);
          p.SetGamma(2, gamma);
          p.SetPower(power);
    
          prob_Sol2 += weight*p.Prob(flvi, flvf, energy);
          }
      }
      //cout << "End bin" << endl;
      output_file << energy << "\t" << cosZ << "\t" << prob_Std << "\t" << prob_Atm << "\t" << prob_Sol1 << "\t" << prob_Sol2 << endl;
          
    }// energy loop
  }// cosZ loop
}