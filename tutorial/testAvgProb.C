#ifndef __CINT__
#include "PMNS_Fast.h"
bool isCINT = false;
#else
bool isCINT = true;
#endif

// Some functions to make nice plots
#include "SetNiceStyle.C"

// Macro to load OscProb library
#include "LoadOscProb.C"

// Check the accuracy of the
// oscillation averaging function
void testAvgProb(){

  // Set nice overall style
  SetNiceStyle();

  // Load OscProb
  if(isCINT) LoadOscProb();
  
  // Get a PMNS object
  OscProb::PMNS_Fast p;

  // Set the baseline through the earth  
  double L = 2*6368 + 18;
  p.SetLength(L);
  
  // Define some fine and coarse binnings
  int navg = 20;
  int nbins = navg * 100;
  
  // Lots of histograms
  TH1D* h1 = new TH1D("","",nbins,0,3);
  TH1D* h2 = new TH1D("","",navg,0,3);
  TH1D* h3 = new TH1D("","",navg,0,3);
  TH1D* h4 = new TH1D("","",navg,0,3);
  TH1D* h5 = new TH1D("","",navg,0,3);

  // Do some fine binning and uniform sampling
  for(int i=1; i<=nbins; i++){
  
    double minE  = pow(10, h1->GetBinLowEdge(i));
    double maxE  = pow(10, h1->GetBinLowEdge(i+1));
    
    double E = 0.5 * (minE + maxE);
    double dE = (maxE - minE);
    
    int abin = h3->FindFixBin(log10(E));
    
    double prob = p.Prob(1,1, E);

    h1->SetBinContent(i, prob);
    h3->AddBinContent(abin, prob * dE);

  }

  // Do the AvgProb sampling and bin center sampling
  for(int i=1; i<=navg; i++){
  
    double minE  = pow(10, h2->GetBinLowEdge(i));
    double maxE  = pow(10, h2->GetBinLowEdge(i+1));
    
    double E  = (maxE + minE) / 2;
    double dE = (maxE - minE);

    h2->SetBinContent(i, p.AvgProb(1,1, E, dE));
    
    h3->SetBinContent(i, h3->GetBinContent(i) / dE);
    h4->SetBinContent(i, p.Prob(1,1, E));
    
    // Store the number of samples per bin
    if(minE < 0.1*E) minE = 0.1 * E;

    double LoE  = (L/minE + L/maxE)/2;
    double dLoE = (L/minE - L/maxE);
    
    int nsamples = (p.GetSamplePoints(LoE, dLoE)).size();
    
    h5->SetBinContent(i, nsamples);

  }

  // Print the sampling statistics
  int sumsample = h5->Integral();
  cout << "Total samples needed: " << sumsample << endl;
  cout << "Mean samples per bin: " << double(sumsample) / navg << endl;

  // Make a long canvas
  MakeLongCanvas();
  
  // Set some nice histograms
  SetHist(h1, kBlack);
  SetHist(h2, kBlue);
  SetHist(h3, kGreen);
  SetHist(h4, kRed);
  
  // Change line styles
  h3->SetLineStyle(7);

  h2->SetLineWidth(3);
  h3->SetLineWidth(3);
  h4->SetLineWidth(5);
  
  // The axis titles
  h1->SetTitle(";Log10[Neutrino Energy (GeV)];P(#nu_{#mu}#rightarrow#nu_{#mu})");
  
  // Draw different samplings
  h1->DrawCopy("curv");
  h4->DrawCopy("hist same ][");
  h2->DrawCopy("hist same ][");
  h3->DrawCopy("hist same ][");

  // Make a new canvas
  gPad->DrawClone();
  
  // Set a nice histogram
  SetHist(h5, kBlue);

  // The axis titles
  h5->SetTitle(";Log10[Neutrino Energy (GeV)];# of samples needed");
  
  // Set the plot scale
  double smax = 1.2 * h5->GetMaximum();
  h5->GetYaxis()->SetRangeUser(0, smax);
  h1->Scale(smax);

  // Draw the number of sample points
  h5->DrawCopy("hist ][");

  // Draw the scaled oscillation
  // probability behind it
  h1->DrawCopy("hist ][ same");
  h5->DrawCopy("hist ][ same");
  
}
