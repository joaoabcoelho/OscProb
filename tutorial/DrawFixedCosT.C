
#include "SetNiceStyle.C"
#include "LoadOscProb.C"

#ifndef __CINT__
#include "PremModel.h"
#include "PMNS_Fast.h"
bool isCINT = false;
#else
bool isCINT = true;
#endif

void DrawFixedCosT(double cosT = -0.7, bool isNuBar = false){

  // Load some nice styles
  SetNiceStyle();

  // Load OscProb library
  if(isCINT) LoadOscProb();

  // Probability Calculator
  OscProb::PMNS_Fast p;
  
  // Set neutrino or antineutrino
  p.SetIsNuBar(isNuBar);
  
  // PREM Model
  OscProb::PremModel prem;
  
  // Fill path for cosT
  prem.FillPath(cosT);
  
  // Give path to calculator
  p.SetPath(prem.GetNuPath());
  
  // Make some histograms
  int nbins = 1000;
  double xmin = 1;
  double xmax = 20;
  
  TH1D* hMuMu_NH = new TH1D("","",nbins,xmin,xmax);
  TH1D* hEMu_NH = new TH1D("","",nbins,xmin,xmax);
  TH1D* hMuMu_IH = new TH1D("","",nbins,xmin,xmax);
  TH1D* hEMu_IH = new TH1D("","",nbins,xmin,xmax);

  double dm31 = 2.5e-3;

  // Fill histograms
  for(int i=1; i<=nbins; i++){
  
    double energy = hMuMu_NH->GetBinCenter(i);
    
    // Set NH
    p.SetDm(3, dm31);
    
    // Fill NH 
    hMuMu_NH->SetBinContent(i, p.Prob(1,1, energy));
    hEMu_NH->SetBinContent(i, p.Prob(0,1, energy));

    // Set IH 
    p.SetDm(3, -dm31 + 7.52e-5);
    
    // Fill IH
    hMuMu_IH->SetBinContent(i, p.Prob(1,1, energy));
    hEMu_IH->SetBinContent(i, p.Prob(0,1, energy));
  
  }
  
  // Set some nice styles
  SetHist(hMuMu_NH, kBlue);
  SetHist(hEMu_NH, kRed);
  SetHist(hMuMu_IH, kBlue);
  SetHist(hEMu_IH, kRed);
  
  // Make IH dashed
  hMuMu_IH->SetLineStyle(7);
  hEMu_IH->SetLineStyle(7);
  
  // Set axis titles
  if(isNuBar) hMuMu_NH->SetTitle(";Neutrino Energy (GeV);P(#bar{#nu}_{x}#rightarrow#bar{#nu}_{#mu})");
  else        hMuMu_NH->SetTitle(";Neutrino Energy (GeV);P(#nu_{x}#rightarrow#nu_{#mu})");
  
  // Set y range
  hMuMu_NH->GetYaxis()->SetRangeUser(0,1);

  // Make a long canvas
  MakeLongCanvas();
  
  // Draw everything
  hMuMu_NH->Draw("curv");
  hEMu_NH->Draw("curv same");
  hMuMu_IH->Draw("curv same");
  hEMu_IH->Draw("curv same");

  // Print cosT in canvas
  MiscText(0.8, 0.88, 0.04, TString::Format("cos#theta_{z} = %0.1f", cosT) );
  
  TLegend* leg = new TLegend(0.7,0.6,0.9,0.8);
  
  if(isNuBar){

    leg->AddEntry(hMuMu_NH, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{#mu}) - NH", "l");
    leg->AddEntry(hMuMu_IH, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{#mu}) - IH", "l");
    leg->AddEntry(hEMu_NH, "P(#bar{#nu}_{e}#rightarrow#bar{#nu}_{#mu}) - NH", "l");
    leg->AddEntry(hEMu_IH, "P(#bar{#nu}_{e}#rightarrow#bar{#nu}_{#mu}) - IH", "l");

  }
  else{
  
    leg->AddEntry(hMuMu_NH, "P(#nu_{#mu}#rightarrow#nu_{#mu}) - NH", "l");
    leg->AddEntry(hMuMu_IH, "P(#nu_{#mu}#rightarrow#nu_{#mu}) - IH", "l");
    leg->AddEntry(hEMu_NH, "P(#nu_{e}#rightarrow#nu_{#mu}) - NH", "l");
    leg->AddEntry(hEMu_IH, "P(#nu_{e}#rightarrow#nu_{#mu}) - IH", "l");
  
  }
  
  SetLeg(leg);
  
  leg->Draw();
  
}
