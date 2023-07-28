
#include "PremModel.h"
#include "PMNS_Fast.h"
#include "PMNS_Decay.h"

#include "SetNiceStyle.C"

void DrawDecayCompare(double cosT = -0.7, bool isNuBar = false, double alpha3=1e-4){

  // Load some nice styles
  SetNiceStyle();

  // Probability Calculator
  OscProb::PMNS_Fast p;
  OscProb::PMNS_Decay d;
  // Set neutrino or antineutrino
  p.SetIsNuBar(isNuBar);
  d.SetIsNuBar(isNuBar);
  // Set alpha3 value
  d.SetAlpha3(alpha3);
  // PREM Model
  OscProb::PremModel prem;
  
  // Fill path for cosT
  prem.FillPath(cosT);
  
  // Give path to calculator
  p.SetPath(prem.GetNuPath());
  d.SetPath(prem.GetNuPath());
  // Make some histograms
  int nbins = 1000;
  double xmin = 1;
  double xmax = 20;
  
  TH1D* hMuMu_Fast = new TH1D("","",nbins,xmin,xmax);
  TH1D* hEMu_Fast = new TH1D("","",nbins,xmin,xmax);
  TH1D* hMuMu_Decay = new TH1D("","",nbins,xmin,xmax);
  TH1D* hEMu_Decay = new TH1D("","",nbins,xmin,xmax);

  double dm31 = 2.5e-3;
  p.SetDm(3, dm31);
  d.SetDm(3,dm31);
  // Fill histograms
  for(int i=1; i<=nbins; i++){
  
    double energy = hMuMu_Fast->GetBinCenter(i);
    
    
    // Fill Fast (No Decay)
    hMuMu_Fast->SetBinContent(i, p.Prob(1,1, energy));
    hEMu_Fast->SetBinContent(i, p.Prob(0,1, energy));

    // Fill Decay
    hMuMu_Decay->SetBinContent(i, d.Prob(1,1, energy));
    hEMu_Decay->SetBinContent(i, d.Prob(0,1, energy));
  
  }
  
  // Set some nice styles
  SetHist(hMuMu_Fast, kBlue);
  SetHist(hEMu_Fast, kRed);
  SetHist(hMuMu_Decay, kBlue);
  SetHist(hEMu_Decay, kRed);
  
  // Make IH dashed
  hMuMu_Decay->SetLineStyle(7);
  hEMu_Decay->SetLineStyle(7);
  
  // Set axis titles
  if(isNuBar) hMuMu_Fast->SetTitle(";Neutrino Energy (GeV);P(#bar{#nu}_{x}#rightarrow#bar{#nu}_{#mu})");
  else        hMuMu_Fast->SetTitle(";Neutrino Energy (GeV);P(#nu_{x}#rightarrow#nu_{#mu})");
  
  // Set y range
  hMuMu_Fast->GetYaxis()->SetRangeUser(0,1);

  // Make a long canvas
  MakeLongCanvas();
  
  // Draw everything
  hMuMu_Fast->Draw("curv");
  hEMu_Fast->Draw("curv same");
  hMuMu_Decay->Draw("curv same");
  hEMu_Decay->Draw("curv same");

  // Print cosT in canvas
  MiscText(0.8, 0.88, 0.04, TString::Format("cos#theta_{z} = %0.1f", cosT) );
  
  TLegend* leg = new TLegend(0.7,0.6,0.9,0.8);
  
  if(isNuBar){

    leg->AddEntry(hMuMu_Fast, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{#mu}) - No Decay", "l");
    leg->AddEntry(hMuMu_Decay, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{#mu}) - Decay", "l");
    leg->AddEntry(hEMu_Fast, "P(#bar{#nu}_{e}#rightarrow#bar{#nu}_{#mu}) - No Decay", "l");
    leg->AddEntry(hEMu_Decay, "P(#bar{#nu}_{e}#rightarrow#bar{#nu}_{#mu}) - Decay", "l");

  }
  else{
  
    leg->AddEntry(hMuMu_Fast, "P(#nu_{#mu}#rightarrow#nu_{#mu}) - No Decay", "l");
    leg->AddEntry(hMuMu_Decay, "P(#nu_{#mu}#rightarrow#nu_{#mu}) - Decay", "l");
    leg->AddEntry(hEMu_Fast, "P(#nu_{e}#rightarrow#nu_{#mu}) - No Decay", "l");
    leg->AddEntry(hEMu_Decay, "P(#nu_{e}#rightarrow#nu_{#mu}) - Decay", "l");
  
  }
  
  SetLeg(leg);
  
  leg->Draw();
  
}
