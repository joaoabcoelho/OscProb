
#include "PremModel.h"
#include "PMNS_Deco.h"

#include "SetNiceStyle.C"
#include "prem_default.hpp"

void DrawFixedBaseline_Deco(double L = 1300, bool isNuBar = false){

  // Load some nice styles
  SetNiceStyle();

  // Probability Calculator
  OscProb::PMNS_Deco p;

  // Set neutrino or antineutrino
  p.SetIsNuBar(isNuBar);
  
  // PREM Model
  OscProb::PremModel prem(PREM_DIR + "/prem_lbl.txt");
  
  // Fill path for L
  prem.FillPath(prem.GetCosT(L));
  
  // Give path to calculator
  p.SetPath(prem.GetNuPath());
  
  // Make some histograms
  int nbins = 1000;
  double xmin = 0.5;
  double xmax = 100;
  vector<double> xbins = GetLogAxis(nbins, xmin, xmax);
  
  TH1D* hMuMu_Std = new TH1D("","",nbins,&xbins[0]);
  TH1D* hMuE_Std = new TH1D("","",nbins,&xbins[0]);
  TH1D* hMuMu_Deco = new TH1D("","",nbins,&xbins[0]);
  TH1D* hMuE_Deco = new TH1D("","",nbins,&xbins[0]);

  // Fill histograms
  for(int i=1; i<=nbins; i++){
  
    double energy = hMuMu_Std->GetBinCenter(i);
    
    // Set gamma value  
    p.SetGamma(3, 0);
  
    // Fill std. osc. 
    hMuMu_Std->SetBinContent(i, p.Prob(1,1, energy));
    hMuE_Std->SetBinContent(i, p.Prob(1,0, energy));

    // Set gamma value
    p.SetGamma(3, 2.3e-23);
    
    // Fill decoherence
    hMuMu_Deco->SetBinContent(i, p.Prob(1,1, energy));
    hMuE_Deco->SetBinContent(i, p.Prob(1,0, energy));
  
  }
  
  // Set some nice styles
  SetHist(hMuMu_Std, kBlue);
  SetHist(hMuE_Std, kGreen+1);
  SetHist(hMuMu_Deco, kRed);
  SetHist(hMuE_Deco, kMagenta+2);
  
  // Make lines thicker
  hMuMu_Std->SetLineWidth(3);
  hMuMu_Deco->SetLineWidth(3);
  hMuE_Std->SetLineWidth(3);
  hMuE_Deco->SetLineWidth(3);
  
  // Make deco dashed
  hMuMu_Deco->SetLineStyle(9);
  hMuE_Deco->SetLineStyle(9);
  
  // Set axis titles
  if(isNuBar) hMuMu_Std->SetTitle(";Neutrino Energy (GeV);P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{x})");
  else        hMuMu_Std->SetTitle(";Neutrino Energy (GeV);P(#nu_{#mu}#rightarrow#nu_{x})");
  
  // Set y range
  hMuMu_Std->GetYaxis()->SetRangeUser(0,1);

  // Make a long canvas
  MakeLongCanvas();
  
  // Draw everything
  hMuMu_Std->Draw("curv");
  hMuE_Std->Draw("curv same");
  hMuMu_Deco->Draw("curv same");
  hMuE_Deco->Draw("curv same");

  // Print L in canvas
  MiscText(0.75, 0.85, 0.04, TString::Format("L = %0.0f km", L) );
  
  TLegend* leg = new TLegend(0.7,0.6,0.9,0.8);
  
  if(isNuBar){

    leg->AddEntry(hMuMu_Std, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{#mu}) - Std.", "l");
    leg->AddEntry(hMuMu_Deco, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{#mu}) - Decoh.", "l");
    leg->AddEntry(hMuE_Std, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{e}) - Std.", "l");
    leg->AddEntry(hMuE_Deco, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{e}) - Decoh.", "l");

  }
  else{
  
    leg->AddEntry(hMuMu_Std, "P(#nu_{#mu}#rightarrow#nu_{#mu}) - Std.", "l");
    leg->AddEntry(hMuMu_Deco, "P(#nu_{#mu}#rightarrow#nu_{#mu}) - Decoh.", "l");
    leg->AddEntry(hMuE_Std, "P(#nu_{#mu}#rightarrow#nu_{e}) - Std.", "l");
    leg->AddEntry(hMuE_Deco, "P(#nu_{#mu}#rightarrow#nu_{e}) - Decoh.", "l");
  
  }
  
  SetLeg(leg);
  
  leg->Draw();
  
  gPad->SetLogx();
  
  hMuMu_Std->GetXaxis()->SetMoreLogLabels();
  
}
