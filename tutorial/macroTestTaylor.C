
#include "PremModel.h"
#include "PMNS_Fast.h"
#include "PMNS_TaylorExp.h"

#include "SetNiceStyle.C"


void macroTestTaylor(double cosT = -0.7, bool isNuBar = false){

  // Load some nice styles
  SetNiceStyle();

  // Probability Calculator
  OscProb::PMNS_Fast f;
  OscProb::PMNS_TaylorExp t;

  // Set neutrino or antineutrino
  f.SetIsNuBar(isNuBar);
  t.SetIsNuBar(isNuBar);

  // PREM Model
  OscProb::PremModel prem;    

  // Fill path for cosT
  prem.FillPath(cosT);

  // Give path to calculator
  f.SetPath(prem.GetNuPath());
  t.SetPath(prem.GetNuPath());

  // Make some histograms
  int nbins = 10000;
  double xmin = 1;
  double xmax = 20;
  double widthBin = (xmax-xmin) / nbins;

  TH1D* hMuMu_fast = new TH1D("","",nbins,xmin,xmax);
  TH1D* hEMu_fast = new TH1D("","",nbins,xmin,xmax);
  TH1D* hMuMu_taylor = new TH1D("","",nbins,xmin,xmax);
  TH1D* hEMu_taylor = new TH1D("","",nbins,xmin,xmax);

  double dm31 = 2.5e-3;

  // Fill histograms
  for(int i=1; i<=nbins; i++){

    double energy = hMuMu_fast->GetBinCenter(i);

    // Set NH
    f.SetDm(3, dm31);
    t.SetDm(3, dm31);

    // Fill NH
    hMuMu_fast->SetBinContent(i, f.Prob(1,1, energy));
    hEMu_fast->SetBinContent(i, f.Prob(0,1, energy));

    // Fill IH
    hMuMu_taylor->SetBinContent(i, t.avgProbTaylor(1,1, energy,widthBin));
    hEMu_taylor->SetBinContent(i, t.avgProbTaylor(0,1, energy,widthBin));

  }

  // Set some nice styles
  SetHist(hMuMu_fast, kPink);
  SetHist(hEMu_fast, kYellow);
  SetHist(hMuMu_taylor, kBlack);
  SetHist(hEMu_taylor, kGreen);

  // Make IH dashed
  hMuMu_taylor->SetLineStyle(7);
  hEMu_taylor->SetLineStyle(7);

  // Set axis titles
  if(isNuBar) hMuMu_fast->SetTitle(";Neutrino Energy (GeV);P(#bar{#nu}_{x}#rightarrow#bar{#nu}_{#mu})");
  else        hMuMu_fast->SetTitle(";Neutrino Energy (GeV);P(#nu_{x}#rightarrow#nu_{#mu})");

  // Set y range
  hMuMu_fast->GetYaxis()->SetRangeUser(0,1);

  // Make a long canvas
  MakeLongCanvas();

  // Draw everything
  hMuMu_fast->Draw("curv");
  hEMu_fast->Draw("curv same");
  hMuMu_taylor->Draw("curv same");
  hEMu_taylor->Draw("curv same");

  // Print cosT in canvas
  MiscText(0.8, 0.88, 0.04, TString::Format("cos#theta_{z} = %0.1f", cosT) );

  TLegend* leg = new TLegend(0.7,0.6,0.9,0.8);

  if(isNuBar){

    leg->AddEntry(hMuMu_fast, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{#mu}) - fast", "l");
    leg->AddEntry(hMuMu_taylor, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{#mu}) - taylor", "l");
    leg->AddEntry(hEMu_fast, "P(#bar{#nu}_{e}#rightarrow#bar{#nu}_{#mu}) - fast", "l");
    leg->AddEntry(hEMu_taylor, "P(#bar{#nu}_{e}#rightarrow#bar{#nu}_{#mu}) - taylor", "l");

  }
  else{

    leg->AddEntry(hMuMu_fast, "P(#nu_{#mu}#rightarrow#nu_{#mu}) - fast", "l");
    leg->AddEntry(hMuMu_taylor, "P(#nu_{#mu}#rightarrow#nu_{#mu}) - taylor", "l");
    leg->AddEntry(hEMu_fast, "P(#nu_{e}#rightarrow#nu_{#mu}) - fast", "l");
    leg->AddEntry(hEMu_taylor, "P(#nu_{e}#rightarrow#nu_{#mu}) - taylor", "l");

  }

  SetLeg(leg);

  leg->Draw();



}
