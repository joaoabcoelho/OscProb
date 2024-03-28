#include "PremModel.h"
#include "PMNS_Fast.h"
#include "PMNS_OQS.h"

#include "SetNiceStyle.C"
#include "prem_default.hpp"

void testOQS(double L = 300, bool isNuBar = false){

  // Load some nice styles
  SetNiceStyle();

  // Probability Calculator
  OscProb::PMNS_Fast p;
  OscProb::PMNS_OQS p2;

  // Set neutrino or antineutrino
  p.SetIsNuBar(isNuBar);
  p2.SetIsNuBar(isNuBar);

  // PREM Model
  OscProb::PremModel prem(PREM_DIR + "/prem_lbl.txt");

  /*
  // Fill path for L
  prem.FillPath(prem.GetCosT(L));

  //Give path to calculator
  p.SetPath(prem.GetNuPath());
  p2.SetPath(prem.GetNuPath());
  */

  
  double dens = 0;
  p.SetDensity(dens); // Set the matter density in g/cm^3       
  p.SetLength(L); // Set baseline in km                         

  p2.SetDensity(dens); // Set the matter density in g/cm^3         
  p2.SetLength(L); // Set baseline in km

  p.SetPath(L, dens);
  p2.SetPath(L, dens);
  

  // Make some histograms
  int nbins = 300;
  double xmin = 1;
  double xmax = 1.1;
  vector<double> xbins = GetLogAxis(nbins, xmin, xmax);

  
  TH1D* hMuMu_Std = new TH1D("","",nbins,&xbins[0]);
  TH1D* hMuE_Std = new TH1D("","",nbins,&xbins[0]);
  TH1D* hMuMu_Deco = new TH1D("","",nbins,&xbins[0]);
  TH1D* hMuE_Deco = new TH1D("","",nbins,&xbins[0]);
  TH1D* hMuMu_OQS = new TH1D("","",nbins,&xbins[0]);
  TH1D* hMuE_OQS = new TH1D("","",nbins,&xbins[0]);
  TH1D* hMuMu_OQS2 = new TH1D("","",nbins,&xbins[0]);
  

  /*
  double Emax = 100;
  
  TH1D* hMuMu_Std = new TH1D("","",nbins,1, Emax);
  TH1D* hMuE_Std = new TH1D("","",nbins,1, Emax);
  TH1D* hMuMu_Deco = new TH1D("","",nbins,1, Emax);
  TH1D* hMuE_Deco = new TH1D("","",nbins,1, Emax);
  TH1D* hMuMu_OQS = new TH1D("","",nbins,1, Emax);
  TH1D* hMuE_OQS = new TH1D("","",nbins,1, Emax);
  */


  int mh=1;
  
  double dm21 = 7.5e-5;
  double dm31 = mh>0 ? 2.457e-3 : -2.449e-3 + dm21;

  p.SetDm(2, dm21);
  p.SetDm(3, dm31);

  p2.SetDm(2, dm21);
  p2.SetDm(3, dm31);

  int flv_initial = 0;
  int flvf = 0;
  
  // Fill histograms
  for(int i=1; i<=nbins; i++){

    //double minE  = pow(10, hMuMu_Std->GetBinLowEdge(i));
    //double maxE  = pow(10, hMuMu_Std->GetBinLowEdge(i+1));

    double minE  = hMuMu_Std->GetBinLowEdge(i);
    double maxE  = hMuMu_Std->GetBinLowEdge(i+1);

    double E  = (maxE + minE) / 2;
    double dE = (maxE - minE);
        
    double energy = hMuMu_Std->GetBinCenter(i);

    // Fill std. osc.
    //    hMuMu_Std->SetBinContent(i, p.AvgProb(1,1, E, dE));
    //hMuE_Std->SetBinContent(i, p.AvgProb(1,0, E, dE));

    hMuMu_Std->SetBinContent(i, p.Prob(flv_initial, flvf, energy));
        
    //hMuMu_OQS->SetBinContent(i, p2.AvgProb(1,1, E, dE));
    //hMuE_OQS->SetBinContent(i, p2.AvgProb(1,0, E, dE));

    hMuMu_OQS->SetBinContent(i, p2.Prob(flv_initial, flvf, energy));
  }
  
  // Set some nice styles
  SetHist(hMuMu_Std, kBlack);
  SetHist(hMuMu_OQS, kRed);
  SetHist(hMuE_OQS, kGreen+2);
  SetHist(hMuMu_OQS2, kTeal+2);
  
  // Make lines thicker
  hMuMu_Std->SetLineWidth(3);
  hMuMu_OQS->SetLineWidth(3);
  hMuMu_OQS2->SetLineWidth(3);
  hMuE_Std->SetLineWidth(3);
  hMuE_OQS->SetLineWidth(3);

  // Make deco dashed
  //  hMuMu_Deco->SetLineStyle(9);
  //hMuE_Deco->SetLineStyle(9);

  // Set axis titles
  if(isNuBar) hMuMu_Std->SetTitle(";Neutrino Energy (GeV);P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{x})");
  else        hMuMu_Std->SetTitle(";Neutrino Energy (GeV);P(#nu_{e}#rightarrow#nu_{e})");
  //  else        hMuMu_Std->SetTitle(";Neutrino Energy (GeV);P(#nu_{#mu}#rightarrow#nu_{x})");

  // Set y range
  hMuMu_Std->GetYaxis()->SetRangeUser(0,1);

  // Make a long canvas
  MakeLongCanvas();

  // Draw everything
  hMuMu_Std->Draw("curv");
  //  hMuE_Std->Draw("curv same");
  hMuMu_OQS->Draw("curv same");
  //  hMuMu_OQS2->Draw("curv same");
  // hMuE_OQS->Draw("curv same");

  // Print L in canvas
  MiscText(0.70, 0.85, 0.05, TString::Format("L = %0.0f km", L) );

  TLegend* leg = new TLegend(0.7,0.6,0.80,0.8);
  //  TLegend* leg = new TLegend(0.7,0.55,0.91,0.82);
  
  if(isNuBar){

    leg->AddEntry(hMuMu_Std, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{#mu}) - Std.", "l");
    leg->AddEntry(hMuMu_OQS, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{#mu}) - OQS", "l");
    leg->AddEntry(hMuE_Std, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{e}) - Std.", "l");
    leg->AddEntry(hMuE_OQS, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{e}) - OQS", "l");
  }
  else{

    leg->AddEntry(hMuMu_Std, "Std.", "l");
    //    leg->AddEntry(hMuMu_Deco, "PQD", "l");
    leg->AddEntry(hMuMu_OQS, "OQS", "l");
  }

  //  SetLeg(leg);

  leg->Draw();

  gPad->SetLogx();

  hMuMu_Std->GetXaxis()->SetMoreLogLabels();

}
