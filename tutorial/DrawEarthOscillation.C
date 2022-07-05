
#include "SetNiceStyle.C"
#include "LoadOscProb.C"

#ifndef __CINT__
#include "PremModel.h"
#include "PMNS_Fast.h"
bool isCINT = false;
#else
bool isCINT = true;
#endif
double tot_length(std::vector<OscProb::NuPath> in); 

void DrawEarthOscillation(double energy = 4,  bool isNuBar = false){

  // Load some nice styles
  SetNiceStyle();

  // Load OscProb library
  if(isCINT) LoadOscProb();

  // Probability Calculator
  OscProb::PMNS_Fast p;
  
  // Set neutrino or antineutrino
  p.SetIsNuBar(isNuBar);
  
  // PREM Model
  // OscProb::PremModel prem = new OscProb::PremModel("prem_425layers.txt");
  OscProb::PremModel prem;
  prem.LoadModel("/Users/lukas/Documents/PhD/OscProb/PremTables/prem_425layers.txt");
  
  // Fill straight-upward path
  prem.FillPath(-1);
  
  std::vector<OscProb::NuPath> fullpath = prem.GetNuPath();

  // Make some histograms
  int nbins = fullpath.size()+1;
  double xmin = 1;
  double xmax = 6368;
  
  TGraph* hMuMu_NH = new TGraph();
  TGraph* hEMu_NH  = new TGraph();
  TGraph* hMuMu_IH = new TGraph();
  TGraph* hEMu_IH  = new TGraph();

  double dm31 = 2.5e-3;


  // Fill histograms
  for(int i=1; i<=fullpath.size(); i++){

    std::vector<OscProb::NuPath> path;
    std::copy( fullpath.begin(), fullpath.begin() + i, std::back_inserter(path) ); 

    double l = tot_length(path);

    // Give path to calculator
    p.SetPath(path);
    
    // Set NH
    p.SetDm(3, dm31);
    
    // Fill NH 
    hMuMu_NH->AddPoint(l, p.Prob(1,1, energy));
    hEMu_NH->AddPoint(l, p.Prob(1,0, energy));

    // Set IH 
    p.SetDm(3, -dm31 + 7.52e-5);
    
    // Fill IH
    hMuMu_IH->AddPoint(l, p.Prob(1,1, energy));
    hEMu_IH->AddPoint(l, p.Prob(1,0, energy));
  
  }
  
  // Set some nice styles
  SetGraph(hMuMu_NH, kBlue);
  SetGraph(hEMu_NH, kRed);
  SetGraph(hMuMu_IH, kBlue);
  SetGraph(hEMu_IH, kRed);
  
  // Make IH dashed
  hMuMu_IH->SetLineStyle(7);
  hEMu_IH->SetLineStyle(7);
  
  // Set axis titles
  if(isNuBar) hMuMu_NH->SetTitle(";Neutrino Energy (GeV);P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{x})");
  else        hMuMu_NH->SetTitle(";Neutrino Energy (GeV);P(#nu_{#mu}#rightarrow#nu_{x})");
  
  // Set x range
  hMuMu_NH->GetXaxis()->SetRangeUser(0,6368*2);
  // Set y range
  hMuMu_NH->GetYaxis()->SetRangeUser(0,1);

  // Make a long canvas
  TCanvas* c1 = MakeLongCanvas();
  // Top margin for legend
  c1->SetTopMargin(0.2);
  
  // Draw everything
  hMuMu_NH->Draw("AL"); // to define the axis
  gPad->SetGrid();

  // Draw Earth layers
  double icb = 1221.5;
  double cmb = 3348;
  double rE = 6368;

  // inner core
  Double_t x1[5] = { rE-icb, rE+icb,  rE+icb,  rE-icb, rE-icb};
  Double_t y1[5] = {   0,    0,   1,  1,   0};
  auto excl1 = new TGraph(5,x1,y1);
  excl1->SetLineColor(kGray);
  excl1->SetFillColor(kGray);
  excl1->SetFillStyle(1001);
  // outer core
  Double_t x2[5] = {  rE-cmb, rE+cmb,  rE+cmb,   rE-cmb, rE-cmb};
  Double_t y2[5] = {  0,   0,  1, 1, 0};
  auto excl2 = new TGraph(5,x2,y2);
  excl2->SetLineColor(kOrange);
  excl2->SetFillColor(kOrange);
  excl2->SetFillStyle(1001);
  // mantle
  Double_t x3[5] = { 0, 2*rE,  2*rE,  0, 0};
  Double_t y3[5] = {   0,    0,   1,  1,   0};
  auto excl3 = new TGraph(5,x3,y3);
  excl3->SetLineColor(kYellow);
  excl3->SetFillColor(kYellow);
  excl3->SetFillStyle(1001);

  excl3->Draw("same F");
  excl2->Draw("same F");
  excl1->Draw("same F");

  double ytext = 0.4;
  TLatex* ic = MiscText(0.55, ytext, 0.04, "inner core" );
  ic->SetTextAngle(90);
  ic->Draw("same");
  TLatex* oc = MiscText(0.4, ytext, 0.04, "outer core" );
  oc->SetTextAngle(90);
  oc->Draw("same");
  TLatex* m = MiscText(0.28, ytext, 0.04, "mantle" );
  m->SetTextAngle(90);
  m->Draw("same");


  // Draw plots
  hMuMu_NH->Draw("same");
  hEMu_NH->Draw ("same");
  hMuMu_IH->Draw("same");
  hEMu_IH->Draw ("same");

  // // Print cosT in canvas
  // MiscText(0.8, 0.88, 0.04, TString::Format("cos#theta_{z} = %0.1f", cosT) );
  
  TLegend* leg = new TLegend(0.3,0.7,0.7,.97);
  leg->SetHeader("cos#theta_{z} = -1", "C");
  leg->SetNColumns(2);
  
  if(isNuBar){

    leg->AddEntry(hMuMu_NH, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{#mu}) - NH", "l");
    leg->AddEntry(hMuMu_IH, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{#mu}) - IH", "l");
    leg->AddEntry(hEMu_NH, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{e}) - NH", "l");
    leg->AddEntry(hEMu_IH, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{e}) - IH", "l");

  }
  else{
  
    leg->AddEntry(hMuMu_NH, "P(#nu_{#mu}#rightarrow#nu_{#mu}) - NH", "l");
    leg->AddEntry(hMuMu_IH, "P(#nu_{#mu}#rightarrow#nu_{#mu}) - IH", "l");
    leg->AddEntry(hEMu_NH, "P(#nu_{#mu}#rightarrow#nu_{e}) - NH", "l");
    leg->AddEntry(hEMu_IH, "P(#nu_{#mu}#rightarrow#nu_{e}) - IH", "l");
  
  }
  
  SetLeg(leg);
  
  leg->Draw();
  gPad->RedrawAxis();
  gPad->RedrawAxis("G");


  // c1->Draw();
  c1->Print(TString::Format("EarthOscillation_E%0.0f.png", energy));
  
}

double tot_length(std::vector<OscProb::NuPath> in){
  double l = 0;
  for( int i = 0; i<in.size(); i++){
    l += in[i].length;
  }
  return l;

}
