
#include "PMNS_Fast.h"
#include "PMNS_TaylorExp.h"
#include "PremModel.h"

// Some functions to make nice plots
#include "SetNiceStyle.C"

// Check the accuracy of the
// oscillation averaging function
void testAvgProbCosT(){

  // Set nice overall style
  SetNiceStyle();

  // Get a PMNS object
  OscProb::PMNS_Fast p;
  OscProb::PMNS_TaylorExp taylor;

  // PREM Model
  OscProb::PremModel prem;
  OscProb::PremModel premTaylor;    

  taylor.GetPremLayers(premTaylor.GetPremLayers());

  // Define some fine and coarse binnings
  int navg = 40;
  int nbins = navg * 100;
  double xmin = -0.9;
  double xmax = -0.1;

  double E = 0.3;

  // Lots of histograms
  TH1D* h1 = new TH1D("","",nbins,xmin,xmax);
  TH1D* h2 = new TH1D("","",navg,xmin,xmax);
  TH1D* h3 = new TH1D("","",navg,xmin,xmax);
  TH1D* h4 = new TH1D("","",navg,xmin,xmax);

  


  // Do some fine binning and uniform sampling
  for(int i=1; i<=nbins; i++){

    double minCosT = -pow(10, h1->GetBinLowEdge(i));
    double maxCosT = -pow(10, h1->GetBinLowEdge(i+1));

    double cosT = 0.5 * (minCosT + maxCosT);
    double dcosT = (minCosT - maxCosT);

    int abin = h3->FindFixBin(log10(-cosT));

    prem.FillPath(cosT);
    p.SetPath(prem.GetNuPath());

    double prob = p.Prob(1,1, E);

    h1->SetBinContent(i, prob);
    h3->AddBinContent(abin, prob * dcosT);

  }

  // Do the AvgProb sampling and bin center sampling
  for(int i=1; i<=navg; i++){

    double minCosT  = -pow(10, h2->GetBinLowEdge(i));
    double maxCosT  = -pow(10, h2->GetBinLowEdge(i+1));

    double cosT = 0.5 * (minCosT + maxCosT);
    double dcosT = (minCosT - maxCosT);

    prem.FillPath(cosT);
    p.SetPath(prem.GetNuPath());
    premTaylor.FillPath(cosT);
    taylor.SetPath(premTaylor.GetNuPath());
  
    double a = taylor.AvgProb(1,1, E, cosT,dcosT);
    double b = h3->GetBinContent(i) / dcosT;
    double ab = abs(a-b);

    cout<<"cosT = "<<cosT<<endl;
    cout<<"dcosT = "<<dcosT<<"    diff = "<<ab;
    

    if (ab < 5 * 1E-4)
      cout<<"     TRUE"<<endl;
    else
      cout<<"     FALSE"<<endl;;

    cout<<"r = "<<dcosT/cosT<<endl<<endl;



    h2->SetBinContent(i, a);

    h3->SetBinContent(i, b);
    h4->SetBinContent(i, p.Prob(1,1, E));



  }





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
  h4->SetLineWidth(3);

  // The axis titles
  h1->SetTitle(";Log10[Neutrino Energy (GeV)];P(#nu_{#mu}#rightarrow#nu_{#mu})");

  // Draw different samplings
  h1->DrawCopy("curv");
  h2->DrawCopy("hist same ][");
  h3->DrawCopy("hist same ][");
  //h4->DrawCopy("hist same ][");

  MiscText(0.75, 0.965, 0.04, TString::Format("nbrBin = %0.1d", navg) );
  MiscText(0.63, 0.965, 0.04, TString::Format("E = %0.1f", E) );

  TLegend* leg = new TLegend(0.2,0.9,0.9,0.6);

  //leg->AddEntry(h4," P cst in every bin ", "l");
  leg->AddEntry(h3, " avg P", "l");
  leg->AddEntry(h2, "avg P with Taylor", "l");

  leg->SetLineColor(kBlack);  // Couleur du cadre
  leg->SetLineWidth(2);       // Épaisseur du cadre
  //leg->AddEntry(h2, "P(#nu_{#mu}#rightarrow#nu_{#mu}) - avg with Taylor", "l");
  
  //leg->AddEntry(h4, "P(#nu_{#mu}#rightarrow#nu_{#mu}) - centre bin", "l");

  SetLeg(leg);

  leg->Draw();

  
}
