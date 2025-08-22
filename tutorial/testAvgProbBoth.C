
#include "PMNS_Fast.h"
#include "PMNS_TaylorExp.h"
#include "PremModel.h"

// Some functions to make nice plots
#include "SetNiceStyle.C"

// Check the accuracy of the
// oscillation averaging function
void testAvgProbBoth(){

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
  int navgE = 30;
  int navgCosT = 30;
  int nbinsE = navgE * 100;
  int nbinsCosT = navgCosT * 100;
  double xminE = 0;
  double xmaxE = 3;
  double yminCosT = -0.9;
  double ymaxCosT = -0.1;

  // Lots of histograms
  TH2D* h1 = new TH2D("","",nbinsE,xminE,xmaxE,nbinsCosT,yminCosT,ymaxCosT);
  TH2D* h2 = new TH2D("","",navgE,xminE,xmaxE,navgCosT,yminCosT,ymaxCosT);
  TH2D* h3 = new TH2D("","",navgE,xminE,xmaxE,navgCosT,yminCosT,ymaxCosT);
  TH2D* h4 = new TH2D("","",navgE,xminE,xmaxE,navgCosT,yminCosT,ymaxCosT);
  TH2D* h5 = new TH2D("","",navgE,xminE,xmaxE,navgCosT,yminCosT,ymaxCosT);

  //h2->SetBinContent(x,y,z);


  // Do some fine binning and uniform sampling
  for(int i=1; i<=nbinsCosT; i++){

    double minCosT = -pow(10, h1->GetYaxis()->GetBinLowEdge(i));
    double maxCosT = -pow(10, h1->GetYaxis()->GetBinLowEdge(i+1));

    double cosT = 0.5 * (minCosT + maxCosT);
    double dcosT = (minCosT - maxCosT);

    int abinCosT = h3->GetYaxis()->FindFixBin(log10(-cosT));

    prem.FillPath(cosT);
    p.SetPath(prem.GetNuPath());

    for(int j=1; j<=nbinsE; j++){

        double minE  = pow(10, h1->GetXaxis()->GetBinLowEdge(j));
        double maxE  = pow(10, h1->GetXaxis()->GetBinLowEdge(j+1));

        double E = 0.5 * (minE + maxE);
        double dE = (maxE - minE);

        int abinE = h3->GetXaxis()->FindFixBin(log10(E));

        double prob = p.Prob(1,1, E);

        h1->SetBinContent(j, i, prob);
        h3->AddBinContent(abinE, abinCosT, prob * dE * dcosT);

    }
  }

  

  // Do the AvgProb sampling and bin center sampling
  for(int i=1; i<=navgCosT; i++){

    double minCosT = -pow(10, h2->GetYaxis()->GetBinLowEdge(i));
    double maxCosT = -pow(10, h2->GetYaxis()->GetBinLowEdge(i+1));

    double cosT = 0.5 * (minCosT + maxCosT);
    double dcosT = (minCosT - maxCosT);

    prem.FillPath(cosT);
    p.SetPath(prem.GetNuPath());
    premTaylor.FillPath(cosT);
    taylor.SetPath(premTaylor.GetNuPath());

    for(int j=0; j<=navgE; j++){

        double minE  = pow(10, h2->GetXaxis()->GetBinLowEdge(j));
        double maxE  = pow(10, h2->GetXaxis()->GetBinLowEdge(j+1));

        double E = 0.5 * (minE + maxE);
        double dE = (maxE - minE);

        double a = taylor.AvgProb(1, 1, E, dE, cosT, dcosT);
        double b = h3->GetBinContent(j,i) / ( dcosT * dE);
        double c = p.Prob(1,1, E);
        double ab = abs(a-b);

        /*cout<<"cosT = "<<cosT<<endl;
        cout<<"dcosT = "<<dcosT<<"    diff = "<<ab;
    

        if (ab < 5 * 1E-4)
            cout<<"     TRUE"<<endl;
        else
            cout<<"     FALSE"<<endl;;

        cout<<"r = "<<dcosT/cosT<<endl<<endl;*/

        h2->SetBinContent(j, i, a);

        h3->SetBinContent(j, i, b);
        h4->SetBinContent(j,i, c);

        h5->SetBinContent(j, i, ab);

    }

  }

  // Set nice histogram
  SetHist(h5);

  // Set titles
  h5->SetTitle(";log(E) (GeV); log(cos#theta_{z});P_{#mu#mu}");

  h5->Draw("colz");

  gPad->SetRightMargin(0.18);


  // Make a long canvas
  /*MakeLongCanvas();

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

  leg->Draw();*/

  
}
