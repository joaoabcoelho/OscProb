#include "PMNS_Fast.h"
#include "PMNS_TaylorExp.h"
#include "PremModel.h"

// Some functions to make nice plots
#include "SetNiceStyle.C"

// Check the accuracy of the
// oscillation averaging function
void VariationTestCosT(){

    // Set nice overall style
    SetNiceStyle();

    // Get a PMNS object
    OscProb::PMNS_Fast p;
    OscProb::PMNS_TaylorExp t;

    // PREM Model
    OscProb::PremModel prem;    
    OscProb::PremModel premTaylor; 

    int nbins = 1000;
    double xmax = 0.1;
    double xmin = -xmax;

    double cosT = -0.9;

    TH1D* h1 = new TH1D("","",nbins,xmin,xmax);
    TH1D* h2 = new TH1D("","",nbins,xmin,xmax);
    TH1D* h3 = new TH1D("","",nbins,xmin,xmax);

    int flavori = 0;
    double E = 0.3;
    t.SetEnergy(E);
    //p.SetEnergy(E);

    premTaylor.FillPath(cosT);
    t.SetPath(premTaylor.GetNuPath());
    t.GetPremLayers(premTaylor.GetPremLayers());

    //double L = 2*6368*abs(cosT);
    //t.SetLength(L);

    for(int i = 1 ; i<=nbins ; i++){

        double varPercentage = h1->GetBinCenter(i);

        prem.FillPath(cosT + varPercentage * cosT);
        p.SetPath(prem.GetNuPath());

        //double Lfast = 2*6368*abs(cosT+ varPercentage * cosT);
        //p.SetLength(Lfast);

        h1->SetBinContent(i, p.Prob(flavori,1,E));
        h2->SetBinContent(i, t.interpolationCosT(flavori,1,cosT,varPercentage * cosT));
    
    }

    // Make a long canvas
    MakeLongCanvas();

    // Set some nice histograms
    SetHist(h1, kBlack);
    SetHist(h2, kRed);
    SetHist(h3, kGreen);

    // Change line styles
    h1->SetLineStyle(1);
    h2->SetLineStyle(7);
    h3->SetLineStyle(7);
    

    // The axis titles
    if(flavori == 0) {
        h1->SetTitle(";#epsilon_{#theta } / #theta _{centedred};P(#nu_{e}#rightarrow#nu_{#mu})");
    }
    if(flavori == 1) {
        h1->SetTitle(";#epsilon_{#theta  } / #theta _{centedred};P(#nu_{#mu}#rightarrow#nu_{#mu})");
    }
    if(flavori == 2) {
        h1->SetTitle(";#epsilon_{#theta  } / #theta _{centedred};P(#nu_{#tau}#rightarrow#nu_{#mu})");
    }
    

    // Draw different samplings
    h1->DrawCopy("curv");
    h2->DrawCopy("hist same ][");
    //h3->DrawCopy("hist same ][");

    MiscText(0.85, 0.965, 0.04, TString::Format("Energy = %0.1f", E) );
    MiscText(0.63, 0.965, 0.04, TString::Format("cosT centered = %0.1f", cosT) );


    // Déterminer les limites de Y pour la ligne
    double y_min = gPad->GetUymin();
    double y_max = gPad->GetUymax();

    // Créer et dessiner la ligne verticale
    TLine *line = new TLine(0, y_min, 0, y_max);
    line->SetLineColor(kRed);
    line->SetLineStyle(1); // Ligne pointillée
    line->Draw("same");

}