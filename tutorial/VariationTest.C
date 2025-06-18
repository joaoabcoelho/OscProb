#include "PMNS_Fast.h"
#include "PMNS_TaylorExp.h"
#include "PremModel.h"

// Some functions to make nice plots
#include "SetNiceStyle.C"

// Check the accuracy of the
// oscillation averaging function
void VariationTest(){

    // Set nice overall style
    SetNiceStyle();

    // Get a PMNS object
    OscProb::PMNS_TaylorExp t;

    // PREM Model
    OscProb::PremModel prem;    

    // Fill path for cosT
    double cosT = -0.7;
    prem.FillPath(cosT);

    // Give path to calculator
    t.SetPath(prem.GetNuPath());

    int nbins = 1000;
    double E = 0.3;
    double xmax = 0.1;
    double xmin = -xmax;

    TH1D* h1 = new TH1D("","",nbins,xmin,xmax);
    TH1D* h2 = new TH1D("","",nbins,xmin,xmax);
    TH1D* h3 = new TH1D("","",nbins,xmin,xmax);

    for(int i = 1 ; i<=nbins ; i++){

        double varPercentage = h1->GetBinCenter(i);
        // EN E OU EN LOG?????

        h1->SetBinContent(i, t.Prob(1,1,E + varPercentage * E));
        h2->SetBinContent(i, t.interpolationEnergy(1,1,E,varPercentage * E));
        //h3->SetBinContent(i, t.Prob(1,1, E + varPercentage * E));
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
    h1->SetTitle(";#epsilon_{E } / E_{centedred};P(#nu_{#mu}#rightarrow#nu_{#mu})");

    // Draw different samplings
    h1->DrawCopy("curv");
    h2->DrawCopy("hist same ][");
    //h3->DrawCopy("hist same ][");

    MiscText(0.75, 0.965, 0.04, TString::Format("Centered Energy = %0.1f", E) );
    MiscText(0.63, 0.965, 0.04, TString::Format("cosT = %0.1f", cosT) );


    // Déterminer les limites de Y pour la ligne
    double y_min = gPad->GetUymin();
    double y_max = gPad->GetUymax();

    // Créer et dessiner la ligne verticale
    TLine *line = new TLine(0, y_min, 0, y_max);
    line->SetLineColor(kRed);
    line->SetLineStyle(1); // Ligne pointillée
    line->Draw("same");

}