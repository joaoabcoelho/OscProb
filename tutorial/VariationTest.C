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
    prem.FillPath(-0.6);

    // Give path to calculator
    t.SetPath(prem.GetNuPath());

    int nbins = 100;
    double xmax = 0.1;
    double xmin = -0.1;
    double xCentre = (xmax + xmin) / 2;

    double E = 0.5;

    TH1D* h1 = new TH1D("","",nbins,xmin,xmax);
    TH1D* h2 = new TH1D("","",nbins,xmin,xmax);

    for(int i = 1 ; i<=nbins ; i++){

        double energy = h1->GetBinCenter(i);
        // EN E OU EN LOG?????

        double ext = xCentre + energy ;

        cout<<ext<<endl;

        h1->SetBinContent(i, t.Prob(1,1,E + ext));
        h2->SetBinContent(i, t.interpolationEnergy(1,1,E,ext));
    }

    // Make a long canvas
    MakeLongCanvas();

    // Set some nice histograms
    SetHist(h1, kBlack);
    SetHist(h2, kBlue);

    // Change line styles
    h1->SetLineStyle(3);
    h2->SetLineStyle(7);

    // The axis titles
    //h1->SetTitle(";Log10[Neutrino Energy (GeV)];P(#nu_{#mu}#rightarrow#nu_{#mu})");

    // Draw different samplings
    h1->DrawCopy("curv");
    h2->DrawCopy("hist same ][");


    // Déterminer les limites de Y pour la ligne
    double y_min = gPad->GetUymin();
    double y_max = gPad->GetUymax();

    // Créer et dessiner la ligne verticale
    TLine *line = new TLine(xCentre, y_min, xCentre, y_max);
    line->SetLineColor(kRed);
    line->SetLineStyle(1); // Ligne pointillée
    line->Draw("same");

}