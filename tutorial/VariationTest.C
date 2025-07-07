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
    OscProb::PMNS_Fast p;
    OscProb::PMNS_TaylorExp t;

    // PREM Model
    OscProb::PremModel prem;    

    // Fill path for cosT
    double cosT = -0.9;
    prem.FillPath(cosT);
    t.SetPath(prem.GetNuPath());
    //double L= 2*6368*abs(cosT);
    //t.SetLength(L);

    //double L = 2*6368 + 18;
    //double L = 2*6368*abs(cosT);
    //p.SetLength(L);

    // Give path to calculator
    

    int nbins = 1000;
    double E = 0.3;
    double xmax = 0.1;
    double xmin = -xmax;

    cout<<t.GetDm(2)<<endl;
    cout<<t.GetDm(3)<<endl;
    //t.SetDm(2, -0.01);
    cout<<t.GetDm(2)<<endl;
    cout<<t.GetDm(3)<<endl;


    TH1D* h1 = new TH1D("","",nbins,xmin,xmax);
    TH1D* h2 = new TH1D("","",nbins,xmin,xmax);
    TH1D* h3 = new TH1D("","",nbins,xmin,xmax);

    TH1D* h4 = new TH1D("","",nbins,xmin,xmax);
    TH1D* h5 = new TH1D("","",nbins,xmin,xmax);

    TH1D* h6 = new TH1D("","",nbins,xmin,xmax);
    TH1D* h7 = new TH1D("","",nbins,xmin,xmax);
    

    int flavori = 1;

    for(int i = 1 ; i<=nbins ; i++){

        double varPercentage = h1->GetBinCenter(i);
        // EN E OU EN LOG?????

        h1->SetBinContent(i, t.Prob(flavori,1,E + varPercentage * E));
        h2->SetBinContent(i, t.interpolationEnergy(flavori,1,E,varPercentage * E));
        h4->SetBinContent(i, t.Prob(0,1,E + varPercentage * E));
        h5->SetBinContent(i, t.interpolationEnergy(0,1,E,varPercentage * E));
        h6->SetBinContent(i, t.Prob(2,1,E + varPercentage * E));
        h7->SetBinContent(i, t.interpolationEnergy(2,1,E,varPercentage * E));

        //h3->SetBinContent(i, p.Prob(flavori,1, E + varPercentage * E));
    }

    // Make a long canvas
    MakeLongCanvas();

    // Set some nice histograms
    SetHist(h1, kBlack);
    SetHist(h2, kRed);
    //SetHist(h3, kGreen);
    SetHist(h4, kBlack);
    SetHist(h5, kBlue);

    SetHist(h6, kYellow);
    SetHist(h7, kGreen);

    // Change line styles
    h1->SetLineStyle(1);
    h2->SetLineStyle(7);
    //h3->SetLineStyle(7);
    h4->SetLineStyle(1);
    h5->SetLineStyle(7);
    h6->SetLineStyle(1);
    h7->SetLineStyle(7);
    

    // The axis titles
    //if(flavori == 0) {
        //h1->SetTitle(";#epsilon_{E } / E_{centedred};P(#nu_{e}#rightarrow#nu_{#mu})");
    //}
    //if(flavori == 1) { h1->SetTitle(";#epsilon_{E } / E_{centedred};P(#nu_{#mu}#rightarrow#nu_{#mu})"); }
    //if(flavori == 2) {h1->SetTitle(";#epsilon_{E } / E_{centedred};P(#nu_{#tau}#rightarrow#nu_{#mu})");}
    h1->SetTitle(";#epsilon_{E } / E_{centedred};P(#nu_{#alpha}#rightarrow#nu_{#mu})");

    // Draw different samplings
    h1->DrawCopy("curv");
    h2->DrawCopy("hist same ][");
    //h3->DrawCopy("curv");
    h4->DrawCopy("hist same ][");
    h5->DrawCopy("hist same ][");
    //h6->DrawCopy("hist same ][");
    //h7->DrawCopy("hist same ][");

    MiscText(0.75, 0.965, 0.04, TString::Format("Centered Energy = %0.1f", E) );
    MiscText(0.63, 0.965, 0.04, TString::Format("cosT = %0.1f", cosT) );

    TLegend* leg = new TLegend(0.7,0.6,0.1,1);
    //leg->AddEntry(h1, " P(#nu_{#alpha}#rightarrow#nu_{#mu}) - exacte", "l");
    leg->AddEntry(h2, " P(#nu_{#mu}#rightarrow#nu_{#mu}) - algorithme", "l");
    //leg->AddEntry(h4, " P(#nu_{#e}#rightarrow#nu_{#mu}) - exacte", "l");
    leg->AddEntry(h5, " P(#nu_{e}#rightarrow#nu_{#mu}) - algorithme", "l");

    leg->SetLineColor(kBlack);  // Couleur du cadre
    leg->SetLineWidth(2);  
    SetLeg(leg);
    leg->Draw();


    // Déterminer les limites de Y pour la ligne
    double y_min = gPad->GetUymin();
    double y_max = gPad->GetUymax();

    // Créer et dessiner la ligne verticale
    TLine *line = new TLine(0, y_min, 0, y_max);
    line->SetLineColor(kRed);
    line->SetLineStyle(1); // Ligne pointillée
    line->Draw("same");

}