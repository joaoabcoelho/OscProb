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

    t.SetAngle(1, 2, asin(sqrt(0.303)));
    t.SetAngle(1, 3, asin(sqrt(0.022)));
    t.SetAngle(2, 3, asin(sqrt(0.572)));
    t.SetDm(2, 7.41e-5);
    t.SetDm(3, 2.51e-3);
    t.SetDelta(1,3,197);

    p.SetAngle(1, 2, asin(sqrt(0.303)));
    p.SetAngle(1, 3, asin(sqrt(0.022)));
    p.SetAngle(2, 3, asin(sqrt(0.572)));
    p.SetDm(2, 7.41e-5);
    p.SetDm(3, 2.51e-3);
    p.SetDelta(1,3,197);

    // PREM Model
    OscProb::PremModel prem;    
    OscProb::PremModel premTaylor; 

    int nbins = 1000;
    double xmax = 0.11;
    double xmin = -xmax;

    double cosT = -0.9;
    double Theta = acos(cosT);
    double E = 0.3;
    t.SetEnergy(E);

    TH1D* NumuNumu_Exact = new TH1D("","",nbins,xmin,xmax);
    TH1D* NumuNumu_Order = new TH1D("","",nbins,xmin,xmax);

    TH1D* NueNumu_Exact = new TH1D("","",nbins,xmin,xmax);
    TH1D* NueNumu_Order = new TH1D("","",nbins,xmin,xmax);

    TH1D* NumuNumuBar_Exact = new TH1D("","",nbins,xmin,xmax);
    TH1D* NumuNumuBar_Order = new TH1D("","",nbins,xmin,xmax);

    TH1D* NueNumuBar_Exact = new TH1D("","",nbins,xmin,xmax);
    TH1D* NueNumuBar_Order = new TH1D("","",nbins,xmin,xmax);

    int flavorf = 1;
    bool barOrNot;

    premTaylor.FillPath(cosT);
    t.SetPath(premTaylor.GetNuPath());
    t.GetPremLayers(premTaylor.GetPremLayers());

    for(int i = 1 ; i<=nbins ; i++){

        //double varPercentage = -NumuNumu_Exact->GetBinCenter(i);

        double ThetaVaria = -NumuNumu_Exact->GetBinCenter(i);

        //prem.FillPath(cosT + varPercentage * cosT);
        prem.FillPath(cos(Theta + ThetaVaria) );
        p.SetPath(prem.GetNuPath());

        t.SetIsNuBar(false);
        p.SetIsNuBar(false);

        NumuNumu_Exact->SetBinContent(i, p.Prob(1,flavorf,E));
        NumuNumu_Order->SetBinContent(i, t.interpolationCosT(1,flavorf,cosT, cos(Theta + ThetaVaria) - cosT ));
        NueNumu_Exact->SetBinContent(i, p.Prob(0,flavorf,E));
        NueNumu_Order->SetBinContent(i, t.interpolationCosT(0,flavorf,cosT, cos(Theta + ThetaVaria) - cosT  ));
        
        t.SetIsNuBar(true);
        p.SetIsNuBar(true);

        NumuNumuBar_Exact->SetBinContent(i, p.Prob(1,flavorf,E));
        NumuNumuBar_Order->SetBinContent(i, t.interpolationCosT(1,flavorf,cosT,cos(Theta + ThetaVaria) - cosT ));
        NueNumuBar_Exact->SetBinContent(i, p.Prob(0,flavorf,E));
        NueNumuBar_Order->SetBinContent(i, t.interpolationCosT(0,flavorf,cosT,cos(Theta + ThetaVaria) - cosT ));
        //cos(Theta + ThetaVaria) - cosT
    }

    // Make a long canvas
    MakeLongCanvas();

    // Set some nice histograms
    if (flavorf == 0){
        SetHist(NumuNumu_Exact, kBlue);
        SetHist(NumuNumu_Order, kBlack);
        SetHist(NueNumu_Exact, kRed);
        SetHist(NueNumu_Order, kBlack);

        SetHist(NumuNumuBar_Exact, kBlue);
        SetHist(NumuNumuBar_Order, kBlack);
        SetHist(NueNumuBar_Exact, kRed);
        SetHist(NueNumuBar_Order, kBlack);
    }
    if (flavorf == 1){
        SetHist(NumuNumu_Exact, kGreen);
        SetHist(NumuNumu_Order, kBlack);
        SetHist(NueNumu_Exact, kOrange);
        SetHist(NueNumu_Order, kBlack);

        SetHist(NumuNumuBar_Exact, kGreen);
        SetHist(NumuNumuBar_Order, kBlack);
        SetHist(NueNumuBar_Exact, kOrange);
        SetHist(NueNumuBar_Order, kBlack);
    }

    // Change line styles
    NumuNumu_Exact->SetLineStyle(1);
    NumuNumu_Exact->SetLineWidth(2);
    NumuNumu_Order->SetLineStyle(1);
    NueNumu_Exact->SetLineStyle(1);
    NueNumu_Exact->SetLineWidth(2);
    NueNumu_Order->SetLineStyle(1);

    NumuNumuBar_Exact->SetLineStyle(7);
    NumuNumuBar_Exact->SetLineWidth(2);
    NumuNumuBar_Order->SetLineStyle(7);
    NueNumuBar_Exact->SetLineStyle(7);
    NueNumuBar_Exact->SetLineWidth(2);
    NueNumuBar_Order->SetLineStyle(7);
    
    
    // The axis titles
    if (flavorf == 0)
        NumuNumuBar_Order->SetTitle(";-#epsilon_{cosT } / cosT_{centedred};P(#nu_{#alpha}#rightarrow#nu_{e})");

    if (flavorf == 1)
        NumuNumuBar_Order->SetTitle(";-#epsilon_{cosT } / cosT_{centedred};P(#nu_{#alpha}#rightarrow#nu_{#mu})");
    
    NumuNumuBar_Order->GetYaxis()->SetRangeUser(0,1.1);

    // Draw different samplings
    NumuNumuBar_Order->DrawCopy("curv");
    NumuNumuBar_Exact->DrawCopy("hist same ][");
    NueNumuBar_Order->DrawCopy("hist same ][");
    NueNumuBar_Exact->DrawCopy("hist same ][");

    NumuNumu_Order->DrawCopy("hist same ][");
    NumuNumu_Exact->DrawCopy("hist same ][");
    NueNumu_Order->DrawCopy("hist same ][");
    NueNumu_Exact->DrawCopy("hist same ][");

    MiscText(0.81, 0.965, 0.04, TString::Format("Energy = %0.1f", E) );
    MiscText(0.59, 0.965, 0.04, TString::Format(" Centered cosT = %0.1f", cosT) );

    TLegend* leg = new TLegend(0.7,0.6,0.1,1.005);
    if(flavorf == 0){
        leg->AddEntry(NumuNumu_Exact, " P(#nu_{#mu}#rightarrow#nu_{e})", "l");
        leg->AddEntry(NueNumu_Exact, " P(#nu_{e}#rightarrow#nu_{e})", "l");
    }
    if(flavorf == 1){
        leg->AddEntry(NumuNumu_Exact, " P(#nu_{#mu}#rightarrow#nu_{#mu})", "l");
        leg->AddEntry(NueNumu_Exact, " P(#nu_{e}#rightarrow#nu_{#mu})", "l");
    }

    TLegend* legBar = new TLegend(0.7,0.6,0.37,1.005);
    if(flavorf == 0){
        legBar->AddEntry(NumuNumuBar_Exact, " P(#bar#nu_{#mu}#rightarrow#bar#nu_{e})", "l");
        legBar->AddEntry(NueNumuBar_Exact, " P(#bar#nu_{e}#rightarrow#bar#nu_{e})", "l");
    }
    if(flavorf == 1){
        legBar->AddEntry(NumuNumuBar_Exact, " P(#bar#nu_{#mu}#rightarrow#bar#nu_{#mu})", "l");
        legBar->AddEntry(NueNumuBar_Exact, " P(#bar#nu_{e}#rightarrow#bar#nu_{#mu})", "l");
    }

    leg->SetLineColor(kBlack);  // Couleur du cadre
    leg->SetLineWidth(2);  
    SetLeg(leg);
    leg->Draw();

    legBar->SetLineColor(kBlack);  // Couleur du cadre
    legBar->SetLineWidth(2);  
    SetLeg(legBar);
    legBar->Draw();

    // Déterminer les limites de Y pour la ligne
    double y_min = gPad->GetUymin();
    double y_max = gPad->GetUymax();

    // Créer et dessiner la ligne verticale
    TLine *line = new TLine(0, y_min, 0, y_max);
    line->SetLineColor(kRed);
    line->SetLineStyle(1); // Ligne pointillée
    line->Draw("same");

}