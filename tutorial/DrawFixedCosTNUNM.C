
#include "PremModel.h"
#include "PMNS_Fast.h"

#include "SetNiceStyle.C"

void Plot(double cosT = -1.0, bool isNuBar = false, int k=0, int l=0){

 for(int n=0; n<2; n++){ // proba
  for(int m=1; m<3; m++){

  // --- which proba ---
  int flv_i = n;
  int flv_f = m;

  // --- std pars ---
  double dm31 = 2.5e-3;
  double th23 = 45.0/180.*3.14159;
  double th14, th24, th34, alpha;
  // --- sterile pars --- 
  if (k == l){
    alpha = -0.051;
    if (k == 0){
      th14 = 18.4349/180.*3.14159;
      th24 = 0;
      th34 = 0;
    }
    if (k == 1){
      th14 = 0;
      th24 = 18.4349/180.*3.14159;
      th34 = 0; 
    }
    if (k == 2){
      th14 = 0;
      th24 = 0;
      th34 = 18.4349/180.*3.14159;
    }
  } 
  else {
    alpha = -0.25;
    if (k == 2 && l == 1){
      th14 = 0; 
      th24 = 30 / 180.*3.14159;
      th34 = 30 / 180.*3.14159;
      //th24 = 18.4349/180.*3.14159;
      //th34 = 18.4349/180.*3.14159;
    }
    if ( k == 1 && l == 0){
      th14 = 30/180.*3.14159;
      th24 = 30/180.*3.14159;
      th34 = 0;
    }
    if ( k == 2 && l == 0){
      alpha = -0.0948;
      th14 = 30/180.*3.14159;
      th24 = 0;
      th34 = 30/180.*3.14159;
    }  
  }

  //alpha = +0.1;
  
  double dm41 = 0.5;

  // alpha matrix
  int ind_i = k;
  int ind_j = l;
  //double alpha = -0.1;//-0.015;//-0.2;//-0.051;

  // path to calculator
  double L = 1300; // 3000 // dune baseline
  double rho = 0;//2.85; //3.35 // aproximate mantle density  
  double rho2 = 50;
  
  // make some histograms
  int nbins = 500; // 1000
  double xmin = 5; // 10
  double xmax = 250; // 250

  // Load some nice styles
  SetNiceStyle();

  // Probability Calculator
  OscProb::PMNS_Fast p;
  OscProb::PMNS_Sterile ps(4);
  OscProb::PMNS_Sterile psn0(4);
  OscProb::PMNS_NUNM pa(0);
  OscProb::PMNS_NUNM pan0(0);
  OscProb::PMNS_NUNM pah(1);
  OscProb::PMNS_NUNM pahn0(1);

  // Set neutrino or antineutrino

  // PREM Model
  OscProb::PremModel prem;

  // Fill path for cosT
  prem.FillPath(cosT);
  
  /*
  p.SetPath(L, rho); 
  ps.SetPath(L, rho);
  psn0.SetPath(L, rho);
  pa.SetPath(L, rho);
  pan0.SetPath(L, rho);
   
  p.AddPath(L, rho);
  ps.AddPath(L, rho);
  psn0.AddPath(L, rho);
  pa.AddPath(L, rho);
  pan0.AddPath(L, rho);
  */

  
  p.SetPath(prem.GetNuPath());
  ps.SetPath(prem.GetNuPath());
  psn0.SetPath(prem.GetNuPath());
  pa.SetPath(prem.GetNuPath());
  pan0.SetPath(prem.GetNuPath());
  pah.SetPath(prem.GetNuPath());
  pahn0.SetPath(prem.GetNuPath()); 
  

  TH1D* hMuMu_NH = new TH1D("","",nbins, xmin, xmax);
  TH1D* hMuMu_NH_s = new TH1D("","",nbins, xmin, xmax);
  TH1D* hMuMu_NH_sn0 = new TH1D("","",nbins, xmin, xmax);
  TH1D* hMuMu_NH_a = new TH1D("","",nbins, xmin, xmax);
  TH1D* hMuMu_NH_an0 = new TH1D("","",nbins, xmin, xmax);
  TH1D* hMuMu_NH_ah = new TH1D("","",nbins, xmin, xmax);
  TH1D* hMuMu_NH_ahn0 = new TH1D("","",nbins, xmin, xmax);
  //TH1D* hEMu_NH = new TH1D("","",nbins,xmin, xmax);
  //TH1D* hMuMu_IH = new TH1D("","",nbins, xmin, xmax);
  //TH1D* hEMu_IH = new TH1D("","",nbins, xmin, xmax);

  
    // Std
    p.SetDm(3, dm31);
    p.SetAngle(2,3, th23);
    
    // Ster. 4 sin(th34)^2 = 0.1
    ps.SetDm(3, dm31);
    ps.SetAngle(2,3, th23);
    ps.SetDm(4, dm41);
    ps.SetAngle(3, 4, th34);
    ps.SetAngle(2, 4, th24);
    ps.SetAngle(1, 4, th14);

    // Ster. 4 sin(th34)^2 = 0.1
    psn0.SetDm(3, dm31);
    psn0.SetAngle(2,3, th23);
    psn0.SetDm(4, dm41);
    psn0.SetAngle(3, 4, th34);
    psn0.SetAngle(2, 4, th24);
    psn0.SetAngle(1, 4, th14);
    psn0.SetFracVnc(0.);
    
    // NUNM alpha33 = -0.051
    pa.SetDm(3, dm31);
    pa.SetAngle(2,3, th23);
    pa.SetAlpha(ind_i, ind_j, alpha, 0 );
    
    // NUNM alpha33 neutron density 0
    pan0.SetDm(3, dm31);
    pan0.SetAngle(2,3, th23);
    pan0.SetAlpha(ind_i, ind_j, alpha, 0 );
    pan0.SetFracVnc(0.);
    
    // NUNM alpha33 = -0.051
    pah.SetDm(3, dm31);
    pah.SetAngle(2,3, th23);
    pah.SetAlpha(ind_i, ind_j, alpha, 0 );

    // NUNM alpha33 neutron density 0
    pahn0.SetDm(3, dm31);
    pahn0.SetAngle(2,3, th23);
    pahn0.SetAlpha(ind_i, ind_j, alpha, 0 );
    pahn0.SetFracVnc(0.);
  
  for(int i=1; i<=nbins; i++){

    double energy = hMuMu_NH->GetBinCenter(i);
    //double energy = 25;
    
    // Fill NH
    hMuMu_NH->SetBinContent(i, p.Prob(flv_i, flv_f, energy));
    hMuMu_NH_s->SetBinContent(i, ps.Prob(flv_i, flv_f, energy));
    hMuMu_NH_sn0->SetBinContent(i, psn0.Prob(flv_i, flv_f, energy));
    hMuMu_NH_a->SetBinContent(i, pa.Prob(flv_i, flv_f, energy));
    hMuMu_NH_an0->SetBinContent(i, pan0.Prob(flv_i, flv_f, energy));
    hMuMu_NH_ah->SetBinContent(i, pah.Prob(flv_i, flv_f, energy));
    hMuMu_NH_ahn0->SetBinContent(i, pahn0.Prob(flv_i, flv_f, energy));

    // Set IH
    //p.SetDm(3, -dm31 + 7.52e-5);

    // Fill IH
    //hMuMu_IH->SetBinContent(i, p.Prob(1,1, energy));
    //hEMu_IH->SetBinContent(i, p.Prob(0,1, energy));

  }

  // Set some nice styles
  SetHist(hMuMu_NH, kBlack);
  //SetHist(hEMu_NH, kRed);
  //SetHist(hMuMu_IH, kBlue);
  //SetHist(hEMu_IH, kRed);
  SetHist(hMuMu_NH_s, kBlue);
  SetHist(hMuMu_NH_sn0, kCyan);
  SetHist(hMuMu_NH_a, kRed);
  SetHist(hMuMu_NH_an0, kOrange);
  SetHist(hMuMu_NH_ah, kGreen);
  SetHist(hMuMu_NH_ahn0, kViolet);


  // Make IH dashed
  hMuMu_NH_sn0->SetLineStyle(7);
  hMuMu_NH_an0->SetLineStyle(3);
  hMuMu_NH_ahn0->SetLineStyle(7);
  hMuMu_NH_s->SetLineStyle(2);
  hMuMu_NH_a->SetLineStyle(2);  
  hMuMu_NH_ah->SetLineStyle(2);
  hMuMu_NH_ahn0->SetLineStyle(7);
  //hMuMu_NH_a->SetLineWidth(2);
  //hMuMu_NH_an0->SetLineWidth(2);
  //hMuMu_NH_a->SetLineStyle(7);

  // Set axis titles
  TString myStrings[3] = {"e", "#mu", "#tau"};
  hMuMu_NH->SetTitle(TString::Format(";Neutrino Energy (GeV);P(#nu_{%s}#rightarrow#nu_{%s})", myStrings[flv_i].Data(), myStrings[flv_f].Data() ));

  // Set y range
  hMuMu_NH->GetYaxis()->SetRangeUser(0,1);

  // Make a long canvas
  TCanvas* c1 = MakeLongCanvas();

  // Draw everything
  hMuMu_NH->Draw("curv");
  hMuMu_NH_s->Draw("curv same");
  hMuMu_NH_sn0->Draw("curv same");
  hMuMu_NH_a->Draw("curv same");
  hMuMu_NH_an0->Draw("curv same");
  hMuMu_NH_ah->Draw("curv same");
  hMuMu_NH_ahn0->Draw("curv same");
  //hEMu_NH->Draw("curv same");
  //hMuMu_IH->Draw("curv same");
  //hEMu_IH->Draw("curv same");

  // Print cosT in canvas
  MiscText(0.8, 0.85, 0.04, TString::Format("cos#theta_{z} = %0.1f", cosT) );
  //MiscText(0.65, 0.85, 0.04, TString::Format("baseline = %0.1f km", L) );

  TLegend* leg = new TLegend(0.53,0.3,0.9,0.6);
  
  //if(isNuBar){
  
    leg->AddEntry(hMuMu_NH, "UNM", "l");
    //leg->AddEntry(hMuMu_NH_s, TString::Format("sin(#Theta_{14})=%0.3f; sin(#Theta_{24})=%0.3f; sin(#Theta_{34})=%0.3f",sin(th14), sin(th24), sin(th34)), "l");
    leg->AddEntry(hMuMu_NH_s, TString::Format("sin(#Theta_{34})^{2}=%0.2f", sin(th34)*sin(th34)), "l");
    //leg->AddEntry(hMuMu_NH_sn0, TString::Format("sin(#Theta_{14})=%0.3f; sin(#Theta_{24})=%0.3f; sin(#Theta_{34})=%0.3f; Nn=0",sin(th14), sin(th24), sin(th34)), "l");
    leg->AddEntry(hMuMu_NH_sn0, TString::Format("sin(#Theta_{34})^{2}=%0.2f; #Delta_{n}=0", sin(th34)*sin(th34)) , "l");
    //leg->AddEntry(hMuMu_NH_a, TString::Format("NUNM low s. #alpha_{%i%i}=%0.2f", ind_i+1, ind_j+1 , alpha), "l"); // -0.051", "l");    
    //leg->AddEntry(hMuMu_NH_an0, TString::Format("NUNM low s. #alpha_{%i%i}=%0.2f; #Delta_{n}=0", ind_i+1, ind_j+1, alpha), "l");// -0.051; Nn=0", "l");
    leg->AddEntry(hMuMu_NH_a, "NUNM #alpha_{22}=-0.001", "l"); // -0.051", "l");    
    leg->AddEntry(hMuMu_NH_an0, "NUNM |#alpha_{22}|=0.001; Nn=0", "l");// -0.051; Nn=0", "l");
    leg->AddEntry(hMuMu_NH_ah, TString::Format("NUNM high s. #alpha_{%i%i}=%0.2f", ind_i+1, ind_j+1 , alpha), "l"); // -0.051", "l"); 
    leg->AddEntry(hMuMu_NH_ahn0, TString::Format("NUNM high s. #alpha_{%i%i}=%0.2f; #Delta_{n}=0", ind_i+1, ind_j+1 , alpha), "l"); // -0.051", "l"); 
    //leg->AddEntry(hMuMu_NH_ah, "P(#nu_{#mu}#rightarrow#nu_{#mu}); NUNM high sc. 2. #alpha_{33}=-0.051 ", "l");// -0.051", "l");
    //leg->AddEntry(hMuMu_NH_ahn0, "P(#nu_{#mu}#rightarrow#nu_{#mu}); NUNM high sc. 2. #alpha_{33}=-0.051; Nn=0 ", "l");// -0.051; Nn=0", "l");
    //leg->AddEntry(hMuMu_IH, "P(#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{#mu}) - IH", "l");
    //leg->AddEntry(hEMu_NH, "P(#bar{#nu}_{e}#rightarrow#bar{#nu}_{#mu}) - NH", "l");
    //leg->AddEntry(hEMu_IH, "P(#bar{#nu}_{e}#rightarrow#bar{#nu}_{#mu}) - IH", "l");

    //leg->AddEntry(hMuMu_NH, "P(#nu_{#mu}#rightarrow#nu_{#mu}) - NH", "l");
    //leg->AddEntry(hMuMu_IH, "P(#nu_{#mu}#rightarrow#nu_{#mu}) - IH", "l");
    //leg->AddEntry(hEMu_NH, "P(#nu_{e}#rightarrow#nu_{#mu}) - NH", "l");
    //leg->AddEntry(hEMu_IH, "P(#nu_{e}#rightarrow#nu_{#mu}) - IH", "l");
  
  c1->SetLogx();
  SetLeg(leg);

  leg->Draw();
  c1->SaveAs(TString::Format("plot/NUNM_v0/earth_v1_ssn0ahahn0_a_%i%i_alpha_%0.3f_flv_%i_flv_%i.pdf", ind_i+1, ind_j+1 , alpha, flv_i, flv_f));

  }
 } 
}

void DrawFixedCosTNUNM(){
    //for(int i=0; i<3; i++){
      //for(int j=0; j<i+1; j++){
        Plot( -1.0,  false, 2, 1);
      //}
    //}
}
