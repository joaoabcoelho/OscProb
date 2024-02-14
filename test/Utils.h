
#include "colormod.h"

#include "../tutorial/SetNiceStyle.C"

//.............................................................................
void SetNominalPars(OscProb::PMNS_Base* p){

  // Set to NuFIT 5.2 values (NO w/ SK)
  p->SetDm(2, 7.41e-5);
  p->SetDm(3, 2.507e-3);
  p->SetAngle(1,2, asin(sqrt(0.303)));
  p->SetAngle(1,3, asin(sqrt(0.02225)));
  p->SetAngle(2,3, asin(sqrt(0.451)));
  p->SetDelta(1,3, 232 * TMath::DegToRad());

}

//.............................................................................
OscProb::PMNS_Fast* GetFast(bool is_nominal){

  OscProb::PMNS_Fast* p = new OscProb::PMNS_Fast();
  SetNominalPars(p);
  return p;

}

//.............................................................................
OscProb::PMNS_Iter* GetIter(bool is_nominal){

  OscProb::PMNS_Iter* p = new OscProb::PMNS_Iter();
  SetNominalPars(p);
  return p;

}

//.............................................................................
OscProb::PMNS_Deco* GetDeco(bool is_nominal){

  OscProb::PMNS_Deco* p = new OscProb::PMNS_Deco();
  SetNominalPars(p);
  if(!is_nominal){
    p->SetGamma(2, 1e-23);
    p->SetGamma(3, 1e-22);
  }

  return p;

}

//.............................................................................
OscProb::PMNS_Sterile* GetSterile(bool is_nominal){

  OscProb::PMNS_Sterile* p = new OscProb::PMNS_Sterile(4);
  SetNominalPars(p);
  if(!is_nominal){
    p->SetDm(4, 0.1);
    p->SetAngle(1,4, 0.1);
    p->SetAngle(2,4, 0.1);
    p->SetAngle(3,4, 0.1);
  }

  return p;

}

//.............................................................................
OscProb::PMNS_Decay* GetDecay(bool is_nominal){

  OscProb::PMNS_Decay* p = new OscProb::PMNS_Decay();
  SetNominalPars(p);
  if(!is_nominal){
    p->SetAlpha3(1e-4);
  }

  return p;

}

//.............................................................................
OscProb::PMNS_NSI* GetNSI(bool is_nominal){

  OscProb::PMNS_NSI* p = new OscProb::PMNS_NSI();
  SetNominalPars(p);
  if(!is_nominal){
    p->SetEps(0,0, 0.1, 0);
    p->SetEps(0,1, 0.2, 0);
    p->SetEps(0,2, 0.3, 0);
    p->SetEps(1,1, 0.4, 0);
    p->SetEps(1,2, 0.5, 0);
    p->SetEps(2,2, 0.6, 0);
  }

  return p;

}

//.............................................................................
OscProb::PMNS_SNSI* GetSNSI(bool is_nominal){

  OscProb::PMNS_SNSI* p = new OscProb::PMNS_SNSI();
  SetNominalPars(p);
  if(!is_nominal){
    p->SetEps(0,0, 0.1, 0);
    p->SetEps(0,1, 0.2, 0);
    p->SetEps(0,2, 0.3, 0);
    p->SetEps(1,1, 0.4, 0);
    p->SetEps(1,2, 0.5, 0);
    p->SetEps(2,2, 0.6, 0);
  }

  return p;

}

//.............................................................................
OscProb::PMNS_LIV* GetLIV(bool is_nominal){

  OscProb::PMNS_LIV* p = new OscProb::PMNS_LIV();
  SetNominalPars(p);
  if(!is_nominal){
    p->SetaT(0,0, 0.1e-22, 0);
    p->SetaT(0,1, 0.2e-22, 0);
    p->SetaT(0,2, 0.3e-22, 0);
    p->SetaT(1,1, 0.4e-22, 0);
    p->SetaT(1,2, 0.5e-22, 0);
    p->SetaT(2,2, 0.6e-22, 0);
    p->SetcT(0,0, 0.1e-22, 0);
    p->SetcT(0,1, 0.2e-22, 0);
    p->SetcT(0,2, 0.3e-22, 0);
    p->SetcT(1,1, 0.4e-22, 0);
    p->SetcT(1,2, 0.5e-22, 0);
    p->SetcT(2,2, 0.6e-22, 0);
  }

  return p;

}

//.............................................................................
OscProb::PMNS_NUNM* GetNUNM(bool is_nominal){

  OscProb::PMNS_NUNM* p = new OscProb::PMNS_NUNM(0);
  SetNominalPars(p);
  if(!is_nominal){
    p->SetAlpha(0,0, 0.1, 0);
    p->SetAlpha(1,0, 0.2, 0);
    p->SetAlpha(2,0, 0.3, 0);
    p->SetAlpha(1,1, 0.4, 0);
    p->SetAlpha(2,1, 0.5, 0);
    p->SetAlpha(2,2, 0.6, 0);
  }

  return p;

}

//.............................................................................
OscProb::PMNS_Base* GetModel(string model, bool is_nominal = false){

  if(model == "Iter")    return GetIter(is_nominal);
  if(model == "Deco")    return GetDeco(is_nominal);
  if(model == "Sterile") return GetSterile(is_nominal);
  if(model == "Decay")   return GetDecay(is_nominal);
  if(model == "NSI")     return GetNSI(is_nominal);
  if(model == "LIV")     return GetLIV(is_nominal);
  if(model == "SNSI")    return GetSNSI(is_nominal);
  if(model == "NUNM")    return GetNUNM(is_nominal);

  return GetFast(is_nominal);

}

//.............................................................................
vector<string> GetListOfModels(){

  //return {"Decay"};

  return {"Fast", "Iter", "Sterile", "NSI",
          "Deco", "Decay", "LIV", "SNSI",
          "NUNM"};

}

//.............................................................................
void SetTestPath(OscProb::PMNS_Base* p){

  p->SetPath(1000, 2);
  p->AddPath(1000, 4);
  p->AddPath(1000, 2);

}

//.............................................................................
void SaveTestFile(OscProb::PMNS_Base* p, TString filename){

  SetTestPath(p);

  int nbins = 100;
  vector<double> xbins = GetLogAxis(nbins, 0.1, 10);
  TH1D* h = 0;

  TFile* f = new TFile("data/"+filename, "recreate");

  for(int flvi=0; flvi<3; flvi++){
  for(int flvf=0; flvf<3; flvf++){
  for(int isnb=0; isnb<2; isnb++){
    p->SetIsNuBar(isnb);
    TString hname = TString::Format("h%d%d%d",flvi,flvf,isnb);
    h = new TH1D(hname, "", nbins, &xbins[0]);
    for(int i=1; i<=nbins; i++){
      double energy = h->GetBinCenter(i);
      double dE = h->GetBinWidth(i);
      h->SetBinContent(i, p->AvgProb(flvi, flvf, energy, dE));
    }
    h->Write();
    delete h;
  }}}

  f->Close();
  cout << "Saved new test file: data/" + filename << endl;

}

//.............................................................................
int CheckProb(OscProb::PMNS_Base* p, TString filename){

  SetTestPath(p);

  TFile* f = new TFile("data/"+filename, "read");

  int ntests = 0;
  int fails = 0;

  TCanvas* c1 = 0;
  TH1D* h0 = 0;
  TH1D* h = 0;

  for(int flvi=0; flvi<3; flvi++){
  for(int flvf=0; flvf<3; flvf++){
  for(int isnb=0; isnb<2; isnb++){
    p->SetIsNuBar(isnb);
    TString hname = TString::Format("h%d%d%d",flvi,flvf,isnb);
    h0 = (TH1D*)f->Get(hname);
    h = (TH1D*)h0->Clone();
    bool plot = false;
    for(int i=1; i<=h0->GetNbinsX(); i++){
      double energy = h->GetBinCenter(i);
      double dE = h->GetBinWidth(i);
      double p0 = h0->GetBinContent(i);
      double p1 = p->AvgProb(flvi, flvf, energy, dE);
      ntests++;
      if(abs(p0-p1)>1e-12){
        plot = true;
        fails++;
      }
      h->SetBinContent(i, p1);
    }
    if(plot){
      c1 = new TCanvas();
      c1->Divide(1,2);
      c1->cd(1);
      SetHist(h0, kBlue);
      SetHist(h, kRed);
      TString nu_lab = isnb ? "#bar{#nu}" : "#nu";
      TString flv_lab[3] = {"e","#mu","#tau"};
      TString ylab = "P(" + nu_lab + "_{" + flv_lab[flvi] +
                     "}#rightarrow" + nu_lab + "_{" + flv_lab[flvf] + "})";
      h0->SetTitle(";Energy [GeV];"+ylab+";");
      h->SetLineStyle(7);
      double ymax = max(h->GetMaximum(),
                        h0->GetMaximum());
      double ymin = min(h->GetMinimum(),
                        h0->GetMinimum());
      h0->GetYaxis()->SetRangeUser(ymin, ymax);
      h0->DrawCopy("hist");
      h->DrawCopy("hist same");
      SetTH1Margin();
      gPad->SetLogx();
      c1->cd(2);
      h->SetTitle(";Energy [GeV];#Delta"+ylab+";");
      h->Add(h0,-1);
      h->SetLineStyle(1);
      h->DrawCopy("hist");
      SetTH1Margin();
      gPad->SetLogx();
      MiscText(0.55,0.8,0.1,filename,kGray,1);
      c1->DrawClone();
      TString pngfile = filename;
      pngfile.ReplaceAll(".root",".png");
      c1->SaveAs("plots/Failed_"+hname+"_"+pngfile);
      delete c1;
    }
    delete h0, h;
  }}}

  if(fails>0){
    printf((Color::FAILED + " Found %d differences in %d tests (%.3g%%) in %s\n").c_str(),
            fails, ntests, 100.*fails/ntests, filename.Data());
  }
  else {
    cout << Color::PASSED << " No differences found in " << filename << endl;
  }

  return fails;

}
