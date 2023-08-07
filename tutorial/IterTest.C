
#include "PremModel.h"
#include "PMNS_Fast.h"
#include "PMNS_Iter.h"

#include "SetNiceStyle.C"

//int CLOCKS_PER_SEC = 1e6;

struct TimeIt {

  TimeIt() { reset(); }

  int start;
  int count;
  double time(){ return double(clock() - start) / CLOCKS_PER_SEC; }
  void reset(){ count = 0; start = clock(); }
  void Print(){
    double tpi = time() / count;
    string scale = "s";
    if(tpi<1){ tpi *= 1e3; scale = "ms"; }
    if(tpi<1){ tpi *= 1e3; scale = "µs"; }
    //if(tpi<1){ tpi *= 1e3; scale = "ns"; }
    cout << "Performance = " << tpi << " " << scale << "/iteration" << endl;
  }

};


void IterTest(){

  // Initialize your objects
  OscProb::PMNS_Fast pStd;
  OscProb::PMNS_Iter pItr;

  pStd.SetUseCache();

  pItr.SetPrec(0.05);

  // Use a PremModel to make paths
  // through the earth
  OscProb::PremModel prem;

  // Chose an angle for the neutrino
  // and fill the paths with cosTheta
  // e.g. cosTheta = -1 (vertical up-going)
  prem.FillPath(-0.8);

  // Give the path to the PMNS object
  // and get the probability
  pStd.SetPath(prem.GetNuPath());
  pItr.SetPath(prem.GetNuPath());

  // Set other quantities
  //pStd.SetPath(OscProb::NuPath(1300,2.5));
  //pItr.SetPath(OscProb::NuPath(1300,2.5));

  int nbins = 500;
  vector<double> xbins = GetLogAxis(nbins, 0.5, 20);
  TH1D* hStd = new TH1D("","", nbins, &xbins[0]);
  TH1D* hItr = new TH1D("","", nbins, &xbins[0]);

  TimeIt time;
  for(int k=0; k<1000; k++){ if(time.time()>1) break;
  for(int i=1; i<=nbins; i++){
    double energy = hStd->GetBinCenter(i);
    hStd->SetBinContent(i,pStd.AvgProb(1,0,energy,0.2*energy));
    time.count++;
  }}
  cout << "========== Standard ==========" << endl;
  time.Print();

  time.reset();
  for(int k=0; k<1000; k++){ if(time.time()>1) break;
  for(int i=1; i<=nbins; i++){
    double energy = hItr->GetBinCenter(i);
    hItr->SetBinContent(i,pItr.AvgProb(1,0,energy,0.2*energy));
    time.count++;
  }}
  cout << "========== Iterative ==========" << endl;
  time.Print();

  SetHist(hStd, kBlack);
  SetHist(hItr, kRed);

  TCanvas* c1 = new TCanvas();

  hStd->DrawCopy("curv");
  hItr->DrawCopy("curv same");

  SetTH1Margin();
  gPad->SetLogx();

  TCanvas* c2 = new TCanvas();

  hItr->Add(hStd, -1);

  hItr->Draw("curv");

  SetTH1Margin();
  gPad->SetLogx();

}
