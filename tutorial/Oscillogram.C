
#include "TF1.h"

#include "PremModel.h"
#include "PMNS_Fast.h"
#include "PMNS_TaylorExp.h"

// Some functions to make nice plots
#include "SetNiceStyle.C"


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
    cout << "Total time = " << time() << endl;
  }

};


// Make oscillogram for given final flavour and MH
TH2D* GetOscHist(int flvf = 1, int mh = 1 , int nbinsx =200 , int nbinsy = 100 , 
                string method = "none" , vector<double> & timeC = *(new vector<double>()));

// Draw energy lines
void DrawEnergyLines(TH2* hNH);

// Make an oscillogram for NH
// nue (0), numu (1) or nutau (2)
vector<double> Oscillogram(int flvf = 1 , int nbinsx = 200 , int nbinsy = 100 , string method = "none"){

  // Set a nice overall style
  SetNiceStyle();

  vector<double> timeC;

  // Make the oscillogram
  TH2D* hNH = GetOscHist(flvf,1, nbinsx , nbinsy , method , timeC);

  // Draw the oscillogram
  hNH->Draw("colz");

  // Add space for the colz bar
  gPad->SetRightMargin(0.18);

  // Draw some lines of constant energy
  DrawEnergyLines(hNH);

  return timeC;

}

// Make oscillogram for given final flavour and MH
TH2D* GetOscHist(int flvf, int mh, int nbinsx , int nbinsy , string method , vector<double> & timeC){

  // Use 200 x bins and 100 y bins
  //int nbinsx = 200;
  //int nbinsy = 100;

  // Set parameters to PDG
  double dm21 = 7.5e-5;
  double dm31 = mh>0 ? 2.457e-3 : -2.449e-3 + dm21;
  double th12 = asin(sqrt(0.304));
  double th13 = asin(sqrt(mh>0 ? 0.0218 : 0.0219));
  double th23 = asin(sqrt(mh>0 ? 0.452 : 0.579));
  double dcp  = (mh>0 ? 306 : 254)*TMath::Pi()/180;

  // Create PMNS object
  OscProb::PMNS_TaylorExp myPMNS;
  

  // Set PMNS parameters
  myPMNS.SetDm(2, dm21);
  myPMNS.SetDm(3, dm31);
  myPMNS.SetAngle(1,2, th12);
  myPMNS.SetAngle(1,3, th13);
  myPMNS.SetAngle(2,3, th23);
  myPMNS.SetDelta(1,3, dcp);

  // Set histogram parameters
  double xmin = 0;
  double xmax = 10000;
  double ymin = -0.6;
  double ymax = 0;
  double widthBinX = (xmax-xmin) / nbinsx;
  double widthBinY = (ymax-ymin) / nbinsy; 

  // The oscillogram histogram
  TH2D* h2 = new TH2D("","",nbinsx,xmin,xmax,nbinsy,ymin,ymax);
  //TH2D* h2 = new TH2D("","",nbinsx,0,50*nbinsx,nbinsy,-1,0);

  // Create default PREM Model
  OscProb::PremModel prem;

  TimeIt time;
  //for(int k=0; k<1000; k++){ if(time.time()>1) break;

  // Loop over cos(theta_z) and L/E
  for(int ct=1; ct<=nbinsy; ct++){

    // Get cos(theta_z) from bin center
    double cosT = h2->GetYaxis()->GetBinCenter(ct);

    // Set total path length L
    double L = prem.GetTotalL(cosT);

    // Skip if cosT is unphysical
    if(cosT < -1 || cosT > 1) continue;

    // Fill paths from PREM model
    prem.FillPath(cosT);

    // Set paths in OscProb
    myPMNS.SetPath(prem.GetNuPath());

    // Loop of L/E
    for(int le=1; le<=nbinsx; le++){

      // Set L/E from bin center
      double loe  = h2->GetXaxis()->GetBinCenter(le);

      double widthBinXforE = L * (1 / (loe - widthBinX / 2) - 1 / (loe + widthBinX / 2));

      // Initialize probability
      double prob = 0;

      // Loop over initial flavour and nu or nubar
      for(int flvi = 0; flvi<2; flvi++){
      for(int nunubar = -1; nunubar<2; nunubar+=2){

        // Define some basic weights for nue/numu and nu/nubar
        double weight = (0.75 + 0.25*nunubar) * (0.5 + 0.5*flvi);

        // Add probabilities from OscProb
        myPMNS.SetIsNuBar(nunubar <= 0);
        //prob += weight*myPMNS.Prob(flvi, flvf, L/loe);

        //PRB VIENT D'IC !!!!!!!!!!!!!!!!!!!!!!!!!!!
        //cout<<myPMNS.avgProbTaylorAngle(flvi, flvf,cosT,widthBinY)<<endl;

        if(method == "fast") {prob += weight*myPMNS.AvgProbLoE(flvi, flvf, loe ,widthBinX);}
        if(method == "taylor") {prob += weight*myPMNS.avgProbTaylorAngle(flvi, flvf,cosT,widthBinY);} //L/loe ,widthBinXforE,

        // ICI L/E ET PAS E
        // ESSAYER DECHG DANS CE CODE LE to E
        // ATTENTION AvgProb utilise la converstion E to LoE MAIS pas avgProbTaylor

        prob -= weight*myPMNS.Prob(flvi, flvf, L/loe);

        time.count++;
      }}

      // Fill probabilities in histogram
      h2->SetBinContent(le,ct,prob);

    }// loe loop
  }// cosT loop

  time.Print();
  double t = time.time();
  timeC.push_back(t);
  timeC.push_back(t / time.count);

  // Set nice histogram
  SetHist(h2);

  // Set titles
  h2->SetTitle(";L/E (km/GeV);cos#theta_{z};P_{#mu#mu} + 0.5#timesP_{#bar{#mu#mu}} + 0.5#timesP_{e#mu} + 0.25#timesP_{#bar{e#mu}}");

  return h2;

}

// Get pad NDC from x value
double GetNDCx(double x) {

  return (x - gPad->GetX1())/(gPad->GetX2()-gPad->GetX1());

}

// Get pad NDC from y value
double GetNDCy(double y) {

  return (y - gPad->GetY1())/(gPad->GetY2()-gPad->GetY1());

}

// Draw energy lines
void DrawEnergyLines(TH2* hNH){

  // Get max value of x-axis
  double xmax = hNH->GetXaxis()->GetXmax();

  // Define earth diameter
  double dEarth = 2*6371;

  // Define constant energy line
  TF1* f= new TF1("f","-x*[0]/[1]",0,xmax);

  // Second parameter is earth's diameter
  f->SetParameter(1,dEarth);

  // Dashed black lines
  f->SetLineStyle(7);
  f->SetLineColor(kBlack);

  // Define min and max energies
  // for drawing lines
  double minE = 0.1 * dEarth/xmax;
  double maxE = 10. * dEarth/xmax;

  // Start from a rounded power of 10 energy
  double sampleE = pow(10, floor(log10(minE)) );

  // Define how to increase energies
  int idx = 0;
  double scale[3] = {2,2.5,2};

  // Make a set of energies to draw
  vector<double> energies;
  while(sampleE<maxE){

    if(sampleE>minE){
      energies.push_back(sampleE);
    }

    sampleE *= scale[idx%3];
    idx++;
  }

  // Update pad to get correct axis ranges
  gPad->Update();

  // Loop over energies and draw lines
  for(int i=0; i<int(energies.size()); i++){

    double energy = energies[i];

    // Set first parameter to energy
    f->SetParameter(0, energy);

    // Draw line
    f->DrawCopy("same");

    // Get position in pad to draw label
    double xval = 2*xmax*6371/(2*6371+xmax*energy);
    double ndcx = GetNDCx(xval);
    double ndcy = GetNDCy(f->Eval(xval));

    // Draw labels
    if(energy>=1) MiscText(ndcx,ndcy,0.03,TString::Format("%d GeV", int(energy)));
    if(energy<1)  MiscText(ndcx,ndcy,0.03,TString::Format("%0.1f GeV", energy));

  }// energies loop

}