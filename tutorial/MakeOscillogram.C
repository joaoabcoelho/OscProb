
#include "TF1.h"

#include "PremModel.h"
#include "PMNS_Fast.h"

// Some functions to make nice plots
#include "SetNiceStyle.C"

// Make oscillogram for given final flavour and MH
TH2D* GetOscHist(int flvf = 1, int mh = 1);

// Draw energy lines
void DrawEnergyLines(TH2* hNH);

// Make an oscillogram for NH
// nue (0), numu (1) or nutau (2)
void MakeOscillogram(int flvf = 1){

  // Set a nice overall style
  SetNiceStyle();

  // Make the oscillogram
  TH2D* hNH = GetOscHist(flvf,1);

  // Draw the oscillogram
  hNH->Draw("colz");

  // Add space for the colz bar
  gPad->SetRightMargin(0.18);

  // Draw some lines of constant energy
  DrawEnergyLines(hNH);

}

// Make oscillogram for given final flavour and MH
TH2D* GetOscHist(int flvf, int mh){

  // Use 200 x bins and 100 y bins
  int nbinsx = 200;
  int nbinsy = 100;

  // Set parameters to PDG
  double dm21 = 7.5e-5;
  double dm31 = mh>0 ? 2.457e-3 : -2.449e-3 + dm21;
  double th12 = asin(sqrt(0.304));
  double th13 = asin(sqrt(mh>0 ? 0.0218 : 0.0219));
  double th23 = asin(sqrt(mh>0 ? 0.452 : 0.579));
  double dcp  = (mh>0 ? 306 : 254)*TMath::Pi()/180;

  // Create PMNS object
  OscProb::PMNS_Fast myPMNS;

  // Set PMNS parameters
  myPMNS.SetDm(2, dm21);
  myPMNS.SetDm(3, dm31);
  myPMNS.SetAngle(1,2, th12);
  myPMNS.SetAngle(1,3, th13);
  myPMNS.SetAngle(2,3, th23);
  myPMNS.SetDelta(1,3, dcp);

  // The oscillogram histogram
  TH2D* h2 = new TH2D("","",nbinsx,0,50*nbinsx,nbinsy,-1,0);

  // Create default PREM Model
  OscProb::PremModel prem;

  // Open the flux data file
  ifstream flux("frj-nu-20-01-000.d");
  
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

    // Loop of L/Es(theta_z) and L/E
    for(int le=1; le<=nbinsx; le++){

      // Set L/E from bin center
      double loe  = h2->GetXaxis()->GetBinCenter(le);

      // Get E from L and L/E
      double E = L/loe;  
      
      // Initialize probability
      double prob = 0;
      
      int n_line = 930 - 103*(1+floor(10*cosT));
      string line;
      string E_line_string;
      double E_line = 0;
      int actualLine = 0;

      while(actualLine != n_line || E>=E_line){
        getline(flux, line);
        if(actualLine == n_line){
          istringstream bbb(line);
          bbb >> E_line_string;
          E_line=stod(E_line_string);
          n_line++;
          cout<<E<<"   "<<E_line<<endl;
        }
        actualLine++;
        //cout<<actualLine2<<endl;
      }
      cout<<line<<endl;
      
      flux.clear();
      flux.seekg(0);

      vector<double> flux;

      // Loop over initial flavour and nu or nubar
      for(int flvi = 1; flvi>=0; flvi--){
      for(int nunubar = -1; nunubar<2; nunubar+=2){

        istringstream aaa(line);
        int n_col = 4.5 - 2*flvi -0.5*nunubar;
        int actualcol{0};
        string value;
        while(actualcol != n_col){
          aaa >> value;
          actualcol++;
        }
        flux.push_back(stod(value));
        cout<<flux.size()<<"   ";

        double weight_flux = 0;
        if(flvi == 1){
          weight_flux = 1;
        }
        else{
          weight_flux = flux[2.5+nunubar*0.5]/flux[0.5+nunubar*0.5];
        }
        cout<<weight_flux<<endl;
        

        // Define some basic weights for nue/numu and nu/nubar
        double weight = (0.75 + 0.25*nunubar) * (0.5 + 0.5*flvi);

        double weight_upgrade = (0.75 + 0.25*nunubar) * weight_flux;

        // Add probabilities from OscProb
        myPMNS.SetIsNuBar(nunubar <= 0);
        //prob += weight_upgrade*myPMNS.Prob(flvi, flvf, L/loe);
        prob += prob += (weight-weight_upgrade) * myPMNS.Prob(flvi, flvf, L/loe);

      }}
      cout<<"-----------------------------------"<<endl;

      flux.clear();

      // Fill probabilities in histogram
      h2->SetBinContent(le,ct,prob);

    }// loe loop
  }// cosT loop

  // Close the flux data file
  flux.close();

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
