
#include "TF1.h"
#include "TFile.h"

#include "PremModel.h"
#include "PMNS_Fast.h"
#include "PMNS_TaylorExp.h"

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
  cout<<"mmmmmmmm"<<endl;

  // Make the oscillogram
  TH2D* hNH = GetOscHist(flvf,1);

  // Draw the oscillogram
  hNH->Draw("colz");

  // Add space for the colz bar
  gPad->SetRightMargin(0.18);

  // Draw some lines of constant energy
  DrawEnergyLines(hNH);
  

}

// force the string to take only two decimal like in the file 
string decimal_precision (double value, double precision)
  {
    ostringstream returned_string;
    returned_string<<fixed<<setprecision(2)<<value;
    return returned_string.str();
  }

// Choose the right sentence  
string testCos (double cosT_min , double cosT_max)
{
  if(cosT_min<0){
    if(cosT_max<0){
      return string("average flux in [cosZ =") + decimal_precision(cosT_min,2) + " -- " + decimal_precision(cosT_max,2) + ", phi_Az =   0 -- 360]";
    }
    else{
      return string("average flux in [cosZ =") + decimal_precision(cosT_min,2) + " --  " + decimal_precision(cosT_max,2) + ", phi_Az =   0 -- 360]";
      
    }
  }
  else{
    return string("average flux in [cosZ = ") + decimal_precision(cosT_min,2) + " --  " + decimal_precision(cosT_max,2) + ", phi_Az =   0 -- 360]";
  }

  return 0;
}

// Copy all the flux data in a map 
void get_flux_data(map<string,map<double,map<int,map<int,double>>>> & flux_data){//PROBLEME PREND UNE LIGNE EN TROP A LA FIN 

  // Open the flux data file
  ifstream flux_file("frj-nu-20-01-000.d");

  // First range value of cosT in the file 
  double cosT_min = 0.90;
  double cosT_max = 1.00;

  string linee;
  string indice_cos;
  
  // Loop until the end of the file
  do{
    // Copy the current line of the file
    getline(flux_file, linee);  

    string test = testCos(cosT_min,cosT_max);
    
    // Test if we change the range of cosT
    if(linee == test){
      // Get the value of cosT
      indice_cos = linee;

      // Change the range value of cosT for later
      cosT_min -=0.1;
      cosT_max -=0.1;

      // Skip the line "Enu(GeV) NuMu  NuMubar NuE NuEbar (m^2 sec sr GeV)^-1"
      getline(flux_file, linee);
    }
    else{
      // Get the value of energy 
      istringstream col(linee);
      double indice_E;
      col >> indice_E;

      // Loop over all the columns except the first one (== energy)
      for(int flvi = 1; flvi>=0; flvi--){
        for(int nunubar = 1; nunubar>-2; nunubar-=2){

          // Get the value of one columns (== flux) 
          double flux_value;
          col >> flux_value;

          // register flux on a map 
          flux_data[indice_cos][indice_E][flvi][nunubar] = flux_value;
        }
      }
    }
  }while(!linee.empty());

  // Close the flux data file
  flux_file.close();
}

// Get all the energy data
void get_energy_flux_data (vector<double> & energy_flux_data){

  // Open the energy data file
  ifstream energy_flux_file("data_energy_flux.txt");

  string E;

  // Loop until the end of the file
  do{
    // Copy the current line of the file
    getline(energy_flux_file, E); 

    // Get the value of energys
    double value;
    istringstream e(E);
    e >> value;

    // register energies on a vector 
    energy_flux_data.push_back(value);

  }while(!E.empty());

  // Close the flux data file
  energy_flux_file.close();

}

// Get the energy index for the map
double get_index_E (double E , vector<double> energy_flux_data){  //problème si energy plus grande que celledes donnes 

  double index_E = 0;

  for(long unsigned int i= 0 ; i<energy_flux_data.size() ; i++){
    index_E = energy_flux_data[i];
    if(E<index_E){break;}
  }
  
  return index_E;
}

// Copy all the CS and energy data in a map 
void get_CS (map<double,map<int,double>> &CS_data , vector<double> &energy_CS_data , int flvf){

  // Open the CS data file
  TFile file_CS("crossSection.root");

  // Select graphs from the CS file 
  TGraph *histo_CS_num;
  TGraph *histo_CS_nbm;

  if(flvf == 1){
    histo_CS_num = (TGraph*)file_CS.Get("single_graphs/gnum_CC_E");
    histo_CS_nbm = (TGraph*)file_CS.Get("single_graphs/gnbm_CC_E");
  }
  else{
    histo_CS_num = (TGraph*)file_CS.Get("single_graphs/gnue_CC_E");
    histo_CS_nbm = (TGraph*)file_CS.Get("single_graphs/gnbe_CC_E");
  }
 
  // Get the number of point of these graphs
  int graph_size = histo_CS_num->GetMaxSize();
  
  // Loop over all these points
  for(int i=0 ; i<graph_size ; i++){

    // Get abscissa values (== energies)
    energy_CS_data.push_back(histo_CS_num->GetPointX(i));

    // Get ordinate values (== CS) 
    CS_data[energy_CS_data[i]][1] = histo_CS_num->GetPointY(i);
    CS_data[energy_CS_data[i]][-1] = histo_CS_nbm->GetPointY(i);
  }

  // Close the CS data file
  file_CS.Close();

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
  OscProb::PMNS_TaylorExp testPMNS;
  
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
  
  // Stock all the flux data
  //unordered map
  map<string,map<double,map<int,map<int,double>>>> flux_data;     //[cosT][E][flv][nunubar]
  get_flux_data(flux_data); 

  // Stock all the energy data for the flux 
  vector<double> energy_flux_data;
  get_energy_flux_data(energy_flux_data);

  // // Stock all the CS data and the energy data for the CS 
  map<double,map<int,double>> CS_data;     //[E][nunubar]
  vector<double> energy_CS_data;
  get_CS(CS_data,energy_CS_data,flvf);

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

    // Get the cosT index
    double cosT_min = floor(10*cosT)/10;
    double cosT_max = ceil(10*cosT)/10;
    if(cosT_max == 0){
      cosT_max = abs(cosT_max);
    }
    string index_cosT = testCos(cosT_min,cosT_max);

    // Loop of L/Es(theta_z) and L/E
    for(int le=1; le<=nbinsx; le++){

      // Set L/E from bin center
      double loe  = h2->GetXaxis()->GetBinCenter(le);

      // Get E from L and L/E
      double E = L/loe;  

      // Get the Energy Flux and CS index
      double index_E_flux = get_index_E(E , energy_flux_data);
      double index_E_CS = get_index_E(E , energy_CS_data);

      // Initialize probability
      double prob = 0;

      // Loop over initial flavour and nu or nubar
      for(int flvi = 1; flvi>=0; flvi--){
      for(int nunubar = -1; nunubar<2; nunubar+=2){

        double weight_flux_part = flux_data[index_cosT][index_E_flux][flvi][nunubar] / flux_data[index_cosT][index_E_flux][1][1];
        double weight_CS_part = CS_data[index_E_CS][nunubar] / CS_data[index_E_CS][1];

        // Define some basic weights for nue/numu and nu/nubar
        double weight = (0.75 + 0.25*nunubar) * (0.5 + 0.5*flvi);
        double weight_flux = (0.75 + 0.25*nunubar) * weight_flux_part;
        double weight_CS = weight_CS_part * (0.5 + 0.5*flvi);
        double weight_flux_CS = weight_CS_part * weight_flux_part;

        // Add probabilities from OscProb
        myPMNS.SetIsNuBar(nunubar <= 0);
        prob +=  weight_flux_CS  * myPMNS.Prob(flvi, flvf, L/loe);
        
      }}
    
      // Fill probabilities in histogram
      h2->SetBinContent(le,ct,prob);

    }// loe loop
  }// cosT loop 
  
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
