
// Some functions to make nice plots
#include "SetNiceStyle.C"

// Macro to load OscProb library
#include "LoadOscProb.C"

// Define the PREM tables path
#include "../prem_default.hpp"

// Function to draw paths fractions as a function of cosT
void DrawPath(double cosT, TString opt = "alp", int col = kBlue, int model = 0);

// Plot a series of paths as a
// function of cosTheta.
void showPaths(){

  // Load the OscProb library.
  LoadOscProb();

  // Set a nice overall style
  SetNiceStyle();

  // Set a color palette
  int p = SetRainbowPalette();

  // Tag first graph to plot
  bool isFirst = true;

  // Define cosT at layer boundaries
  vector<double> lys;
  lys.push_back(-1);      // Vertical (sees inner core)
  lys.push_back(-0.98);   // Outer-core boundary
  lys.push_back(-0.837);  // core-mantle boundary
  lys.push_back(-0.4455); // lower-mantle boundary
  lys.push_back(-0.0819); // mantle-crust boundary
  lys.push_back(0);       // Horizontal (Ocean and Air only)

  // Get number of graphs
  int ngr = lys.size();

  // Loop over layer boundaries
  for(int j=0; j<ngr; j++){ 
  
    // Get cosT
    double cosT = lys[j];
    
    // Set color from palette
    int col = p + 255*j/ngr;
    
    // Set drawing option
    TString opt = "lp";
    
    // If first, draw axes
    if(isFirst) opt = "a" + opt;
    
    // Draw the two models
    DrawPath(cosT, opt, col,0);
    DrawPath(cosT, "lp", col,1);
  
    // Tag as not first
    isFirst = false;
  
  }

}


// Function to draw paths fractions as a function of cosT
void DrawPath(double cosT, TString opt, int col, int model){

  // Get the model table paths
  string filename;

  if(model==0) filename = PREM_DIR + "/prem_425layers.txt";
  else         filename = PREM_DIR + "/prem_425layers.txt";

  // Set the PREM model from a table
  OscProb::PremModel prem(filename);
  
  // Fill the paths
  prem.FillPath(cosT);
  
  // Get the lengths and densities
  vector<OscProb::NuPath> p;
  if(model==0) p = prem.GetNuPath();
  else         p = prem.GetMergedPaths();

  // Number of paths
  int nsteps = p.size();
  
  // Create a TGraph to plot
  TGraph* gr = new TGraph(2*nsteps);
  
  // Variables to keep track of path
  double sumL = 0;

  // Total length for normalizing
  double totL = prem.GetTotalL(cosT);
  
  // Variable for setting text position
  double txtpos = 0;
  
  // Loop over paths
  for(int i=0; i<nsteps; i++){
    // Set graph point in path start
    gr->SetPoint(2*i, sumL/totL, p[i].density);
    // Set graph point in path end
    sumL += p[i].length;
    gr->SetPoint(2*i+1, sumL/totL, p[i].density);
    
    // Set text position in maximum point
    if(p[i].density > txtpos) txtpos = p[i].density;
  }
  
  // Setup a nice TGraph
  SetGraph(gr,col);

  // Make 44 layer dashed
  if(model==0){
   gr->SetLineStyle(7);
   gr->SetMarkerStyle(24);
  }
  
  // Set the limits
  gr->GetYaxis()->SetRangeUser(-0.1,15);
  gr->GetXaxis()->SetLimits(0,1);
  
  // Set the axis titles
  gr->SetTitle(";Path Fraction;Density (g/cm^{3})");
  
  // Draw a clone
  gr->DrawClone(opt);
  
  // Draw the cosT labels
  if(model){
  
    // Get NDC from y value
    txtpos = YtoNDC(txtpos) + 0.02;
    
    // Place text below for inner core boundary
    if(fabs(cosT+0.98) < 0.01) txtpos -= 0.09;
  
    // Draw the label
    MiscText(0.43, txtpos, 0.05, TString::Format("cos#theta_{z} = %0.2f", cosT), col);
    
  }

}

