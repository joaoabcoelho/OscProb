
#ifndef __CINT__
#include "PremModel.h"
bool isCINT = false;
#else
bool isCINT = true;
#endif

// Some functions to make nice plots
#include "SetNiceStyle.C"

// Macro to load OscProb library
#include "LoadOscProb.C"

// Define the PREM tables path
#include "prem_default.hpp"

// Function to draw all layers from upgoing path
void DrawModel(TString opt = "alp", int col = kBlue, int model = 0);

// Plot the PREM model
void showModel(){

  // Load the library.
  // You can also include this in your rootrc file
  // The library must be in your path
  if(isCINT) LoadOscProb();

  // Set a nice overall style
  SetNiceStyle();

  // Draw two different models (15 and 44 layers)
  DrawModel("alp", kBlue, 0);
  DrawModel("lp" , kRed, 1);
  
  SetTH1Margin();

}

// Function to draw all layers from upgoing path
void DrawModel(TString opt, int col, int model){

  // Set as up-going
  double cosT = -1;

  // Get the model table paths
  string filename;

  if(model==0) filename = PREM_DIR + "/prem_44layers.txt";
  else         filename = PREM_DIR + "/prem_15layers.txt";
  
  // Set the PREM model from a table
  OscProb::PremModel prem(filename);
  
  // Get the prem layers
  vector<OscProb::PremLayer> pl = prem.GetPremLayers();

  // Number of layers
  int nlayers = pl.size();
  
  // Create a TGraph to plot
  TGraph* gr = new TGraph(2*nlayers);
  
  // Loop over paths
  for(int i=0; i<nlayers; i++){
    // Set graph point in start of layer
    if(i==0) gr->SetPoint(2*i, 0, pl[i].density);
    else     gr->SetPoint(2*i, pl[i-1].radius, pl[i].density);
    // Set graph point in end of layer
    gr->SetPoint(2*i+1, pl[i].radius, pl[i].density);
  }
  
  // Setup a nice TGraph
  SetGraph(gr,col);

  // Make 44 layer dashed
  if(model==0){
   gr->SetLineStyle(7);
   gr->SetMarkerStyle(24);
  }
  
  // Set the limits stopping at the earth's center
  gr->GetYaxis()->SetRangeUser(-0.1,15);
  gr->GetXaxis()->SetLimits(10,pl[nlayers-1].radius);
  
  // Set the axis titles
  gr->SetTitle(";Radius (km);Density (g/cm^{3})");
  
  // Draw a clone
  gr->DrawClone(opt);

}
