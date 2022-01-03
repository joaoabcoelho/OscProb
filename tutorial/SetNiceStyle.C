#include <iostream>
#include <fstream>
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TColor.h"
#include "TROOT.h"

using namespace std;

#ifndef SETNICESTYLE_H
#define SETNICESTYLE_H

// Define nice style for plots
void SetNiceStyle()
{

  // Defaults to classic style, but that's OK, we can fix it
  TStyle* niceStyle = new TStyle("niceStyle", "Nice Style");

  ////////////////////////////////
  /// These are my additions   ///
  ////////////////////////////////

  // Set the size of the default canvas
  niceStyle->SetCanvasDefH(600);
  niceStyle->SetCanvasDefW(730);
  niceStyle->SetCanvasDefX(10);
  niceStyle->SetCanvasDefY(10);

  //set marker style
  niceStyle->SetMarkerStyle(20);
  niceStyle->SetMarkerSize(1);

  // Set margins -- I like to shift the plot a little up and to the
  // right to make more room for axis labels
  niceStyle->SetPadTopMargin(0.05);
  niceStyle->SetPadBottomMargin(0.12);
  niceStyle->SetPadLeftMargin(0.12);
  //niceStyle->SetPadRightMargin(0.05); // Default without colz
  niceStyle->SetPadRightMargin(0.18); // Extra room for colz

  // Set Data/Stat/... and other options
  niceStyle->SetOptDate(0);
  //  niceStyle->SetDateX(0.1);
  //  niceStyle->SetDateY(0.1);
  niceStyle->SetOptFile(0);
  niceStyle->SetStatFormat("6.2f");
  niceStyle->SetFitFormat("8.4f");
  niceStyle->SetOptFit(1);
  niceStyle->SetStatH(0.20);
  niceStyle->SetStatStyle(0);
  niceStyle->SetStatW(0.30);
  niceStyle->SetStatX(0.845);
  niceStyle->SetStatY(0.845);
  niceStyle->SetOptTitle(0);
  niceStyle->SetTitleW(0.75);

  // Set paper size for life in the US
  niceStyle->SetPaperSize(TStyle::kUSLetter);

  niceStyle->SetTimeOffset(0);

  niceStyle->SetEndErrorSize(0);
  niceStyle->SetStripDecimals(false);

  ////////////////////////////////
  /// Begin default niceStyle  ///
  ////////////////////////////////

  // Centre title
  niceStyle->SetTitleAlign(22);
  niceStyle->SetTitleX(.5);
  niceStyle->SetTitleY(.95);
  niceStyle->SetTitleBorderSize(0);

  // No info box
  niceStyle->SetOptStat(0);

  //set the background color to white
  niceStyle->SetFillColor(10);
  niceStyle->SetFrameFillColor(10);
  niceStyle->SetCanvasColor(10);
  niceStyle->SetPadColor(10);
  niceStyle->SetTitleFillColor(0);
  niceStyle->SetStatColor(10);

  // Don't put a colored frame around the plots
  niceStyle->SetFrameBorderMode(0);
  niceStyle->SetCanvasBorderMode(0);
  niceStyle->SetPadBorderMode(0);

  // Set the default line color for a fit function to be red
  niceStyle->SetFuncColor(kRed);

  // Marker settings
  //  niceStyle->SetMarkerStyle(kFullCircle);

  // No border on legends
  niceStyle->SetLegendBorderSize(0);

  // Scientific notation on axes
  //  TGaxis::SetMaxDigits(3);

  // Axis titles
  niceStyle->SetTitleSize(.055, "xyz");
  niceStyle->SetTitleOffset(.8, "xyz");
  // More space for y-axis to avoid clashing with big numbers
  niceStyle->SetTitleOffset(.9, "y");
  // This applies the same settings to the overall plot title
  niceStyle->SetTitleSize(.055, "");
  niceStyle->SetTitleOffset(.8, "");
  // Axis labels (numbering)
  niceStyle->SetLabelSize(.04, "xyz");
  niceStyle->SetLabelOffset(.005, "xyz");

  // Thicker lines
  niceStyle->SetHistLineWidth(2);
  niceStyle->SetFrameLineWidth(2);
  niceStyle->SetFuncWidth(2);

  // Set the number of tick marks to show
  niceStyle->SetNdivisions(506, "xyz");

  // Set the tick mark style
  niceStyle->SetPadTickX(1);
  niceStyle->SetPadTickY(1);

  // Fonts
  const int kNiceFont = 42;
  niceStyle->SetStatFont(kNiceFont);
  niceStyle->SetLabelFont(kNiceFont, "xyz");
  niceStyle->SetTitleFont(kNiceFont, "xyz");
  niceStyle->SetTitleFont(kNiceFont, ""); // Apply same setting to plot titles
  niceStyle->SetTextFont(kNiceFont);
  niceStyle->SetLegendFont(kNiceFont);

  // Get moodier colours for colz
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  niceStyle->SetNumberContours(NCont);

  gROOT->SetStyle("niceStyle");

  // Uncomment this line if you want to force all plots loaded from files
  // to use this same style
  //gROOT->ForceStyle();
}

void SetTH1Margin(){

  gPad->SetRightMargin(0.05);

}

// Make a long canvas
TCanvas* MakeLongCanvas(){

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->SetBottomMargin(0.12);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.05);
  c1->SetTopMargin(0.05);

  return c1;

}

// Draw some text in position x,y
TLatex *MiscText(float x, float y, float size, TString text, int col = kBlack, bool centered=false)
{
   TLatex *l = new TLatex(x,y,text);
   l->SetTextAlign(centered ? 21 : 11);
   l->SetNDC();
   l->SetTextSize(size);
   l->SetTextColor(col);
   l->Draw();
   return l;
}

// Divide canvas nicely
void DivideCanvas(TCanvas *c1, int col, int row){

  c1->Divide(col,row);

  for(int i=1; i<=col*row; i++){
    c1->GetPad(i)->SetBottomMargin(0.19);
    c1->GetPad(i)->SetLeftMargin(0.19);
    c1->GetPad(i)->SetRightMargin(0.03);
    c1->GetPad(i)->SetTopMargin(0.035);
  }

  return;

}

// Get a log-scale axis
vector<double> GetLogAxis(int nbins, double xmin, double xmax){
  
  vector<double> xbins(nbins+1);

  for(int i=0; i<=nbins; i++){
    xbins[i] = pow(10, log10(xmin) + i*(log10(xmax) - log10(xmin))/nbins);
  }
  
  return xbins;

}

// Set some nice histogram style
void SetHist(TH1 *hist,int col=1,bool fill=false){

  hist->SetDirectory(0);

  hist->SetLineColor(col);  
  hist->SetLineWidth(4);  
  hist->SetFillColor(fill*col);  
  hist->SetMarkerColor(col);  
  hist->SetMarkerSize(0);  
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetZaxis()->CenterTitle();
  hist->SetTitleOffset(1.05,"XY");
  hist->SetTitleOffset(1.2,"Z");
  hist->SetTitleSize(0.05,"XYZ"); 
  hist->SetLabelSize(0.04,"XYZ");

  return;

}

// Set some nice THStack style
void SetStack(THStack *stack){

  stack->GetXaxis()->CenterTitle();
  stack->GetYaxis()->CenterTitle();
  stack->GetXaxis()->SetTitleOffset(1.05);
  stack->GetYaxis()->SetTitleOffset(1.05);
  stack->GetXaxis()->SetTitleSize(0.05); 
  stack->GetYaxis()->SetTitleSize(0.05); 
  stack->GetXaxis()->SetLabelSize(0.04);
  stack->GetYaxis()->SetLabelSize(0.04);

  return;

}

// Set some nice TH2 style
void SetHist(TH2 *hist){

  hist->SetDirectory(0);

  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetZaxis()->CenterTitle();
  hist->SetTitleOffset(1.05,"XY");
  hist->SetTitleOffset(1.2,"Z");
  hist->SetTitleSize(0.05,"XYZ"); 
  hist->SetLabelSize(0.04,"XYZ");
/*
  hist->SetContour(2);
  hist->SetContourLevel(0,2.296);
  hist->SetContourLevel(1,4.605);
  hist->SetLineWidth(2);
*/
  return;

}

// Set nice graph style
void SetGraph(TGraph *gr,int col=1){

  gr->SetLineColor(col);  
  gr->SetLineWidth(2);  
  gr->SetMarkerColor(col);  
  gr->SetMarkerStyle(20);  
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleOffset(1);
  gr->GetYaxis()->SetTitleOffset(1);
  gr->GetXaxis()->SetTitleSize(0.05); 
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetLabelSize(0.05);

  return;

}

// Make bins proportional to area
void SetDensity(TH1 *hist){

  int Nbins = hist->GetNbinsX();

  for(int i=1;i<=Nbins;i++){
    hist->SetBinContent(i,hist->GetBinContent(i)/hist->GetBinWidth(i));
    hist->SetBinError(i,hist->GetBinError(i)/hist->GetBinWidth(i));
  }

  return;

}

// Make bins proportional to area
void SetDensity(TH2 *hist){

  int nbinsx = hist->GetNbinsX();
  int nbinsy = hist->GetNbinsY();

  for(int i=1;i<=nbinsx;i++){
  for(int j=1;j<=nbinsy;j++){
    double wx = hist->GetXaxis()->GetBinWidth(i);
    double wy = hist->GetYaxis()->GetBinWidth(j);
    hist->SetBinContent(i, j, hist->GetBinContent(i, j)/wx/wy);
    hist->SetBinError(i,hist->GetBinError(i,j)/wx/wy);
  }}

  return;

}

// Set nice legend style
void SetLeg(TLegend *leg){

  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);

  leg->SetY1(leg->GetY2()-leg->GetNRows()*0.05);

  return;

}

// Poisson 68% C.L. for asymmetric error bars
void GetPoissonError(int n,double &eh,double &el){

  Double_t jimeu[51];
  Double_t jimed[51];
  jimeu[0]=  1.841;  jimed[0]= 0.000 ;
//  jimeu[0]=  1.148;  jimed[0]= 0.000 ;
  jimeu[1]=  3.3;    jimed[1]= 0.173 ;
  jimeu[2]=  4.638;  jimed[2]= 0.708 ;
  jimeu[3]=  5.918;  jimed[3]= 1.367 ;
  jimeu[4]=  7.163;  jimed[4]= 2.086 ;
  jimeu[5]=  8.382;  jimed[5]= 2.84  ;
  jimeu[6]=  9.584;  jimed[6]= 3.62  ;
  jimeu[7]=  10.77;  jimed[7]= 4.419 ;
  jimeu[8]=  11.95;  jimed[8]= 5.232 ;
  jimeu[9]=  13.11;  jimed[9]= 6.057 ;
  jimeu[10]= 14.27;  jimed[10]=6.891 ;
  jimeu[11]= 15.42;  jimed[11]=7.734 ;
  jimeu[12]= 16.56;  jimed[12]=8.585 ;
  jimeu[13]= 17.7;   jimed[13]=9.441;
  jimeu[14]= 18.83;  jimed[14]=10.3  ;
  jimeu[15]= 19.96;  jimed[15]=11.17 ;
  jimeu[16]= 21.08;  jimed[16]=12.04 ;
  jimeu[17]= 22.2 ;  jimed[17]=12.92 ;
  jimeu[18]= 23.32;  jimed[18]=13.8  ;
  jimeu[19]= 24.44;  jimed[19]=14.68 ;
  jimeu[20]= 25.55;  jimed[20]=15.57 ;
  jimeu[21]= 26.66;  jimed[21]=16.45 ;
  jimeu[22]= 27.76;  jimed[22]=17.35 ;
  jimeu[23]= 28.87;  jimed[23]=18.24 ;
  jimeu[24]= 29.97;  jimed[24]=19.14 ;
  jimeu[25]= 31.07;  jimed[25]=20.03 ;  
  jimeu[26]= 32.16;  jimed[26]=20.93 ;
  jimeu[27]= 33.26;  jimed[27]=21.84 ;
  jimeu[28]= 34.35;  jimed[28]=22.74 ;
  jimeu[29]= 35.45;  jimed[29]=23.65 ;
  jimeu[30]= 36.54;  jimed[30]=24.55 ;
  jimeu[31]= 37.63;  jimed[31]=25.46;
  jimeu[32]= 38.72;  jimed[32]=26.37;
  jimeu[33]= 39.80;  jimed[33]=27.28;
  jimeu[34]= 40.89;  jimed[34]=28.20;
  jimeu[35]= 41.97;  jimed[35]=29.11;
  jimeu[36]= 43.06;  jimed[36]=30.03;
  jimeu[37]= 44.14;  jimed[37]=30.94;
  jimeu[38]= 45.22;  jimed[38]=31.86;
  jimeu[39]= 46.30;  jimed[39]=32.78;
  jimeu[40]= 47.32;  jimed[40]=33.70;
  jimeu[41]= 48.36;  jimed[41]=34.62;
  jimeu[42]= 49.53;  jimed[42]=35.55;
  jimeu[43]= 50.61;  jimed[43]=36.47;
  jimeu[44]= 51.68;  jimed[44]=37.39;
  jimeu[45]= 52.76;  jimed[45]=38.32;
  jimeu[46]= 53.83;  jimed[46]=39.24;
  jimeu[47]= 54.90;  jimed[47]=40.17;
  jimeu[48]= 55.98;  jimed[48]=41.10;
  jimeu[49]= 57.05;  jimed[49]=42.02;
  jimeu[50]= 58.12;  jimed[50]=42.95;

  if(n<51){ 
    eh = jimeu[n]-n;
    el = n-jimed[n];
  }
  else{
    eh = sqrt(n + 0.25) + 0.5;
    el = sqrt(n + 0.25) - 0.5;
  }

  return;

}

// Set TGraphAsymmErrors with poisson error bars 
void SetGraphErrors(TGraphAsymmErrors *gr){

  for(int i=0;i<gr->GetN();i++){
    double x,y;
    gr->GetPoint(i,x,y);
    double ey = gr->GetErrorYhigh(i);
    double ex = gr->GetErrorXhigh(i);
    double ph,pl;
    double a = ey*ey/y;
    int n = pow(y/ey,2)+0.5;
    GetPoissonError(n,ph,pl);
    ph *= a;
    pl *= a;
    gr->SetPointEYhigh(i,ph);
    gr->SetPointEYlow(i,pl);
  }

  return;

}

int SetNicePalette(char* filename){

  ifstream ifs(filename);
  
  double r, g, b;
  
  vector<double> stops;
  vector<double> red;
  vector<double> green;
  vector<double> blue;
  
  while(ifs >> r >> g >> b){
     red.push_back(double(r)/255);
     green.push_back(double(g)/255);
     blue.push_back(double(b)/255);
  }
  
  int ncols = red.size();
  for(int i=0; i<ncols; i++){
    stops.push_back(double(i)/ncols);
  }

  int NCont = 255;
  int p = TColor::CreateGradientColorTable(ncols, &stops[0], &red[0], &green[0], &blue[0], NCont);
  gStyle->SetNumberContours(NCont);

  return p;

}

int SetNicePalette(int ncols, int* cols){

  if(ncols<2) return 0;

  int NCont = 255;
  
  vector<double> stops(ncols);
  vector<double> red(ncols);
  vector<double> green(ncols);
  vector<double> blue(ncols);
  
  for(int i=0; i<ncols; i++){
    TColor* c = gROOT->GetColor(cols[i]);
    stops[i] = double(i)/(ncols-1);
    red[i] = c->GetRed();
    green[i] = c->GetGreen();
    blue[i] = c->GetBlue();
  }
  
  int p = TColor::CreateGradientColorTable(ncols, &stops[0], &red[0], &green[0], &blue[0], NCont);
  gStyle->SetNumberContours(NCont);

  return p;

}

// A nice rainbow color palette
int SetRainbowPalette(){

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  int p = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  return p;

}

// Set ROOT6 kLightTemperature palette
int SetLightTemperature(){

  const Int_t NRGBs = 9;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
  Double_t red[NRGBs]   = {  31./255.,  71./255., 123./255., 160./255., 210./255., 222./255., 214./255., 199./255., 183./255.};
  Double_t green[NRGBs] = {  40./255., 117./255., 171./255., 211./255., 231./255., 220./255., 190./255., 132./255.,  65./255.};
  Double_t blue[NRGBs]  = { 234./255., 214./255., 228./255., 222./255., 210./255., 160./255., 105./255.,  60./255.,  34./255.};
  int p = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  return p;

}

// Set ROOT6 kTemperatureMap palette
int SetTemperatureMap(){

  const Int_t NRGBs = 9;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
  Double_t red[NRGBs]   = {  34./255.,  70./255., 129./255., 187./255., 225./225., 226./255., 216./255., 193./255., 179./255.};
  Double_t green[NRGBs] = {  48./255.,  91./255., 147./255., 194./255., 226./226., 229./255., 196./255., 110./255.,  12./255.};
  Double_t blue[NRGBs]  = { 234./255., 212./255., 216./255., 224./255., 206./206., 110./255.,  53./255.,  40./255.,  29./255.};
  int p = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  return p;

}

// Set ROOT6 kBird palette
int SetBirdPalette(){

  const Int_t NRGBs = 9;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
  Double_t red[NRGBs]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  Double_t green[NRGBs] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  Double_t blue[NRGBs]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
  int p = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  return p;

}


// A nice blue white red color palette
int SetBlueRedPalette(){

  const Int_t NRGBs = 3;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.50, 1.00 };
  Double_t red[NRGBs]   = { 0.20, 1.00, 1.00};
  Double_t green[NRGBs] = { 0.20, 1.00, 0.20};
  Double_t blue[NRGBs]  = { 1.00, 1.00, 0.20};
  int p = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  return p;

  /*
  int cols[] = {kBlue, kWhite, kRed};
  
  return SetNicePalette(3, cols);
  */

}

// A nice temperature palette
int SetNiceTempPalette(){

  int cols[] = {kBlue+3, kBlue, kAzure-4, kWhite, kOrange, kRed, kRed+3};

  return SetNicePalette(7, cols);

}

// A nice white red color palette
int SetWhiteRedPalette(){

  int cols[] = {kWhite, kOrange, kRed};

  return SetNicePalette(3, cols);

}

// A nice blue white color palette
int SetBlueWhitePalette(){

  int cols[] = {kBlue, kAzure-4, kWhite};

  return SetNicePalette(3, cols);

}

// A nice blue white color palette
int SetWhiteBluePalette(){

  int cols[] = {kWhite, kAzure-4, kBlue};

  return SetNicePalette(3, cols);

}

// Draw a progress bar to follow loops
void ProgBar(int k, int numk, double perc = 5){

    int ck = numk*0.01*perc + 0.5;
    if(ck==0) ck= 1;

    if(k%ck==0){

      int npts = numk/ck;
      //if(numk%ck>0) npts++;

      cout << "\r[";
      for(int i=0; i<npts; i++){
        if(i <= k/ck) cout << "=";
        else cout << " ";
      }
      cout << "]" << flush;

    }
    if(k==numk-1) cout << endl;

}

// Get NDC from x value
double XtoNDC(double x){

  if(!gPad) return 0;
  
  gPad->SetBatch(true);
  gPad->Update();

  double dpx  = gPad->GetX2() - gPad->GetX1();
  double xp1  = gPad->GetX1();

  gPad->SetBatch(false);

  return (x-xp1)/dpx;

}

// Get NDC from y value
double YtoNDC(double y){

  if(!gPad) return 0;

  gPad->SetBatch(true);
  gPad->Update();

  double dpy  = gPad->GetY2() - gPad->GetY1();
  double yp1  = gPad->GetY1();

  gPad->SetBatch(false);

  return (y-yp1)/dpy;

}

#endif
