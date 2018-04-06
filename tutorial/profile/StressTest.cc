#include <ctime>
#include <iostream>
#include "PremModel.h"
#include "PMNS_Fast.h"

using namespace std;

int main(int argc, char **argv){

  int minM = 0;
  int maxM = 2;

  if(argc>1){
    string test = argv[1];
    if(test=="new") maxM = 1;
    if(test=="old") minM = 1;
  }
  
  int ntries = 1e4;
  if(maxM-minM<2) ntries = 1e3;

  OscProb::PMNS_Fast p;
  
  // PREM Model
  OscProb::PremModel prem;

  // Fill path for cosT
  prem.FillPath(-1);

  // Give path to calculator
  p.SetPath(prem.GetNuPath());

  double oldTime = 0;
  double newTime = 0;
  
  for(int m=minM; m<maxM; m++){
    
    int myTime = clock();
  
    for(int i=0; i<ntries; i++){
  
      //p.SetDm(3, (1 -2*(i%2)) * 2.5e-3);
      //p.SetEnergy(2 + (i%2));
      p.SetOldProp(m);
    
      p.Prob(1,0);
  
    }
  
    double timeSec = double(clock() - myTime) / CLOCKS_PER_SEC;
    if(m) oldTime = timeSec;
    else  newTime = timeSec;

  }

  cout << "Old time = " << oldTime << endl;
  cout << "New time = " << newTime << endl;
  if(maxM-minM==2) cout << "Speedup  = " << oldTime/newTime << " times" << endl;
  
}
