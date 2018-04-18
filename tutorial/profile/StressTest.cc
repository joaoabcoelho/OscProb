#include <ctime>
#include <iostream>
#include "PremModel.h"
#include "PMNS_Fast.h"
//#include "gperftools/profiler.h"

using namespace std;

int main(int argc, char **argv){

  //ProfilerStart("cpu.prof");

  int minM = 0;
  int maxM = 2;

  if(argc>1){
    string test = argv[1];
    if(test=="new") maxM = 1;
    if(test=="old") minM = 1;
  }
  
  OscProb::PMNS_Fast p;
  
  // PREM Model
  OscProb::PremModel prem;

  // Fill path for cosT
  prem.FillPath(-1);

  // Give path to calculator
  //p.SetPath(prem.GetNuPath());
  
  cout << "Using " << p.GetPath().size() << " layers..." << endl;

  double oldTime = 0;
  double newTime = 0;
  
  for(int m=minM; m<maxM; m++){
    
    int myTime = clock();
    int count = 0;
  
    //for(int i=0; i<ntries; i++){
    while(clock() - myTime < 5*CLOCKS_PER_SEC){
  
      //p.SetDm(3, (1 -2*(count%2)) * 2.5e-3);
      //p.SetEnergy(2 + (count%2));
      p.SetOldProp(m);
    
      p.Prob(1,0);
      
      count++;
  
    }
  
    double timeSec = double(clock() - myTime) / CLOCKS_PER_SEC;
    if(m){
      oldTime = 1e6*timeSec / count;
      cout << "Old time = " << oldTime << " µs/prob" << endl;
    }
    else {
      newTime = 1e6*timeSec / count;
      cout << "New time = " << newTime << " µs/prob" << endl;
    }

  }

  if(maxM-minM==2) cout << "Speedup  = " << oldTime/newTime << " times" << endl;
  
  //ProfilerStop();

}
