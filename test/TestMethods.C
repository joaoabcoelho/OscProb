#include "Utils.h"

bool TestAvgProb(OscProb::PMNS_Base* p, double energy, double &max_diff){
 
  double dE = 0.1 * energy;
  double prec = 1e-4;
  p->SetAvgProbPrec(prec);
  
  for(int flvi=0; flvi<3; flvi++){
  for(int flvf=0; flvf<3; flvf++){

    double avg_prob = p->AvgProb(flvi, flvf, energy, dE);

    int npts = 1000;
    double avg_test = 0;
    for(int i=0; i<npts; i++){
      double energy_i = energy + dE*((i+0.5)/npts - 0.5);
      avg_test += p->Prob(flvi, flvf, energy_i) / npts;
    }
    
    if(abs(avg_prob - avg_test) > max_diff){
      max_diff = abs(avg_prob - avg_test);
    }

    if(abs(avg_prob - avg_test) > 5*prec){
      cerr << "Found mismatch in AvgProb(" 
           << flvi << ", " << flvf << ", "
           << energy << ", " << dE << "): "
           << avg_prob << " - " << avg_test
           << " = " << avg_prob - avg_test
           << "; with IsNuBar = " << p->GetIsNuBar()
           << endl;
      return false;
    }

  }}
  
  return true;

}

bool TestParallel(OscProb::PMNS_Base* p, double energy){
 
  p->SetEnergy(energy);

  OscProb::matrixD pmatrix = p->ProbMatrix(3,3);
 
  for(int flvi=0; flvi<3; flvi++){

    OscProb::vectorD pvector = p->ProbVector(flvi);

    for(int flvf=0; flvf<3; flvf++){

      double prob = p->Prob(flvi, flvf);

      if(abs(prob - pmatrix[flvi][flvf]) > 1e-12){
        cerr << "Found mismatch in ProbMatrix(" 
             << flvi << ", " << flvf << "): "
             << prob << " - " << pmatrix[flvi][flvf]
             << " = " << prob - pmatrix[flvi][flvf]
             << "; with IsNuBar = " << p->GetIsNuBar()
             << " and E = " << p->GetEnergy()
             << endl;
        return false;
      }

      if(abs(prob - pvector[flvf]) > 1e-12){
        cerr << "Found mismatch in ProbVector(" 
             << flvi << ", " << flvf << "): "
             << prob << " - " << pvector[flvf]
             << " = " << prob - pvector[flvf]
             << "; with IsNuBar = " << p->GetIsNuBar()
             << " and E = " << p->GetEnergy()
             << endl;
        return false;
      }

    }

  }
  
  return true;

}

bool TestMethodsModel(OscProb::PMNS_Base* p, string model){

  SetTestPath(p);

  int nbins = 30;
  vector<double> xbins = GetLogAxis(nbins, 0.1, 100);
  double max_diff = 0;
  for(int isnb=0; isnb<2; isnb++){

    p->SetIsNuBar(isnb);
      
    for(int i=0; i<=nbins; i++){

      double energy = xbins[i];
      
      if(!TestAvgProb(p, energy, max_diff)){
        cerr << "FAILED: Model PMNS_" << model 
             << " AvgProb" << endl;
        return false;
      }

      if(!TestParallel(p, energy)){
        cerr << "FAILED: Model PMNS_" << model 
             << " parallel methods" << endl;
        return false;
      }

    }

  }
  
  cout << "PASSED: PMNS_" << model 
       << " with maximum AvgProb error: "
       << max_diff << endl;
  
  return true;

}

int TestMethods(){

  OscProb::PMNS_Fast p_fast = GetFast();
  if(!TestMethodsModel(&p_fast, "Fast")) return 1;

  OscProb::PMNS_Iter p_iter = GetIter();
  if(!TestMethodsModel(&p_iter, "Iter")) return 1;

  OscProb::PMNS_Sterile p_sterile = GetSterile();
  if(!TestMethodsModel(&p_sterile, "Sterile")) return 1;

  OscProb::PMNS_NSI p_nsi = GetNSI();
  if(!TestMethodsModel(&p_nsi, "NSI")) return 1;

  OscProb::PMNS_Deco p_deco = GetDeco();
  if(!TestMethodsModel(&p_deco, "Deco")) return 1;

  OscProb::PMNS_Decay p_decay = GetDecay();
  if(!TestMethodsModel(&p_decay, "Decay")) return 1;

  OscProb::PMNS_LIV p_liv = GetLIV();
  if(!TestMethodsModel(&p_liv, "LIV")) return 1;

  return 0;

}
