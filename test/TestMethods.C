#include "Utils.h"

//.............................................................................
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
           << energy << ", " << dE << ") "
           << " with IsNuBar = " << p->GetIsNuBar() << endl
           << "Difference is "
           << "abs(" << avg_prob << " - " << avg_test << ")"
           << " = " << abs(avg_prob - avg_test)
           << " > " << 5*prec << " threshold"
           << endl;
      return false;
    }

  }}

  return true;

}

//.............................................................................
bool TestParallel(OscProb::PMNS_Base* p, double energy){

  p->SetEnergy(energy);

  OscProb::matrixD pmatrix = p->ProbMatrix(3,3);

  for(int flvi=0; flvi<3; flvi++){

    OscProb::vectorD pvector = p->ProbVector(flvi);

    for(int flvf=0; flvf<3; flvf++){

      double prob = p->Prob(flvi, flvf);

      if(abs(prob - pmatrix[flvi][flvf]) > 1e-12){
        cerr << "Found mismatch in ProbMatrix("
             << flvi << ", " << flvf << ") "
             << " with IsNuBar = " << p->GetIsNuBar()
             << " and E = " << p->GetEnergy() << endl
             << "Difference is "
             << "abs(" << prob << " - " << pmatrix[flvi][flvf] << ")"
             << " = " << abs(prob - pmatrix[flvi][flvf])
             << endl;
        return false;
      }

      if(abs(prob - pvector[flvf]) > 1e-12){
        cerr << "Found mismatch in ProbVector("
             << flvi << ", " << flvf << "): "
             << " with IsNuBar = " << p->GetIsNuBar()
             << " and E = " << p->GetEnergy() << endl
             << "Difference is "
             << "abs(" << prob << " - " << pvector[flvf] << ")"
             << " = " << abs(prob - pvector[flvf])
             << endl;
        return false;
      }

    }

  }

  return true;

}

//.............................................................................
bool TestNominal(string model, double &max_diff){

  max_diff = 0;

  double prec = 1e-12;

  if(model == "Fast") return true;

  OscProb::PMNS_Fast p0;
  OscProb::PMNS_Base* p = GetModel(model, true);

  if(model == "Iter"){
    prec = 1e-3;
    static_cast<OscProb::PMNS_Iter*> (p)->SetPrec(prec);
  }

  SetTestPath(&p0);
  SetTestPath(p);

  int nbins = 100;
  vector<double> xbins = GetLogAxis(nbins, 0.1, 100);
  for(int isnb=0; isnb<2; isnb++){

    p->SetIsNuBar(isnb);
    p0.SetIsNuBar(isnb);

    for(int i=0; i<=nbins; i++){

      double energy = xbins[i];

      for(int flvi=0; flvi<3; flvi++){
      for(int flvf=0; flvf<3; flvf++){

        double nom_prob = p0.Prob(flvi, flvf, energy);
        double alt_prob = p->Prob(flvi, flvf, energy);

        if(abs(nom_prob - alt_prob) > max_diff){
          max_diff = abs(nom_prob - alt_prob);
        }

        if(abs(nom_prob - alt_prob) > prec){

          cerr << "Found mismatch in nominal Prob("
               << flvi << ", " << flvf << ", " << energy << ") "
               << " with IsNuBar = " << p->GetIsNuBar() << endl
               << "Difference is "
               << "abs(" << nom_prob << " - " << alt_prob << ")"
               << " = " << abs(nom_prob - alt_prob)
               << " > " << prec << " threshold"
               << endl;

          cerr << Color::FAILED << " Model PMNS_" << model
               << " comparison to nominal" << endl;

          return false;

        }

      }}

    }

  }

  delete p;

  return true;

}

//.............................................................................
bool TestMethodsModel(OscProb::PMNS_Base* p, string model){

  double nom_max_diff;
  if(!TestNominal(model, nom_max_diff)) return false;

  SetTestPath(p);

  int nbins = 30;
  vector<double> xbins = GetLogAxis(nbins, 0.1, 100);
  double avg_max_diff = 0;
  for(int isnb=0; isnb<2; isnb++){

    p->SetIsNuBar(isnb);

    for(int i=0; i<=nbins; i++){

      double energy = xbins[i];

      if(!TestParallel(p, energy)){
        cerr << Color::FAILED << " Model PMNS_" << model
             << " parallel methods" << endl;
        return false;
      }

      if(!TestAvgProb(p, energy, avg_max_diff)){
        cerr << Color::FAILED << " Model PMNS_" << model
             << " AvgProb" << endl;
        return false;
      }

    }

  }

  cout << Color::PASSED << " PMNS_" << model
       << " with max nominal error: "
       << nom_max_diff
       << " and max AvgProb error: "
       << avg_max_diff << endl;

  return true;

}

//.............................................................................
int TestMethods(){

  vector<string> models = GetListOfModels();

  bool pass_all = true;

  cout << endl << string(28,'=') << " Test All Prob Methods " << string(29,'=') << endl;

  for(int i=0; i<models.size(); i++){

    OscProb::PMNS_Base* p = GetModel(models[i]);
    pass_all *= TestMethodsModel(p, models[i]);

    delete p;

  }
  cout << string(80,'=') << endl;

  cout << endl << string(28,'=') << " Test Benchmark Values " << string(29,'=') << endl;

  for(int i=0; i<models.size(); i++){

    OscProb::PMNS_Base* p = GetModel(models[i]);
    pass_all *= !CheckProb(p, "PMNS_"+models[i]+"_test_values.root");

    delete p;

  }
  cout << string(80,'=') << endl;

  return !pass_all;

}
