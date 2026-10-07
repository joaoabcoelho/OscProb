
#include "Utils.h"

void MakeTestSamples(vector<string> models = {}){

  if(!models.size()) models = GetListOfModels();

  for(int i=0; i<models.size(); i++){

    TString filename = "PMNS_"+models[i]+"_test_values.root";
    if(!gSystem->AccessPathName("data/"+filename)) continue;

    OscProb::PMNS_Base* p = GetModel(models[i]);
    SaveTestFile(p, filename);
    delete p;

  }

}
