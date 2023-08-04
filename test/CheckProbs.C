#include "Utils.h"

int CheckProbs(){

  vector<string> models = GetListOfModels();

  for(int i=0; i<models.size(); i++){

    OscProb::PMNS_Base* p = GetModel(models[i]);

    if(CheckProb(p, "PMNS_"+models[i]+"_test_values.root")) return 1;
  }

  return 0;

}
