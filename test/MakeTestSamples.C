
#include "Utils.h"

void MakeTestSamples(){

  OscProb::PMNS_Fast p_fast = GetFast();
  SaveTestFile(&p_fast, "PMNS_Fast_test_values.root");

  OscProb::PMNS_Iter p_iter = GetIter();
  SaveTestFile(&p_iter, "PMNS_Iter_test_values.root");

  OscProb::PMNS_Sterile p_sterile = GetSterile();
  SaveTestFile(&p_sterile, "PMNS_Sterile_test_values.root");

  OscProb::PMNS_NSI p_nsi = GetNSI();
  SaveTestFile(&p_nsi, "PMNS_NSI_test_values.root");

  OscProb::PMNS_Deco p_deco = GetDeco();
  SaveTestFile(&p_deco, "PMNS_Deco_test_values.root");

  OscProb::PMNS_Decay p_decay = GetDecay();
  SaveTestFile(&p_decay, "PMNS_Decay_test_values.root");

  OscProb::PMNS_LIV p_liv = GetLIV();
  SaveTestFile(&p_liv, "PMNS_LIV_test_values.root");

}
