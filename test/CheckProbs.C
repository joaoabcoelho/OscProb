#include "Utils.h"

int CheckProbs(){

  OscProb::PMNS_Fast p_fast = GetFast();
  if(CheckProb(&p_fast, "PMNS_Fast_test_values.root")) return 1;

  OscProb::PMNS_Iter p_iter = GetIter();
  if(CheckProb(&p_iter, "PMNS_Iter_test_values.root")) return 1;

  OscProb::PMNS_Sterile p_sterile = GetSterile();
  if(CheckProb(&p_sterile, "PMNS_Sterile_test_values.root")) return 1;

  OscProb::PMNS_NSI p_nsi = GetNSI();
  if(CheckProb(&p_nsi, "PMNS_NSI_test_values.root")) return 1;

  OscProb::PMNS_Deco p_deco = GetDeco();
  if(CheckProb(&p_deco, "PMNS_Deco_test_values.root")) return 1;

  OscProb::PMNS_Decay p_decay = GetDecay();
  if(CheckProb(&p_decay, "PMNS_Decay_test_values.root")) return 1;

  OscProb::PMNS_LIV p_liv = GetLIV();
  if(CheckProb(&p_liv, "PMNS_LIV_test_values.root")) return 1;

  return 0;
  
}
