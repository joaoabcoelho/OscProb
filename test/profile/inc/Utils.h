
#include "PMNS_LIV.h"
#include "PMNS_SNSI.h"
#include "PMNS_Deco.h"
#include "PMNS_Iter.h"
#include "PMNS_Decay.h"
#include "PMNS_Sterile.h"
#include "PMNS_TaylorExp.h"

//.............................................................................
OscProb::PMNS_Fast* GetFast(){

  return new OscProb::PMNS_Fast();

}

//.............................................................................
OscProb::PMNS_Iter* GetIter(){

  return new OscProb::PMNS_Iter();

}

//.............................................................................
OscProb::PMNS_Deco* GetDeco(){

  OscProb::PMNS_Deco* p = new OscProb::PMNS_Deco();
  p->SetGamma(2, 1e-23);
  p->SetGamma(3, 1e-22);

  return p;

}

//.............................................................................
OscProb::PMNS_Sterile* GetSterile(){

  OscProb::PMNS_Sterile* p = new OscProb::PMNS_Sterile(4);
  p->SetDm(4, 0.1);
  p->SetAngle(1,4, 0.1);
  p->SetAngle(2,4, 0.1);
  p->SetAngle(3,4, 0.1);

  return p;

}

//.............................................................................
OscProb::PMNS_Decay* GetDecay(){

  OscProb::PMNS_Decay* p = new OscProb::PMNS_Decay();
  p->SetAlpha3(1e-4);

  return p;

}

//.............................................................................
OscProb::PMNS_NSI* GetNSI(){

  OscProb::PMNS_NSI* p = new OscProb::PMNS_NSI();
  p->SetEps(0,0, 0.1, 0);
  p->SetEps(0,1, 0.2, 0);
  p->SetEps(0,2, 0.3, 0);
  p->SetEps(1,1, 0.4, 0);
  p->SetEps(1,2, 0.5, 0);
  p->SetEps(2,2, 0.6, 0);

  return p;

}

//.............................................................................
OscProb::PMNS_SNSI* GetSNSI(){

  OscProb::PMNS_SNSI* p = new OscProb::PMNS_SNSI();
  p->SetEps(0,0, 0.1, 0);
  p->SetEps(0,1, 0.2, 0);
  p->SetEps(0,2, 0.3, 0);
  p->SetEps(1,1, 0.4, 0);
  p->SetEps(1,2, 0.5, 0);
  p->SetEps(2,2, 0.6, 0);

  return p;

}

//.............................................................................
OscProb::PMNS_LIV* GetLIV(){

  OscProb::PMNS_LIV* p = new OscProb::PMNS_LIV();
  p->SetaT(0,0, 0.1e-22, 0);
  p->SetaT(0,1, 0.2e-22, 0);
  p->SetaT(0,2, 0.3e-22, 0);
  p->SetaT(1,1, 0.4e-22, 0);
  p->SetaT(1,2, 0.5e-22, 0);
  p->SetaT(2,2, 0.6e-22, 0);
  p->SetcT(0,0, 0.1e-22, 0);
  p->SetcT(0,1, 0.2e-22, 0);
  p->SetcT(0,2, 0.3e-22, 0);
  p->SetcT(1,1, 0.4e-22, 0);
  p->SetcT(1,2, 0.5e-22, 0);
  p->SetcT(2,2, 0.6e-22, 0);

  return p;

}

//.............................................................................
OscProb::PMNS_TaylorExp* GetTaylorExp(){

  return new OscProb::PMNS_TaylorExp();

}

//.............................................................................
OscProb::PMNS_Base* GetModel(std::string model){

  if(model == "Iter")    return GetIter();
  if(model == "Deco")    return GetDeco();
  if(model == "Sterile") return GetSterile();
  if(model == "Decay")   return GetDecay();
  if(model == "NSI")     return GetNSI();
  if(model == "LIV")     return GetLIV();
  if(model == "SNSI")    return GetSNSI();
  if(model == "TaylorExp")    return GetTaylorExp();

  return GetFast();

}



