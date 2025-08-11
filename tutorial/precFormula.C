#include "PMNS_Fast.h"
#include "PMNS_TaylorExp.h"
#include "PremModel.h"

// Some functions to make nice plots
#include "SetNiceStyle.C"


void precFormula(){

  // Set nice overall style
  SetNiceStyle();

  // Get a PMNS object
  OscProb::PMNS_Fast p;
  OscProb::PMNS_TaylorExp taylor;

  // Set the baseline through the earth
  //double L = 2*6368 + 18;
  //p.SetLength(L);
  //taylor.SetLength(L);

  // PREM Model
  OscProb::PremModel prem;    

  double cosT = -0.6;
  double L = 6368 * abs(cosT); //L max = 6368 et non 2*6368

  // Fill path for cosT
  prem.FillPath(cosT);

  // Give path to calculator
  p.SetPath(prem.GetNuPath());
  taylor.SetPath(prem.GetNuPath());


  

  double dE = 0.1 * energy;
  double prec = 1e-1;
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

}

