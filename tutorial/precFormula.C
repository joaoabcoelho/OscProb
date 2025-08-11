#include "PMNS_Fast.h"
#include "PMNS_TaylorExp.h"
#include "PremModel.h"

// Some functions to make nice plots
#include "SetNiceStyle.C"


bool precFormula(){

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

  //double cosT = -0.6;
  //double L = 6368 * abs(cosT); //L max = 6368 et non 2*6368

  // Fill path for cosT
  //prem.FillPath(cosT);

  // Give path to calculator
  //p.SetPath(prem.GetNuPath());
  //taylor.SetPath(prem.GetNuPath());




  //double dE = 0.1 * energy;
  double prec = 1e-4;
  p.SetAvgProbPrec(prec);
  taylor.SetAvgProbPrec(prec);

  //for (double b = 0 ; b<1 ; b+=0.2){

    for(int flvi=0; flvi<3; flvi++){
      for(int flvf=0; flvf<3; flvf++){

        for(double cosZ=-0.95; cosZ<1; cosZ+=0.1){

          prem.FillPath(cosZ);
          p.SetPath(prem.GetNuPath());
          taylor.SetPath(prem.GetNuPath());

          for(double E=1; E<100; E++){

            double dE = E * 0.1;

            int npts = 1000;
            double avg_test = 0;
            for(int i=0; i<npts; i++){
              double energy_i = E + dE*((i+0.5)/npts - 0.5);
              avg_test += p.Prob(flvi, flvf, energy_i) / npts;
            }

            double avgTaylor = taylor.AvgProb(flvi,flvf, E, dE,0.8, 0.1 ,0.5,10);

            if(abs(avgTaylor - avg_test) > 5*prec){
              cerr << "Found mismatch in AvgProb("
                   << flvi << ", " << flvf << ", "
                   << E << ", " << dE << ") "
                   << " with IsNuBar = " << taylor.GetIsNuBar() << endl
                   << "Difference is "
                   << "abs(" << avgTaylor << " - " << avg_test << ")"
                   << " = " << abs(avgTaylor - avg_test)
                   << " > " << 5*prec << " threshold"
                   //<< endl << "Value b = " << b
                   << endl;

              return false;
            }

          }
        }
    
      }
    }

  //}

  return true;

}

