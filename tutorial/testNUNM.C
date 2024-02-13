
#include "PremModel.h"
#include "PMNS_Fast.h"
#include "PMNS_NSI.h"
#include "PMNS_Sterile.h"
#include "PMNS_NUNM.h"

void testNUNM(){

  // Initialize your objects
  OscProb::PMNS_Fast myPMNS;
  
  // Other examples
  OscProb::PMNS_NSI  myNSI; // NSI extension
  OscProb::PMNS_Sterile mySterile(4); // 3+2 Model
  OscProb::PMNS_NUNM  myNUNM;
  
  OscProb::PremModel prem;
  prem.FillPath(-1);
  
  myPMNS.SetEnergy(25.0);
  myNUNM.SetEnergy(25.0);
  mySterile.SetEnergy(25.0);
  
  int fNumNus = 3;
  // Set oscillation paramaters
  // By default, they are set to
  // PDG values with SetStdPars()
  // in the constructor
//  myPMNS.SetDm(3, -2.4e-3);   // Set Dm31 in eV^2
  //myPMNS.SetAngle(1,2, 45./180.*3.14159); // Set Theta13 in radians
  //myNUNM.SetAngle(1,2, 45./180.*3.14159);
  //mySterile.SetAngle(1,2, 45./180.*3.14159);
  //  myPMNS.SetAngle(1,3, 0);
//  myPMNS.SetAngle(2,3, 0);
  //myPMNS.SetDelta(1,3, 1.57); // Set Delta13 in radians

  // Set other quantities
  //myPMNS.SetDensity(2.5); // Set the matter density in g/cm^3
  myPMNS.SetPath(prem.GetNuPath());
  myNUNM.SetPath(prem.GetNuPath());
  mySterile.SetPath(prem.GetNuPath());

  // Say whether you want antineutrinos
  // Default is false, i.e. neutrinos
  //myPMNS.SetIsNuBar(true);

  // Calculate probability of
  cerr << " ----------------------------- PMNS ----------------------------- " << endl;
  auto start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < fNumNus; i++) {
    for (int j = 0; j < fNumNus; j++) {
        cerr << "P( "<< i << "-->" << j <<" ) = " << myPMNS.Prob(i,j) << endl;
    }
    cerr << "sum P( x --> "<< i << " ) = " << myPMNS.Prob(i,0)+ myPMNS.Prob(i,1)+ myPMNS.Prob(i,2) << endl;
  }  
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> elapsed = end - start;
  cerr << "elapsed time: " << elapsed.count() << " ms\n"<< endl;  

  cerr << " ----------------------------- Sterile 0 ----------------------------- " << endl;
  start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < fNumNus; i++) {
    for (int j = 0; j < fNumNus; j++) {
        cerr << "P( "<< i << "-->" << j <<" ) = " << mySterile.Prob(i,j) << endl;
    }
    cerr << "sum P( x --> "<< i << " ) = " << mySterile.Prob(i,0)+ mySterile.Prob(i,1)+ mySterile.Prob(i,2) << endl;
  }
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cerr << "elapsed time: " << elapsed.count() << " ms\n"<< endl;

  cerr << " ----------------------------- Sterile4 10°  ----------------------------- " << endl;
  mySterile.SetAngle(1, 4, 10/180.*3.14159);  
  mySterile.SetDm(4, 100.0);
  start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < fNumNus; i++) {
    for (int j = 0; j < fNumNus; j++) {
        cerr << "P( "<< i << "-->" << j <<" ) = " << mySterile.Prob(i,j) << endl;
    }
    cerr << "sum P( x --> "<< i << " ) = " << mySterile.Prob(i,0)+ mySterile.Prob(i,1)+ mySterile.Prob(i,2) << endl;
  }
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cerr << "elapsed time: " << elapsed.count() << " ms\n"<< endl;
  
  cerr << " ----------------------------- NUNM alpha = 0 fracVnc = 0 ----------------------------- " << endl;
  myNUNM.SetFracVnc(0.);
  start = std::chrono::high_resolution_clock::now();
                  for (int i = 0; i < fNumNus; i++) {
                          for (int j = 0; j < fNumNus; j++) {
                                  cerr << "P( "<< i << "-->" << j <<" ) = " << myNUNM.Prob(i,j) << endl;
                          }
                  }
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cerr << "elapsed time: " << elapsed.count() << " ms\n"<< endl;
  
  cerr << " ----------------------------- NUNM alpha = 0 fracVnc = 1 ----------------------------- " << endl;
  myNUNM.SetFracVnc(1.);
  start = std::chrono::high_resolution_clock::now();
                  for (int i = 0; i < fNumNus; i++) {
                          for (int j = 0; j < fNumNus; j++) {
                                  cerr << "P( "<< i << "-->" << j <<" ) = " << myNUNM.Prob(i,j) << endl;
                          }
			  cerr << "sum P( x --> "<< i << " ) = " << myNUNM.Prob(i,0)+ myNUNM.Prob(i,1)+ myNUNM.Prob(i,2) << endl;
                  }
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cerr << "elapsed time: " << elapsed.count() << " ms\n"<< endl;
  
  /*
  cerr << " ----------------------------- NUNM one alpha =! 0 fracVnc = 1 ----------------------------- " << endl;
  myNUNM.SetFracVnc(1.);
  for (int i = 0; i < fNumNus; i++) {
	  for (int j = 0; j < i+1; j++) {
		  myNUNM.SetAlpha(i, j, 0.1, 0 );
		  start = std::chrono::high_resolution_clock::now();
		  for (int k = 0; k < fNumNus; k++) {
			  for (int l = 0; l < fNumNus; l++) {
				  cerr << "P( "<< k << "-->" << l <<" ) = " << myNUNM.Prob(k,l) << endl;
			  }
			  cerr << "sum_x P( "<< k << " --> x ) = " << myNUNM.Prob(k,0)+ myNUNM.Prob(k,1)+ myNUNM.Prob(k,2) << endl;
		  }
		  end = std::chrono::high_resolution_clock::now();
		  elapsed = end - start;
		  cerr << "elapsed time: " << elapsed.count() << " ms\n"<< endl;
		  myNUNM.SetAlpha(i, j, 0, 0 );
	  }		
  }
  */

  /*  
  cerr << " ----------------------------- NUNM one alpha =! 0 fracVnc = 0 ----------------------------- " << endl;
  myNUNM.SetFracVnc(0.);
  for (int i = 0; i < fNumNus; i++) {
          for (int j = 0; j < i+1; j++) {
                  myNUNM.SetAlpha(i, j, 0.1, 0 );
                  start = std::chrono::high_resolution_clock::now();
		  for (int k = 0; k < fNumNus; k++) {
			  for (int l = 0; l < fNumNus; l++) {
                                  cerr << "P( "<< k << "-->" << l <<" ) = " << myNUNM.Prob(k,l) << endl;
                          }
			  cerr << "sum_x P( "<< k << " --> x ) = " << myNUNM.Prob(k,0)+ myNUNM.Prob(k,1)+ myNUNM.Prob(k,2) << endl;
                  }
                  end = std::chrono::high_resolution_clock::now();
                  elapsed = end - start;
                  cerr << "elapsed time: " << elapsed.count() << " ms\n"<< endl;
		  myNUNM.SetAlpha(i, j, 0, 0 );
          }
  }
  */

  //myNUNM.SetAlpha(0, 0, 1, 0 );
  //myNUNM.SetAlpha(1, 1, 1, 0 );
  //myNUNM.SetAlpha(2, 2, 0.173648, 0 );
  //myNUNM.SetAlpha(1, 0, 1, 0 );
  //myNUNM.SetAlpha(2, 0, 1, 0 );
  myNUNM.SetAlpha(2, 1, 0.173648, 0 );
  
  cerr << " ----------------------------- NUNM alpha32 = 0.17 fracVnc = 0 ----------------------------- " << endl;
  myNUNM.SetFracVnc(0.);
  //myNUNM.Prob(0,0);// commented
  start = std::chrono::high_resolution_clock::now();  
                for (int i = 0; i < fNumNus; i++) {
                          for (int j = 0; j < fNumNus; j++) {
                                  cerr << "P( "<< i << "-->" << j <<" ) = " << myNUNM.Prob(i,j) << endl;
                          }
			  cerr << "sum_x P( " << i << " --> x ) = " << myNUNM.Prob(i,0)+ myNUNM.Prob(i,1)+ myNUNM.Prob(i,2) << endl;
                }
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cerr << "elapsed time: " << elapsed.count() << " ms\n"<< endl;
  
  cerr << " ----------------------------- NUNM alpha33 = 0.17 fracVnc = 1 ----------------------------- " << endl;
  myNUNM.SetFracVnc(1.);
  //myNUNM.Prob(0,0);// commented
  start = std::chrono::high_resolution_clock::now();
                for (int i = 0; i < fNumNus; i++) {
                          for (int j = 0; j < fNumNus; j++) {
                                  cerr << "P( "<< i << "-->" << j <<" ) = " << myNUNM.Prob(i,j) << endl;
                          }
                          cerr << "sum_x P( " << i << " --> x ) = " << myNUNM.Prob(i,0)+ myNUNM.Prob(i,1)+ myNUNM.Prob(i,2) << endl;
                }
  end = std::chrono::high_resolution_clock::now();
  elapsed = end - start;
  cerr << "elapsed time: " << elapsed.count() << " ms\n"<< endl;

  // Set 3 paths in order
  //myPMNS.SetPath(1000, 2.50); // Set initial path with L = 1000 km and rho = 2.50 g/cm^3
  //myPMNS.AddPath(2000, 4.25); // Add a second path with L = 2000 km and rho = 4.25 g/cm^3
  //myPMNS.AddPath(1500, 6.00); // Add a third path with L = 1500 km and rho = 6.00 g/cm^3

  // Get numubar to nuebar probability
  //Pme = myPMNS.Prob(1,0);



  // Use a PremModel to make paths
  // through the earth
  //OscProb::PremModel prem;

  // Chose an angle for the neutrino
  // and fill the paths with cosTheta
  // e.g. cosTheta = -1 (vertical up-going)
  //prem.FillPath(-1);

  // Give the path to the PMNS object
  // and get the probability
  //myPMNS.SetPath(prem.GetNuPath());
  //Pme = myPMNS.Prob(1,0);


}
