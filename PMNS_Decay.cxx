////////////////////////////////////////////////////////////////////////
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework with Neutrino Decay in the 3rd state. 
//
// This class expands the PMNS_Fast class including the decay of the third mass// state of the neutrino through a decay constant alpha3=m3/tau3 (eV^2), where 
// m3 is the mass in the restframe and tau3 is the lifetime in the restframe.
//
// Diagonalization of the matrix is done using the Eigen library. 
//
// \author Joao Coelho - coelho\@lal.in2p3.fr and
// \colaborator Victor Carretero - vcarretero\@km3net.de  
////////////////////////////////////////////////////////////////////////

#include "PMNS_Decay.h"
#include <iostream>
#include <cassert>
#include <stdlib.h>
#include "complexsolver.h"
using namespace OscProb;

//......................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_Decay::PMNS_Decay() : PMNS_Base(), fHam() {
	falpha    = vector<double>(fNumNus, 0);
	fEvalI = vector<double>(fNumNus, 0); /////
	fEvecinv = vector< vector<complexD> >(fNumNus, vector<complexD>(fNumNus,zero)); /////
	 
	}

//......................................................................
///
/// Nothing to clean.
///
PMNS_Decay::~PMNS_Decay(){}

//......................................................................
///
/// Set all mixing parameters at once.
/// @param th12    - The value of the mixing angle theta_12
/// @param th23    - The value of the mixing angle theta_23
/// @param th13    - The value of the mixing angle theta_13
/// @param deltacp - The value of the CP phase delta_13
///
void PMNS_Decay::SetMix(double th12, double th23, double th13, double deltacp) 
{

  SetAngle(1,2, th12);
  SetAngle(1,3, th13);
  SetAngle(2,3, th23);
  SetDelta(1,3, deltacp);

}
//......................................................................                                                       
///     Setting Alpha3 parameter, it must be possitive, and units are eV^2 .
///     Alpha3=m3/tau3, mass and lifetime of the third state in the restframe 
void PMNS_Decay::SetAlpha3(double alpha3) 
{
if(alpha3<0){
	cout << "Alpha3 must be positive" << endl;
	return;
}

fBuiltHms *= (falpha[2] == alpha3);


 falpha[0]=0;
 falpha[1]=0;
 falpha[2] = alpha3;

}

void PMNS_Decay::SetIsNuBar(bool isNuBar)
{

  // Check if value is actually changing
  fGotES *= (fIsNuBar == isNuBar);
  fBuiltHms *= (fIsNuBar == isNuBar);
  fIsNuBar = isNuBar;
}
//......................................................................
///
/// Set both mass-splittings at once.
///
/// These are Dm_21 and Dm_32 in eV^2.\n
/// The corresponding Dm_31 is set in the class attributes.
///
/// @param dm21 - The solar mass-splitting Dm_21 
/// @param dm32 - The atmospheric mass-splitting Dm_32 
///
void PMNS_Decay::SetDeltaMsqrs(double dm21, double dm32) 
{

  SetDm(2, dm21);
  SetDm(3, dm32 + dm21);

}

//
/// Get alpha3
double PMNS_Decay::GetAlpha3() 
{
return falpha[2];
}
///Build the Hamiltonian and rotate it
void PMNS_Decay::RotateH(int i,int j){

  // Do nothing if angle is zero
  if(fTheta[i][j]==0) return;

  double fSinBuffer = sin(fTheta[i][j]);
  double fCosBuffer = cos(fTheta[i][j]);

  double  fHmsBufferD;
  complexD fHmsBufferC;

  // With Delta
  if(i+1<j){
    complexD fExpBuffer;
    if(fIsNuBar==true){
      fExpBuffer = complexD(cos(-fDelta[i][j]), -sin(-fDelta[i][j]));
    }
    else{
      fExpBuffer = complexD(cos(fDelta[i][j]), -sin(fDelta[i][j]));
    }

    // General case
    if(i>0){
      // Top columns
      for(int k=0; k<i; k++){
        fHmsBufferC = fHms[k][i];

        fHms[k][i] *= fCosBuffer;
        fHms[k][i] += fHms[k][j] * fSinBuffer * conj(fExpBuffer);

        fHms[k][j] *= fCosBuffer;
        fHms[k][j] -= fHmsBufferC * fSinBuffer * fExpBuffer;
      }

      // Middle row and column
      for(int k=i+1; k<j; k++){
        fHmsBufferC = fHms[k][j];

        fHms[k][j] *= fCosBuffer;
        fHms[k][j] -= conj(fHms[i][k]) * fSinBuffer * fExpBuffer;

        fHms[i][k] *= fCosBuffer;
        fHms[i][k] += fSinBuffer * fExpBuffer * conj(fHmsBufferC);
      }

      // Nodes ij
      fHmsBufferC = fHms[i][i];
      fHmsBufferD = real(fHms[j][j]);

      fHms[i][i] *= fCosBuffer * fCosBuffer;
      fHms[i][i] += 2 * fSinBuffer * fCosBuffer * real(fHms[i][j] * conj(fExpBuffer));
      fHms[i][i] += fSinBuffer * fHms[j][j] * fSinBuffer;

      fHms[j][j] *= fCosBuffer * fCosBuffer;
      fHms[j][j] += fSinBuffer * fHmsBufferC * fSinBuffer;
      fHms[j][j] -= 2 * fSinBuffer * fCosBuffer * real(fHms[i][j] * conj(fExpBuffer));

      fHms[i][j] -= 2 * fSinBuffer * real(fHms[i][j] * conj(fExpBuffer)) * fSinBuffer * fExpBuffer;
      fHms[i][j] += fSinBuffer * fCosBuffer * (fHmsBufferD - fHmsBufferC) * fExpBuffer;

    }
    // First rotation on j (No top columns)
    else{
      // Middle rows and columns
      for(int k=i+1; k<j; k++){
        fHms[k][j] = -conj(fHms[i][k]) * fSinBuffer * fExpBuffer;

        fHms[i][k] *= fCosBuffer;
      }

      // Nodes ij
      fHmsBufferD = real(fHms[i][i]);

      fHms[i][j] = fSinBuffer * fCosBuffer * (fHms[j][j] - fHmsBufferD) * fExpBuffer;

      fHms[i][i] *= fCosBuffer * fCosBuffer;
      fHms[i][i] += fSinBuffer * fHms[j][j] * fSinBuffer;

      fHms[j][j] *= fCosBuffer * fCosBuffer;
      fHms[j][j] += fSinBuffer * fHmsBufferD * fSinBuffer;
    }

  }
  // Without Delta (No middle rows or columns: j = i+1)
  else{
    // General case
    if(i>0){
      // Top columns
      for(int k=0; k<i; k++){
        fHmsBufferC = fHms[k][i];

        fHms[k][i] *= fCosBuffer;
        fHms[k][i] += fHms[k][j] * fSinBuffer;

        fHms[k][j] *= fCosBuffer;
        fHms[k][j] -= fHmsBufferC * fSinBuffer;
      }

      // Nodes ij
      fHmsBufferC = fHms[i][i];
      fHmsBufferD = real(fHms[j][j]);

      fHms[i][i] *= fCosBuffer * fCosBuffer;
      fHms[i][i] += 2 * fSinBuffer * fCosBuffer * real(fHms[i][j]);
      fHms[i][i] += fSinBuffer * fHms[j][j] * fSinBuffer;

      fHms[j][j] *= fCosBuffer * fCosBuffer;
      fHms[j][j] += fSinBuffer * fHmsBufferC * fSinBuffer;
      fHms[j][j] -= 2 * fSinBuffer * fCosBuffer * real(fHms[i][j]);

      fHms[i][j] -= 2 * fSinBuffer * real(fHms[i][j]) * fSinBuffer;
      fHms[i][j] += fSinBuffer * fCosBuffer * (fHmsBufferD - fHmsBufferC);

    }
    // First rotation (theta12)
    else{

      fHms[i][j] = fSinBuffer * fCosBuffer * fHms[j][j];

      fHms[i][i] = fSinBuffer * fHms[j][j] * fSinBuffer;

      fHms[j][j] *= fCosBuffer * fCosBuffer;

    }
  }

}


void PMNS_Decay::RotateHd(int i,int j){

  // Do nothing if angle is zero
  if(fTheta[i][j]==0) return;

  double fSinBuffer = sin(fTheta[i][j]);
  double fCosBuffer = cos(fTheta[i][j]);

  double  fHdBufferD;
  complexD fHdBufferC;

  // With Delta
  if(i+1<j){
    complexD fExpBuffer;
    if(fIsNuBar==true){
      fExpBuffer = complexD(cos(-fDelta[i][j]), -sin(-fDelta[i][j]));
}
    else{
      fExpBuffer = complexD(cos(fDelta[i][j]), -sin(fDelta[i][j]));
    
}

    // General case
    if(i>0){
      // Top columns
      for(int k=0; k<i; k++){
        fHdBufferC = fHd[k][i];

        fHd[k][i] *= fCosBuffer;
        fHd[k][i] += fHd[k][j] * fSinBuffer * conj(fExpBuffer);

        fHd[k][j] *= fCosBuffer;
        fHd[k][j] -= fHdBufferC * fSinBuffer * fExpBuffer;
      }

      // Middle row and column
      for(int k=i+1; k<j; k++){
        fHdBufferC = fHd[k][j];

        fHd[k][j] *= fCosBuffer;
        fHd[k][j] -= conj(fHd[i][k]) * fSinBuffer * fExpBuffer;

        fHd[i][k] *= fCosBuffer;
        fHd[i][k] += fSinBuffer * fExpBuffer * conj(fHdBufferC);
      }

      // Nodes ij
      fHdBufferC = fHd[i][i];
      fHdBufferD = real(fHd[j][j]);

      fHd[i][i] *= fCosBuffer * fCosBuffer;
      fHd[i][i] += 2 * fSinBuffer * fCosBuffer * real(fHd[i][j] * conj(fExpBuffer));
      fHd[i][i] += fSinBuffer * fHd[j][j] * fSinBuffer;

      fHd[j][j] *= fCosBuffer * fCosBuffer;
      fHd[j][j] += fSinBuffer * fHdBufferC * fSinBuffer;
      fHd[j][j] -= 2 * fSinBuffer * fCosBuffer * real(fHd[i][j] * conj(fExpBuffer));

      fHd[i][j] -= 2 * fSinBuffer * real(fHd[i][j] * conj(fExpBuffer)) * fSinBuffer * fExpBuffer;
      fHd[i][j] += fSinBuffer * fCosBuffer * (fHdBufferD - fHdBufferC) * fExpBuffer;

    }
    // First rotation on j (No top columns)
    else{
      // Middle rows and columns
      for(int k=i+1; k<j; k++){
        fHd[k][j] = -conj(fHd[i][k]) * fSinBuffer * fExpBuffer;

        fHd[i][k] *= fCosBuffer;
      }

      // Nodes ij
      fHdBufferD = real(fHd[i][i]);

      fHd[i][j] = fSinBuffer * fCosBuffer * (fHd[j][j] - fHdBufferD) * fExpBuffer;

      fHd[i][i] *= fCosBuffer * fCosBuffer;
      fHd[i][i] += fSinBuffer * fHd[j][j] * fSinBuffer;

      fHd[j][j] *= fCosBuffer * fCosBuffer;
      fHd[j][j] += fSinBuffer * fHdBufferD * fSinBuffer;
    }

  }
  // Without Delta (No middle rows or columns: j = i+1)
  else{
    // General case
    if(i>0){
      // Top columns
      for(int k=0; k<i; k++){
        fHdBufferC = fHd[k][i];

        fHd[k][i] *= fCosBuffer;
        fHd[k][i] += fHd[k][j] * fSinBuffer;

        fHd[k][j] *= fCosBuffer;
        fHd[k][j] -= fHdBufferC * fSinBuffer;
      }

      // Nodes ij
      fHdBufferC = fHd[i][i];
      fHdBufferD = real(fHd[j][j]);

      fHd[i][i] *= fCosBuffer * fCosBuffer;
      fHd[i][i] += 2 * fSinBuffer * fCosBuffer * real(fHd[i][j]);
      fHd[i][i] += fSinBuffer * fHd[j][j] * fSinBuffer;

      fHd[j][j] *= fCosBuffer * fCosBuffer;
      fHd[j][j] += fSinBuffer * fHdBufferC * fSinBuffer;
      fHd[j][j] -= 2 * fSinBuffer * fCosBuffer * real(fHd[i][j]);

      fHd[i][j] -= 2 * fSinBuffer * real(fHd[i][j]) * fSinBuffer;
      fHd[i][j] += fSinBuffer * fCosBuffer * (fHdBufferD - fHdBufferC);

    }
    // First rotation (theta12)
    else{

      fHd[i][j] = fSinBuffer * fCosBuffer * fHd[j][j];

      fHd[i][i] = fSinBuffer * fHd[j][j] * fSinBuffer;

      fHd[j][j] *= fCosBuffer * fCosBuffer;

    }
  }

}



void PMNS_Decay::BuildHam()
{

  // Check if anything changed
  if(fBuiltHms) return;
   
  // Tag to recompute eigensystem
  fGotES = false;
  
  ///Reset everything
   for(int i=0; i<fNumNus; i++){
   for(int j=0; j<fNumNus; j++){
	   fHms[i][j]= 0;
	   }
   }
  for(int j=0; j<fNumNus; j++){
    // Set mass splitting
    fHms[j][j] = fDm[j];
    // Reset off-diagonal elements
    for(int i=0; i<j; i++){
      fHms[i][j] = 0;
    }
 
    //Rotate j neutrinos
    for(int i=0; i<j; i++){
      RotateH(i,j);
    }
  }
  
   for(int i=0; i<fNumNus; i++){
   for(int j=i+1; j<fNumNus; j++){
	   fHms[j][i]= conj(fHms[i][j]);
	   }
   }
  

  ///Introduction of the alpha3
  ///Reset everything
   for(int i=0; i<fNumNus; i++){
   for(int j=0; j<fNumNus; j++){
	   fHd[i][j]= 0;
	   }
   }
  
   for(int j=0; j<fNumNus; j++){
	   
    // Set alpha
    fHd[j][j] = falpha[j];
    
    // Reset off-diagonal elements
    for(int i=0; i<j; i++){
      fHd[i][j] = 0;
    }
 
    // Rotate j neutrinos
    for(int i=0; i<j; i++){
      RotateHd(i,j);
    }
  }
  
   for(int i=0; i<fNumNus; i++){
   for(int j=i+1; j<fNumNus; j++){
	   fHd[j][i]= conj(fHd[i][j]);
	   }
   }
  

  const   complex<double> numi(0.0,1.0);  
  ///Construct the total Hms+Hd
  for(int j=0; j<fNumNus; j++){
	for(int i=0; i<fNumNus; i++){
		fHt[i][j]=fHms[i][j]-numi*fHd[i][j];
  }
}  
 
 
 
  //ClearCache();

  // Tag as built
  fBuiltHms = true;
  

 

}


//......................................................................
///
/// Build the full Hamiltonian in matter.
///
/// Here we divide the mass squared matrix Hms by the 2E
/// to obtain the vacuum Hamiltonian in eV. Then, the matter
/// potential is added to the electron component.
///


void PMNS_Decay::UpdateHam()
{

  double lv = 2 * kGeV2eV*fEnergy;     // 2E in eV 

  double kr2GNe = kK2*M_SQRT2*kGf;
  kr2GNe *= fPath.density * fPath.zoa; // Matter potential in eV

  // Finish building Hamiltonian in matter with dimension of eV
  for(int i=0;i<fNumNus;i++){
	  for(int j=0; j<fNumNus; j++){
    fHam[i][j] = fHt[i][j]/lv;
    
    
    }
  }
 

  if(!fIsNuBar) fHam[0][0] += kr2GNe;
  else          fHam[0][0] -= kr2GNe;
  
}
//......................................................................
///
/// Solve the full Hamiltonian for eigenvectors and eigenvalues.
///
/// This is using a method from the Eigen library
///
void PMNS_Decay::SolveHam()
{


  // Build Hamiltonian
  BuildHam();

  // Check if anything changed  
  if(fGotES) return;

  // Try caching if activated
  //if(TryCache()) return;

  UpdateHam();


  complexD fEvalEigen[3];
  complexD fEvecEigen[3][3];
  complexD fEvecEigeninv[3][3];
  // Solve Hamiltonian for eigensystem using the Eigen library method 
  complexsolver(fHam,fEvecEigen,fEvecEigeninv,fEvalEigen);



  // Fill fEval and fEvec vectors from Eigen arrays  
  for(int i=0;i<fNumNus;i++){
    fEval[i] = fEvalEigen[i].real();
    fEvalI[i] = fEvalEigen[i].imag();
    for(int j=0;j<fNumNus;j++){
      fEvec[i][j] = fEvecEigen[i][j];
      fEvecinv[i][j] = fEvecEigeninv[i][j];
    }
  }
  
  fGotES = true;

  // Fill cache if activated
  //FillCache();

}

bool PMNS_Decay::CheckProb(int flvi)
{
	double controler=0;
for( int i=0; i<fNumNus; i++){
	controler+=Prob(flvi,i);
}

if(controler <= 1.0005)
{
	return true;
	}
else{
	return false;}
}
//.....................................................................
///
/// Set the eigensystem to the analytic solution in vacuum.
///
/// We know the vacuum eigensystem, so just write it explicitly
///

void PMNS_Decay::PropagatePath(NuPath p)
{

  // Set the neutrino path
  SetCurPath(p);

  // Solve for eigensystem
  SolveHam();

  double LengthIneV = kKm2eV * p.length;
  for(int i=0; i<fNumNus; i++){

    double arg = fEval[i] * LengthIneV;
    double argI = fEvalI[i] * LengthIneV;
    std::complex<double> numi=std::complex<double>(0,1.0);
   // fPhases[i] = complexD(cos(arg)*exp(argI), -sin(arg)*exp(argI));
   fPhases[i] = exp(-numi*(arg+numi*argI));
   
   
}
  for(int i=0;i<fNumNus;i++){
    fBuffer[i] = 0;
    for(int j=0;j<fNumNus;j++){
		//if(falpha[2]==0){
      fBuffer[i] += fEvecinv[i][j] * fNuState[j];
        //                }
     // else{
	//	  fBuffer[i] += conj(fEvec[j][i]) * fNuState[j];
    }
    fBuffer[i] *= fPhases[i];
  }

  // Propagate neutrino state
  for(int i=0;i<fNumNus;i++){
    fNuState[i] = 0;
    for(int j=0;j<fNumNus;j++){
      fNuState[i] +=  fEvec[i][j] * fBuffer[j];
      
    }
    
	
		
		
	
  }

}


////////////////////////////////////////////////////////////////////////
