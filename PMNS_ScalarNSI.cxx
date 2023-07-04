////////////////////////////////////////////////////////////////////////
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework. 
//
// coelho@lal.in2p3.fr
////////////////////////////////////////////////////////////////////////

#include "PMNS_ScalarNSI.h"
#include "MatrixDecomp/zheevh3.h"

#include <iostream>
#include <cassert>
#include <stdlib.h>

using namespace OscProb;

using namespace std;
//......................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_ScalarNSI::PMNS_ScalarNSI() : PMNS_Fast(), fEta()
{
  SetStdPath();
  SetScalarNSI(0.,0.,0.,0.,0.,0.,0.,0.,0.);
  SetM(0);
}

//......................................................................
///
/// Nothing to clean.
///
PMNS_ScalarNSI::~PMNS_ScalarNSI(){}

//......................................................................
///
/// Set all NSIs parameters at once.
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param eta_ee      - The real parameter eta_ee
/// 
/// @param eta_mumu    - The real parameter eta_mumu
/// 
/// @param eta_tautau  - The real parameter eta_tautau
/// 
///non-diagonal elements are zero
void PMNS_ScalarNSI::SetScalarNSI(double eta_ee,      double eta_emu,      double eta_etau,
                double eta_mumu,    double eta_mutau,    double eta_tautau,
                double delta_emu, double delta_etau, double delta_mutau
                       ) 
{

  SetEta(0,0, eta_ee,     0);
  SetEta(1,1, eta_mumu,   0);
  SetEta(2,2, eta_tautau, 0);

  SetEta(0,1, eta_emu,  delta_emu );
  SetEta(0,2,eta_etau ,  delta_etau);
  SetEta(1,2, eta_mutau, delta_mutau);

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
void PMNS_ScalarNSI::SetEta(int flvi, int flvj, double val, double phase){

  if(flvi > flvj){
    cout << "First argument should be smaller or equal to second argument" << endl;
    cout << "Setting reverse order (Eta_" << flvj << flvi << "). " << endl;
    int temp = flvi;
    flvi = flvj;
    flvj = temp;
  }
  if(flvi<0 || flvi>2 || flvj < flvi || flvj > 2){
    cout << "Eta_" << flvi << flvj << " not valid for " << fNumNus;
    cout << " neutrinos. Doing nothing." << endl;
    return;
  }

  complexD h = val;  

  if(flvi != flvj) h *= complexD(cos(phase), sin(phase)); 

  bool isSame = (fEta[flvi][flvj] == h);
  
  if(!isSame) ClearCache();

  fGotES *= isSame;
  
  fEta[flvi][flvj] = h;
  
}

//......................................................................
///
/// Get any given NSI parameter.
///
/// Flavours are:\n
/// - 0 = nue, 1 = numu, 2 = nutau
///
/// Requires that flvi < flvj. Will notify you if input is wrong.
/// If flvi > flvj, will assume reverse order and swap flvi and flvj.
///
/// @param flvi  - The first flavour index
/// @param flvj  - The second flavour index
///
complexD PMNS_ScalarNSI::GetEta(int flvi, int flvj){

  if(flvi > flvj){
    cout << "First argument should be smaller or equal to second argument" << endl;
    cout << "Setting reverse order (Eta_" << flvj << flvi << "). " << endl;
    int temp = flvi;
    flvi = flvj;
    flvj = temp;
  }
  if(flvi<0 || flvi>2 || flvj < flvi || flvj > 2){
    cout << "Eta_" << flvi << flvj << " not valid for " << fNumNus;
    cout << " neutrinos. Returning 0." << endl;
    return zero;
  }

  return fEta[flvi][flvj];

}
/// Set lightest neutrino mass
void PMNS_ScalarNSI::SetM(double m){

fM=m;

}
//......................................................................
///
/// Set eta_ee parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a - The value of the real parameter eta_ee
///
void PMNS_ScalarNSI::SetEta_ee(double a){ SetEta(0,0, a, 0); }

//......................................................................
///
/// Set eta_mumu parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a - The value of the real parameter eta_mumu
///
void PMNS_ScalarNSI::SetEta_mumu(double a){ SetEta(1,1, a, 0); }

//......................................................................
///
/// Set eta_tautau parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a - The value of the real parameter eta_tautau
///
void PMNS_ScalarNSI::SetEta_tautau(double a){ SetEta(2,2, a, 0); }

//......................................................................

void PMNS_ScalarNSI::SetEta_emu(double a, double phi){ SetEta(0,1, a, phi); }

//......................................................................
///
/// Set eta_mumu parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a - The value of the real parameter eta_mumu
///
void PMNS_ScalarNSI::SetEta_mutau(double a, double phi){ SetEta(1,2, a, phi); }

//......................................................................
///
/// Set eta_tautau parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a - The value of the real parameter eta_tautau
///
void PMNS_ScalarNSI::SetEta_etau(double a, double phi){ SetEta(0,2, a, phi); }


////Get lightest neutrino mass


double PMNS_ScalarNSI::GetM(){ return fM;}
/// Set eps_mutau parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a   - The absolute value of the parameter eps_mutau
///
/// Build the full Hamiltonian in matter.
///
/// Here we divide the mass squared matrix Hms by the 2E
/// to obtain the vacuum Hamiltonian in eV. Then, the matter
/// potential is added to the electron component.
///
void PMNS_ScalarNSI::UpdateHam()
{

  double lv = 2 * kGeV2eV*fEnergy;     // 2E in eV 

  double kr2GNe = kK2*M_SQRT2*kGf;
  kr2GNe *= fPath.density * fPath.zoa; // Matter potential in eV

  // PMNS matrix U
  
  complexD U[3][3];
  int i,j,k;
  for( i=0;i<3;i++)
  {
    
     for( j=0;j<3;j++)
     {
        U[i][j]=fEvec[i][j];
     }
  }

  



complexD M[3][3], eta[3][3];

/// Defining mass matrix M=Diag{m1,m2,m3}, and scalar NSI matrix for NO
  if (fDm[2]>0)
  {
  
  double m1=fM;
  double m2=sqrt(fDm[1]+m1*m1);
  double m3=sqrt(fDm[2]+m1*m1);
  
   
  
  M[0][0]=m1;
  M[1][1]=m2;
  M[2][2]=m3;
  
  
  
    for( i=0;i<3;i++)
    {
      for( j=i+1;j<3;j++)
      {
      M[i][j]=0;
      }
    }



eta[0][0]=fEta[0][0]*sqrt(fDm[2]);
eta[1][1]=fEta[1][1]*sqrt(fDm[2]);
eta[2][2]=fEta[2][2]*sqrt(fDm[2]);

 eta[0][1]=fEta[0][1]*sqrt(fDm[2]);
eta[0][2]=fEta[0][2]*sqrt(fDm[2]);
eta[1][2]=fEta[1][2]*sqrt(fDm[2]);


eta[1][0]=conj(eta[0][1]);
eta[2][0]=conj(eta[0][2]);
eta[2][1]=conj(eta[1][2]);
 } 
/// Defining mass matrix M=Diag{m1,m2,m3}, and scalar NSI matrix Eta for IO
if(fDm[2]<0)
{
  double m3=fM;
  double m1=sqrt(m3*m3-fDm[2]);
  double m2=sqrt(fDm[1]+m1*m1);
  
   
  
  M[0][0]=m1;
  M[1][1]=m2;
  M[2][2]=m3;
  
  
  
    for( i=0;i<3;i++)
    {
      for( j=i+1;j<3;j++)
      {
      M[i][j]=0;
      }
    }


eta[0][0]=fEta[0][0]*sqrt(-fDm[2]);
eta[1][1]=fEta[1][1]*sqrt(-fDm[2]);
eta[2][2]=fEta[2][2]*sqrt(-fDm[2]);

 eta[0][1]=fEta[0][1]*sqrt(-fDm[2]);
eta[0][2]=fEta[0][2]*sqrt(-fDm[2]);
eta[1][2]=fEta[1][2]*sqrt(-fDm[2]);


eta[1][0]=conj(eta[0][1]);
eta[2][0]=conj(eta[0][2]);
eta[2][1]=conj(eta[1][2]);


}
 /// Calculating U^{\dagger}
  complexD Udag[3][3];
  
  for(int i=0;i<3;i++)
  {
    
    for(int j=0;j<3;j++)
    {
      Udag[i][j]=conj(U[j][i]);
    }
  }
  
/// Calculating MU^\dagger
   complexD MUdag[3][3];
  
  for(i=0;i<3;i++)    
  {    
    for(j=0;j<3;j++)    
    {    
    MUdag[i][j]=0;    
      for(k=0;k<3;k++)    
      {    
      MUdag[i][j]+=M[i][k]*Udag[k][j];    
      }    
    }    
  } 
/// Calculating UMU^\dagger
 complexD UMUdag[3][3];

 for(i=0;i<3;i++)    
 {    
   for(j=0;j<3;j++)    
   {    
   UMUdag[i][j]=0;    
     for(k=0;k<3;k++)    
     {    
     UMUdag[i][j]+=U[i][k]*MUdag[k][j];    
     }    
   }    
 } 

 
/// Calculating UmUdagger+Eta
  
   complexD MpdM[3][3];
  
  for( i=0;i<3;i++)
  {
    for( j=0;j<3;j++)
    {
      MpdM[i][j]=UMUdag[i][j]+eta[i][j];
    }
  }
/// Calculating (UmUdagger+Eta)^\dagger 
  complexD MpdMdag[3][3];
  
  
  for(int i=0;i<3;i++)
  {
    for(int j=0;j<3;j++)
    {
      MpdMdag[i][j]=conj(MpdM[j][i]);
    }
  }
/// Building the vacuum Hamiltonian H=1/(2E)*((UmUdagger+Eta)(UmUdagger+Eta)^\dagger)  
  for(i=0;i<3;i++)    
  {    
    for(j=0;j<3;j++)    
    {    
    fHam[i][j]=0;    
      for(k=0;k<3;k++)    
      {    
      fHam[i][j]+=MpdM[i][k]*MpdMdag[k][j]/lv;    
      }    
    }    
  } 
  
  /// Adding matter potential
  if(!fIsNuBar) fHam[0][0] += kr2GNe;
  else          fHam[0][0] -= kr2GNe;

}


////////////////////////////////////////////////////////////////////////
