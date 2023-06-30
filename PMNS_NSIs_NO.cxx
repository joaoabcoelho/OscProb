#include "PMNS_NSIs_NO.h"
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
PMNS_NSIs_NO::PMNS_NSIs_NO() : PMNS_Fast(), fEta()
{
  SetStdPath();
  SetNSIs(0.,0.,0.,0.,0.,0.,0.,0.,0.);
  SetM1(0);
}

//......................................................................
///
/// Nothing to clean.
///
PMNS_NSIs_NO::~PMNS_NSIs_NO(){}

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
void PMNS_NSIs_NO::SetNSIs(double eta_ee,      double eta_emu,      double eta_etau,
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
void PMNS_NSIs_NO::SetEta(int flvi, int flvj, double val, double phase){

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
complexD PMNS_NSIs_NO::GetEta(int flvi, int flvj){

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

void PMNS_NSIs_NO::SetLowestMass(double m1){

double h=m1;
fM1=h;

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
void PMNS_NSIs_NO::SetEta_ee(double a){ SetEta(0,0, a, 0); }

//......................................................................
///
/// Set eta_mumu parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a - The value of the real parameter eta_mumu
///
void PMNS_NSIs_NO::SetEta_mumu(double a){ SetEta(1,1, a, 0); }

//......................................................................
///
/// Set eta_tautau parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a - The value of the real parameter eta_tautau
///
void PMNS_NSIs_NO::SetEta_tautau(double a){ SetEta(2,2, a, 0); }

//......................................................................

void PMNS_NSIs_NO::SetEta_emu(double a, double phi){ SetEta(0,1, a, phi); }

//......................................................................
///
/// Set eta_mumu parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a - The value of the real parameter eta_mumu
///
void PMNS_NSIs_NO::SetEta_mutau(double a, double phi){ SetEta(1,2, a, phi); }

//......................................................................
///
/// Set eta_tautau parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a - The value of the real parameter eta_tautau
///
void PMNS_NSIs_NO::SetEta_etau(double a, double phi){ SetEta(0,2, a, phi); }


////Set lightest neutrino mass
void PMNS_NSIs_NO::SetM1(double m1){ SetLowestMass(m1);}

double PMNS_NSIs_NO::GetM1(){ return fM1;}
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
void PMNS_NSIs_NO::UpdateHam()
{

  double lv = 2 * kGeV2eV*fEnergy;     // 2E in eV 

  double kr2GNe = kK2*M_SQRT2*kGf;
  kr2GNe *= fPath.density * fPath.zoa; // Matter potential in eV

  // Finish building Hamiltonian in matter with dimension of eV
  double  s12, s23, s13, c12, c23, c13;
  complexD idelta(0.0, fDelta[0][2]);
  if(fIsNuBar) idelta = conj(idelta);

  s12 = sin(fTheta[0][1]);  s23 = sin(fTheta[1][2]);  s13 = sin(fTheta[0][2]);
  c12 = cos(fTheta[0][1]);  c23 = cos(fTheta[1][2]);  c13 = cos(fTheta[0][2]);
  complexD U[3][3];
  U[0][0] =  c12*c13;
  U[0][1] =  s12*c13;
  U[0][2] =  s13*exp(-idelta);

  U[1][0] = -s12*c23-c12*s23*s13*exp(idelta);
  U[1][1] =  c12*c23-s12*s23*s13*exp(idelta);
  U[1][2] =  s23*c13;

  U[2][0] =  s12*s23-c12*c23*s13*exp(idelta);
  U[2][1] = -c12*s23-s12*c23*s13*exp(idelta);
  U[2][2] =  c23*c13;
  
  double dm21=fDm[1];
  double dm31=fDm[2];
  double m1=fM1;
  double m2=sqrt(dm21+m1*m1);
  double m3=sqrt(dm31+m1*m1);
  
   complexD M[3][3];
  
  M[0][0]=m1;
  M[1][1]=m2;
  M[2][2]=m3;
  
  int i,j,k;
  
  for( i=0;i<3;i++){
    
    for( j=i+1;j<3;j++){
      M[i][j]=0;
    }
  }
  
   complexD Udag[3][3];
  
  for(int i=0;i<3;i++){
    
    for(int j=0;j<3;j++){
      Udag[i][j]=conj(U[j][i]);
    }
  }
  
  
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

 complexD eta[3][3];

eta[0][0]=fEta[0][0]*sqrt(dm31);
eta[1][1]=fEta[1][1]*sqrt(dm31);
eta[2][2]=fEta[2][2]*sqrt(dm31);

 eta[0][1]=fEta[0][1]*sqrt(dm31);
eta[0][2]=fEta[0][2]*sqrt(dm31);
eta[1][2]=fEta[1][2]*sqrt(dm31);


eta[1][0]=conj(eta[0][1]);
eta[2][0]=conj(eta[0][2]);
eta[2][1]=conj(eta[1][2]);

  
   complexD MpdM[3][3];
  
  for( i=0;i<3;i++){
    
    for( j=0;j<3;j++){
      MpdM[i][j]=UMUdag[i][j]+eta[i][j];
    }
  }
  
  complexD MpdMdag[3][3];
  
  
  for(int i=0;i<3;i++){
    
    for(int j=0;j<3;j++){
      MpdMdag[i][j]=conj(MpdM[j][i]);
    }
  }
  
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
  
  
  if(!fIsNuBar) fHam[0][0] += kr2GNe;
  else          fHam[0][0] -= kr2GNe;

}


////////////////////////////////////////////////////////////////////////

