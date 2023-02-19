///////////////////////////////////////////////////////////////////////////////                   
/// class OscProb::PMNS_LIV                                                                   
///                                                                                                     
/// Implementation of neutrino oscillations in matter in a                                                
/// three-neutrino framework with LIV as modelled by the SME.                                               
/// The SME coefficients are included up to the 8th order.                                            
///                                                                                                        
/// This developement is part of the QGRANT project with                                                   
/// ID: 101068013,                                                                                        
/// founded by the HORIZON-MSCA-2021-PF-01-01 programme.                                               
///                                                                                                     
/// \author Nafis R. K. Chowdhury - nrkhanchowdhury\@km3net.de                                             
/// \author Joao Coelho - jcoelho\@apc.in2p3.fr                                                           
/// \author Alba Domi - alba.domi\@fau.de                       
///////////////////////////////////////////////////////////////////////////////  

#include "PMNS_LIV.h"
#include "MatrixDecomp/zheevh3.h"

#include <iostream>
#include <cassert>
#include <stdlib.h>

using namespace OscProb;

using namespace std;

//...................................................................... 
/// Constructor. 
/// This class is restricted to 3 neutrino flavours.
/// By default, all LIV coefficients are set to zero.
///
PMNS_LIV::PMNS_LIV() : PMNS_Fast()
{
  SetStdPath();

  for(int dim = 3; dim <= 7; dim += 2){
    SetaT(0,0, dim, 0, 0);
    SetaT(1,1, dim, 0, 0);
    SetaT(2,2, dim, 0, 0);
    SetaT(0,1, dim, 0, 0);
    SetaT(0,2, dim, 0, 0);
    SetaT(1,2, dim, 0, 0);
  }

  for(int dim = 4; dim <= 8; dim += 2){
    SetcT(0,0, 4, 0, 0);
    SetcT(1,1, 4, 0, 0);
    SetcT(2,2, 4, 0, 0);
    SetcT(0,1, 4, 0, 0);
    SetcT(0,2, 4, 0, 0);
    SetcT(1,2, 4, 0, 0);
  }
}

//......................................................................
///
/// Distructor.
///
PMNS_LIV::~PMNS_LIV(){}

//......................................................................
/// Set any given aT parameter of a chosen dimension.
///
/// Flavours are:\n
/// - 0 = nue, 1 = numu, 2 = nutau
///
/// Requires that flvi < flvj. Will notify you if input is wrong.
/// If flvi > flvj, will assume reverse order and swap flvi and flvj.
///
/// @param flvi  - First flavour index
/// @param flvj  - Second flavour index
/// @param dim   - Dimension of the coefficient: can be 3,5,7.
/// @param val   - Absolute value of the parameter
/// @param phase - Complex phase of the parameter in radians
///
void PMNS_LIV::SetaT(int flvi, int flvj, int dim, double val, double phase){

  if(flvi > flvj){
    cout << "First argument should be smaller or equal to second argument" << endl;
    cout << "Setting reverse order (aT_" << flvj << flvi << "). " << endl;
    int temp = flvi;
    flvi = flvj;
    flvj = temp;
  }
  if(flvi<0 || flvi>2 || flvj < flvi || flvj > 2){
    cout << "aT_" << flvi << flvj << " not valid for " << fNumNus;
    cout << " neutrinos. Doing nothing." << endl;
    return;
  }
  if(dim != 3 && dim != 5 && dim != 7){
    cout << "aT coefficient dimension should be an odd value, starting from 3\n";
    return;
  }

  int pos = dim == 3 ? 0: dim == 5 ? 1 : dim == 7 ? 2 : 4;
  complexD h[3];
  h[pos] = val;

  if(flvi != flvj) h[pos] *= complexD(cos(phase), sin(phase));

  bool isSame = (faT[flvi][flvj][pos] == h[pos]);

  if(!isSame) ClearCache();

  fGotES *= isSame;

  faT[flvi][flvj][pos] = h[pos];

}

//......................................................................                                  
/// Set any given cT parameter of a chosen dimension.                         
///                                                                                                       
/// Flavours are:\n                                                                                            
/// - 0 = nue, 1 = numu, 2 = nutau                                                                         
///                                                                                                        
/// Requires that flvi < flvj. Will notify you if input is wrong.                                           
/// If flvi > flvj, will assume reverse order and swap flvi and flvj.                                      
///                                                                                                         
/// @param flvi  - First flavour index                                                                          
/// @param flvj  - Second flavour index                                                                      
/// @param dim   - Dimension of the coefficient: can be 4,6,8.      
/// @param val   - Absolute value of the parameter                                                        
/// @param phase - Complex phase of the parameter in radians                                         
///                                                                    
void PMNS_LIV::SetcT(int flvi, int flvj, int dim, double val, double phase){

  if(flvi > flvj){
    cout << "First argument should be smaller or equal to second argument" << endl;
    cout << "Setting reverse order (cT_" << flvj << flvi << "). " << endl;
    int temp = flvi;
    flvi = flvj;
    flvj = temp;
  }
  if(flvi<0 || flvi>2 || flvj < flvi || flvj > 2){
    cout << "cT_" << flvi << flvj << " not valid for " << fNumNus;
    cout << " neutrinos. Doing nothing." << endl;
    return;
  }
  if(dim != 4 && dim != 6 && dim != 8){
    cout << "cT coefficient dimension should be an even value, starting from 4\n";
    return;
  }

  int pos = dim == 4 ? 0: dim == 6 ? 1 : dim == 8 ? 2 : 4;
  complexD h[3];
  h[pos] = val;

  if(flvi != flvj) h[pos] *= complexD(cos(phase), sin(phase));

  bool isSame = (fcT[flvi][flvj][pos] == h[pos]);

  if(!isSame) ClearCache();

  fGotES *= isSame;

  fcT[flvi][flvj][pos] = h[pos];

}

//......................................................................
///
/// Get any given aT parameter of a chosen dimension.  
///
/// Flavours are:\n
/// - 0 = nue, 1 = numu, 2 = nutau
///
/// Requires that flvi < flvj. Will notify you if input is wrong.
/// If flvi > flvj, will assume reverse order and swap flvi and flvj.
///
/// @param flvi  - The first flavour index
/// @param flvj  - The second flavour index
/// @param dim   - Dimension of the coefficient: can be 3,5,7.  
///
complex<double> PMNS_LIV::GetaT(int flvi, int flvj, int dim){

  if(flvi > flvj){
    cout << "First argument should be smaller or equal to second argument" << endl;
    cout << "Setting reverse order (aT_" << flvj << flvi << "). " << endl;
    int temp = flvi;
    flvi = flvj;
    flvj = temp;
  }
  if(flvi<0 || flvi>2 || flvj < flvi || flvj > 2){
    cout << "aT_" << flvi << flvj << " not valid for " << fNumNus;
    cout << " neutrinos. Returning 0." << endl;
    return zero;
  }
  if(dim != 3 && dim != 5 && dim != 7){
    cout << "aT coefficient dimension should be an odd value, starting from 3\n";
    return zero;
  }

  int pos = dim == 3 ? 0: dim == 5 ? 1 : dim == 7 ? 2 : 4;
  
  return faT[flvi][flvj][pos];

}

//......................................................................
///
/// Get any given cT parameter of a chosen dimension.
///
/// Flavours are:\n
/// - 0 = nue, 1 = numu, 2 = nutau
///
/// Requires that flvi < flvj. Will notify you if input is wrong.
/// If flvi > flvj, will assume reverse order and swap flvi and flvj.
///
/// @param flvi  - The first flavour index
/// @param flvj  - The second flavour index
/// @param dim   - Dimension of the coefficient: can be 4,6,8. 
///
complex<double> PMNS_LIV::GetcT(int flvi, int flvj, int dim){

  if(flvi > flvj){
    cout << "First argument should be smaller or equal to second argument" << endl;
    cout << "Setting reverse order (cT_" << flvj << flvi << "). " << endl;
    int temp = flvi;
    flvi = flvj;
    flvj = temp;
  }
  if(flvi<0 || flvi>2 || flvj < flvi || flvj > 2){
    cout << "cT_" << flvi << flvj << " not valid for " << fNumNus;
    cout << " neutrinos. Returning 0." << endl;
    return zero;
  }
  if(dim != 4 && dim != 6 && dim != 8){
    cout << "cT coefficient dimension should be an even value, starting from 4\n";
    return zero;
  }
  
  int pos = dim == 4 ? 0: dim == 6 ? 1 : dim == 8 ? 2 : 4;
  
  return fcT[flvi][flvj][pos];

}


//......................................................................
///
/// Build the full Hamiltonian in matter
///
/// Here we divide the mass squared matrix Hms by the 2E
/// to obtain the vacuum Hamiltonian in eV. Then, the matter
/// potential is added to each flavour pair.
///
void PMNS_LIV::UpdateHam()
{

  double lv = 2 * kGeV2eV*fEnergy;           // 2*E in eV

  double kr2GNe = kK2*M_SQRT2*kGf*fPath.density;
  kr2GNe *=  fPath.zoa; // Matter potential in eV

  // Finish build Hamiltonian in matter with dimension of eV
  for(int i=0;i<fNumNus;i++){
    for(int j=i;j<fNumNus;j++){
      if(!fIsNuBar){
	fHam[i][j] = fHms[i][j]/lv + kGeV2eV*(faT[i][j][0]
					      -4./3. * fEnergy * fcT[i][j][0]
					      +std::pow(fEnergy, 2) * faT[i][j][1]
					      -4./3. * std::pow(fEnergy, 3) * fcT[i][j][1]
					      +std::pow(fEnergy, 4) * faT[i][j][2]
					      -4./3. * std::pow(fEnergy, 5) * fcT[i][j][2]);

      } else {
	fHam[i][j] = conj(fHms[i][j]/lv + kGeV2eV*(-faT[i][j][0]
						   -4./3. * fEnergy * fcT[i][j][0]
						   -std::pow(fEnergy, 2) * faT[i][j][1]
						   -4./3. * std::pow(fEnergy, 3) * fcT[i][j][1]
						   -std::pow(fEnergy, 4) * faT[i][j][2]
						   -4./3. * std::pow(fEnergy, 5) * fcT[i][j][2]));
      }
    }
  }

  if(!fIsNuBar) fHam[0][0] += kr2GNe;
  else          fHam[0][0] -= kr2GNe;

}

////////////////////////////////////////////////////////////////////////
