////////////////////////////////////////////////////////////////////////
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework with LIV.
//
// This  class inherits from the PMNS_Fast class
//
// N. R. Khan Chowdhury
////////////////////////////////////////////////////////////////////////

#include "PMNS_LIV.h"
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
//PMNS_LIV::PMNS_LIV() : PMNS_Fast(), faT(), fcTT()
PMNS_LIV::PMNS_LIV() : PMNS_Fast()
{
  SetStdPath();
  //SetLIV(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
  SetLIV(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
}

//......................................................................
///
/// Nothing to clean.
///
PMNS_LIV::~PMNS_LIV(){}

//......................................................................
///
/// Set all LIV parameters at once.
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param aT_ij       - The absolute value of the parameter aT_ij
/// @param delta_aT_ij - The phase of the complex parameter aT_ij in radians
/// @param cT_ij       - The absolute value of the complex parameter cT_ij
/// @param delta_cT_ij - The phase of the complex parameter cT_ij in radians
///
void PMNS_LIV::SetLIV(double aT_ee,     double aT_mumu,      double aT_tautau,
		                double aT_emu,     double aT_etau,      double aT_mutau,
		                double cT_ee,     double cT_mumu,      double cT_tautau,
                                double cT_emu,     double cT_etau,      double cT_mutau,
			        double delta_aT_emu, double delta_aT_etau, double delta_aT_mutau,
			        double delta_cT_emu, double delta_cT_etau, double delta_cT_mutau)
{

  SetaT(0,0, aT_ee, 0);
  SetaT(1,1, aT_mumu, 0);
  SetaT(2,2, aT_tautau, 0);

  SetaT(0,1, aT_emu, delta_aT_emu);
  SetaT(0,2, aT_etau, delta_aT_etau);
  SetaT(1,2, aT_mutau, delta_aT_mutau);

  SetcT(0,0, cT_ee, 0);
  SetcT(1,1, cT_mumu, 0);
  SetcT(2,2, cT_tautau, 0);

  SetcT(0,1, cT_emu, delta_cT_emu);
  SetcT(0,2, cT_etau, delta_cT_etau);
  SetcT(1,2, cT_mutau, delta_cT_mutau);


}

//......................................................................
///
/// Set any given LIV1 parameter (aT & cT).
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// Flavours are:\n
/// - 0 = nue, 1 = numu, 2 = nutau
///
/// Requires that flvi < flvj. Will notify you if input is wrong.
/// If flvi > flvj, will assume reverse order and swap flvi and flvj.
///
/// @param flvi  - The first flavour index
/// @param flvj  - The second flavour index
/// @param val   - The absolute value of the parameter
/// @param phase - The complex phase of the parameter in radians
///
void PMNS_LIV::SetaT(int flvi, int flvj, double val, double phase){

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

  complexD h = val;

  if(flvi != flvj) h *= complexD(cos(phase), sin(phase));

  //bool isSame = (faT[flvi][flvj] == h);

  //if(!isSame) ClearCache();

  fGotES *= (faT[flvi][flvj] == h);

  faT[flvi][flvj] = h;

}

void PMNS_LIV::SetcT(int flvi, int flvj, double val, double phase){

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

  complexD h = val;

  if(flvi != flvj) h *= complexD(cos(phase), sin(phase));

  //bool isSame = (fcTT[flvi][flvj] == h);

//  if(!isSame) ClearCache();

  //fGotES *= isSame;
  fGotES *= (faT[flvi][flvj] == h);

  fcT[flvi][flvj] = h;

}

//......................................................................
///
/// Get any given LIV1 parameter (aT & cTT).
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
complex<double> PMNS_LIV::GetaT(int flvi, int flvj){

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

  return faT[flvi][flvj];

}

complex<double> PMNS_LIV::GetcT(int flvi, int flvj){

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

  return fcT[flvi][flvj];

}

//......................................................................
///
/// Set aT_ij parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a   - The absolute value of the parameter
/// @param phi - The complex phase of the parameter in radians
///
void PMNS_LIV::SetaT_ee(double a){ SetaT(0,0, a, 0); }
void PMNS_LIV::SetaT_mumu(double a){ SetaT(1,1, a, 0); }
void PMNS_LIV::SetaT_tautau(double a){ SetaT(2,2, a, 0); }

void PMNS_LIV::SetaT_emu(double a, double phi){ SetaT(0,1, a, phi); }
void PMNS_LIV::SetaT_etau(double a, double phi){ SetaT(0,2, a, phi); }
void PMNS_LIV::SetaT_mutau(double a, double phi){ SetaT(1,2, a, phi); }

//......................................................................
///
/// Set cTT_ij parameter
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param a   - The absolute value of the parameter
/// @param phi - The complex phase of the parameter in radians
///
void PMNS_LIV::SetcT_ee(double a){ SetcT(0,0, a, 0); }
void PMNS_LIV::SetcT_mumu(double a){ SetcT(1,1, a, 0); }
void PMNS_LIV::SetcT_tautau(double a){ SetcT(2,2, a, 0); }

void PMNS_LIV::SetcT_emu(double a, double phi){ SetcT(0,1, a, phi); }
void PMNS_LIV::SetcT_etau(double a, double phi){ SetcT(0,2, a, phi); }
void PMNS_LIV::SetcT_mutau(double a, double phi){ SetcT(1,2, a, phi); }

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
  double fEnergy_ev =     kGeV2eV*fEnergy;           // E in eV

  double kr2GNe   = kK2*M_SQRT2*kGf*fPath.density;
  kr2GNe *=  fPath.zoa; // Matter potential in eV

  // Finish build Hamiltonian in matter with dimension of eV
  for(int i=0;i<fNumNus;i++){
    for(int j=i;j<fNumNus;j++){
      if(!fIsNuBar) fHam[i][j] = fHms[i][j]/lv + 1e9*faT[i][j] - (4.*fEnergy_ev / 3.)*1e9*fcT[i][j];
      else          fHam[i][j] = conj(fHms[i][j]/lv - 1e9*faT[i][j] - (4.*fEnergy_ev / 3.)*1e9*fcT[i][j]);
    }
  }

  if(!fIsNuBar) fHam[0][0] += kr2GNe;
  else          fHam[0][0] -= kr2GNe;

}

////////////////////////////////////////////////////////////////////////
