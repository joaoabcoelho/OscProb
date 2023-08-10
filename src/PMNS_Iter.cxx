///////////////////////////////////////////////////////////////////////////////
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework.
//
// jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#include "PMNS_Iter.h"

using namespace OscProb;

//.............................................................................
///
/// Constructor. \sa PMNS_Fast::PMNS_Fast
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_Iter::PMNS_Iter() : PMNS_Fast(), fPrec(1e-3) {}

//.............................................................................
///
/// Nothing to clean.
///
PMNS_Iter::~PMNS_Iter(){}


//.............................................................................
///
/// Propagate neutrino through matter component.
///
void PMNS_Iter::SetExpVL(NuPath p)
{

  double kr2GNe = kK2*M_SQRT2*kGf;
  kr2GNe *= p.density * p.zoa; // Matter potential in eV

  fVL = kr2GNe * kKm2eV * p.length;

  fExpVL = complexD(cos(fVL), -sin(fVL));
  if(fIsNuBar) fExpVL = conj(fExpVL);

}

//.............................................................................
///
/// Propagate neutrino through matter component.
///
void PMNS_Iter::PropMatter()
{

  fNuState[0] *= fExpVL;

}

//.............................................................................
///
/// Just use the vacuum Hamiltonian to start.
///
void PMNS_Iter::SolveHam()
{

  if(fGotES) return;

  // Do vacuum oscillation
  if(!fBuiltHms){

    fPath.density = 0;
    PMNS_Fast::SolveHam();
    fBuiltHms = true;

  }
  else {

    for(int i=1; i<fNumNus; i++){
      fEval[i] *= fPrevEnergy / fEnergy;
    }

  }

  fPrevEnergy = fEnergy;

  fGotES = true;
  return;

}

//.............................................................................
///
/// Set iterator precision
///
void PMNS_Iter::SetPrec(double prec)
{
  fPrec = prec;
  if(fPrec<=0) fPrec = 1e-3;
}

//.............................................................................
///
/// Propagate neutrino state through split path
///
void PMNS_Iter::PropagatePath(NuPath p)
{

  SetExpVL(p);

  double dm = 0;
  for(int i=0; i<fNumNus; i++){
    if(dm<fabs(fDm[i])) dm = fabs(fDm[i]);
  }
  dm *= kKm2eV * p.length / (2 * kGeV2eV * fEnergy);

  int nsplit = sqrt(0.065 * dm * fVL / fPrec);
  if(nsplit>1){
    p.length /= nsplit;
    SetExpVL(p);
  }
  else nsplit = 1;

  for(int i=0; i<nsplit; i++){

    PMNS_Base::PropagatePath(p);
    PropMatter();

  }

}

///////////////////////////////////////////////////////////////////////////////
