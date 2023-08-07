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
PMNS_Iter::PMNS_Iter() : PMNS_Fast(), fPrec(0.1) {}

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

  // Do vacuum oscillation in low density
  SetVacuumEigensystem();

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

  int nsplit = sqrt(dm * fVL / fPrec);
  nsplit++;

  p.length /= nsplit;
  SetExpVL(p);

  for(int i=0; i<nsplit; i++){

    PMNS_Base::PropagatePath(p);
    PropMatter();

  }

}

//.............................................................................
///
/// Set the eigensystem to the analytic solution in vacuum.
///
/// We know the vacuum eigensystem, so just write it explicitly
///
void PMNS_Iter::SetVacuumEigensystem()
{

  if(!fBuiltHms){

    double  s12, s23, s13, c12, c23, c13;
    complexD idelta(0.0, fDelta[0][2]);
    if(fIsNuBar) idelta = conj(idelta);

    s12 = sin(fTheta[0][1]);  s23 = sin(fTheta[1][2]);  s13 = sin(fTheta[0][2]);
    c12 = cos(fTheta[0][1]);  c23 = cos(fTheta[1][2]);  c13 = cos(fTheta[0][2]);

    fEvec[0][0] =  c12*c13;
    fEvec[0][1] =  s12*c13;
    fEvec[0][2] =  s13*exp(-idelta);

    fEvec[1][0] = -s12*c23-c12*s23*s13*exp(idelta);
    fEvec[1][1] =  c12*c23-s12*s23*s13*exp(idelta);
    fEvec[1][2] =  s23*c13;

    fEvec[2][0] =  s12*s23-c12*c23*s13*exp(idelta);
    fEvec[2][1] = -c12*s23-s12*c23*s13*exp(idelta);
    fEvec[2][2] =  c23*c13;

  }

  fBuiltHms = true;

  fEval[0] = 0;
  fEval[1] = fDm[1] / (2 * kGeV2eV*fEnergy);
  fEval[2] = fDm[2] / (2 * kGeV2eV*fEnergy);

}

///////////////////////////////////////////////////////////////////////////////
