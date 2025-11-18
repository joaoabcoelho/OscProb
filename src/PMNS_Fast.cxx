///////////////////////////////////////////////////////////////////////////////
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework.
//
// jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#include "MatrixDecomp/zheevh3.h"

#include "PMNS_Fast.h"

using namespace OscProb;

//.............................................................................
///
/// Constructor. \sa PMNS_Maltoni::PMNS_Maltoni
///
/// This class is restricted to 3 neutrino flavours.
///
// PMNS_Fast::PMNS_Fast() : PMNS_Base(), fHam() {}
PMNS_Fast::PMNS_Fast() : PMNS_Maltoni(3), fHam() {}

//.............................................................................
///
/// Nothing to clean.
///
PMNS_Fast::~PMNS_Fast() {}

//.............................................................................
///
/// Set all mixing parameters at once.
/// @param th12    - The value of the mixing angle theta_12
/// @param th23    - The value of the mixing angle theta_23
/// @param th13    - The value of the mixing angle theta_13
/// @param deltacp - The value of the CP phase delta_13
///
void PMNS_Fast::SetMix(double th12, double th23, double th13, double deltacp)
{
  SetAngle(1, 2, th12);
  SetAngle(1, 3, th13);
  SetAngle(2, 3, th23);
  SetDelta(1, 3, deltacp);
}

//.............................................................................
///
/// Set both mass-splittings at once.
///
/// These are Dm_21 and Dm_32 in eV^2.\n
/// The corresponding Dm_31 is set in the class attributes.
///
/// @param dm21 - The solar mass-splitting Dm_21
/// @param dm32 - The atmospheric mass-splitting Dm_32
///
void PMNS_Fast::SetDeltaMsqrs(double dm21, double dm32)
{
  SetDm(2, dm21);
  SetDm(3, dm32 + dm21);
}

//.............................................................................
///
/// Build the full Hamiltonian in matter.
///
/// Here we divide the mass squared matrix Hms by the 2E
/// to obtain the vacuum Hamiltonian in eV. Then, the matter
/// potential is added to the electron component.
///
void PMNS_Fast::UpdateHam()
{
  double lv = 2 * kGeV2eV * fEnergy; // 2E in eV

  double kr2GNe = kK2 * M_SQRT2 * kGf;
  kr2GNe *= fPath.density * fPath.zoa; // Matter potential in eV

  // Finish building Hamiltonian in matter with dimension of eV
  for (int i = 0; i < fNumNus; i++) {
    fHam[i][i] = fHms[i][i] / lv;
    for (int j = i + 1; j < fNumNus; j++) {
      if (!fIsNuBar)
        fHam[i][j] = fHms[i][j] / lv;
      else
        fHam[i][j] = conj(fHms[i][j]) / lv;
    }
  }
  if (!fIsNuBar)
    fHam[0][0] += kr2GNe;
  else
    fHam[0][0] -= kr2GNe;
}

//.............................................................................
///
/// Solve the full Hamiltonian for eigenvectors and eigenvalues.
///
/// If vacuum, just use the PMNS matrix, otherwise solve in Matter using GLoBES.
///
void PMNS_Fast::SolveHam()
{
  // Do vacuum oscillation in low density
  if (fPath.density < 1.0e-6) {
    SetVacuumEigensystem();
    return;
  }

  SolveHamMatter();
}

//.............................................................................
///
/// Solve the full Hamiltonian for eigenvectors and eigenvalues.
///
/// This is using a method from the GLoBES software available at
/// http://www.mpi-hd.mpg.de/personalhomes/globes/3x3/ \n
/// We should cite them accordingly
///
void PMNS_Fast::SolveHamMatter()
{
  // Build Hamiltonian
  BuildHms();

  // Check if anything changed
  if (fGotES) return;

  // Try caching if activated
  if (TryCache()) return;

  UpdateHam();

  double   fEvalGLoBES[3];
  complexD fEvecGLoBES[3][3];

  // Solve Hamiltonian for eigensystem using the GLoBES method
  zheevh3(fHam, fEvecGLoBES, fEvalGLoBES);

  // Fill fEval and fEvec vectors from GLoBES arrays
  for (int i = 0; i < fNumNus; i++) {
    fEval[i] = fEvalGLoBES[i];
    for (int j = 0; j < fNumNus; j++) { fEvec[i][j] = fEvecGLoBES[i][j]; }
  }

  fGotES = true;

  // Fill cache if activated
  FillCache();
}

//.............................................................................
///
/// Set the eigensystem to the analytic solution in vacuum.
///
/// We know the vacuum eigensystem, so just write it explicitly
///
void PMNS_Fast::SetVacuumEigensystem()
{
  double   s12, s23, s13, c12, c23, c13;
  complexD idelta(0.0, fDelta[0][2]);
  if (fIsNuBar) idelta = conj(idelta);

  s12 = sin(fTheta[0][1]);
  s23 = sin(fTheta[1][2]);
  s13 = sin(fTheta[0][2]);
  c12 = cos(fTheta[0][1]);
  c23 = cos(fTheta[1][2]);
  c13 = cos(fTheta[0][2]);

  fEvec[0][0] = c12 * c13;
  fEvec[0][1] = s12 * c13;
  fEvec[0][2] = s13 * exp(-idelta);

  fEvec[1][0] = -s12 * c23 - c12 * s23 * s13 * exp(idelta);
  fEvec[1][1] = c12 * c23 - s12 * s23 * s13 * exp(idelta);
  fEvec[1][2] = s23 * c13;

  fEvec[2][0] = s12 * s23 - c12 * c23 * s13 * exp(idelta);
  fEvec[2][1] = -c12 * s23 - s12 * c23 * s13 * exp(idelta);
  fEvec[2][2] = c23 * c13;

  fEval[0] = 0;
  fEval[1] = fDm[1] / (2 * kGeV2eV * fEnergy);
  fEval[2] = fDm[2] / (2 * kGeV2eV * fEnergy);
}

///////////////////////////////////////////////////////////////////////////////
