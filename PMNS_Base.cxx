////////////////////////////////////////////////////////////////////////
//
// Base class for oscillations of neutrinos in matter in a
// n-neutrino framework.
//
//......................................................................
//
// coelho@lal.in2p3.fr
//
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cassert>
#include <stdlib.h>
#include <algorithm>
#include "PMNS_Base.h"

using namespace std;

using namespace OscProb;

// Some usefule complex numbers
const complexD PMNS_Base::zero(0,0);
const complexD PMNS_Base::one(1,0);

// Define some constants from PDG 2015
const double PMNS_Base::kGeV2eV = 1.0e+09;                    // GeV to eV conversion
const double PMNS_Base::kKm2eV  = 1.0 / 1.973269788e-10;      // (hbar.c [eV.km])^-1
const double PMNS_Base::kNA     = 6.022140857e23;             // Avogadro constant (N_A)
const double PMNS_Base::kK2     = 1e-3 * kNA / pow(kKm2eV,3); // N_A * (hbar*c [GeV.cm])^3 * kGeV2eV
const double PMNS_Base::kGf     = 1.1663787e-05;              // G_F/(hbar*c)^3 [GeV^-2]


//......................................................................
///
/// Constructor.
///
/// Sets the number of neutrinos and initializes attributes
///
/// Default starts with a 2 GeV muon neutrino.
///
/// Path is set to the default 1000 km in crust density.
///
/// Oscillation parameters are from PDG for NH by default.
///
/// @param numNus - the number of neutrino flavours
///
PMNS_Base::PMNS_Base(int numNus) :
fGotES(false), fBuiltHms(false), fMaxCache(1e6), fProbe(numNus)
{

  SetUseCache(false);  // Don't cache eigensystems

  fNumNus = numNus;    // Set the number of neutrinos

  SetStdPath();        // Set some default path
  SetEnergy(2);        // Set default energy to 2 GeV
  SetIsNuBar(false);   // Neutrino by default

  InitializeVectors(); // Initialize all vectors

  SetStdPars();        // Set PDG parameters

  ResetToFlavour(1);   // Numu by default
  
}

//......................................................................
///
/// Nothing to clean.
///
PMNS_Base::~PMNS_Base(){}

//......................................................................
///
/// Set vector sizes and initialize elements to zero.
///
void PMNS_Base::InitializeVectors()
{

  fDm    = vector<double>(fNumNus, 0);
  fTheta = vector< vector<double> >(fNumNus, vector<double>(fNumNus,0));
  fDelta = vector< vector<double> >(fNumNus, vector<double>(fNumNus,0));

  fNuState = vector<complexD>(fNumNus, zero);
  fHms     = vector< vector<complexD> >(fNumNus, vector<complexD>(fNumNus,zero));

  fPhases = vector<complexD>(fNumNus, zero);
  fBuffer = vector<complexD>(fNumNus, zero);

  fEval = vector<double>(fNumNus, 0);
  fEvec = vector< vector<complexD> >(fNumNus, vector<complexD>(fNumNus,zero));

}

//......................................................................
///
/// Turn on/off caching of eigensystems.
/// This can save a lot of CPU time by avoiding recomputing eigensystems
/// if we've already seen them recently.
/// Especially useful when running over multiple earth layers and even more
/// if multiple baselines will be computed, e.g. for atmospheric neutrinos.
///
/// @param u - flag to set caching on (default: true)
///
void PMNS_Base::SetUseCache(bool u)
{
  fUseCache = u; 
}

//......................................................................
///
/// Clear the cache
///
void PMNS_Base::ClearCache()
{
  fMixCache.clear();
}
        
//......................................................................
///
/// Set maximum number of cached eigensystems.
/// Finding eigensystems can become slow and take up memory.
/// This protects the cache from becoming too large.
///
/// @param mc - Max cache size (default: 1e6)
///
void PMNS_Base::SetMaxCache(int mc)
{
  fMaxCache = mc;
}

//......................................................................
///
/// Try to find a cached version of this eigensystem.
///
bool PMNS_Base::TryCache()
{

  if(fUseCache && !fMixCache.empty()){

    fProbe.SetVars(fEnergy, fPath, fIsNuBar);
    
    std::set<EigenPoint>::iterator it = fMixCache.find(fProbe);

    if(it != fMixCache.end()){
      for(int i=0; i<fNumNus; i++){
        fEval[i] = (*it).fEval[i] * (*it).fEnergy / fEnergy;
        for(int j=0; j<fNumNus; j++){
          fEvec[i][j] = (*it).fEvec[i][j];
        }
      }
      return true;
    }

  }
  
  return false;

}

//......................................................................
///
/// If using caching, save the eigensystem in memory
///
void PMNS_Base::FillCache()
{

  if(fUseCache){
    if(fMixCache.size()>fMaxCache){
      fMixCache.erase(fMixCache.begin());
      fMixCache.erase(--fMixCache.end());
    }
    for(int i=0; i<fNumNus; i++){
      fProbe.fEval[i] = fEval[i];
      for(int j=0; j<fNumNus; j++){
        fProbe.fEvec[i][j] = fEvec[i][j];
      }
    }
    fMixCache.insert(fProbe);
  }

}

//......................................................................
///
/// Set standard oscillation parameters from PDG 2015.
///
/// For two neutrinos, Dm is set to the muon disappearance
/// effective mass-splitting and mixing angle.
///
void PMNS_Base::SetStdPars()
{

  if(fNumNus>2){
    // PDG values for 3 neutrinos
    // Also applicable for 3+N neutrinos
    SetAngle(1,2, asin(sqrt(0.304)));
    SetAngle(1,3, asin(sqrt(0.0219)));
    SetAngle(2,3, asin(sqrt(0.514)));
    SetDm(2, 7.53e-5);
    SetDm(3, 2.52e-3);
  }
  else if(fNumNus==2){
    // Effective muon disappearance values
    // for two-flavour approximation
    SetAngle(1,2, 0.788);
    SetDm(2, 2.47e-3);
  }

}

//......................................................................
///
/// Set standard single path.
///
/// Length is 1000 km, so ~2 GeV peak energy.
///
/// Density is approximate from CRUST2.0 (~2.8 g/cm^3).
/// Z/A is set to a round 0.5.
///
void PMNS_Base::SetStdPath(){

  NuPath p;

  p.length  = 1000; // 1000 km default
  p.density = 2.8;  // Crust density
  p.zoa     = 0.5;  // Crust Z/A
  p.layer   = 0;    // Single layer

  SetPath(p);

}

//......................................................................
///
/// Set neutrino energy in GeV.
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param E - The neutrino energy in GeV
///
void PMNS_Base::SetEnergy(double E)
{

  // Check if value is actually changing
  fGotES *= (fEnergy == E);

  fEnergy = E;

}

//......................................................................
///
/// Set anti-neutrino flag.
///
/// This will check if value is changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param isNuBar - Set to true for anti-neutrino and false for neutrino.
///
void PMNS_Base::SetIsNuBar(bool isNuBar)
{

  // Check if value is actually changing
  fGotES *= (fIsNuBar == isNuBar);

  fIsNuBar = isNuBar;

}

//......................................................................
///
/// Get the neutrino energy in GeV.
///
double PMNS_Base::GetEnergy() {

  return fEnergy;

}

//......................................................................
///
/// Get the anti-neutrino flag.
///
bool PMNS_Base::GetIsNuBar() {

  return fIsNuBar;

}

//......................................................................
///
/// Set the path currentlyin use by the class.
///
/// This will be used to know what path to propagate through next.
///
/// It will also check if values are changing to keep track of whether
/// the eigensystem needs to be recomputed.
///
/// @param p - A neutrino path segment
///
void PMNS_Base::SetCurPath(NuPath p)
{

  // Check if relevant value are actually changing
  fGotES *= (fPath.density == p.density);
  fGotES *= (fPath.zoa == p.zoa);

  fPath = p;

}

//......................................................................
///
/// Clear the path vector.
///
void PMNS_Base::ClearPath(){

  fNuPaths.clear();

}

//......................................................................
///
/// Set vector of neutrino paths.
/// @param paths - A sequence of neutrino paths
///
void PMNS_Base::SetPath(std::vector<NuPath> paths){

  fNuPaths=paths;

}

//......................................................................
///
/// Get the vector of neutrino paths.
///
vector<NuPath> PMNS_Base::GetPath(){

  return fNuPaths;

}

//......................................................................
///
/// Add a path to the sequence.
/// @param p - A neutrino path segment
///
void PMNS_Base::AddPath(NuPath p){

  fNuPaths.push_back(p);

}

//......................................................................
///
/// Add a path to the sequence defining attributes directly.
/// @param length  - The length of the path segment in km
/// @param density - The density of the path segment in g/cm^3
/// @param zoa     - The effective Z/A of the path segment
/// @param layer   - An index to identify the layer type (e.g. earth inner core)
///
void PMNS_Base::AddPath(double length, double density, double zoa, int layer){

  AddPath(NuPath(length, density, zoa, layer));

}

//......................................................................
///
/// Set a single path.
///
/// This destroys the current path sequence and creates a new first path.
///
/// @param p - A neutrino path segment
///
void PMNS_Base::SetPath(NuPath p){

  ClearPath();
  AddPath(p);

}

//......................................................................
///
/// Set a single path defining attributes directly.
///
/// This destroys the current path sequence and creates a new first path.
///
/// @param length  - The length of the path segment in km
/// @param density - The density of the path segment in g/cm^3
/// @param zoa     - The effective Z/A of the path segment
/// @param layer   - An index to identify the layer type (e.g. earth inner core)
///
void PMNS_Base::SetPath(double length, double density, double zoa, int layer){

  SetPath(NuPath(length, density, zoa, layer));

}

//......................................................................
///
/// Set some single path attribute.
///
/// An auxiliary function to set individual attributes in a single path.
///
/// If the path sequence is not a single path, a new single path will
/// be created and the previous sequence will be lost. Use with care.
///
/// @param att - The value of the attribute
/// @param idx - The index of the attribute (0,1,2,3) = (L, Rho, Z/A, Layer)
///
void PMNS_Base::SetAtt(double att, int idx){

  if(fNuPaths.size() != 1){

    cout << "Warning: Clearing path vector and starting new single path." << endl;
    cout << "To avoid possible issues, use the SetPath function." << endl;

    SetStdPath();

  }

  switch(idx){
    case 0: fNuPaths[0].length  = att; break;
    case 1: fNuPaths[0].density = att; break;
    case 2: fNuPaths[0].zoa     = att; break;
    case 3: fNuPaths[0].layer   = att; break;
  }

}

//......................................................................
///
/// Set the length for a single path.
///
/// If the path sequence is not a single path, a new single path will
/// be created and the previous sequence will be lost. Use with care.
///
/// @param L - The length of the path segment in km
///
void PMNS_Base::SetLength(double L){

  SetAtt(L, 0);

}

//......................................................................
///
/// Set single path density.
///
/// If the path sequence is not a single path, a new single path will
/// be created and the previous sequence will be lost. Use with care.
///
/// @param rho - The density of the path segment in g/cm^3
///
void PMNS_Base::SetDensity(double rho){

  SetAtt(rho, 1);

}

//......................................................................
///
/// Set single path Z/A.
///
/// If the path sequence is not a single path, a new single path will
/// be created and the previous sequence will be lost. Use with care.
///
/// @param zoa - The effective Z/A of the path segment
///
void PMNS_Base::SetZoA(double zoa){

  SetAtt(zoa, 2);

}

//......................................................................
///
/// Set all values of a path attribute.
///
/// An auxiliary function to set individual attributes in a path sequence.
///
/// If the path sequence is of a different size, a new path sequence will
/// be created and the previous sequence will be lost. Use with care.
///
/// @param att - The values of the attribute
/// @param idx - The index of the attribute (0,1,2,3) = (L, Rho, Z/A, Layer)
///
void PMNS_Base::SetAtt(std::vector<double> att, int idx){

  // Get the sizes of the attribute and
  // path sequence vectors
  int nA = att.size();
  int nP = fNuPaths.size();

  // If the vector sizes are equal, update this attribute
  if(nA == nP){

    for(int i=0; i<nP; i++){

      switch(idx){
        case 0: fNuPaths[i].length  = att[i]; break;
        case 1: fNuPaths[i].density = att[i]; break;
        case 2: fNuPaths[i].zoa     = att[i]; break;
        case 3: fNuPaths[i].layer   = att[i]; break;
      }

    }

  }
  // If the vector sizes differ, create a new path sequence
  // and set value for this attribute. Other attributes will
  // be taken from default single path.
  else{

    cout << "Warning: New vector size. Starting new path vector." << endl;
    cout << "To avoid possible issues, use the SetPath function." << endl;

    // Start a new standard path just
    // to set default values
    SetStdPath();

    // Create a path segment with default values
    NuPath p = fNuPaths[0];

    // Clear the path sequence
    ClearPath();

    // Set this particular attribute's value
    // and add the path segment to the sequence
    for(int i=0; i<nA; i++){

      switch(idx){
        case 0: p.length  = att[i]; break;
        case 1: p.density = att[i]; break;
        case 2: p.zoa     = att[i]; break;
        case 3: p.layer   = att[i]; break;
      }

      AddPath(p);

    }

  }

}

//......................................................................
///
/// Set multiple path lengths.
///
/// If the path sequence is of a different size, a new path sequence will
/// be created and the previous sequence will be lost. Use with care.
///
/// @param L - The lengths of the path segments in km
///
void PMNS_Base::SetLength(std::vector<double> L){

  SetAtt(L, 0);

}

//......................................................................
///
/// Set multiple path densities.
///
/// If the path sequence is of a different size, a new path sequence will
/// be created and the previous sequence will be lost. Use with care.
///
/// @param rho - The densities of the path segments in g/cm^3
///
void PMNS_Base::SetDensity(std::vector<double> rho){

  SetAtt(rho, 1);

}

//......................................................................
///
/// Set multiple path Z/A values.
///
/// If the path sequence is of a different size, a new path sequence will
/// be created and the previous sequence will be lost. Use with care.
///
/// @param zoa - The effective Z/A of the path segments
///
void PMNS_Base::SetZoA(std::vector<double> zoa){

  SetAtt(zoa, 2);

}

//......................................................................
///
/// Set multiple path layer indices.
///
/// If the path sequence is of a different size, a new path sequence will
/// be created and the previous sequence will be lost. Use with care.
///
/// @param lay - Indices to identify the layer types (e.g. earth inner core)
///
void PMNS_Base::SetLayers(std::vector<int> lay){

  vector<double> lay_double(lay.begin(), lay.end());

  SetAtt(lay_double, 3);

}

//......................................................................
///
/// Set the mixing angle theta_ij in radians.
///
/// Requires that i<j. Will notify you if input is wrong.
/// If i>j, will assume reverse order and swap i and j.
///
/// This will check if value is changing to keep track of whether
/// the hamiltonian needs to be rebuilt.
///
/// @param i,j - the indices of theta_ij
/// @param th  - the value of theta_ij
///
void PMNS_Base::SetAngle(int i, int j, double th)
{

  if(i>j){
    cout << "Warning: First argument should be smaller than second argument" << endl;
    cout << "         Setting reverse order (Theta" << j << i << "). " << endl;
    int temp = i;
    i = j;
    j = temp;
  }
  if(i<1 || i>fNumNus-1 || j<2 || j>fNumNus){
    cout << "ERROR: Theta" << i << j << " not valid for " << fNumNus;
    cout << " neutrinos. Doing nothing." << endl;
    return;
  }

  // Check if value is actually changing
  fBuiltHms *= (fTheta[i-1][j-1] == th);

  fTheta[i-1][j-1] = th;

}

//......................................................................
///
/// Get the mixing angle theta_ij in radians.
///
/// Requires that i<j. Will notify you if input is wrong.
/// If i>j, will assume reverse order and swap i and j.
///
/// @param i,j - the indices of theta_ij
///
double PMNS_Base::GetAngle(int i, int j)
{

  if(i>j){
    cout << "Warning: First argument should be smaller than second argument" << endl;
    cout << "         Setting reverse order (Theta" << j << i << "). " << endl;
    int temp = i;
    i = j;
    j = temp;
  }
  if(i<1 || i>fNumNus-1 || j<2 || j>fNumNus){
    cout << "ERROR: Theta" << i << j << " not valid for " << fNumNus;
    cout << " neutrinos. Returning zero." << endl;
    return 0;
  }

  return fTheta[i-1][j-1];

}


//......................................................................
///
/// Set the CP phase delta_ij in radians.
///
/// Requires that i+1<j. Will notify you if input is wrong.
/// If i>j, will assume reverse order and swap i and j.
///
/// This will check if value is changing to keep track of whether
/// the hamiltonian needs to be rebuilt.
///
/// @param i,j    - the indices of delta_ij
/// @param delta  - the value of delta_ij
///
void PMNS_Base::SetDelta(int i, int j, double delta)
{

  if(i>j){
    cout << "Warning: First argument should be smaller than second argument" << endl;
    cout << "         Setting reverse order (Delta" << j << i << "). " << endl;
    int temp = i;
    i = j;
    j = temp;
  }
  if(i<1 || i>fNumNus-1 || j<2 || j>fNumNus){
    cout << "ERROR: Delta" << i << j << " not valid for " << fNumNus;
    cout << " neutrinos. Doing nothing." << endl;
    return;
  }
  if(i+1==j){
    cout << "Warning: Rotation " << i << j << " is real. Doing nothing." << endl;
    return;
  }

  // Check if value is actually changing
  fBuiltHms *= (fDelta[i-1][j-1] == delta);

  fDelta[i-1][j-1] = delta;

}

//......................................................................
///
/// Get the CP phase delta_ij in radians.
///
/// Requires that i+1<j. Will notify you if input is wrong.
/// If i>j, will assume reverse order and swap i and j.
///
/// @param i,j    - the indices of delta_ij
///
double PMNS_Base::GetDelta(int i, int j)
{

  if(i>j){
    cout << "Warning: First argument should be smaller than second argument" << endl;
    cout << "         Setting reverse order (Delta" << j << i << "). " << endl;
    int temp = i;
    i = j;
    j = temp;
  }
  if(i<1 || i>fNumNus-1 || j<2 || j>fNumNus){
    cout << "ERROR: Delta" << i << j << " not valid for " << fNumNus;
    cout << " neutrinos. Returning zero." << endl;
    return 0;
  }
  if(i+1==j){
    cout << "Warning: Rotation " << i << j << " is real. Returning zero." << endl;
    return 0;
  }

  return fDelta[i-1][j-1];

}

//......................................................................
///
/// Set the mass-splitting dm_j1 = (m_j^2 - m_1^2) in eV^2
///
/// Requires that j>1. Will notify you if input is wrong.
///
/// This will check if value is changing to keep track of whether
/// the hamiltonian needs to be rebuilt.
///
/// @param j    - the index of dm_j1
/// @param dm   - the value of dm_j1
///
void PMNS_Base::SetDm(int j, double dm)
{

  if(j<2 || j>fNumNus){
    cout << "ERROR: Dm" << j << "1 not valid for " << fNumNus;
    cout << " neutrinos. Doing nothing." << endl;
    return;
  }

  // Check if value is actually changing
  fBuiltHms *= (fDm[j-1] == dm);

  fDm[j-1] = dm;

}

//......................................................................
///
/// Get the mass-splitting dm_j1 = (m_j^2 - m_1^2) in eV^2
///
/// Requires that j>1. Will notify you if input is wrong.
///
/// @param j    - the index of dm_j1
///
double PMNS_Base::GetDm(int j)
{

  if(j<2 || j>fNumNus){
    cout << "ERROR: Dm" << j << "1 not valid for " << fNumNus;
    cout << " neutrinos. Returning zero." << endl;
    return 0;
  }

  return fDm[j-1];

}

//......................................................................
///
/// Get the effective mass-splitting dm_j1 in matter in eV^2
///
/// Requires that j>1. Will notify you if input is wrong.
///
/// @param j    - the index of dm_j1
///
double PMNS_Base::GetDmEff(int j)
{

  if(j<2 || j>fNumNus){
    cout << "ERROR: Dm" << j << "1 not valid for " << fNumNus;
    cout << " neutrinos. Returning zero." << endl;
    return 0;
  }

  // Solve the Hamiltonian to update eigenvalues
  SolveHam();
  
  // Sort eigenvalues in same order as vacuum Dm^2
  vector<int> TrueIdx(fNumNus, 0);
  vector<double> TrueVals(fNumNus, 0);
  vector<int> EffIdx(fNumNus, 0);
  for(int i=0; i<fNumNus; i++){
    TrueIdx[i] = i;
    EffIdx[i] = i;
  }
  sort(TrueIdx.begin(), TrueIdx.end(), IdxCompare(fDm));
  for(int i=0; i<fNumNus; i++) TrueVals[i] = TrueIdx[i];
  sort(TrueIdx.begin(), TrueIdx.end(), IdxCompare(TrueVals));
  sort(EffIdx.begin(), EffIdx.end(), IdxCompare(fEval));

  // Return eigenvalues * 2E
  return (fEval[EffIdx[TrueIdx[j-1]]] - fEval[EffIdx[TrueIdx[0]]]) * fEnergy * 2e9;

}


//......................................................................
///
/// Rotate the Hamiltonian by the angle theta_ij and phase delta_ij.
///
/// The rotations assume all off-diagonal elements with i > j are zero.
/// This is correct if the order of rotations is chosen appropriately
/// and it speeds up computation by skipping null terms
///
/// @param i,j - the indices of the rotation ij
///
void PMNS_Base::RotateH(int i,int j){

  // Do nothing if angle is zero
  if(fTheta[i][j]==0) return;

  double fSinBuffer = sin(fTheta[i][j]);
  double fCosBuffer = cos(fTheta[i][j]);

  double  fHmsBufferD;
  complexD fHmsBufferC;

  // With Delta
  if(i+1<j){
    complexD fExpBuffer = complexD(cos(fDelta[i][j]), -sin(fDelta[i][j]));

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

//......................................................................
///
/// Build Hms = H*2E, where H is the Hamiltonian in vacuum on flavour basis
/// and E is the neutrino energy in eV. Hms is effectively the matrix of
/// masses squared.
///
/// This is a hermitian matrix, so only the
/// upper triangular part needs to be filled
///
/// The construction of the Hamiltonian avoids computing terms that
/// are simply zero. This has a big impact in the computation time.
///
void PMNS_Base::BuildHms()
{

  // Check if anything changed
  if(fBuiltHms) return;
  
  // Tag to recompute eigensystem
  fGotES = false;

  for(int j=0; j<fNumNus; j++){
    // Set mass splitting
    fHms[j][j] = fDm[j];
    // Reset off-diagonal elements
    for(int i=0; i<j; i++){
      fHms[i][j] = 0;
    }
    // Rotate j neutrinos
    for(int i=0; i<j; i++){
      RotateH(i,j);
    }
  }

  ClearCache();

  // Tag as built
  fBuiltHms = true;

}

//.....................................................................
///
/// Propagate the current neutrino state through a given path
/// @param p - A neutrino path segment
///
void PMNS_Base::PropagatePath(NuPath p)
{

  // Set the neutrino path
  SetCurPath(p);

  // Solve for eigensystem
  SolveHam();

  double LengthIneV = kKm2eV * p.length;
  for(int i=0; i<fNumNus; i++){
    double arg = fEval[i] * LengthIneV;
    fPhases[i] = complexD(cos(arg), -sin(arg));
  }
  
  for(int i=0;i<fNumNus;i++){
    fBuffer[i] = 0;
    for(int j=0;j<fNumNus;j++){
      fBuffer[i] += conj(fEvec[j][i]) * fNuState[j];
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

//......................................................................
///
/// Propagate neutrino state through full path
///
void PMNS_Base::Propagate()
{

  for(int i=0; i<int(fNuPaths.size()); i++){

    PropagatePath(fNuPaths[i]);

  }

}

//......................................................................
///
/// Reset the neutrino state back to a pure flavour where it starts
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
/// @param flv - The neutrino starting flavour.
///
void PMNS_Base::ResetToFlavour(int flv)
{
  assert(flv>=0 && flv<fNumNus);
  for (int i=0; i<fNumNus; ++i){
    if (i==flv) fNuState[i] = one;
    else        fNuState[i] = zero;
  }
}

//......................................................................
///
/// Compute oscillation probability of flavour flv from current state
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
/// @param flv - The neutrino final flavour.
///
/// @return Neutrino oscillation probability
///
double PMNS_Base::P(int flv)
{
  assert(flv>=0 && flv<fNumNus);
  return norm(fNuState[flv]);
}

//.....................................................................
///
/// Compute the probability of flvi going to flvf.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
///
/// @return Neutrino oscillation probability
///
double PMNS_Base::Prob(int flvi, int flvf)
{

  ResetToFlavour(flvi);

  return Prob(fNuState, flvf);

}

//.....................................................................
///
/// Compute the probability of nu_in going to flvf.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
/// @param nu_in - The neutrino initial state in flavour basis.
/// @param flvf  - The neutrino final flavour.
///
/// @return Neutrino oscillation probability
///
double PMNS_Base::Prob(vector<complexD> nu_in, int flvf)
{

  assert(nu_in.size() == fNumNus);
  assert(flvf >= 0 && flvf < fNumNus);

  fNuState = nu_in;

  Propagate();

  return P(flvf);

}

//.....................................................................
///
/// Compute the probability of nu_in going to flvf for a given energy in GeV.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
/// @param nu_in - The neutrino initial state in flavour basis.
/// @param flvf  - The neutrino final flavour.
/// @param E     - The neutrino energy in GeV
///
/// @return Neutrino oscillation probability
///
double PMNS_Base::Prob(vector<complexD> nu_in, int flvf, double E)
{

  SetEnergy(E);

  return Prob(nu_in, flvf);

}

//.....................................................................
///
/// Compute the probability of flvi going to flvf for a given energy in GeV.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param E    - The neutrino energy in GeV
///
/// @return Neutrino oscillation probability
///
double PMNS_Base::Prob(int flvi, int flvf, double E)
{

  SetEnergy(E);

  return Prob(flvi, flvf);

}

//.....................................................................
///
/// Compute the probability of nu_in going to flvf for a given
/// energy in GeV and distance in km in a single path.
///
/// If the path sequence is not a single path, a new single path will
/// be created and the previous sequence will be lost.
///
/// Don't use this if you want to propagate over multiple path segments.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
/// @param nu_in - The neutrino initial state in flavour basis.
/// @param flvf  - The neutrino final flavour.
/// @param E     - The neutrino energy in GeV
/// @param L     - The neutrino path length in km
///
/// @return Neutrino oscillation probability
///
double PMNS_Base::Prob(vector<complexD> nu_in, int flvf, double E, double L)
{

  SetEnergy(E);
  SetLength(L);

  return Prob(nu_in, flvf);

}

//.....................................................................
///
/// Compute the probability of flvi going to flvf for a given
/// energy in GeV and distance in km in a single path.
///
/// If the path sequence is not a single path, a new single path will
/// be created and the previous sequence will be lost.
///
/// Don't use this if you want to propagate over multiple path segments.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param E    - The neutrino energy in GeV
/// @param L    - The neutrino path length in km
///
/// @return Neutrino oscillation probability
///
double PMNS_Base::Prob(int flvi, int flvf, double E, double L)
{

  SetEnergy(E);
  SetLength(L);

  return Prob(flvi, flvf);

}

//.....................................................................
///
/// Compute the average probability of flvi going to flvf
/// over a bin of energy E with width dE.
///
/// This gets transformed into L/E, since the oscillation terms
/// have arguments linear in L/E and not E.
///
/// This function works best for single paths.
/// In multiple paths the accuracy may be somewhat worse.
/// If needed, average over smaller energy ranges.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param E    - The neutrino energy in the bin center in GeV
/// @param dE   - The energy bin width in GeV
///
/// @return Average neutrino oscillation probability
///
double PMNS_Base::AvgProb(int flvi, int flvf, double E, double dE)
{

  ResetToFlavour(flvi);

  return AvgProb(fNuState, flvf, E, dE);

}

//.....................................................................
///
/// Compute the average probability of nu_in going to flvf
/// over a bin of energy E with width dE.
///
/// This gets transformed into L/E, since the oscillation terms
/// have arguments linear in L/E and not E.
///
/// This function works best for single paths.
/// In multiple paths the accuracy may be somewhat worse.
/// If needed, average over smaller energy ranges.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
/// @param nu_in - The neutrino initial state in flavour.
/// @param flvf  - The neutrino final flavour.
/// @param E     - The neutrino energy in the bin center in GeV
/// @param dE    - The energy bin width in GeV
///
/// @return Average neutrino oscillation probability
///
double PMNS_Base::AvgProb(vector<complexD> nu_in, int flvf, double E, double dE)
{

  // Do nothing if energy is not positive
  if(E<=0) return 0;

  if(fNuPaths.empty()) return 0;

  // Don't average zero width
  if(dE<=0) return Prob(nu_in, flvf, E);

  // Make sure fPath is set
  // Use average if multiple paths
  SetCurPath(AvgPath(fNuPaths));

  // Define L/E variables
  double LoE = 0;
  double dLoE = 0;

  // Set a minimum energy
  double minE = 0.1 * E;

  // Transform range to L/E
  // Full range if low edge > minE
  if(E-dE/2 > minE){
    LoE = 0.5 * (fPath.length/(E-dE/2) + fPath.length/(E+dE/2));
    dLoE = fPath.length/(E-dE/2) - fPath.length/(E+dE/2);
  }
  // Else start at minE
  else{
    LoE = 0.5 * (fPath.length/minE + fPath.length/(E+dE/2));
    dLoE = fPath.length/minE - fPath.length/(E+dE/2);
  }

  // Compute average in LoE
  return AvgProbLoE(nu_in, flvf, LoE, dLoE);

}

//.....................................................................
///
/// Compute the average probability of flvi going to flvf
/// over a bin of L/E with width dLoE.
///
/// The probabilities are weighted by (L/E)^-2 so that event
/// density is flat in energy. This avoids giving too much
/// weight to low energies. Better approximations would be
/// achieved if we used an interpolated event density.
///
/// This function works best for single paths.
/// In multiple paths the accuracy may be somewhat worse.
/// If needed, average over smaller L/E ranges.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
/// @param flvi - The neutrino starting flavour.
/// @param flvf - The neutrino final flavour.
/// @param LoE  - The neutrino  L/E value in the bin center in km/GeV
/// @param dLoE   - The L/E bin width in km/GeV
///
/// @return Average neutrino oscillation probability
///
double PMNS_Base::AvgProbLoE(int flvi, int flvf, double LoE, double dLoE)
{

  ResetToFlavour(flvi);

  return AvgProbLoE(fNuState, flvf, LoE, dLoE);

}

//.....................................................................
///
/// Compute the average probability of nu_in going to flvf
/// over a bin of L/E with width dLoE.
///
/// The probabilities are weighted by (L/E)^-2 so that event
/// density is flat in energy. This avoids giving too much
/// weight to low energies. Better approximations would be
/// achieved if we used an interpolated event density.
///
/// This function works best for single paths.
/// In multiple paths the accuracy may be somewhat worse.
/// If needed, average over smaller L/E ranges.
///
/// Flavours are:
/// <pre>
///   0 = nue, 1 = numu, 2 = nutau
///   3 = sterile_1, 4 = sterile_2, etc.
/// </pre>
/// @param nu_in - The neutrino intial state in flavour basis.
/// @param flvf  - The neutrino final flavour.
/// @param LoE   - The neutrino  L/E value in the bin center in km/GeV
/// @param dLoE  - The L/E bin width in km/GeV
///
/// @return Average neutrino oscillation probability
///
double PMNS_Base::AvgProbLoE(vector<complexD> nu_in, int flvf, double LoE, double dLoE)
{

  // Do nothing if L/E is not positive
  if(LoE<=0) return 0;

  if(fNuPaths.empty()) return 0;

  // Make sure fPath is set
  // Use average if multiple paths
  SetCurPath(AvgPath(fNuPaths));

  // Set the energy at bin center
  SetEnergy(fPath.length/LoE);

  // Don't average zero width
  if(dLoE<=0) return Prob(nu_in, flvf);

  // Get sample points for this bin
  vector<double> samples = GetSamplePoints(LoE, dLoE);

  // Variables to fill sample
  // probabilities and weights
  double sumw = 0;
  double prob = 0;
  double length = fPath.length;

  // Loop over all sample points
  for(int j=0; j<int(samples.size()); j++){

    // Set (L/E)^-2 weights
    double w = 1./pow(samples[j],2);

    // Add weighted probability
    prob += w * Prob(nu_in, flvf, length / samples[j]);

    // Increment sum of weights
    sumw += w;

  }

  // Return weighted average of probabilities
  return prob / sumw;

}

//.....................................................................
///
/// Compute the sample points for a bin of L/E with width dLoE
///
/// This is used for averaging the probability over a bin of L/E.
/// It should be a private function, but I'm keeping it public for
/// now for debugging purposes. The number of sample points seems
/// too high for most purposes. The number of subdivisions needs
/// to be optimized.
///
/// @param LoE  - The neutrino  L/E value in the bin center in km/GeV
/// @param dLoE   - The L/E bin width in km/GeV
///
vector<double> PMNS_Base::GetSamplePoints(double LoE, double dLoE)
{

  // Solve Hamiltonian to get eigenvalues
  SolveHam();

  // Define conversion factor [km/GeV -> 1/(4 eV^2)]
  const double k1267 = kKm2eV / (4 * kGeV2eV);

  // Get list of all effective Dm^2
  vector<double> effDm;

  for(int i=0; i<fNumNus-1; i++){
    for(int j=i+1; j<fNumNus; j++){
      effDm.push_back( 2 * kGeV2eV * fEnergy * fabs(fEval[j] - fEval[i]) );
    }
  }

  int numDm = effDm.size();

  // Sort the effective Dm^2 list
  sort(effDm.begin(), effDm.end());

  // Set a number of sub-divisions to achieve "good" accuracy
  // This needs to be studied better
  int n_div = ceil( 20 * pow(dLoE/LoE,0.8) );
  //int n_div = 1;

  // A vector to store sample points
  vector<double> allSamples;

  // Loop over sub-divisions
  for(int k=0; k<n_div; k++){

    // Define sub-division center and width
    double bctr = LoE - dLoE/2 + (k+0.5)*dLoE/n_div;
    double bwdt = dLoE/n_div;

    // Make a vector of L/E sample values
    // Initialized in the sub-division center
    vector<double> samples;
    samples.push_back(bctr);

    // Loop over all Dm^2 to average each frequency
    // This will recursively sample points in smaller
    // bins so that all relevant frequencies are used
    for(int i=0; i<numDm; i++){

      // Copy the list of sample L/E values
      vector<double> prev = samples;

      // Redefine bin width to lie within full sub-division
      double Width = 2*min(prev[0] - (bctr - bwdt/2), (bctr + bwdt/2) - prev[0]);

      // Compute oscillation argument sorted from lowest  to highest
      const double arg = k1267 * effDm[i] * Width;

      // Skip small oscillation values.
      // If it's the last one, lower the tolerance
      if(i < numDm-1){
        if(arg<0.9) continue;
      }
      else{
        if(arg<0.1) continue;
      }

      // Reset samples to redefine them
      samples.clear();

      // Loop over previous samples
      for(int j=0; j<int(prev.size()); j++){

        // Compute new sample points around old samples
        // This is based on a oscillatory quadrature rule
        double sample = (1/sqrt(3)) * (Width/2);
        if(arg!=0) sample = acos(sin(arg)/arg)/arg * (Width/2);

        // Add samples above and below center
        samples.push_back(prev[j]-sample);
        samples.push_back(prev[j]+sample);

      }

    }// End of loop over Dm^2

    // Add sub-division samples to the end of allSamples vector
    allSamples.insert(allSamples.end(), samples.begin(), samples.end());

  }// End of loop over sub-divisions

  // Return all sample points
  return allSamples;

}

////////////////////////////////////////////////////////////////////////
