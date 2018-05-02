////////////////////////////////////////////////////////////////////////
//
// Implements an earth model with spherical shells
//
// coelho@lal.in2p3.fr
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <math.h>

#include "PremModel.h"
#include "prem_default.hpp"

using namespace std;

using namespace OscProb;

///
/// Define the tolerance of the detector around
/// a boundary when considering whether to add
/// a new layer for the detector position.
///
const double PremModel::DET_TOL = 0.2; // Tolerance in km

//......................................................................
///
/// Constructor.
///
/// By default this implements the model stored in PremTables/prem_default.txt
/// with the detector at the bottom of the ocean layer (radius = 6368 km).
///
/// @param filename - The txt file containing a table of earth layers
///
PremModel::PremModel(string filename) :
fDetLayer(0), fRemoveSmallPaths(false)
{

  SetDetPos(6368);
  LoadModel(filename);

}

//......................................................................
///
/// Nothing to clean.
///
PremModel::~PremModel(){}

//......................................................................
///
/// Set the detector position in km.
///
/// If the position is within 200m of a layer boundary,
/// the detector is considered to be on the boundary.
/// If not, an extra boundary is inserted in the detector
/// position to distinguish what parts of the earth are above
/// and below the detector.
///
/// This must be done before loading the earth model file.
///
/// @param pos - The radius where the detector is in km
///
void PremModel::SetDetPos(double pos){ fDetPos = pos; }

//......................................................................
///
/// Get the current neutrino path sequence
///
/// The path needs to be filled for a given cosTheta before
/// calling this function.
///
vector<NuPath> PremModel::GetNuPath(){ return fNuPath; }

//......................................................................
///
/// Get the set of earth layers
///
/// This returns the set of PremLayer's for this earth model
/// and detector position.
///
vector<PremLayer> PremModel::GetPremLayers(){ return fPremLayers; }

//......................................................................
///
/// Clear the earth model information.
///
void PremModel::ClearModel()
{

  fDetLayer = 0;
  fPremLayers.clear();

}

//......................................................................
///
/// Add a layer to the earth model.
///
/// @param radius  - The outer radius of the layer in km
/// @param density - The density of the layer in g/cm^3
/// @param zoa     - The effective Z/A value of the layer
/// @param layer   - An index to identify the matter type (e.g. earth inner core)
///
void PremModel::AddLayer(double radius, double density,
                         double zoa,    double layer)
{

  fPremLayers.push_back( PremLayer(radius, density, zoa, layer) );

}

//......................................................................
///
/// Load an earth model from a file.
///
/// By default it loads the model stored in PremTables/prem_default.txt
///
/// The decision of whether to create a layer for the detector position
/// is made in this function, so the detector position must be set before
/// calling this function.
///
/// @param filename - The txt file containing a table of earth layers
///
void PremModel::LoadModel(string filename)
{

  // Clear the current model
  ClearModel();

  // Use default if no file provided
  if(filename == ""){
    filename = PREM_DEFAULT;
  }

  // Open the file
  ifstream fin;
  fin.open(filename.c_str());

  if(!fin){
    cout << "ERROR: File " << filename << " not found!" << endl;
    return;
  }

  // Variables for storing table rows
  float radius, density, zoa, layer;

  // Flag to mark we've passed the detector
  bool crossed_det = false;

  // Keep track of previous radius
  double rprev = 0;

  // Loop over table rows
  while(fin >> radius >> density >> zoa >> layer){

    // Radii must be ordered in model file
    if(radius <= rprev){
      cout << "ERROR: Radii are not sorted in increasing order in the model file" << endl;
      ClearModel();
      return;
    }

    // See if we passed the detector and decide whether
    // to create a special layer for it
    if(radius > fDetPos - DET_TOL && !crossed_det){

      crossed_det = true;

      fDetLayer = fPremLayers.size();

      // If detector is not near boundary, add a special layer
      if(radius > fDetPos + DET_TOL){
        AddLayer(fDetPos, density, zoa, layer);
      }

    }

    // Add this layer to the model
    AddLayer(radius, density, zoa, layer);

  }

}

//......................................................................
///
/// Set the effective Z/A value for all layers of a given type.
///
/// Use this to change the Z/A of indexed layer,
/// e.g. all outer-core layers
///
/// @param layer - The index of the layer type
/// @param zoa   - The effective Z/A value to use
///
void PremModel::SetLayerZoA(int layer, double zoa)
{

  int nlayers = fPremLayers.size();

  // Loop over all layers and change the ones
  // with the given index
  for(int i=0; i<nlayers; i++){

    if(fPremLayers[i].layer != layer) continue;

    fPremLayers[i].zoa = zoa;

  }

}

//......................................................................
///
/// Get the effective Z/A value for all layers of a given type, e.g. all
/// outer-core layers.
///
/// @param layer - The index of the layer type
///
double PremModel::GetLayerZoA(int layer)
{

  int nlayers = fPremLayers.size();

  // This check assumes the layer types are in increasing order
  if(layer > fPremLayers.back().layer){
    cout << "ERROR: Not that many layer types" << endl;
    cout << "Returning 0" << endl;
    return 0;
  }

  for(int i=0; i<nlayers; i++){

    if(fPremLayers[i].layer != layer) continue;

    else return fPremLayers[i].zoa;

  }

  // End of vector reached without finding input type
  cout << "ERROR: layer type not found" << endl;
  cout << "Returning 0" << endl;
  return 0;
}

//......................................................................
///
/// Set the radius of the outermost layer of the model.
///
/// This usually corresponds to the atmosphere and is useful for
/// computing oscillations with variable neutrino production height.
///
/// @param thick - The thickness of the outer layer in km
///
void PremModel::SetTopLayerSize(double thick)
{

  if(thick <= 0){
    cout << "Layer thickness should be positive. Do nothing." << endl;
    return;
  }

  int nlayers = fPremLayers.size();

  if(nlayers < 1){
    cout << "PremModel has no layers. Do nothing." << endl;
    return;
  }

  double bottomRadius;
  if(nlayers>1) bottomRadius = fPremLayers[nlayers-2].radius;
  else          bottomRadius = 0;

  double topRadius = bottomRadius + thick;

  fPremLayers[nlayers-1].radius = topRadius;

}

//......................................................................
///
/// Get the total baseline for a given cosTheta.
///
/// The total baseline contains both from above and below the detector.
///
/// @param cosT - The cosine of the neutrino direction
///
double PremModel::GetTotalL(double cosT)
{

  if(fabs(cosT) > 1) return 0;

  double rAbove = fPremLayers.back().radius;     // Radius above detector
  double rBelow = fPremLayers[fDetLayer].radius; // Radius below detector

  double sinsqrT = 1 - cosT*cosT;

  return -rBelow*cosT + sqrt(rAbove*rAbove - rBelow*rBelow*sinsqrT);

}

//......................................................................
///
/// Get the cosTheta for a given total baseline.
///
/// Given a baseline, find the direction of the neutrino.
/// This could be useful for experiments with fixed baselines for example.
///
/// The baseline must be within the range of possible values in this
/// earth model. Will return vertical neutrinos otherwise.
///
/// @param L - The total baseline of the neutrino
///
double PremModel::GetCosT(double L)
{

  double rAbove = fPremLayers.back().radius;     // Radius above detector
  double rBelow = fPremLayers[fDetLayer].radius; // Radius below detector

  if(L < rAbove - rBelow) return  1;
  if(L > rAbove + rBelow) return -1;

  return (rAbove*rAbove - rBelow*rBelow - L*L) / (2*rBelow*L);

}

//......................................................................
///
/// Add a path segment to the sequence.
///
/// For a given PremLayer, adds a path of a given length in that material
///
/// @param length  - The length of the path segment in km
/// @param pl      - The layer we are crossing
///
void PremModel::AddPath(double length, PremLayer pl)
{

  fNuPath.push_back( NuPath(length, pl.density, pl.zoa, pl.layer) );

}

//......................................................................
///
/// Fill the path sequence in a vector.
///
/// This will start at the upper-most layer and find straight paths
/// to the boundary of the next layer down. When reaching the inner-most
/// layer in the given direction, the paths move back to outer layers
/// until hitting the detector.
///
/// The path sequence is stored as an attribute and can be retrieved with
/// the function GetNuPath.
///
/// @param cosT - The cosine of the neutrino direction
/// @return The number of path segments in the sequence
///
int PremModel::FillPath(double cosT)
{

  // Clear current path sequence
  fNuPath.clear();

  // Do nothing if cosine is unphysical
  if(fabs(cosT) > 1) return 0;

  // Define the minimum path radius
  double minR = fPremLayers[fDetLayer].radius * sqrt(1 - cosT*cosT);

  // Set the top layer index
  int toplayer = fPremLayers.size() - 1;

  // Find the inner-most crossed layer
  int minlayer = 0;
  while(fPremLayers[minlayer].radius < minR) minlayer++;

  // Compute the number of path segments needed
  int nsteps = toplayer - fDetLayer;
  if(cosT < 0) nsteps += 2*(fDetLayer - minlayer) + 1;

  // Start at the top layer and go down
  int layer = toplayer;
  int dl = -1;

  // Loop over all path segments
  for(int i=0; i<nsteps; i++){

    // Get square of the path length between this layer's
    // outer radius and  inner-most radius
    double L1 = pow(fPremLayers[layer].radius,2) - minR*minR;

    // If L1 is negative, outer radius is not crossed.
    // This only happens if detector is at the top layer and the
    // neutrino is coming from above.
    if(L1 < 0) return true;

    // Get square of the path length between this layer's
    // inner radius and  inner-most radius
    double L2 = -minR*minR;
    if(layer>0) L2 += pow(fPremLayers[layer-1].radius,2);

    // If L2 is negative, inner radius is not crossed,
    // so set this as the minimum layer.
    bool ismin = (L2<=0 && cosT < 0);

    // Store the path segment length
    double dL;

    // If it's the minimum layer, connect two outer radius
    // crossing points. If not, compute difference between
    // inner and outer radius crossings.
    if(ismin)      dL = 2 * sqrt(L1);
    else if(L2>=0) dL = sqrt(L1) - sqrt(L2);
    else           dL = sqrt(L1); // This should never actually happen,
                                  // but protect in case L2 is slightly
                                  // negative when it should be zero.
                                  // e.g. arriving at the detector from above.

    // Add this path segment to the sequence
    AddPath(dL, fPremLayers[layer]);

    // If we reached the inner-most layer,
    // start moving up again.
    if(ismin) dl = 1;

    // Move to next layer
    layer += dl;

  }

  // Return the number of path segments
  return fNuPath.size();

}

//......................................................................
///
/// Merge similar paths to reduce number of steps
///
/// This method will merge consecutive paths and take their averages
/// until it finds a large enough gap to start a new merged path.
///
/// The merged paths will be returned, and the original detailed path
/// will not be changed and will stay stored as an attribute.
///
/// @param prec - The precision to merge paths in g/cm^3
/// @return The vector of merged path segments
///
vector<NuPath> PremModel::GetMergedPaths(double prec){

  // The output vector
  vector<NuPath> mergedPath;

  // Start with the first path
  OscProb::NuPath path = fNuPath[0];

  // Track the total length
  double totL = 0;

  // Loop over all paths starting from second
  for(int i=1; i<fNuPath.size(); i++){

    // If this path electron density is beyond the tolerance
    if( fabs(path.density*path.zoa - fNuPath[i].density*fNuPath[i].zoa) > prec*path.zoa ){

      // Add merged path to vector
      mergedPath.push_back(path);

      // Set this path as new merged path
      path = fNuPath[i];

    }
    // If path is within tolerance
    else{

      // Merge the path with current merged path
      path = AvgPath(path, fNuPath[i]);

    }

    // Increment total length
    totL += fNuPath[i].length;

  }// End of loop over paths

  // Add the final merged path to vector
  mergedPath.push_back(path);

  // If tag is true, remove small paths
  if(fRemoveSmallPaths){

    // Start at first path
    int k = 0;

    // While not at the end of vector
    while(k + 1 < mergedPath.size()){

      // If length is less than 1% of total
      if(mergedPath[k].length < 0.01*totL){

        // Merge path with following path
        mergedPath = MergePaths(mergedPath, k, k+1);

      }
      // If path is long enough skip it
      else k++;

    }// End of while loop

  }// End of if statement

  // return the merged vector
  return mergedPath;

}

//......................................................................
///
/// Get the merged average of two paths
///
/// This method will merge two paths and take their average density
/// weighted by Z/A and path length.
///
/// The Z/A of the first path will be kept in the merged path
///
/// @param p1 - The first path to merge
/// @param p2 - The second path to merge
/// @return The merged path
///
NuPath PremModel::AvgPath(NuPath p1, NuPath p2){

  // Start with the first path
  NuPath mergedPath = p1;

  // Add the second length
  mergedPath.length += p2.length;

  // Compute weighted average of density
  mergedPath.density = (p1.density*p1.zoa*p1.length + p2.density*p2.zoa*p2.length) / (p1.zoa*p1.length + p2.zoa*p2.length);

  // return merged path
  return mergedPath;

}

//......................................................................
///
/// Merge two specific paths by their indices in a path vector
///
/// @param inputPath - The original vecotr of paths to merge
/// @param j,k - The indices of the two paths to merge
/// @return The merged vector of paths
///
vector<NuPath> PremModel::MergePaths(vector<NuPath> inputPath, int j, int k){

  // Output vector
  vector<NuPath> mergedPath;

  // Loop over input paths
  for(int i=0; i<inputPath.size(); i++){

    // If first index, merge the paths j and k
    if(i==j) mergedPath.push_back(AvgPath(inputPath[j], inputPath[k]));
    // If not second index add the path as is
    else if(i!=k) mergedPath.push_back(inputPath[i]);

  }// End of loop

  // return merged vector
  return mergedPath;

}

//......................................................................
///
/// Set the boolean to tag whether to remove small paths when merging
///
/// @param rp - Boolean value to set
///
void PremModel::SetRemoveSmallPaths(bool rp){ fRemoveSmallPaths = rp; }
