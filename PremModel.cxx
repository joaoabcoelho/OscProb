////////////////////////////////////////////////////////////////////////
//
// Implements an earth model with spherical shells
//
// jcoelho\@apc.in2p3.fr
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
fDetLayer(0)
{

  SetRemoveSmallPaths(false);
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
/// Set the coordinates of the detector:
///   radius in km, latitude in degrees, longitude in degrees
///
/// @param rad - The distance from the detector to the Earth's center in km 
/// @param lat - The latitude of the detector in deg N (between -90 and 90)
/// @param lon - The longitude of the detector in deg E (between 0 and 360)
///
void PremModel::SetDetPos(double rad, double lat, double lon)
{

  SetDetectorCoordinates(rad, lat, lon);

}

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
    if(radius > fDetRadius - DET_TOL && !crossed_det){

      crossed_det = true;

      fDetLayer = fPremLayers.size();

      // If detector is not near boundary, add a special layer
      if(radius > fDetRadius + DET_TOL){
        AddLayer(fDetRadius, density, zoa, layer);
      } else if(radius != fDetRadius) {
        //update detector radius
        cout << "WARNING: Adjusting detector radius from " << fDetRadius << " km to " << radius << " km." << endl;
        fDetRadius = radius;
      }

    }

    // Add this layer to the model
    AddLayer(radius, density, zoa, layer);

  }

  //Set the maximum radius in the model
  fRadiusMax = fPremLayers.back().radius;

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
/// @return Z/A corresponding to layer
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

  //Update max radius
  fRadiusMax = fPremLayers.back().radius;

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

  AddPathSegment(length, pl.density, pl.zoa, pl.layer);

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
/// @param cosT - The cosine of the zenith angle of the neutrino direction
/// @param phi - The azimuthal angle of the neutrino direction (default = 0; not used)
/// @return The number of path segments in the sequence
///
int PremModel::FillPath(double cosT, double phi)
{

  // Clear current path sequence
  fNuPath.clear();

  // Do nothing if cosine is unphysical
  if(fabs(cosT) > 1) return 0;

  // Define the minimum path radius
  double minR = fDetRadius * sqrt(1 - cosT*cosT);
  double minRsq = minR*minR;

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
    double L1 = pow(fPremLayers[layer].radius,2) - minRsq;

    // If L1 is negative, outer radius is not crossed.
    // This only happens if detector is at the top layer and the
    // neutrino is coming from above.
    if(L1 < 0) return true;

    // Get square of the path length between this layer's
    // inner radius and  inner-most radius
    double L2 = -minRsq;
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