/*************************************************************************
 * EarthModel.cxx
 * Implements a 3D earth model (binned in depth, latitude, and longitude)
 * for the programming library OscProb
 * Modified from PREMmodel.cxx (by coelho@lal.in2p3.fr)
 * by Rebekah Pestes (last modified: 2/6/23)
 *************************************************************************/

#include <iostream>
#include <fstream>
#include <math.h>

#include "EarthModel.h"
#include "EarthTables/model_defaults.hpp"

using namespace std;

using namespace OscProb;

///
/// Define the tolerance of the detector around
/// a boundary when considering whether to add
/// a new layer for the detector position.
///
const double EarthModel::DET_TOL = 0.2; // Tolerance in km

//......................................................................
///
/// Constructor.
///
/// By default this implements the model stored in EarthTables/prem_default.txt
/// with the detector at the bottom of the ocean layer (radius = 6368 km or depth = 3 km).
///
/// @param filename - The txt file containing a table of earth layers
///
EarthModel::EarthModel(string filename) :
fDetBin(0), fRemoveSmallPaths(false)
{

  SetDetPos(3,0,0);
  LoadModel(filename);

}

//......................................................................
///
/// Nothing to clean.
///
EarthModel::~EarthModel(){}

//......................................................................
///
/// Set the coordinates of the detector:
///   depth in km, latitude in degrees, longitude in degrees
///
/// This must be done before loading the earth model file.
///
/// @param dep - The depth of the detector in km
/// @param lat - The latitude of the detector in degrees
/// @param lon - The longitude of the detector in degrees
///
void EarthModel::SetDetPos(double dep, double lat, double lon)
{
  fDetDepth = dep;
  fDetLat = lat;
  fDetLon = lon;
}

//......................................................................
///
/// Get the current neutrino path sequence
///
/// The path needs to be filled for a given cosTheta before
/// calling this function.
///
vector<NuPath> EarthModel::GetNuPath(){ return fNuPath; }

//......................................................................
///
/// Get the set of Earth bins
///
/// This returns the set of bins for this Earth model.
///
vector<EarthBin> EarthModel::GetEarthBins(){ return fEarthBins; }

//......................................................................
///
/// Clear the earth model information.
///
void EarthModel::ClearModel()
{

  fDetBin = 0;
  fnDepthBins = 0;
  fnLatBins = 0;
  fEarthBins.clear();

}

//......................................................................
///
/// Add a layer to the earth model.
///
/// @param depth_out  - The outer depth of the bin in km
/// @param depth_in  - The inner depth of the bin in km
/// @param latitude - The latitude of the bin center in degrees
/// @param longitude - The longitude of the bin center in degrees
/// @param density - The density of the matter in the bin in g/cm^3
/// @param zoa     - The effective Z/A value of the matter in the bin
/// @param index   - An index to identify the matter type (e.g. earth inner core)
///
void EarthModel::AddBin(double depth_out, double depth_in, double latitude, double longitude,
                         double density, double zoa, double index)
{

  fEarthBins.push_back( EarthBin(depth_out, depth_in, latitude, longitude, density, zoa, index) );

}

//......................................................................
///
/// Load an earth model from a file.
/// ***Need to check what depth will be in the Earth table!!!***
///
/// By default it loads the model stored in EarthTables/prem_default.txt
///
/// The row format for the model table needs to be:
/// 	inner depth (km)	latitude (deg)	longitude (deg)	density (g/cm^3?)	Z/A	layer_index
/// The data in the table must be in order of increasing depth, followed
/// by increasing latitude, and then increasing longitude.
/// Latitude bin widths and longitude bin widths must each be constant.
/// In each bin of one coordinate, the bins for the other coordinates
/// must be the same.
///
/// The bin location of the detector is calculated in this function,
/// so the detector position must be set before calling this function.
///
/// @param filename - The txt file containing a table of earth layers
///
void EarthModel::LoadModel(string filename)
{

  // Clear the current model
  ClearModel();

  // Use default if no file provided
  if(filename == ""){
    filename = PREM3D_DEFAULT;
  }
  cout << "Loading " << filename << " Earth model table..." << endl;

  // Open the file
  ifstream fin;
  fin.open(filename.c_str());

  if(!fin){
    cout << "ERROR: File " << filename << " not found!" << endl;
    return;
  }

  // Variables for storing table rows
  float depth, latitude, longitude, density, zoa, index;

  // Flag to mark we've passed the detector
  bool crossed_det_depth = false;

  // Flag to alert that we need to add a layer for the detector
//  bool add_layer = false;

  // Keep track of previous depth/latitude/longitude
  double depth_prev = 0;
  double lat_prev = -90.0;
  double lon_prev = 0;

  // Hold onto previous different depth, so we have depth minimum for bin
  double depth_min = 0;

  // Loop over table rows
  while(fin >> depth >> latitude >> longitude >> density >> zoa >> index){

    // Move minimum depth (and reset latitude/longitude) if depth has changed from previous one
    if(depth > depth_prev){
      depth_min = depth_prev;
      lat_prev = -90.0;
      lon_prev = 0.0;
      fnDepthBins++;
    } else if(depth < depth_prev){ // Depths must be ordered in model file
      cout << "ERROR: Depths are not sorted in increasing order in the model file" << endl;
      ClearModel();
      return;
    } else if(latitude > lat_prev){ // Reset longitude if latitude has changed from previous one
      lon_prev = 0.0;
      fnLatBins++;
    } else if(latitude < lat_prev){ // Latitudes must be ordered in model file (after depths)
      cout << "ERROR: Latitudes are not sorted in increasing order (after depths) in the model file" << endl;
      ClearModel();
      return;
    } else if(longitude < lon_prev){ // Longitudes must be ordered in model file (after latitudes)
      cout << "ERROR: Longitudes are not sorted in increasing order (after latitudes) in the model file" << endl;
      ClearModel();
      return;
    }

    // See if we passed the detector in depth
    if(!crossed_det_depth && depth > fDetDepth - DET_TOL){
        crossed_det_depth = true;
        fDetBin = fEarthBins.size(); //Acutally, first bin at detector depth
        // If detector is not near boundary, flag to add layer later
  //      if(depth > fDetPos + DET_TOL){
  //        add_layer = true;
  //      }

    }

    // Add this bin to the model
    AddBin(depth_min, depth, latitude, longitude, density, zoa, index);

    // Set previous coordinates for next bin
    depth_prev = depth;
    lat_prev = latitude;
    lon_prev = longitude;

  }

  //Set Earth Radius
  fEarthRadius = depth_prev;

  //Modify number of latitude bins to how many there actually are (instead of how many times we changed the latitude bin without changing the depth bin)
  fnLatBins = fnLatBins/fnDepthBins + 1;
  double LatBinWidth = 180/fnLatBins;

  //Find detector bin index
  int nLonBins = fEarthBins.size()/(fnLatBins*fnDepthBins);
  fDetDepthBin = fDetBin/(fnLatBins*nLonBins); //Detector Depth Bin
  double LonBinWidth = 360/nLonBins;
  fDetBin += floor((fDetLat+90)/LatBinWidth)*nLonBins + floor(fDetLon/LonBinWidth);  //Detector Bin
  if(abs(fEarthBins[fDetBin].latitude - fDetLat) > LatBinWidth) {
    cout << "ERROR: Detector latitude (" << fDetLat << ") is not within Earth bin " << fDetBin << " (bin center latitude = " << fEarthBins[fDetBin].latitude << ", latitude bin width = " << LatBinWidth << "). " << endl;
    ClearModel();
    return;
  } else if(abs(fEarthBins[fDetBin].longitude - fDetLon) > LonBinWidth) {
    cout << "ERROR: Detector longitude (" << fDetLon << ") is not within Earth bin " << fDetBin << " (bin center longitude = " << fEarthBins[fDetBin].longitude << ", longitude bin width = " << LonBinWidth << "). " << endl;
    ClearModel();
    return;
  }

  //Error if first bin center deviates from about half of the calculated bin widths
  if(abs(LatBinWidth-2*(fEarthBins[0].latitude+90))>LatBinWidth*0.1 || abs(LonBinWidth-2*fEarthBins[0].longitude)>LonBinWidth*0.1) { //Warn if first bin center isn't about half of the bin width
    cout << "ERROR: Latitude/longitude of 1st bin center (" << fEarthBins[0].latitude << ", "<< fEarthBins[0].longitude << ") is not in the middle of the first bin, whose sizes should be " << LatBinWidth << " deg and " << LonBinWidth << " deg, as calculated from having " << fnLatBins << " latitude bins and " << nLonBins << " longitude bins." << endl;
    ClearModel();
    return;
  }

  cout << "\t...done (" << fEarthBins.size() << " bins: " << fnDepthBins << " depth, " << fnLatBins << " latitude, and " << nLonBins << " longitude)." << endl;
  cout << "Earth Radius using: " << fEarthRadius << " km" << endl;
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
void EarthModel::SetLayerZoA(int layer, double zoa)
{

  int nbins = fEarthBins.size();

  // Loop over all bins and change the ones
  // with the given index
  for(int i=0; i<nbins; i++){

    if(fEarthBins[i].index != layer) continue;

    fEarthBins[i].zoa = zoa;

  }

}

//......................................................................
///
/// Get the effective Z/A value for all layers of a given type, e.g. all
/// outer-core layers.
/// (Assumes that all bins of given layer type have same Z/A value)
///
/// @param layer - The index of the layer type
///
double EarthModel::GetLayerZoA(int layer)
{

  int nbins = fEarthBins.size();

  // This check assumes the layer types are in increasing order
  if(layer > fEarthBins.back().index){
    cout << "ERROR: Not that many layer types" << endl;
    cout << "Returning 0" << endl;
    return 0;
  }

  for(int i=0; i<nbins; i++){

    if(fEarthBins[i].index != layer) continue;

    else return fEarthBins[i].zoa;

  }

  // End of vector reached without finding input type
  cout << "ERROR: layer type not found" << endl;
  cout << "Returning 0" << endl;
  return 0;
}

//......................................................................
///
/// Get the total baseline for a given cosTheta.
///
/// @param cosT - The cosine of the neutrino direction
///
double EarthModel::GetTotalL(double cosT)
{

  if(fabs(cosT) > 1) return 0;

  double rDet = fEarthRadius - fDetDepth;

  double sinsqrT = 1 - cosT*cosT;

  return -rDet*cosT + sqrt(fEarthRadius*fEarthRadius - rDet*rDet*sinsqrT);

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
double EarthModel::GetCosT(double L)
{

  double rDet = fEarthRadius - fDetDepth;

  if(L < fDetDepth) return  1;
  if(L > fEarthRadius + rDet) return -1;

  return (fEarthRadius*fEarthRadius - rDet*rDet - L*L) / (2*rDet*L);

}

//......................................................................
///
/// Add a path segment to the sequence.
///
/// For a given EarthBin, adds a path of a given length in that material
///
/// @param length  - The length of the path segment in km
/// @param bin     - The bin we are crossing
///
void EarthModel::AddPath(double length, EarthBin bin)
{

  fNuPath.push_back( NuPath(length, bin.density, bin.zoa, bin.index) );

}

//......................................................................
///
/// Fill the path sequence in a vector.
/// /***NEEDS UPDATING!!!***/
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
int EarthModel::FillPath(double cosT)
{

  // Clear current path sequence
  fNuPath.clear();

  // Do nothing if cosine is unphysical
  if(fabs(cosT) > 1) return 0;

  // Define the maximum depth along the path (detector depth, if cosT is non-negative)
  double maxDepth = fDetDepth;
  double rDet = fEarthRadius - fDetDepth;
  double sinSqT = 1 - cosT*cosT;
  if(cosT < 0) {maxDepth = fEarthRadius - rDet * sqrt(sinSqT);}

  // Get initial Detector Distance
  double rDetCosT = rDet*cosT;
  double rDetSqSinSqT = rDet*rDet*sinSqT;
  double prev_DetDistance = -rDetCosT + sqrt(fEarthRadius*fEarthRadius - rDetSqSinSqT);

  // Get initial Earth bin and # of bins per depth layer
  int index = 0;
  int binsPerDepth = fEarthBins.size()/fnDepthBins;
  cout << "\t..." << binsPerDepth << " Earth bins per depth layer..." << endl;

  // Start at the top layer (when neutrino enters Earth) and go down to minimum, looping over depth bins
  int dBin = 0;
  for(dBin = 0; dBin<fnDepthBins; dBin++) {
    //Check if max depth is contained in depth bin
    double d_in = fEarthBins[index].depth_in; //inner-most depth
    if(d_in > maxDepth) break;
    double r_in = fEarthRadius - d_in; //inner-most radius

    //Distance from detector at inner-most radius
    double DetDistance = -rDetCosT + sqrt(r_in*r_in - rDetSqSinSqT);

    //Check for latitude/longitude bin crossings
    /***WRITE CODE!!!***/

    //Add Segment in depth bin to path
    AddPath(prev_DetDistance - DetDistance, fEarthBins[index]);

    //Reset variables for "previous" depth bin
    prev_DetDistance = DetDistance;
    index += binsPerDepth;
  }

  //Do minimum depth bin and then go up to the detector
  for(int i=dBin; i>fDetDepthBin; i--) {
    double r_out = fEarthRadius - fEarthBins[index].depth_out; //outer-most radius

    //Distance from detector at outer-most radius
    double DetDistance = -rDetCosT - sqrt(r_out*r_out - rDetSqSinSqT);

    //Check for latitude/longitude bin crossings
    /***WRITE CODE!!!***/

    //Add Segment in depth bin to path
    AddPath(prev_DetDistance - DetDistance, fEarthBins[index]);

    //Reset variables for "previous" depth bin
    prev_DetDistance = DetDistance;
    index -= binsPerDepth;
  }

  //Do path to detector
    //Check for latitude/longitude bin crossings
    /***WRITE CODE!!!***/
  AddPath(prev_DetDistance, fEarthBins[index]);

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
vector<NuPath> EarthModel::GetMergedPaths(double prec){

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

    } else { // If path is within tolerance

      // Merge the path with current merged path
      path = OscProb::AvgPath(path, fNuPath[i]);

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
/// Set the boolean to tag whether to remove small paths when merging
/// Small is defined as <1% of the total baseline
///
/// @param rp - Boolean value to set
///
void EarthModel::SetRemoveSmallPaths(bool rp){ fRemoveSmallPaths = rp; }
