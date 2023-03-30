/*************************************************************************
 * EarthModel.cxx
 * Implements a 3D earth model (binned in depth, latitude, and longitude)
 * for the programming library OscProb
 * Modified from PREMmodel.cxx (by coelho@lal.in2p3.fr)
 * by Rebekah Pestes (last modified: 2/6/23)
 *************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include "EarthModel.h"
#include "EarthTables/model_defaults.hpp"

using namespace std;

using namespace OscProb;

//......................................................................
///
/// Constructor.
///
/// By default this implements the model stored in EarthTables/prem_default.txt
/// with the detector 3 km below sea level, having a latitude of 0 deg N
/// and a longitude of 0 deg E.
///
/// @param filename - The txt file containing a table of earth layers
///
EarthModel::EarthModel(string filename) :
fRemoveSmallPaths(false), fLonError(0)
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
/// Assumes the Earth's sea level is 6371 km from its center.
/// Note: Not only is a detector in the center of the Earth impractical,
/// it messes up the equation for the distance from the detector as a
/// function of latitude in the neutrino path calculation, so this
/// restricts the detector depth to being strictly less than 6371 km.
///
/// @param dep - The depth of the detector below sea level in km
/// @param lat - The latitude of the detector in deg N (between -90 and 90)
/// @param lon - The longitude of the detector in deg E (between 0 and 360)
///
void EarthModel::SetDetPos(double dep, double lat, double lon)
{
  fEarthRadius = 6371.0; //km
  if(dep >= fEarthRadius) {
    cout << "ERROR: Really??? Are you sure you want the detector in the Earth's core? Don't... just don't. (Inputted depth: " << dep << " km; Earth's radius: " << fEarthRadius << " km)" << endl;
    return;
  }
  fDetRadius = fEarthRadius-dep;

  //Force latitude to be between -90 and 90 deg
  lat -= floor((lat+90.0)/360.0)*360.0;
  if (lat > 90)
    lat = 180.0-lat;
  fDetLat = lat/180.0*M_PI; //convert to radians

  //Force longitude to be between 0 and 360 deg
  lon -= floor(lon/360.0)*360.0;
  fDetLon = lon/180.0*M_PI; //convert to radians
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

  fnDepthBins = 0;
  fnLonBins = 0;
  fnLatBins = 0;
  fEarthBins.clear();

}

//......................................................................
///
/// Add a bin to the earth model.
///
/// @param radius_out  - The outer depth of the bin in km
/// @param radius_in  - The inner depth of the bin in km
/// @param latitude - The latitude of the bin center in degrees
/// @param longitude - The longitude of the bin center in degrees
/// @param density - The density of the matter in the bin in g/cm^3
/// @param zoa     - The effective Z/A value of the matter in the bin
/// @param index   - An index to identify the matter type (e.g. earth inner core)
///
void EarthModel::AddBin(double radius_out, double radius_in, double latitude, double longitude,
                         double density, double zoa, double index)
{

  fEarthBins.push_back( EarthBin(radius_out, radius_in, latitude, longitude, density, zoa, index) );

}

//......................................................................
///
/// Load an earth model from a file.
///
/// By default it loads the model stored in EarthTables/prem_default.txt
///
/// The row format for the model table needs to be:
/// 	longitude (deg)	latitude (deg)	outer depth (km)	density (g/cm^3?)	Z/A	layer_index
/// The data in the table must be in order of decreasing depth, followed
/// by increasing longitude, and then increasing latitude.
/// Latitude bin widths and longitude bin widths must each be constant.
/// In each bin of one coordinate, the bins centers for the other
/// coordinates must be the same.
/// This assumes that the center of the Earth is at a depth of 6371 km.
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
  cout << "Loading Earth model table from " << filename << "..." << endl;

  // Open the file
  ifstream fin;
  fin.open(filename.c_str());
  if(!fin){
    cout << "ERROR: File " << filename << " not found!" << endl;
    return;
  }

  // Variables for storing table rows
  float depth, latitude, longitude, density, zoa, index;

  // Keep track of previous depth/latitude/longitude
  double radius_prev = 0;
  double lat_prev = -90.0;
  double lon_prev = 0;

  // Hold onto previous different depth, so we have depth maximum for bin
  double radius_min = radius_prev;

  // Loop over table rows
  while(fin >> longitude  >> latitude >> depth >> density >> zoa >> index){
    double outer_radius = fEarthRadius - depth;

    // Move minimum depth (and reset latitude/longitude) if depth has changed from previous one
    if(outer_radius > radius_prev){
      radius_min = radius_prev;
      lat_prev = -90.0;
      lon_prev = 0.0;
      fnDepthBins++;
    } else if(outer_radius < radius_prev){ // Depths must be ordered in model file
      cout << "ERROR: Depths are not sorted in decreasing order in the model file (or depth is greater than " << fEarthRadius << " km)." << endl;
      ClearModel();
      return;
    } else if(longitude > lon_prev){ // Reset latitude if longitude has changed from previous one
      lat_prev = -90.0;
      fnLonBins++;
    } else if(longitude < lon_prev){ // Longitudes must be ordered in model file (after depths)
      cout << "ERROR: Longitudes are not sorted in increasing order (after depths) in the model file." << endl;
      ClearModel();
      return;
    } else if(latitude < lat_prev){ // Latitudes must be ordered in model file (after longitudes)
      cout << "ERROR: Latitudes are not sorted in increasing order (after longitudes) in the model file." << endl;
      ClearModel();
      return;
    }

    // Add this bin to the model
    AddBin(outer_radius, radius_min, latitude, longitude, density, zoa, index);

    // Set previous coordinates for next bin
    radius_prev = outer_radius;
    lat_prev = latitude;
    lon_prev = longitude;

  }

  //Set Outermost Radius
  fRadiusMax = radius_prev;

  //Modify number of longitude bins to how many there actually are (instead of how many times we changed the longitude bin without changing the depth bin = nLonBins-1 times per depth bin)
  fnLonBins = fnLonBins/fnDepthBins + 1;
  fInvLonBinWidth = 0.5*fnLonBins/M_PI; //in radians^(-1)
  fHalfLonBinWidth = M_PI/fnLonBins; //in radians

  //Calculate number of latitude bins
  fnLatBins = fEarthBins.size()/(fnLonBins*fnDepthBins);
  fInvLatBinWidth = fnLatBins/M_PI; //in radians^(-1)
  fHalfLatBinWidth = 0.5*M_PI/fnLatBins; //in radians

  //Error if first bin center deviates from about half of the calculated bin widths
  if(abs(fHalfLatBinWidth-(fEarthBins[0].latitude+0.5*M_PI))>fHalfLatBinWidth*0.1 || abs(fHalfLonBinWidth-fEarthBins[0].longitude)>fHalfLonBinWidth*0.1) { //Warn if first bin center isn't about half of the bin width
    cout << "ERROR: Latitude/longitude of 1st bin center (" << fEarthBins[0].latitude << ", "<< fEarthBins[0].longitude << ") is not in the middle of the first bin, whose sizes should be " << 2*fHalfLatBinWidth << " rad and " << 2*fHalfLonBinWidth << " rad, as calculated from having " << fnLatBins << " latitude bins and " << fnLonBins << " longitude bins." << endl;
    ClearModel();
    return;
  }

  cout << "\t...done (" << fEarthBins.size() << " bins: " << fnDepthBins << " depth, " << fnLatBins << " latitude, and " << fnLonBins << " longitude)." << endl;
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

  double sinsqrT = 1 - cosT*cosT;

  return -fDetRadius*cosT + sqrt(fRadiusMax*fRadiusMax - fDetRadius*fDetRadius*sinsqrT);

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

  if(L < fRadiusMax - fDetRadius) return  1;
  if(L > fRadiusMax + fDetRadius) return -1;

  return (fRadiusMax*fRadiusMax - fDetRadius*fDetRadius - L*L) / (2*fDetRadius*L);

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
/// Calculate the distance to the detector along a neutrino trajectory,
/// specified by cosT and the azimuthal angle, at the edge of
/// the current latitude bin in the direction specified by dir.
///
/// All the variables not listed below are the same as defined inside the
/// FillPath function.  In an attempt to make this calculation more robust,
/// this function turns the latitude bin iteration around when sin^2(lat)
/// gets within 10^(-14) of its extreme limit.
///
/// @param cur_index  - Index of current Earth bin
/// @param dir        - "+1" ("-1") if latitude is increasing (decreasing)
/// @param dLat       - Half of the latitude bin width
/// @param max_reached - Indicates whether the maximum lat has been reached (1), or not (0)
///
double EarthModel::DetDistForNextLatBin(int cur_index, int &dir, double dLat, double maxSinSqLat, int &max_reached, double beta, double gamma, double gammaSq, double rDetGammaSinDetLat, double rDetCosAcosDetLat, double rDetSinDetLat, double rDetCosT, double rDetSinT)
{

  //Find next lat bin boundary
  double NewLat = fEarthBins[cur_index].latitude + dir*dLat;
  //Turn around if lat reaches +/- 90 degrees
  if (abs(NewLat)-0.5*M_PI > -1.0e-14) {
    if (max_reached > 0)
      return -1; //Only turn around once
    max_reached = 1;
    dir = -dir;
    NewLat = fEarthBins[cur_index].latitude + dir*dLat;
    if (abs(NewLat)-0.5*M_PI > -1.0e-14)
      return -1; //This would happen if there's only one latitude bin...
  }

  //Calculate parts of answer
  double sinNewLat = sin(NewLat);
  double NsinSqNewLat = -sinNewLat*sinNewLat;
  double radicand = maxSinSqLat+NsinSqNewLat;
  //Turn around if the radicand is 0 (close enough) or less
  if (radicand < 1.0e-14) {
    if (max_reached > 0)
      return -1; //Only turn around once
    max_reached = 1;
    dir = -dir;
    NewLat = fEarthBins[cur_index].latitude + dir*dLat;
    sinNewLat = sin(NewLat);
    NsinSqNewLat = -sinNewLat*sinNewLat;
    radicand = maxSinSqLat+NsinSqNewLat;
    if (radicand < 1.0e-14)
      return -1; //The rest of the track is contained in this lon bin
  }
  double denominator = gammaSq+NsinSqNewLat;

  //Calculate distance from the detector at next lat bin boundary
  double answer = 0;
  if (denominator == 0) {
     answer = 0.5*(rDetCosAcosDetLat/beta-rDetSinDetLat/gamma);
  } else {
    answer = -(rDetCosT*NsinSqNewLat+rDetGammaSinDetLat+dir*rDetSinT*sinNewLat*sqrt(radicand))/denominator;
  }
  //Indicates that going to the edge of the bin jumped to an invalid part of the function
  if (max_reached > 0 && answer >= rDetCosAcosDetLat/beta) {
    return -1; //The rest of the track is contained in this lon bin
  }

  return answer;

}

//......................................................................
///
/// Calculate the distance to the detector at the edge of the current
/// longitude bin along a neutrino trajectory specified by cosT and the
/// azimuthal angle, incrementing lon_bin in the process. The maximum/
/// minimum longitude are expressed around the 0-2pi boundary. So, if
/// the track spans that boundary, the actual allowed range for lon is
/// [min_lon,2pi)U[0,max_lon].
///
/// All the variables not listed below are the same as defined inside the
/// FillPath function.
///
/// @param prev_lon - Value at center of initial longitude bin
/// @param dLon     - Difference between prev_lon and edge of next lon bin
/// @param new_bin  - Index of longitude bin to be incremented
/// @param min_lon - Minimum longitude allowed in neutrino's trajectory
/// @param max_lon - Maximum longiutde allowed in neutrino's trajectory
///
double EarthModel::DetDistForNextLonBin(double prev_lon, double dLon, int &lon_bin, double min_lon, double max_lon, double sinDetLon, double cosDetLon, double alpha, double sinTsinAsinDetLon, double sinTsinAcosDetLon, double rDetCosDetLat)
{

  //Increment to next bin in direction dLon
  double lon = prev_lon + dLon;
  if (dLon > 0) {
    lon_bin++;
    if(lon_bin >= fnLonBins) {
      //cross from 2pi to 0
      lon_bin = 0;
      lon = 0.0;
    }
  } else {
    lon_bin--;
    if(lon_bin < 0) {
      //cross from 0 to 2pi
      lon_bin = fnLonBins - 1;
      lon = 0.0;
    }
  }

  //Check if lon is outside of range covered by neutrino's trajectory
  if (max_lon < min_lon) {
    if (min_lon > lon && max_lon < lon)
      return -1;
  } else {
    if (min_lon > lon || max_lon < lon)
      return -1;
  } 

  //Parts of the answer
  double cosNewLon = cos(lon);
  double sinNewLon = sin(lon);
  double expression = sinDetLon*cosNewLon-sinNewLon*cosDetLon;
  double denominator = alpha*expression+sinTsinAcosDetLon*cosNewLon+sinTsinAsinDetLon*sinNewLon;

  //Calculate & return detector distance
  if (denominator == 0) { //Is there ever a case when this is true???
    //Error message will be printed at the end of the FillPath function
    fErrorMessage_LonInfo.assign("Longitude along neutrino trajectory = "+std::to_string(lon));
    fLonError = -1;
    return -1;
  }
  return rDetCosDetLat*expression/denominator;

}

//......................................................................
///
/// ...
///
/// All the variables not listed below are the same as defined inside the
/// FillPath function.
///
/// @param detDist_nextDbin  - ...
///
void EarthModel::RecordLatLonBinCrossings(double detDist_nextDbin, double &DetDist, int &index, int &latBin, int &nextLatBin, double &detDist_nextLatBin, int &sign, double dLat, double maxSinSqLat, int &maxlatreached, int &lonBin, int &nextLonBin, double &detDist_nextLonBin, double dLon, double min_lon, double max_lon, double alpha, double beta, double gamma, double gammaSq, double rDetGammaSinDetLat, double rDetCosAcosDetLat, double rDetSinDetLat, double rDetCosDetLat, double rDetCosT, double rDetSinT, double sinTsinAsinDetLon, double sinTsinAcosDetLon, double sinDetLon, double cosDetLon)
{

  while(detDist_nextDbin < detDist_nextLatBin || detDist_nextDbin < detDist_nextLonBin) {

    if (detDist_nextLatBin > detDist_nextLonBin) {
//    cout << "\tChanging lat bins (to " << nextLatBin << ") at x=" << detDist_nextLatBin << "km" << endl;

      //Add segment to next lat bin to path
      AddPath(DetDist - detDist_nextLatBin, fEarthBins[index]);

      //Move lat bins
      DetDist = detDist_nextLatBin;
      index += nextLatBin - latBin;
      latBin = nextLatBin;

      //Find next lat bin
      if (dLat == 0) {
        detDist_nextLatBin = -1;
      } else {
        detDist_nextLatBin = DetDistForNextLatBin(index, sign, dLat, maxSinSqLat, maxlatreached, beta, gamma, gammaSq, rDetGammaSinDetLat, rDetCosAcosDetLat, rDetSinDetLat, rDetCosT, rDetSinT);
        nextLatBin = latBin + sign;
      }
    } else {
//    cout << "\tChanging lon bins (to " << nextLonBin << ") at x=" << detDist_nextLonBin << "km" << endl;

      //Add segment to next lon bin to path
      AddPath(DetDist - detDist_nextLonBin, fEarthBins[index]);

      //Move lon bins
      DetDist = detDist_nextLonBin;
      index += (nextLonBin - lonBin)*fnLatBins;
      lonBin = nextLonBin;

      //Find next lon bin
      if (dLon == 0) {
        detDist_nextLonBin = -1;
      } else {
        detDist_nextLonBin = DetDistForNextLonBin(fEarthBins[index].longitude, dLon, nextLonBin, min_lon, max_lon, sinDetLon, cosDetLon, alpha, sinTsinAsinDetLon, sinTsinAcosDetLon, rDetCosDetLat);
      }
    }

  }

}

//......................................................................
///
/// Fill the path sequence in a vector.
///
/// This will start at the upper-most layer (where the neutrino's path
/// started) and find path lengths to the next bin crossed. When reaching
/// the inner-most depth in the given direction, the paths move back to
/// outer layers until hitting the detector.  The neutrino direction is
/// specified as the zenith angle and azimuthal angle of the vector pointing
/// from the detector to where the neutrino entered the Earth.
///
/// The path sequence is stored as an attribute and can be retrieved with
/// the function GetNuPath.
///
/// @param cosT - The cosine of the zenith angle of the neutrino direction
/// @param phi - The azimuthal angle (in degrees from North) of the neutrino direction
/// @return The number of path segments in the sequence
///
int EarthModel::FillPath(double cosT, double phi)
{

  // Clear current path sequence
  fNuPath.clear();

  // Do nothing if cosine is unphysical
  if(fabs(cosT) > 1) return 0;

  // Calculate sines/cosines of angles
  double Az = phi/180.0*M_PI; //convert from degrees to radians
  double cosA = cos(Az);
  double sinA = sin(Az);
  double cosDetLat = cos(fDetLat); //could calculate earlier
  double sinDetLat = sin(fDetLat); //could calculate earlier
  double cosDetLon = cos(fDetLon); //could calculate earlier
  double sinDetLon = sin(fDetLon); //could calculate earlier

  // Convert cosT to sinT
  double sinSqT = 1 - cosT*cosT;
  double sinT = sqrt(sinSqT);

  // Handy combinations of variables
  double rDetSinT = fDetRadius*sinT;
  double rDetCosT = fDetRadius*cosT;
  double rDetSqSinSqT = rDetSinT*rDetSinT;
  double rDetCosDetLat = fDetRadius*cosDetLat; //could calculate earlier
  double rDetSinDetLat = fDetRadius*sinDetLat; //could calculate earlier
  double cosAcosDetLat = cosA*cosDetLat;
  double rDetCosAcosDetLat = cosA*rDetCosDetLat;
  double cosTcosDetLat = cosT*cosDetLat;
  double sinTsinDetLat = sinT*sinDetLat;
  double sinTcosAsinDetLat = cosA*sinTsinDetLat;
  double cosTcosAcosDetLat = cosA*cosTcosDetLat;
  double sinTsinA = sinT*sinA;
  double sinTsinAcosDetLon = sinTsinA*cosDetLon;
  double sinTsinAsinDetLon = sinTsinA*sinDetLon;
  double alpha = sinTcosAsinDetLat-cosTcosDetLat; //better name?
  double beta = sinTsinDetLat-cosTcosAcosDetLat; //better name?
  double gamma = sinT*cosAcosDetLat+cosT*sinDetLat;
  double gammaSq = gamma*gamma;
  double rDetGammaSinDetLat = rDetSinDetLat*gamma;
  double maxSinSqLat = 1-pow(sinA*cosDetLat,2);
  double distEdgeToMin = sqrt(fRadiusMax*fRadiusMax-rDetSqSinSqT); //distance from point where neutrino enters to minimum depth (neglecting stopping at the detector)
  double baseline = distEdgeToMin-rDetCosT; //total baseline divided by the total Earth model radius

  // Define the maximum depth (minimum radius) along the path (detector depth/radius, if cosT is non-negative)
  double minRadius = fDetRadius;
  if(cosT < 0) {minRadius = rDetSinT;}

  // Starting latitude
  double init_lat = asin((rDetSinT*beta+gamma*distEdgeToMin)/fRadiusMax); //in radians

  // Starting longitude
  double init_lon = atan2(sinDetLon*rDetCosDetLat-baseline*(alpha*sinDetLon+sinTsinAcosDetLon), cosDetLon*rDetCosDetLat+baseline*(sinTsinAsinDetLon-alpha*cosDetLon)); //in radians
  double max_lon = init_lon+M_PI;
  if (init_lon < 0) {
    init_lon += 2*M_PI;
  }
  double min_lon = init_lon;

  // Get initial Earth bin and # of bins per depth layer
  int binsPerDepth = fEarthBins.size()/fnDepthBins;
  int init_lonBin = floor(init_lon*fInvLonBinWidth);
  int init_latBin = floor((init_lat+M_PI/2)*fInvLatBinWidth);
  int index = fEarthBins.size()-binsPerDepth + init_lonBin*fnLatBins + init_latBin;
//  cout << "\t...Initial (Lat,Long) = (" << init_lat << "," << init_lon << ") => Initial bin: " << index << "..." << endl;

  /* Find detector distances for 1st latitude/longitude bin changes */
  double detDist_nextLatBin;
  int nextLatBin = init_latBin;
  double detDist_nextLonBin;
  int nextLonBin = init_lonBin;
  double dLat = fHalfLatBinWidth; //change in lat from bin center to next bin (0 if no more than 1 lat bin change)
  //Condition calculation for starting sign in x(lat) equation
  int sign = 1; //+ => lat incr., - => lat decr. (with decr. x)
  int maxlatreached = 0; //this changes to 1 if past the lat func transition
  if (beta < 0)
    sign = -1;
  if (beta == 0) { if (cosA > 0) sign = -1;
  } else if (baseline < rDetCosAcosDetLat/beta) {
    sign = -sign;
    maxlatreached++;
  }
  //Condition calculation for incr./decr. Longitude (with decr. x)
  double dLon = fHalfLonBinWidth; //change in lon from bin center to next bin (includes direction; 0 if no more than 1 lon bin change)
  if (sinA < 0) {
    dLon = -dLon;
    min_lon = max_lon;
    max_lon = init_lon;
  }
  //Remainder of detDist calculations
  if (sinT == 0) {
    //both latitude and logitude flip at center of the Earth
    detDist_nextLonBin = -rDetCosT; //formula is divided by cosT, but cosT = +-1
    nextLonBin = floor(fDetLon*fInvLonBinWidth);
    dLon = 0;
    detDist_nextLatBin = detDist_nextLonBin;
    nextLatBin = floor((fDetLat+0.5*M_PI)*fInvLatBinWidth);
    dLat = 0;
  } else {
    if (sinA == 0) {
      if (alpha == 0) { //crosses "at" infinity
        detDist_nextLonBin = -1;
      } else { //longitude flips at center of the Earth
        detDist_nextLonBin = rDetCosDetLat/alpha;
        nextLonBin = floor(fDetLon*fInvLonBinWidth);
      }
      dLon = 0;
    } else if (cosDetLat == 0) {
      detDist_nextLonBin = -1;
    } else {
      //x(lon) calculation
      detDist_nextLonBin = DetDistForNextLonBin(fEarthBins[index].longitude, dLon, nextLonBin, min_lon, max_lon, sinDetLon, cosDetLon, alpha, sinTsinAsinDetLon, sinTsinAcosDetLon, rDetCosDetLat);
    }
    //x(lat) calculation
    detDist_nextLatBin = DetDistForNextLatBin(index, sign, dLat, maxSinSqLat, maxlatreached, beta, gamma, gammaSq, rDetGammaSinDetLat, rDetCosAcosDetLat, rDetSinDetLat, rDetCosT, rDetSinT);
    nextLatBin = init_latBin + sign;
  }

  // Get initial Detector Distance
  double prev_DetDistance = baseline;

  /* Start at the top layer of Earth model and go down to minimum (not including minimum), looping over depth bins */
  int dBin = 0;
  int latBin = init_latBin;
  int lonBin = init_lonBin;
  for(dBin = 0; dBin<fnDepthBins-1; dBin++) {
    //Check if minimum radius is contained in depth bin
    double r_in = fEarthBins[index].radius_in; //inner-most radius
    if(r_in < minRadius) break;

    //Distance from detector at inner-most radius
    double DetDistance = -rDetCosT + sqrt(r_in*r_in - rDetSqSinSqT);

    //Check for latitude/longitude bin crossings
    RecordLatLonBinCrossings(DetDistance, prev_DetDistance, index, latBin, nextLatBin, detDist_nextLatBin, sign, dLat, maxSinSqLat, maxlatreached, lonBin, nextLonBin, detDist_nextLonBin, dLon, min_lon, max_lon, alpha, beta, gamma, gammaSq, rDetGammaSinDetLat, rDetCosAcosDetLat, rDetSinDetLat, rDetCosDetLat, rDetCosT, rDetSinT, sinTsinAsinDetLon, sinTsinAcosDetLon, sinDetLon, cosDetLon);

    //Add Segment to next depth bin to path
    AddPath(prev_DetDistance - DetDistance, fEarthBins[index]);

    //Reset variables for next depth bin
    prev_DetDistance = DetDistance;
    index -= binsPerDepth;
  }

  /* Do minimum depth bin and then go up to the detector */
  int index2 = index;
  for(index2 = index; fEarthBins[index2].radius_out < fDetRadius; index2+=binsPerDepth) {
    double r_out = fEarthBins[index2].radius_out; //outer-most radius

    //Distance from detector at outer-most radius
    double DetDistance = -rDetCosT - sqrt(r_out*r_out - rDetSqSinSqT);

    //Check for latitude/longitude bin crossings
    RecordLatLonBinCrossings(DetDistance, prev_DetDistance, index2, latBin, nextLatBin, detDist_nextLatBin, sign, dLat, maxSinSqLat, maxlatreached, lonBin, nextLonBin, detDist_nextLonBin, dLon, min_lon, max_lon, alpha, beta, gamma, gammaSq, rDetGammaSinDetLat, rDetCosAcosDetLat, rDetSinDetLat, rDetCosDetLat, rDetCosT, rDetSinT, sinTsinAsinDetLon, sinTsinAcosDetLon, sinDetLon, cosDetLon);

    //Add Segment in depth bin to path
    AddPath(prev_DetDistance - DetDistance, fEarthBins[index2]);

    //Reset variables for "previous" depth bin
    prev_DetDistance = DetDistance;
  }

  /* Do path to detector */
  //Check for latitude/longitude bin crossings
  RecordLatLonBinCrossings(0, prev_DetDistance, index2, latBin, nextLatBin, detDist_nextLatBin, sign, dLat, maxSinSqLat, maxlatreached, lonBin, nextLonBin, detDist_nextLonBin, dLon, min_lon, max_lon, alpha, beta, gamma, gammaSq, rDetGammaSinDetLat, rDetCosAcosDetLat, rDetSinDetLat, rDetCosDetLat, rDetCosT, rDetSinT, sinTsinAsinDetLon, sinTsinAcosDetLon, sinDetLon, cosDetLon);
  //Add the path segment within the detector bin
  AddPath(prev_DetDistance, fEarthBins[index2]);

  //Longitude calculation error message, if needed
  if (fLonError < 0) {
    cout << "ERROR:  Oops... It turns out that I was wrong about the denominator of the x(long) equation never being 0 apart from the cases where sin(Zenith) = 0, sin(Azimuthal) = 0, or cos(Detector Longitude) = 0.  Send the following message to rpestes@apc.in2p3.fr:" << endl << endl;
    cout << "You were wrong about the denominator of x(long)... Here are the values of the variables that broke it:" << endl;
    cout << "\tDetector Coordinates (r,lat,lon) = (" << fDetRadius << ", " << fDetLat << ", " << fDetLon << ")" << endl;
    cout << "\tcos(Zenith Angle) = " << cosT << endl;
    cout << "\tAzimuthal Angle = " << Az << endl;
    cout << "\t" << fErrorMessage_LonInfo << endl;
    cout << "Please fix this ASAP!" << endl;
    return -1;
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
