/*************************************************************************
 * EarthModelBinned.cxx
 * Implements a 3D earth model (binned in depth, latitude, and longitude)
 * for the programming library OscProb
 * Created by Rebekah Pestes (last modified: 4/4/23)
 *************************************************************************/

#include <cmath>
#include <fstream>
#include <iostream>

#include "EarthModelBinned.h"
#include "prem_default.hpp"

using namespace std;

using namespace OscProb;

//.............................................................................
///
/// Update values of zenith angle and azimuthal angle for neutrino
/// trajectory.
///
/// @param cosTheta - Cosine of the zenith angle
/// @param phi      - The azimuthal angle in radians
///
void TrajConstants::UpdateNuAngles(double cosTheta, double phi)
{
  cosT   = cosTheta;
  sinSqT = 1 - cosT * cosT;
  sinT   = sqrt(sinSqT);

  cosA = cos(phi);
  sinA = sin(phi);

  sinTsinA = sinT * sinA;
  sinTcosA = sinT * cosA;
}

//.............................................................................
///
/// Update detector position for neutrino trajectory calculations.
///
/// @param rDet   - The distance from the center of the Earth to the detector
/// (in km)
/// @param DetLat - The latitude of the detector (in rad)
/// @param DetLon - The longitude of the detector (in rad)
///
void TrajConstants::UpdateDetPos(double rDet, double DetLat, double DetLon)
{
  sinDetLat = sin(DetLat);
  cosDetLon = cos(DetLon);
  sinDetLon = sin(DetLon);
  DetRadius = rDet;

  cosDetLat     = 1 - sinDetLat * sinDetLat;
  rDetSinDetLat = DetRadius * sinDetLat;
  rDetCosDetLat = DetRadius * cosDetLat;
}

//.............................................................................
///
/// Recalculate constants that use combinations of detector position
/// variables and neutrino direction angles.
///
void TrajConstants::Recalculate()
{
  cosTcosDetLat      = cosT * cosDetLat;
  sinTsinAsinDetLon  = sinTsinA * sinDetLon;
  sinTsinAcosDetLon  = sinTsinA * cosDetLon;
  rDetCosT           = DetRadius * cosT;
  rDetSinT           = DetRadius * sinT;
  rDetCosAcosDetLat  = cosA * rDetCosDetLat;
  alpha              = sinTcosA * sinDetLat - cosTcosDetLat;
  beta               = sinT * sinDetLat - cosA * cosTcosDetLat;
  gamma              = sinTcosA * cosDetLat + cosT * sinDetLat;
  gammaSq            = gamma * gamma;
  rDetGammaSinDetLat = rDetSinDetLat * gamma;
  maxSinSqLat        = 1 - pow(sinA * cosDetLat, 2);
}

//.............................................................................
///
/// Constructor.
///
/// By default this implements the model stored in
/// EarthTables/earth_binned_default.txt
/// with the detector 6368 km from the center of the Earth, having a
/// latitude of 0 deg N and a longitude of 0 deg E.
///
/// @param filename - The txt file containing a table of earth layers
///
EarthModelBinned::EarthModelBinned(string filename)
{
  SetDetPos(6368, 0, 0);
  LoadModel(filename);
  SetRemoveSmallPaths(false);
  fLonError = 0;
}

//.............................................................................
///
/// Nothing to clean.
///
EarthModelBinned::~EarthModelBinned() {}

//.............................................................................
///
/// Get the set of Earth bins
///
/// This returns the set of bins for this Earth model.
///
vector<EarthBin> EarthModelBinned::GetEarthBins() { return fEarthBins; }

//.............................................................................
///
/// Clear the earth model information.
///
void EarthModelBinned::ClearModel()
{
  fnDepthBins = 0;
  fnLonBins   = 0;
  fnLatBins   = 0;
  fEarthBins.clear();
}

//.............................................................................
///
/// Set the coordinates of the detector:
///   radius in km, latitude in degrees, longitude in degrees
/// Also, updates the detector-coordinate-dependent parts of the trajectory
/// constants.
///
/// @param rad - The distance from the detector to the Earth's center in km
/// @param lat - The latitude of the detector in deg N (between -90 and 90)
/// @param lon - The longitude of the detector in deg E (between 0 and 360)
///
void EarthModelBinned::SetDetPos(double rad, double lat, double lon)
{
  SetDetectorCoordinates(rad, lat, lon);
  fC.UpdateDetPos(fDetRadius, fDetLat, fDetLon);
}

//.............................................................................
///
/// Add a bin to the earth model.
///
/// @param radius_out  - The outer depth of the bin in km
/// @param radius_in  - The inner depth of the bin in km
/// @param latitude - The latitude of the bin center in degrees
/// @param longitude - The longitude of the bin center in degrees
/// @param density - The density of the matter in the bin in g/cm^3
/// @param zoa     - The effective Z/A value of the matter in the bin
/// @param index   - Region index (all bins with same index have same Z/A)
///
void EarthModelBinned::AddBin(double radius_out, double radius_in,
                              double latitude, double longitude, double density,
                              double zoa, double index)
{
  fEarthBins.push_back(EarthBin(radius_out, radius_in, latitude, longitude,
                                density, zoa, index));
}

//.............................................................................
///
/// Load an earth model from a file.
///
/// By default it loads the model stored in
/// EarthTables/earth_binned_default.txt
///
/// The row format for the model table needs to be:
/// 	longitude (deg)	latitude (deg)	outer depth (km)	density
/// (g/cm^3?)	Z/A	region_index The data in the table must be in order of
/// decreasing depth, followed by increasing longitude, and then increasing
/// latitude. Latitude bin widths and longitude bin widths must each be
/// constant. In each bin of one coordinate, the bins centers for the other
/// coordinates must be the same.
/// This assumes that the center of the Earth is at a depth of 6371 km.
///
/// @param filename - The txt file containing a table of earth layers
///
void EarthModelBinned::LoadModel(string filename)
{
  // Clear the current model
  ClearModel();
  double EarthRadius = 6371.0;

  // Use default if no file provided
  if (filename == "") { filename = PREM3D_DEFAULT; }
  cout << "Loading Earth model table from " << filename << "..." << endl;

  // Open the file
  ifstream fin;
  fin.open(filename.c_str());
  if (!fin) {
    cerr << "ERROR: File " << filename << " not found!" << endl;
    return;
  }

  // Variables for storing table rows
  float depth, latitude, longitude, density, zoa, index;

  // Keep track of previous depth/latitude/longitude
  double radius_prev = 0;
  double lat_prev    = -90.0;
  double lon_prev    = 0;

  // Hold onto previous different depth, so we have depth maximum for bin
  double radius_min = radius_prev;

  // Initialize max region index
  fmaxRegIndex = -1;

  // Loop over table rows
  while (fin >> longitude >> latitude >> depth >> density >> zoa >> index) {
    double outer_radius = EarthRadius - depth;

    // Move minimum depth (and reset latitude/longitude) if depth has changed
    // from previous one
    if (outer_radius > radius_prev) {
      radius_min = radius_prev;
      lat_prev   = -90.0;
      lon_prev   = 0.0;
      fnDepthBins++;
    }
    else if (outer_radius <
             radius_prev) { // Depths must be ordered in model file
      cerr << "ERROR: Depths are not sorted in decreasing order in the model "
              "file "
           << "(or depth is greater than " << EarthRadius << " km)." << endl;
      ClearModel();
      return;
    }
    else if (longitude > lon_prev) { // Reset latitude if longitude has changed
                                     // from previous one
      lat_prev = -90.0;
      fnLonBins++;
    }
    else if (longitude < lon_prev) { // Longitudes must be ordered in model file
                                     // (after depths)
      cerr << "ERROR: Longitudes are not sorted in increasing order (after "
              "depths) in the model file."
           << endl;
      ClearModel();
      return;
    }
    else if (latitude < lat_prev) { // Latitudes must be ordered in model file
                                    // (after longitudes)
      cerr << "ERROR: Latitudes are not sorted in increasing order (after "
              "longitudes) in the model file."
           << endl;
      ClearModel();
      return;
    }

    // Add this bin to the model
    AddBin(outer_radius, radius_min, latitude, longitude, density, zoa, index);

    // Check region index
    if (index > fmaxRegIndex) { fmaxRegIndex = index; }

    // Set previous coordinates for next bin
    radius_prev = outer_radius;
    lat_prev    = latitude;
    lon_prev    = longitude;
  }

  // Set Outermost Radius
  fRadiusMax = radius_prev;

  // Modify number of longitude bins to how many there actually are (instead of
  // how many times we changed the longitude bin without changing the depth bin
  // = nLonBins-1 times per depth bin)
  fnLonBins        = fnLonBins / fnDepthBins + 1;
  fInvLonBinWidth  = 0.5 * fnLonBins / M_PI; // in radians^(-1)
  fHalfLonBinWidth = M_PI / fnLonBins;       // in radians

  // Calculate number of latitude bins
  fnLatBins        = fEarthBins.size() / (fnLonBins * fnDepthBins);
  fInvLatBinWidth  = fnLatBins / M_PI;       // in radians^(-1)
  fHalfLatBinWidth = 0.5 * M_PI / fnLatBins; // in radians

  // Error if first latitude bin center deviates from about half of the
  // calculated latitude bin width
  if (abs(fHalfLatBinWidth - (fEarthBins[0].latitude + 0.5 * M_PI)) >
      fHalfLatBinWidth * 0.1) { // Warn if first latitude bin center isn't about
                                // half of the latitude bin width
    cerr << "ERROR: Latitude of 1st bin center (" << fEarthBins[0].latitude
         << ") is not in the middle of the first bin, whose size should be "
         << 2 * fHalfLatBinWidth << " rad, as calculated from having "
         << fnLatBins << " latitude bins." << endl;
    ClearModel();
    return;
  }

  cout << "\t...done (" << fEarthBins.size() << " bins: " << fnDepthBins
       << " depth, " << fnLatBins << " latitude, and " << fnLonBins
       << " longitude)." << endl;
}

//.............................................................................
///
/// Set the effective Z/A value for all bins with given region index.
///
/// Use this to change the Z/A of indexed bin,
/// e.g. all outer-core layers
///
/// @param index - The region index
/// @param zoa   - The effective Z/A value to use
///
void EarthModelBinned::SetRegionZoA(int index, double zoa)
{
  if (index > fmaxRegIndex) {
    cerr << "WARNING: Index " << index << " does not exist in model "
         << "(max index = " << fmaxRegIndex << ")." << endl
         << "Doing nothing." << endl;
    return;
  }

  int nbins = fEarthBins.size();

  // Loop over all bins and change the ones
  // with the given index
  for (int i = 0; i < nbins; i++) {
    if (fEarthBins[i].index != index) continue;

    fEarthBins[i].zoa = zoa;
  }
}

//.............................................................................
///
/// Get the effective Z/A value for all bins in region specified by
/// index, e.g. all outer-core layers.
/// (Assumes that all bins of given index have same Z/A value)
///
/// @param index - The region index
/// @return Z/A corresponding to index
///
double EarthModelBinned::GetRegionZoA(int index)
{
  if (index > fmaxRegIndex) {
    cerr << "ERROR: Index " << index << " does not exist in model "
         << "(max index = " << fmaxRegIndex << "). Returning 0" << endl;
    return 0;
  }

  int nbins = fEarthBins.size();

  for (int i = 0; i < nbins; i++) {
    if (fEarthBins[i].index != index)
      continue;

    else
      return fEarthBins[i].zoa;
  }

  // End of vector reached without finding index
  cerr << "ERROR: Index " << index << " not found! Returning 0" << endl;
  return 0;
}

//.............................................................................
///
/// Scale density by scaling factor for all bins in region specified by
/// given index.
///
/// @param index  - The region index
/// @param factor - The scaling factor for changing the density
///
void EarthModelBinned::ScaleRegionDensity(int index, double factor)
{
  if (index > fmaxRegIndex) {
    cerr << "WARNING: Index " << index << " does not exist in model "
         << "(max index = " << fmaxRegIndex << ")." << endl
         << "Doing nothing." << endl;
    return;
  }

  int nbins = fEarthBins.size();

  // Loop over all bins and change the ones
  // with the given index
  for (int i = 0; i < nbins; i++) {
    if (fEarthBins[i].index != index) continue;

    double density_init   = fEarthBins[i].density;
    fEarthBins[i].density = factor * density_init;
  }
}

//.............................................................................
///
/// Add a path segment to the sequence.
///
/// For a given EarthBin, adds a path of a given length in that material
///
/// @param length  - The length of the path segment in km
/// @param bin     - The bin we are crossing
///
void EarthModelBinned::AddPath(double length, EarthBin bin)
{
  AddPathSegment(length, bin.density, bin.zoa, bin.index);
}

//.............................................................................
///
/// Calculate the distance to the detector along a neutrino trajectory,
/// specified by cosT and the azimuthal angle, at the edge of
/// the current latitude bin in the direction specified by dir.
///
/// In an attempt to make this calculation more robust, this function turns
/// the latitude bin iteration around when sin^2(lat) gets within 10^(-14)
/// of its extreme limit.
///
/// @param cur_index - Index of current Earth bin
/// @param L         - Information for the lat bin crossings
/// @return Calculated detector distance
///
double EarthModelBinned::DetDistForNextLatBin(int cur_index, LatBinInfo& L)
{
  // Find next lat bin boundary
  double NewLat = fEarthBins[cur_index].latitude + L.sign * L.dLat;
  // Turn around if lat reaches +/- 90 degrees
  if (abs(NewLat) - 0.5 * M_PI > -1.0e-14) {
    if (L.maxreached > 0) return -1; // Only turn around once
    L.maxreached = 1;
    L.sign       = -L.sign;
    NewLat       = fEarthBins[cur_index].latitude + L.sign * L.dLat;
    if (abs(NewLat) - 0.5 * M_PI > -1.0e-14)
      return -1; // This would happen if there's only one latitude bin...
  }

  // Calculate parts of answer
  double sinNewLat    = sin(NewLat);
  double NsinSqNewLat = -sinNewLat * sinNewLat;
  double radicand     = fC.maxSinSqLat + NsinSqNewLat;
  // Turn around if the radicand is 0 (close enough) or less
  if (radicand < 1.0e-14) {
    if (L.maxreached > 0) return -1; // Only turn around once
    L.maxreached = 1;
    L.sign       = -L.sign;
    NewLat       = fEarthBins[cur_index].latitude + L.sign * L.dLat;
    sinNewLat    = sin(NewLat);
    NsinSqNewLat = -sinNewLat * sinNewLat;
    radicand     = fC.maxSinSqLat + NsinSqNewLat;
    if (radicand < 1.0e-14)
      return -1; // The rest of the track is contained in this lon bin
  }
  double denominator = fC.gammaSq + NsinSqNewLat;

  // Calculate distance from the detector at next lat bin boundary
  double answer = 0;
  if (denominator == 0) {
    answer =
        0.5 * (fC.rDetCosAcosDetLat / fC.beta - fC.rDetSinDetLat / fC.gamma);
  }
  else {
    answer = -(fC.rDetCosT * NsinSqNewLat + fC.rDetGammaSinDetLat +
               L.sign * fC.rDetSinT * sinNewLat * sqrt(radicand)) /
             denominator;
  }
  // Indicates that going to the edge of the bin jumped to an invalid part of
  // the function
  if (L.maxreached > 0 && answer >= fC.rDetCosAcosDetLat / fC.beta) {
    return -1; // The rest of the track is contained in this lon bin
  }

  return answer;
}

//.............................................................................
///
/// Calculate the distance to the detector at the edge of the current
/// longitude bin along a neutrino trajectory specified by cosT and the
/// azimuthal angle, setting the index for the next longitude bin in the
/// process. The maximum/minimum longitude are expressed around the 0-2pi
/// boundary. So, if the track spans that boundary, the actual allowed
/// range for lon is [min_lon,2pi)U[0,max_lon].
///
/// @param prev_lon - Value at center of initial longitude bin
/// @param L        - Information about longitude bin crossings
/// @return Calculated detector distance
///
double EarthModelBinned::DetDistForNextLonBin(double prev_lon, LonBinInfo& L)
{
  double twoPi = 2 * M_PI;
  // Increment to next bin in direction dLon
  double lon = prev_lon + L.dLon;
  if (lon < 0 || lon >= twoPi) lon -= floor(lon / twoPi) * twoPi;
  if (L.dLon > 0) {
    L.nextBin = L.bin + 1;
    if (L.nextBin >= fnLonBins) {
      // cross from 2pi to 0
      L.nextBin = 0;
    }
  }
  else {
    L.nextBin = L.bin - 1;
    if (L.nextBin < 0) {
      // cross from 0 to 2pi
      L.nextBin = fnLonBins - 1;
    }
  }

  // Check if lon is outside of range covered by neutrino's trajectory
  if (L.max < L.min) {
    if (0 > lon - L.min && L.max - lon < 0) return -1;
  }
  else {
    if (0 > lon - L.min || L.max - lon < 0) return -1;
  }

  // Parts of the answer
  double cosNewLon   = cos(lon);
  double sinNewLon   = sin(lon);
  double expression  = fC.sinDetLon * cosNewLon - sinNewLon * fC.cosDetLon;
  double denominator = fC.alpha * expression +
                       fC.sinTsinAcosDetLon * cosNewLon +
                       fC.sinTsinAsinDetLon * sinNewLon;

  // Calculate & return detector distance
  if (denominator == 0) { // Is there ever a case when this is true???
    // Error message will be printed at the end of the FillPath function
    L.err_message.assign("Longitude along neutrino trajectory = " +
                         std::to_string(lon));
    L.error = -1;
    return -1;
  }
  return fC.rDetCosDetLat * expression / denominator;
}

//.............................................................................
///
/// Record all the path segments as the neutrino crosses into new latitude
/// and longitude bins until it reaches the next depth bin.
///
/// All the variables not listed below are the same as the arguments for
/// the DetDistForNextLonBin and DetDistForNextLatBin functions.
///
/// @param detDist_nextDbin - Distance from the detector at which neutrino
/// passes into next depth bin
/// @param DetDist          - Current distance from the detector
/// @param index            - Index of current Earth bin
/// @param latI             - Information about lat bin crossings
/// @param lonI             - Information about lon bin crossings
///
void EarthModelBinned::RecordLatLonBinCrossings(double  detDist_nextDbin,
                                                double& DetDist, int& index,
                                                LatBinInfo& latI,
                                                LonBinInfo& lonI)
{
  while (detDist_nextDbin < latI.detDist_nextBin ||
         detDist_nextDbin < lonI.detDist_nextBin) {
    if (latI.detDist_nextBin > lonI.detDist_nextBin) {
      // Add segment to next lat bin to path
      AddPath(DetDist - latI.detDist_nextBin, fEarthBins[index]);

      // Move lat bins
      DetDist = latI.detDist_nextBin;
      index += latI.nextBin - latI.bin;
      latI.bin = latI.nextBin;

      // Find next lat bin
      if (latI.dLat == 0) { latI.detDist_nextBin = -1; }
      else {
        latI.detDist_nextBin = DetDistForNextLatBin(index, latI);
        latI.nextBin         = latI.bin + latI.sign;
      }
    }
    else {
      // Add segment to next lon bin to path
      AddPath(DetDist - lonI.detDist_nextBin, fEarthBins[index]);

      // Move lon bins
      DetDist = lonI.detDist_nextBin;
      index += (lonI.nextBin - lonI.bin) * fnLatBins;
      lonI.bin = lonI.nextBin;

      // Find next lon bin
      if (lonI.dLon == 0) { lonI.detDist_nextBin = -1; }
      else {
        lonI.detDist_nextBin =
            DetDistForNextLonBin(fEarthBins[index].longitude, lonI);
      }
    }
  }
}

//.............................................................................
///
/// Find and return longitude bin index that contains longitude.
///
/// @param longitude - longitude for which bin is found
///
int EarthModelBinned::LonBinIndex(double longitude)
{
  // Find difference between longitude and beginning of first lon bin
  double londiff = longitude - (fEarthBins[0].longitude - fHalfLonBinWidth);

  // Shift londiff to be between 0 and 2pi
  if (londiff < 0 || londiff >= 2 * M_PI)
    londiff -= floor(londiff / (2 * M_PI));

  // Number of binwidths after from first lon bin
  return floor(londiff * fInvLonBinWidth);
}

//.............................................................................
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
/// @param phi - The azimuthal angle (in degrees from North) of the neutrino
/// direction
/// @return The number of path segments in the sequence
///
int EarthModelBinned::FillPath(double cosT, double phi)
{
  // Clear current path sequence
  fNuPath.clear();

  // Do nothing if cosine is unphysical
  if (fabs(cosT) > 1) return 0;

  // Calculate useful combinations of variables
  double Az = phi / 180.0 * M_PI; // convert from degrees to radians
  fC.UpdateNuAngles(cosT, Az);
  fC.Recalculate();
  double rDetSqSinSqT = pow(fC.rDetSinT, 2);
  double distEdgeToMin =
      sqrt(fRadiusMax * fRadiusMax -
           rDetSqSinSqT); // distance from point where neutrino enters to
                          // minimum depth (neglecting stopping at the detector)
  double baseline = distEdgeToMin - fC.rDetCosT; // total baseline

  // Define the maximum depth (minimum radius) along the path (detector
  // depth/radius, if cosT is non-negative)
  double minRadius = fDetRadius;
  if (cosT < 0) { minRadius = fC.rDetSinT; }

  // Starting latitude
  double init_lat = asin((fC.rDetSinT * fC.beta + fC.gamma * distEdgeToMin) /
                         fRadiusMax); // in radians

  // Starting longitude
  double init_lon = atan2(
      fC.sinDetLon * fC.rDetCosDetLat -
          baseline * (fC.alpha * fC.sinDetLon + fC.sinTsinAcosDetLon),
      fC.cosDetLon * fC.rDetCosDetLat +
          baseline *
              (fC.sinTsinAsinDetLon -
               fC.alpha * fC.cosDetLon)); // in radians (between -M_PI and M_PI)
  if (init_lon < 0) { init_lon += 2 * M_PI; }

  // Get initial Earth bin and # of bins per depth layer
  LatBinInfo latinfo;
  LonBinInfo loninfo;
  int        binsPerDepth = fEarthBins.size() / fnDepthBins;
  int        init_lonBin  = LonBinIndex(init_lon);
  int        init_latBin  = floor((init_lat + M_PI / 2) * fInvLatBinWidth);
  int        index =
      fEarthBins.size() - binsPerDepth + init_lonBin * fnLatBins + init_latBin;

  latinfo.bin = init_latBin;
  loninfo.bin = init_lonBin;

  /* Find detector distances for 1st latitude/longitude bin changes */
  latinfo.dLat = fHalfLatBinWidth; // change in lat from bin center to next bin
                                   // (0 if no more than 1 lat bin change)
  // Condition calculation for starting sign in x(lat) equation
  latinfo.sign       = 1; //+ => lat incr., - => lat decr. (with decr. x)
  latinfo.maxreached = 0; // this changes to 1 if past the lat func transition
  if (fC.beta < 0) latinfo.sign = -1;
  if (fC.beta == 0) {
    if (fC.cosA > 0) latinfo.sign = -1;
  }
  else if (baseline < fC.rDetCosAcosDetLat / fC.beta) {
    latinfo.sign = -(latinfo.sign);
    latinfo.maxreached++;
  }
  // Condition calculation for incr./decr. Longitude (with decr. x)
  loninfo.dLon =
      fHalfLonBinWidth;   // change in lon from bin center to next bin (includes
                          // direction; 0 if no more than 1 lon bin change)
  loninfo.min = init_lon; // initial lon if lon is inc.
  loninfo.max = fDetLon;  // detector lon if lon is inc.
  if (fC.sinA < 0) {
    loninfo.dLon = -(loninfo.dLon);
    loninfo.min  = loninfo.max;
    loninfo.max  = init_lon;
    // Switch max/min if dec., instead of inc.
    loninfo.max = init_lon;
    loninfo.min = fDetLon;
  }
  // Remainder of detDist (and nextBin) calculations
  if (fC.sinT == 0) {
    // both latitude and logitude flip at center of the Earth
    loninfo.detDist_nextBin =
        -fC.rDetCosT; // formula is divided by cosT, but cosT = +-1
    loninfo.nextBin         = LonBinIndex(fDetLon);
    loninfo.dLon            = 0;
    latinfo.detDist_nextBin = loninfo.detDist_nextBin;
    latinfo.nextBin         = floor((fDetLat + 0.5 * M_PI) * fInvLatBinWidth);
    latinfo.dLat            = 0;
  }
  else {
    if (fC.sinA == 0) {
      if (fC.alpha == 0) { // crosses "at" infinity
        loninfo.detDist_nextBin = -1;
      }
      else { // longitude flips at center of the Earth
        loninfo.detDist_nextBin = fC.rDetCosDetLat / fC.alpha;
        loninfo.nextBin         = LonBinIndex(fDetLon);
      }
      loninfo.dLon = 0;
    }
    else if (fC.cosDetLat == 0) {
      loninfo.detDist_nextBin = -1;
    }
    else {
      // x(lon) calculation
      loninfo.detDist_nextBin = DetDistForNextLonBin(
          fEarthBins[index].longitude, loninfo); // also sets loninfo.nextBin
    }
    // x(lat) calculation
    latinfo.detDist_nextBin = DetDistForNextLatBin(index, latinfo);
    latinfo.nextBin         = init_latBin + latinfo.sign;
  }

  // Get initial Detector Distance
  double prev_DetDistance = baseline;

  /* Start at the top layer of Earth model and go down to minimum (not including
   * minimum), looping over depth bins */
  int dBin = 0; // depth bin index
  for (dBin = 0; dBin < fnDepthBins - 1; dBin++) {
    // Check if minimum radius is contained in depth bin
    double r_in = fEarthBins[index].radius_in; // inner-most radius
    if (r_in < minRadius) break;

    // Distance from detector at inner-most radius
    double DetDistance = -fC.rDetCosT + sqrt(r_in * r_in - rDetSqSinSqT);

    // Check for latitude/longitude bin crossings
    RecordLatLonBinCrossings(DetDistance, prev_DetDistance, index, latinfo,
                             loninfo);

    // Add Segment to next depth bin to path
    AddPath(prev_DetDistance - DetDistance, fEarthBins[index]);

    // Reset variables for next depth bin
    prev_DetDistance = DetDistance;
    index -= binsPerDepth;
  }

  /* Do minimum depth bin and then go up to the detector */
  int index2 = index;
  for (index2 = index; fEarthBins[index2].radius_out < fDetRadius;
       index2 += binsPerDepth) {
    double r_out = fEarthBins[index2].radius_out; // outer-most radius

    // Distance from detector at outer-most radius
    double DetDistance = -fC.rDetCosT - sqrt(r_out * r_out - rDetSqSinSqT);

    // Check for latitude/longitude bin crossings
    RecordLatLonBinCrossings(DetDistance, prev_DetDistance, index2, latinfo,
                             loninfo);

    // Add Segment in depth bin to path
    AddPath(prev_DetDistance - DetDistance, fEarthBins[index2]);

    // Reset variables for "previous" depth bin
    prev_DetDistance = DetDistance;
  }

  /* Do path to detector */
  // Check for latitude/longitude bin crossings
  RecordLatLonBinCrossings(0, prev_DetDistance, index2, latinfo, loninfo);
  // Add the path segment within the detector bin
  AddPath(prev_DetDistance, fEarthBins[index2]);

  // Longitude calculation error message, if needed
  if (loninfo.error < 0) {
    cerr
        << "ERROR:  Oops... It turns out that I was wrong about the denominator"
        << " of the x(long) equation never being 0 apart from the cases where "
        << "sin(Zenith) = 0, sin(Azimuthal) = 0, or cos(Detector Longitude) = "
           "0."
        << "  Send the following message to rpestes@apc.in2p3.fr:" << endl
        << endl
        << "You were wrong about the denominator of x(long)... "
        << "Here are the values of the variables that broke it:" << endl
        << "\tDetector Coordinates (r,lat,lon) = (" << fDetRadius << ", "
        << fDetLat << ", " << fDetLon << ")" << endl
        << "\tcos(Zenith Angle) = " << cosT << endl
        << "\tAzimuthal Angle = " << phi << " deg" << endl
        << "\t" << loninfo.err_message << endl
        << "Please fix this ASAP!" << endl;
    return -1;
  }

  // Return the number of path segments
  return fNuPath.size();
}
