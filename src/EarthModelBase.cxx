/*************************************************************************
 * EarthModelBase.cxx
 * Base class for all Earth model clases in OscProb
 * (created by Rebekah Pestes based on PremModel.cxx)
 *************************************************************************/

#include <cmath>
#include <iostream>

#include "EarthModelBase.h"

using namespace std;

using namespace OscProb;

//......................................................................
///
/// Set the coordinates of the detector:
///   radius in km, latitude in degrees, longitude in degrees
///
/// @param rad - The distance from the detector to the Earth's center in km
/// @param lat - The latitude of the detector in deg N (between -90 and 90)
/// @param lon - The longitude of the detector in deg E (between 0 and 360)
///
void EarthModelBase::SetDetectorCoordinates(double rad, double lat, double lon)
{

  //Force radius to be non-negative
  if(rad < 0) {
    cerr << "WARNING: Negative radius detected. Setting to absolute value." << endl;
    rad = -rad;
  }
  fDetRadius = rad;

  //Force latitude to be between -90 and 90 deg
  lat -= floor((lat+90.0)/360.0)*360.0;
  if (lat > 90) {
    lat = 180.0-lat;
  }
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
vector<NuPath> EarthModelBase::GetNuPath(){ return fNuPath; }

//......................................................................
///
/// Add a path segment to the sequence.
///
/// For a given EarthBin, adds a path of a given length in that material
///
/// @param length  - The length of the path segment in km
/// @param density - The density along the path segment
/// @param zoa     - Z/A along the path segment
/// @param index   - Index for the matter along the path segment
///
void EarthModelBase::AddPathSegment(double length, double density, double zoa, int index)
{

  fNuPath.push_back( NuPath(length, density, zoa, index) );

}

//......................................................................
///
/// Get the total baseline for a given cosTheta.
///
/// @param cosT - The cosine of the neutrino direction
///
double EarthModelBase::GetTotalL(double cosT)
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
double EarthModelBase::GetCosT(double L)
{

  if(L < fRadiusMax - fDetRadius) return  1;
  if(L > fRadiusMax + fDetRadius) return -1;

  return (fRadiusMax*fRadiusMax - fDetRadius*fDetRadius - L*L) / (2*fDetRadius*L);

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
vector<NuPath> EarthModelBase::GetMergedPaths(double prec){

  // The output vector
  vector<NuPath> mergedPath;

  // Start with the first path
  NuPath path = fNuPath[0];

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
/// Set the boolean to tag whether to remove small paths when merging
/// Small is defined as <1% of the total baseline
///
/// @param rp - Boolean value to set
///
void EarthModelBase::SetRemoveSmallPaths(bool rp){ fRemoveSmallPaths = rp; }
