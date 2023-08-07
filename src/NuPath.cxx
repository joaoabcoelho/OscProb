////////////////////////////////////////////////////////////////////////
//
// Implements functions to average neutrino paths
//
// jcoelho\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

#include "NuPath.h"

using namespace std;

using OscProb::NuPath;

//......................................................................
///
/// Get the merged average of two paths
///
/// This method will merge two paths and take their average density
/// weighted by Z/A and path length.
///
/// The Z/A will be the average weighted by path length
///
/// @param p1 - The first path to merge
/// @param p2 - The second path to merge
/// @return The merged path
///
NuPath OscProb::AvgPath(NuPath& p1, NuPath& p2){

  // Start with the first path
  NuPath mergedPath = p1;

  // Add the second length
  mergedPath.length += p2.length;

  // Compute weighted average of Z/A
  mergedPath.zoa = (p1.zoa*p1.length + p2.zoa*p2.length) / (p1.length + p2.length);

  // Compute weighted average of density
  mergedPath.density = (p1.density*p1.zoa*p1.length + p2.density*p2.zoa*p2.length) / (p1.zoa*p1.length + p2.zoa*p2.length);

  // return merged path
  return mergedPath;

}

//......................................................................
///
/// Get the merged average of a vector of paths
///
/// This method will merge a set of paths and take their average density
/// weighted by Z/A and path length.
///
/// The Z/A will be the average weighted by path length
///
/// @param pv - vector of paths to merge
/// @return The merged path
///
NuPath OscProb::AvgPath(std::vector<NuPath>& pv){

  // Get size of vector
  int np = pv.size();

  // Start with the first path
  NuPath mergedPath;

  // If vector is not empty, start on first path
  if(np>0) mergedPath = pv[0];
  else return mergedPath;

  // Merge each of the following paths
  for(int i=1; i<np; i++){
    mergedPath = AvgPath(mergedPath, pv[i]);
  }

  // return merged path
  return mergedPath;

}

//......................................................................
///
/// Merge two specific paths by their indices in a path vector
///
/// @param inputPath - The original vector of paths to merge
/// @param j,k - The indices of the two paths to merge
/// @return The merged vector of paths
///
vector<NuPath> OscProb::MergePaths(std::vector<NuPath>& inputPath, int j, int k){

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

