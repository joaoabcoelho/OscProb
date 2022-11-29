////////////////////////////////////////////////////////////////////////
/// \struct OscProb::NuPath
///
/// \brief A struct representing a neutrino path segment
///
/// This struct stores the properties of a neutrino path segment
/// so that the neutrino propagation through a path is done
/// consistently in the PMNS classes.
///
/// \author jcoelho\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
/// \struct OscProb::PremLayer
///
/// \brief A struct representing a spherical shell of matter 
/// for earth models
///
/// This struct stores the properties of a spherical shell
/// to be used by the PremModel class in order to build an
/// earth model. Only the outer radius of the shell is stored,
/// so PremLayer's need to be assembled in order inside a vector.
///
/// \author jcoelho\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

#ifndef NUPATH_H
#define NUPATH_H

#include <vector>

namespace OscProb {

  struct NuPath
  {
  
    ///
    /// \brief Constructor.
    ///
    /// Constructor.
    ///
    /// By default it creates a path of zero length and zero density.
    /// The effective Z/A value is set to 0.5 by default. 
    ///
    /// The properties of the path can be given directly in the construction.
    ///
    /// @param l  - The length of the path segment in km
    /// @param d  - The density of the path segment in g/cm^3
    /// @param z  - The effective Z/A value of the path segment
    /// @param ly - An index to identify the matter type (e.g. earth inner core)
    ///
    NuPath(double l=0, double d=0, double z=0.5, int ly=0){
      SetPath(l,d,z,ly);
    }

    ///
    /// \brief Set the properties of the neutrino path.
    ///
    /// Set the properties of the neutrino path.
    ///
    /// By default it sets the path to zero length and zero density.
    /// The effective Z/A value is set to 0.5 by default. 
    ///
    /// @param l  - The length of the path segment in km
    /// @param d  - The density of the path segment in g/cm^3
    /// @param z  - The effective Z/A value of the path segment
    /// @param ly - An index to identify the matter type (e.g. earth inner core)
    ///
    void SetPath(double l=0, double d=0, double z=0.5, int ly=0){
      length = l;
      density = d;
      zoa = z;
      layer = ly;
    }

    double length;  ///< The length of the path segment in km
    double density; ///< The density of the path segment in g/cm^3
    double zoa;     ///< The effective Z/A value of the path segment
    int layer;      ///< An index to identify the matter type

  };

  struct PremLayer
  {

    ///
    /// \brief Constructor.
    ///
    /// Constructor.
    ///
    /// By default it creates a layer of zero radius and zero density.
    /// The effective Z/A value is set to 0.5 by default.
    ///
    /// The properties of the layer can be given directly in the construction.
    ///
    /// @param r  - The outer radius of the layer in km
    /// @param d  - The density of the layer in g/cm^3
    /// @param z  - The effective Z/A value of the layer
    /// @param ly - An index to identify the matter type (e.g. earth inner core)
    ///
    PremLayer(double r=0, double d=0, double z=0.5, int ly=0){
      SetLayer(r,d,z,ly);
    }

    ///
    /// \brief Set the properties of the layer.
    ///
    /// Set the properties of the layer.
    ///
    /// By default it sets the layer to zero radius and zero density.
    /// The effective Z/A value is set to 0.5 by default.
    ///
    /// @param r  - The outer radius of the layer in km
    /// @param d  - The density of the layer in g/cm^3
    /// @param z  - The effective Z/A value of the layer
    /// @param ly - An index to identify the matter type (e.g. earth inner core)
    ///
    void SetLayer(double r=0, double d=0, double z=0.5, int ly=0){
      radius = r;
      density = d;
      zoa = z;
      layer = ly;
    }

    double radius;  ///< The outer radius of the layer in km
    double density; ///< The density of the layer in g/cm^3
    double zoa;     ///< The effective Z/A value of the layer
    int layer;      ///< An index to identify the matter type

  };
  
  NuPath AvgPath(NuPath& p1, NuPath& p2); ///< Get the average of two paths
  NuPath AvgPath(std::vector<NuPath>& pv); ///< Get the average of a vector of paths
  std::vector<NuPath> MergePaths(std::vector<NuPath> &inputPath, int j, int k); ///< Merge paths j and k in vector

}

#endif
