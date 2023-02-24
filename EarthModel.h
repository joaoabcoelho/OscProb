////////////////////////////////////////////////////////////////////////
/// \class OscProb::EarthModel
///
/// \brief Implements an earth model with depth/latitude/longitude bins
///
/// This class implements a 3D (asymmetric) model of the earth
/// using EarthBin's to store bins with different properties.
///
/// The class is then able to produce path sequences through the Earth
/// as a function of the azimuthal angle and cosine of the zenith angle
/// with respect to the detector.
///
/// The detector can be positioned at any point within the Earth, and
/// the path sequences will take into account the fact that some layers
/// are above the detector.
///
/// By default this implements the model stored in EarthTables/prem_default.txt
/// with the detector at the bottom of the ocean layer (radius = 6368 km) at
/// the south pole.
///
/// This class inherits from TObject and can be saved in ROOT files.
///
/// \author rpestes\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
/// \struct OscProb::EarthBins
///
/// \brief A struct representing a 3D bin (depth, latitude, and longitude)
/// of matter for asymmetric Earth models
///
/// This struct stores the properties of a 3D bin (part of a spherical shell)
/// to be used by the EarthModel class in order to build an
/// earth model.
///
/// \author rpestes\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

#ifndef EARTHMODEL_H
#define EARTHMODEL_H

#include <string>
#include <vector>

#include <TObject.h>

#include <NuPath.h>

namespace OscProb {

  struct EarthBin
  {

    ///
    /// \brief Constructor.
    ///
    /// Constructor.
    ///
    /// By default it creates a bin of zero radius, zero latitude, zero longitude, and zero density.
    /// The effective Z/A value is set to 0.5 by default.
    ///
    /// The properties of the layer can be given directly in the construction.
    ///
    /// @param d_out  - The outer depth of the bin in km
    /// @param d_in  - The inner depth of the bin in km
    /// @param lat - Latitude of bin center in deg
    /// @param lon - Longitude of bin center in deg
    /// @param den  - The density of the matter in the bin in g/cm^3
    /// @param z  - The effective Z/A value of the matter in the bin
    /// @param n - An index to identify the matter type (e.g. earth inner core)
    ///
    EarthBin(double d_out=0, double d_in=0, double lat=0, double lon=0, double den=0, double z=0.5, int n=0){
      SetLayer(d_out, d_in, lat, lon, den, z, n);
    }

    ///
    /// \brief Set the properties of the bin.
    ///
    /// Set the properties of the bin.
    ///
    /// By default it creates a bin of zero radius, zero latitude, zero longitude, and zero density.
    /// The effective Z/A value is set to 0.5 by default.
    ///
    /// @param d_out  - The outer depth of the bin in km
    /// @param d_in  - The inner depth of the bin in km
    /// @param lat - Latitude of bin center in deg
    /// @param lon - Longitude of bin center in deg
    /// @param den  - The density of the matter in the bin in g/cm^3
    /// @param z  - The effective Z/A value of the matter in the bin
    /// @param n - An index to identify the matter type (e.g. earth inner core)
    ///
    void SetLayer(double d_out=0, double d_in=0, double lat=0, double lon=0, double den=0, double z=0.5, int n=0){
      depth_out = d_out;
      depth_in = d_in;
      latitude = lat;
      longitude = lon;
      density = den;
      zoa = z;
      index = n;
    }

    double depth_out;  ///< The outer depth of the bin in km
    double depth_in;  ///< The inner depth of the bin in km
    double latitude; ///< The latitude of the bin center in degrees
    double longitude; ///< The longitude of the bin center in degrees
    double density; ///< The density of the matter in the bin in g/cm^3
    double zoa;     ///< The effective Z/A value of the matter in the bin
    int index;      ///< An index to identify the matter type

  };

  class EarthModel : public TObject {

    public:

      EarthModel(std::string filename=""); ///< Constructor
      virtual ~EarthModel();          ///< Destructor

      virtual int FillPath(double cosT); ///< Fill the path sequence in a vector

      virtual std::vector<OscProb::NuPath> GetNuPath(); ///< Get the current neutrino path sequence

      virtual std::vector<OscProb::NuPath> GetMergedPaths(double prec = 0.25); ///< Get merged path sequence in a vector

      virtual double GetTotalL(double cosT); ///< Get the total baseline for a given cosTheta
      virtual double GetCosT(double L); ///< Get the cosTheta for a given total baseline

      virtual void SetLayerZoA(int layer, double zoa); ///< Set Z/A of all layers of a given type
      virtual double GetLayerZoA(int layer); ///< Get Z/A of all layers of a given type

//      virtual void SetTopLayerSize(double thick);      ///< Set the outermost layer thickness in km ***THIS NEEDS HELP!!!***

      virtual void LoadModel(std::string filename); ///< Load an earth model from a file

      virtual void SetDetPos(double dep, double lat, double lon); ///< Set the detector position (dep = depth in km, lat/lon in deg)

      virtual std::vector<OscProb::EarthBin> GetEarthBins(); ///< Get the set of earth layers

      virtual void SetRemoveSmallPaths(bool rp = true); ///< Set tag to remove small paths

    protected:

      virtual void ClearModel(); ///< Clear the earth model information

      virtual void AddBin(double depth_out, double depth_in, double latitude, double longitude, double density,
                            double zoa,    double layer); ///< Add a bin to the model

//      virtual void AddLayer(double radius, double density,
//                            double zoa,    double layer); ///< Add a layer to the model  ***THIS NEEDS HELP!!!***

      virtual void AddPath(double length, EarthBin bin);  ///< Add a path segment to the sequence

      std::vector<OscProb::EarthBin> fEarthBins; ///< The layers in the earth model
      int fnDepthBins; ///< Total number of depth bins
      int fnLatBins; ///< Total number of latitude bins

      std::vector<OscProb::NuPath> fNuPath; ///< The current neutrino path sequence

      double fEarthRadius; ///< Biggest depth value in Earth model

      int fDetBin;  ///< The index of bin containing the detector
      int fDetDepthBin;  ///< The index of depth bin containing the detector
//      double fDetPos; ///< The radius where the detector lives
      double fDetDepth; ///< The depth where the detector lives
      double fDetLat; ///< The latitude where the detector lives
      double fDetLon; ///< The longitude where the detector lives

      bool fRemoveSmallPaths; ///< Tag whether to merge small paths

      static const double DET_TOL; ///< The detector position tolerance near boundaries

      // Required for saving in ROOT files
      ClassDef(EarthModel, 1);

  };

}
#endif
