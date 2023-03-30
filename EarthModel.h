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
#include <math.h>

#include <TObject.h>

#include "NuPath.h"

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
    /// @param r_out  - The outer radius of the bin in km
    /// @param r_in  - The inner radius of the bin in km
    /// @param lat - Latitude of bin center in deg
    /// @param lon - Longitude of bin center in deg
    /// @param den  - The density of the matter in the bin in g/cm^3
    /// @param z  - The effective Z/A value of the matter in the bin
    /// @param n - An index to identify the matter type (e.g. earth inner core)
    ///
    EarthBin(double r_out=0, double r_in=0, double lat=0, double lon=0, double den=0, double z=0.5, int n=0){
      SetLayer(r_out, r_in, lat, lon, den, z, n);
    }

    ///
    /// \brief Set the properties of the bin.
    ///
    /// Set the properties of the bin.
    ///
    /// By default it creates a bin of zero radius, zero latitude, zero longitude, and zero density.
    /// The effective Z/A value is set to 0.5 by default.
    ///
    /// @param r_out  - The outer radius of the bin in km
    /// @param r_in  - The inner radius of the bin in km
    /// @param lat - Latitude of bin center in deg
    /// @param lon - Longitude of bin center in deg
    /// @param den  - The density of the matter in the bin in g/cm^3
    /// @param z  - The effective Z/A value of the matter in the bin
    /// @param n - An index to identify the matter type (e.g. earth inner core)
    ///
    void SetLayer(double r_out=0, double r_in=0, double lat=0, double lon=0, double den=0, double z=0.5, int n=0){
      radius_out = r_out;
      radius_in = r_in;
      latitude = lat*M_PI/180.0; //convert to radians
      longitude = lon*M_PI/180.0; //convert to radians
      density = den;
      zoa = z;
      index = n;
    }

    double radius_out;  ///< The outer radius of the bin in km
    double radius_in;  ///< The inner radius of the bin in km
    double latitude; ///< The latitude of the bin center in radians
    double longitude; ///< The longitude of the bin center in radians
    double density; ///< The density of the matter in the bin in g/cm^3
    double zoa;     ///< The effective Z/A value of the matter in the bin
    int index;      ///< An index to identify the matter type

  };

  class EarthModel : public TObject {

    public:

      EarthModel(std::string filename=""); ///< Constructor
      virtual ~EarthModel();          ///< Destructor

      virtual int FillPath(double cosT, double phi); ///< Fill the path sequence in a vector (phi in degrees)

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

      virtual void AddBin(double radius_out, double radius_in, double latitude, double longitude, double density,
                            double zoa,    double layer); ///< Add a bin to the model (angles in degrees)

      virtual void AddPath(double length, EarthBin bin);  ///< Add a path segment to the sequence

      virtual double DetDistForNextLatBin(int cur_index, int &dir, double dLat, double maxSinSqLat, int &max_reached, double beta, double gamma, double gammaSq, double rDetGammaSinDetLat, double rDetCosAcosDetLat, double rDetSinDetLat, double rDetCosT, double rDetSinT); ///< Calculate the detector distance at the edge of the current lat bin along the neutrino's trajectory

      virtual double DetDistForNextLonBin(double prev_lon, double dLon, int &lon_bin, double min_lon, double max_lon, double sinDetLon, double cosDetLon, double alpha, double sinTsinAsinDetLon, double sinTsinAcosDetLon, double rDetCosDetLat); ///< Calculate the detector distance at the edge of the current lon bin along the neutrino's trajectory and increment lon_bin

      virtual void RecordLatLonBinCrossings(double detDist_nextDbin, double &DetDist, int &index, int &latBin, int &nextLatBin, double &detDist_nextLatBin, int &sign, double dLat, double maxSinSqLat, int &maxlatreached, int &lonBin, int &nextLonBin, double &detDist_nextLonBin, double dLon, double min_lon, double max_lon, double alpha, double beta, double gamma, double gammaSq, double rDetGammaSinDetLat, double rDetCosAcosDetLat, double rDetSinDetLat, double rDetCosDetLat, double rDetCosT, double rDetSinT, double sinTsinAsinDetLon, double sinTsinAcosDetLon, double sinDetLon, double cosDetLon); ///< Record path segments for each latitude/longitude bin crossed before reaching detDist_nextDbin

      std::vector<OscProb::EarthBin> fEarthBins; ///< The layers in the earth model
      int fnDepthBins; ///< Total number of depth bins
      int fnLonBins; ///< Total number of longitude bins
      int fnLatBins; ///< Total number of latitude bins
      double fInvLonBinWidth; ///< 1/binwidth for each longitude bin
      double fInvLatBinWidth; ///< 1/binwidth for each latitude bin
      double fHalfLonBinWidth; ///< Half-width of each longitude bin
      double fHalfLatBinWidth; ///< Half-width of each latitude bin

      std::vector<OscProb::NuPath> fNuPath; ///< The current neutrino path sequence

      double fEarthRadius; ///< Biggest depth value in Earth model
      double fRadiusMax; ///< Minimum depth value (converted to radius) in Earth model

      double fDetRadius; ///< The radius where the detector lives
      double fDetLat; ///< The latitude (in rad) where the detector lives
      double fDetLon; ///< The longitude (in rad) where the detector lives

      bool fRemoveSmallPaths; ///< Tag whether to merge small paths

      static const double DET_TOL; ///< The detector position tolerance near boundaries

      //For Latitude Calculation Error Message
      std::string fErrorMessage_LonInfo; ///< Part of error message containing the longitude information
      int fLonError;

      // Required for saving in ROOT files
      ClassDef(EarthModel, 1);

  };

}
#endif
