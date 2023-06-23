////////////////////////////////////////////////////////////////////////
/// \class OscProb::EarthModelBinned
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
/// Using indices to specify various regions (which must have a constant
/// Z/A value), the Z/A value for a region can be changed, and the density
/// throughout the region can be scaled by a constant factor.
///
/// By default this implements the model stored in
/// EarthTables/earth_binned_default.txt
/// with the detector at the bottom of the ocean layer (radius = 6368 km)
/// where the prime meridian intersects the equator.
///
/// This class inherits from EarthModelBase and can be saved in ROOT files.
///
/// \author rpestes\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
/// \struct OscProb::EarthBins
///
/// \brief A struct representing a 3D bin (depth, latitude, and longitude)
/// of matter for asymmetric Earth models
///
/// This struct stores the properties of a 3D bin (part of a spherical
/// shell) to be used by the EarthModel class in order to build an
/// earth model.
///
/// \author rpestes\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
/// \struct OscProb::TrajConstants
///
/// \brief A struct holding useful combinations of trajectory variables
///
/// This struct holds combinations of variables that are useful when
/// doing the calculations needed for FindPath within EarthModelBinned.
///
/// \author rpestes\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
/// \struct OscProb::EarthModelBinned::LatBinInfo
///
/// \brief A struct holding information about upcoming latitude bin
/// crossings along the neutrino's trajectory
///
/// \author rpestes\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
/// \struct OscProb::EarthModelBinned::LonBinInfo
///
/// \brief A struct holding information about upcoming longitude bin
/// crossings along the neutrino's trajectory
///
/// \author rpestes\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

#ifndef EARTHMODELBINNED_H
#define EARTHMODELBINNED_H

#include <string>
#include <vector>
#include <math.h>

#include "EarthModelBase.h"
//#include "TObject.h"
#include "NuPath.h"

namespace OscProb {

  struct TrajConstants
  {
    ///
    /// \brief Constructor.
    ///
    /// Constructor.
    ///
    /// By default, it sets cosT, Az, latD, and lonD equal to 0, and it uses
    /// 6368km for rD.
    ///
    /// @param cosT - Cosine of the zenith angle for the neutrino trajectory
    /// @param phi - The azimuthal angle for the neutrino trajectory (in rad)
    /// @param DetLat - The latitude of the detector (in rad)
    /// @param DetLon - The longitude of the detector (in rad)
    /// @param rDet - The distance from the center of the Earth to the detector (in km)
    ///
    TrajConstants(double cosT=0, double phi=0, double DetLat=0, double DetLon=0, double rDet=6368) {
      UpdateNuAngles(cosT, phi);
      UpdateDetPos(DetLat, DetLon, rDet);
    }

    void UpdateNuAngles(double cosTheta, double phi); ///< Update values of zenith angle and azimuthal angle for neutrino trajectory
    void UpdateDetPos(double rDet, double DetLat, double DetLon); ///< Update detector position for neutrino trajectory calculations
    void Recalculate(); ///< Calculate constants that use combinations of detector position variables and neutrino direction angles

    //Variables defining trajectory
    double cosT; ///< cosT
    double cosA; ///< cos(phi)
    double sinA; ///< sin(phi)
    double sinDetLat; ///< sin(DetLat)
    double cosDetLon; ///< cos(DetLon)
    double sinDetLon; ///< sin(DetLon)
    double DetRadius; ///< rDet

    //Derived Variables Calculated by UpdateNuAngles()
    double sinSqT; ///< sin^2(T) = 1 - (cosT)^2
    double sinT; ///< sin(T) = sqrt(sinSqT)
    double sinTsinA; ///< sin(T)*sin(phi)
    double sinTcosA; ///< sin(T)*cos(phi)

    //Derived Variables Calculated by UpdateDetPos()
    double cosDetLat; ///< cos(DetLat)
    double rDetSinDetLat; ///< rDet*sin(DetLat)
    double rDetCosDetLat; ///< rDet*cos(DetLat)

    //Derived Variables Calculated by Recalculate()
    double sinTsinAsinDetLon; ///< sin(T)*sin(phi)*sin(DetLon)
    double sinTsinAcosDetLon; ///< sin(T)*sin(phi)*cos(DetLon)
    double cosTcosDetLat; ///< cosT*cos(DetLat)
    double rDetCosT; ///< rDet*cosT
    double rDetSinT; ///< rDet*sin(T)
    double rDetCosAcosDetLat; ///< rDet*cos(phi)*cos(DetLat)
    double alpha; ///< sin(T)*cos(phi)*sin(DetLat)-cosT*cos(DetLat)
    double beta; ///< sin(T)*sin(DetLat)-cos(phi)*cosT*cos(DetLat)
    double gamma; ///< sin(T)*cos(phi)*cos(DetLat)+cosT*sin(DetLat)
    double gammaSq; ///< [sin(T)*cos(phi)*cos(DetLat)+cosT*sin(DetLat)]^2
    double rDetGammaSinDetLat; ///< rDet*sin(DetLat)*[sin(T)*cos(phi)*cos(DetLat)+cosT*sin(DetLat)]
    double maxSinSqLat; ///< 1 - [sin(phi)*cos(DetLat)]^2

  };

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
    /// @param n - Region index
    ///
    EarthBin(double r_out=0, double r_in=0, double lat=0, double lon=0, double den=0, double z=0.5, int n=0){
      SetBin(r_out, r_in, lat, lon, den, z, n);
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
    /// @param n - Region index
    ///
    void SetBin(double r_out=0, double r_in=0, double lat=0, double lon=0, double den=0, double z=0.5, int n=0){
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
    int index;      ///< Region index

  };

  class EarthModelBinned : public EarthModelBase {

    public:

      EarthModelBinned(std::string filename=""); ///< Constructor
      virtual ~EarthModelBinned();          ///< Destructor

      void SetDetPos(double rad, double lat=0, double lon=0); ///< Set the detector position (rad = radius in km, lat/lon in deg)

      int FillPath(double cosT, double phi=0); ///< Fill the path sequence in a vector (phi in degrees)

      virtual void LoadModel(std::string filename); ///< Load an earth model from a file

      virtual std::vector<OscProb::EarthBin> GetEarthBins(); ///< Get the set of earth layers

      virtual void SetRegionZoA(int index, double zoa); ///< Set Z/A of all bins with specified region index
      virtual double GetRegionZoA(int index); ///< Get Z/A of all bins with specified region index
      virtual void ScaleRegionDensity(int index, double scalingfactor); ///< Set Z/A of all bins with specified region index

    protected:

      struct LatBinInfo
      {
        int bin; ///< Index of current latitude bin
        int nextBin; ///< Index of next latitude bin
        double detDist_nextBin; ///< Distance along the neutrino trajectory from the edge of next latitude bin to the detector
        int sign; ///< Indicates whether latitude is increasing or decreasing with respect to decreasing distance from the detector (+1 => inc, -1 => dec)
        double dLat; ///< Change in lat from bin center to next bin (excludes direction; 0 if no more than 1 lat bin change)
        int maxreached; ///< Indicates whether the latitude function transition is still to come (0 => yes, 1 => no)
      };

      struct LonBinInfo
      {
        int bin; ///< Index of current longitude bin
        int nextBin; ///< Index of next longitude bin
        double detDist_nextBin; ///< Distance along the neutrino trajectory from the edge of next longitude bin to the detector
        double dLon; ///< Change in lon from bin center to next bin (includes direction; 0 if no more than 1 lon bin change)
        double min; ///< "Minimum" longitude
        double max; ///< "Maximum" longitude
        int error = 0; ///< Indicates if an error has been detected (0 => no, -1 => yes)
        std::string err_message; ///< Part of error message specific to piece of path
      };

      virtual void ClearModel(); ///< Clear the earth model information

      virtual void AddBin(double radius_out, double radius_in, double latitude, double longitude, double density, double zoa, double layer); ///< Add a bin to the model (angles in degrees)

      virtual void AddPath(double length, EarthBin bin);  ///< Add a path segment to the sequence

      virtual double DetDistForNextLatBin(int cur_index, LatBinInfo &L); ///< Calculate the detector distance at the edge of the current lat bin along the neutrino's trajectory

      virtual double DetDistForNextLonBin(double prev_lon, LonBinInfo &L); ///< Calculate the detector distance at the edge of the current lon bin along the neutrino's trajectory and increment L.bin

      virtual void RecordLatLonBinCrossings(double detDist_nextDbin, double &DetDist, int &index, LatBinInfo &latI, LonBinInfo &lonI); ///< Record path segments for each latitude/longitude bin crossed before reaching detDist_nextDbin

      virtual int LonBinIndex(double longitude); ///< Find lon bin index containing longitude

      OscProb::TrajConstants fC; ///< Useful constants for the trajectory

      std::vector<OscProb::EarthBin> fEarthBins; ///< The bins in the earth model
      int fnDepthBins; ///< Total number of depth bins
      int fnLonBins; ///< Total number of longitude bins
      int fnLatBins; ///< Total number of latitude bins
      double fInvLonBinWidth; ///< 1/binwidth for each longitude bin
      double fInvLatBinWidth; ///< 1/binwidth for each latitude bin
      double fHalfLonBinWidth; ///< Half-width of each longitude bin
      double fHalfLatBinWidth; ///< Half-width of each latitude bin

      int fmaxRegIndex; ///< Largest region index in model

      //For Latitude Calculation Error Message
      std::string fErrorMessage_LonInfo; ///< Part of error message containing the longitude information
      int fLonError;

      // Required for saving in ROOT files
      ClassDef(EarthModelBinned, 1);

  };

}
#endif
