////////////////////////////////////////////////////////////////////////
/// \class OscProb::EarthModelBase
///
/// \brief Base class for implementing an earth model
///
/// This abstract class provides the structure to able to produce path
/// sequences through the Earth as a function of the azimuthal angle
/// and cosine of the zenith angle with respect to the detector.
///
/// The detector can be positioned at any point within the Earth (except
/// at the exact center), and the path sequences will take into account
/// the fact that some layers are above the detector.
///
/// This class inherits from TObject and can be saved in ROOT files.
///
/// \author rpestes\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

#ifndef EARTHMODELBASE_H
#define EARTHMODELBASE_H

#include "TObject.h"
#include "NuPath.h"

namespace OscProb {

  class EarthModelBase : public TObject {

    public:

      ///
      /// Set the coordinates of the detector:
      ///   radius in km, latitude in degrees, longitude in degrees
      ///
      /// Not implemented in base class.
      ///
      virtual void SetDetPos(double rad, double lat=0, double lon=0) = 0; ///< Set the detector position (rad = radius in km, lat/lon in deg)

      ///
      /// Construct the neutrino path through the Earth having
      /// zenith angle defined by cosT and azimuthal angle phi.
      ///
      /// Not implemented in base class.
      ///
      virtual int FillPath(double cosT, double phi) = 0; ///< Fill the path sequence in a vector (phi in degrees)

      virtual std::vector<NuPath> GetNuPath(); ///< Get the current neutrino path sequence

      virtual std::vector<NuPath> GetMergedPaths(double prec = 0.25); ///< Get merged path sequence in a vector

      virtual double GetTotalL(double cosT); ///< Get the total baseline for a given cosTheta
      virtual double GetCosT(double L); ///< Get the cosTheta for a given total baseline

      virtual void SetRemoveSmallPaths(bool rp = true); ///< Set tag to remove small paths

    protected:

      virtual void SetDetectorCoordinates(double rad, double lat, double lon); ///< Set the coordinates of the detector (rad = radius in km, lat/lon in deg)

      virtual void AddPathSegment(double length, double density, double zoa, int index);  ///< Add a path segment to the sequence

      std::vector<NuPath> fNuPath; ///< The current neutrino path sequence

      double fRadiusMax; ///< Maximum radius in Earth model (in km)
      double fDetRadius; ///< The radius where the detector sits (in km)
      double fDetLat; ///< The latitude (in rad) where the detector sits
      double fDetLon; ///< The longitude (in rad) where the detector sits

      bool fRemoveSmallPaths; ///< Tag whether to merge small paths

      // Required for saving in ROOT files
      ClassDef(EarthModelBase, 1);

  };

}
#endif
