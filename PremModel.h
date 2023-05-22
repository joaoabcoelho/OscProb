////////////////////////////////////////////////////////////////////////
/// \class OscProb::PremModel
///
/// \brief Implements an earth model with spherical shells
///
/// This class implements a spherically symmetric model of the earth
/// using PremLayer's to store spherical shells with different properties.
///
/// The class is then able to produce path sequences through the different
/// earth layers as a function of the cosine of the zenith angle with
/// respect to the detector.
///
/// The detector can be positioned at any radius within the model and
/// the path sequences will take into account the fact that some layers
/// are above the detector.
///
/// By default this implements the model stored in PremTables/prem_default.txt
/// with the detector at the bottom of the ocean layer (radius = 6368 km).
///
/// This class inherits from EarthModelBase and can be saved in ROOT files.
///
/// \author jcoelho\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

#ifndef PREMMODEL_H
#define PREMMODEL_H

#include <string>

#include "EarthModelBase.h"
//#include "TObject.h"

#include "Definitions.h"
#include "NuPath.h"

namespace OscProb {

  class PremModel : public EarthModelBase /*: public TObject*/ {

    public:

      PremModel(std::string filename=""); ///< Constructor
      virtual ~PremModel();          ///< Destructor

      void SetDetPos(double rad, double lat=0, double lon=0); ///< Set the detector position (rad = radius in km, lat/lon in deg)

      int FillPath(double cosT, double phi=0); ///< Fill the path sequence in a vector

      virtual void LoadModel(std::string filename); ///< Load an earth model from a file

      virtual std::vector<PremLayer> GetPremLayers(); ///< Get the set of earth layers

      virtual void SetLayerZoA(int layer, double zoa); ///< Set Z/A of all layers of a given type
      virtual double GetLayerZoA(int layer); ///< Get Z/A of all layers of a given type

      virtual void SetTopLayerSize(double thick);      ///< Set the outermost layer thickness in km

    protected:

      virtual void ClearModel(); ///< Clear the earth model information

      virtual void AddLayer(double radius, double density,
                            double zoa,    double layer); ///< Add a layer to the model

      virtual void AddPath(double length, PremLayer pl);  ///< Add a path segment to the sequence

      std::vector<PremLayer> fPremLayers; ///< The layers in the earth model

      int fDetLayer;  ///< The layer index of the detector

      static const double DET_TOL; ///< The detector position tolerance near boundaries

      // Required for saving in ROOT files
      ClassDef(PremModel, 1);

  };

}
#endif
