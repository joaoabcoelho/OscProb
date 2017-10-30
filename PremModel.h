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
/// This class inherits from TObject and can be saved in ROOT files.
///
/// \author jcoelho\@apc.in2p3.fr
////////////////////////////////////////////////////////////////////////

#ifndef PREMMODEL_H
#define PREMMODEL_H

#include <string>
#include <vector>

#include "TObject.h"

#include "NuPath.h"

namespace OscProb {

  class PremModel : public TObject {

    public:
 
      PremModel(std::string filename=""); ///< Constructor
      virtual ~PremModel();          ///< Destructor
      
      virtual int FillPath(double cosT); ///< Fill the path sequence in a vector

      virtual std::vector<OscProb::NuPath> GetNuPath(); ///< Get the current neutrino path sequence

      virtual std::vector<OscProb::NuPath> GetMergedPaths(double prec = 0.25); ///< Get merged path sequence in a vector

      virtual double GetTotalL(double cosT); ///< Get the total baseline for a given cosTheta
      virtual double GetCosT(double L); ///< Get the cosTheta for a given total baseline
      
      virtual void SetLayerZoA(int layer, double zoa); ///< Set Z/A of all layers of a given type

      virtual void LoadModel(std::string filename); ///< Load an earth model from a file

      virtual void SetDetPos(double pos); ///< Set the detector position in km
      
      virtual std::vector<OscProb::PremLayer> GetPremLayers(); ///< Get the set of earth layers

      virtual OscProb::NuPath AvgPath(OscProb::NuPath p1, OscProb::NuPath p2); ///< Get the average of two paths
      
      virtual void SetRemoveSmallPaths(bool rp = true); ///< Set tag to remove small paths
      
    protected:

      virtual void ClearModel(); ///< Clear the earth model information

      virtual void AddLayer(double radius, double density,
                            double zoa,    double layer); ///< Add a layer to the model

      virtual void AddPath(double length, PremLayer pl);  ///< Add a path segment to the sequence
      
      virtual std::vector<OscProb::NuPath> MergePaths(std::vector<OscProb::NuPath> inputPath, int j, int k); ///< Merge paths j and k in vector

      std::vector<OscProb::PremLayer> fPremLayers; ///< The layers in the earth model

      std::vector<OscProb::NuPath> fNuPath; ///< The current neutrino path sequence

      int fDetLayer;  ///< The layer index of the detector
      double fDetPos; ///< The radius where the detector lives
      
      bool fRemoveSmallPaths; ///< Tag whether to merge small paths
      
      static const double DET_TOL; ///< The detector position tolerance near boundaries 
      
      // Required for saving in ROOT files
      ClassDef(PremModel, 1);

  };

}
#endif

