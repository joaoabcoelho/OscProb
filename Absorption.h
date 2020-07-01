//#ifndef Absorption
//#define Absorption

#include "NuPath.h"
#include <vector>
#include <numeric>
#include <iterator>
#include <iostream>

namespace OscProb {

  class Absorption{
    public:
    Absorption();
    virtual ~Absorption();

    //virtual double Trans(int flvi, double E);
    virtual double Trans(double xsec);
    
    virtual void SetPath(std::vector<OscProb::NuPath> paths);  ///< Set a path sequence
    protected:
    std::vector<OscProb::NuPath> fNuPaths; ///< Vector of neutrino paths
    // Set the neutrino path
    //virtual void SetPath(OscProb::NuPath p);              ///< Set a single path
    //virtual void SetPath(double length, double density,
    //                          double zoa=0.5,    int layer=0); ///< Set a single path
    
    //virtual void SetPath(std::vector<OscProb::NuPath> paths);  ///< Set a path sequence

    static const double kNA; //Avogardo constant
  };

}
