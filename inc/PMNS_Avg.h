///////////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_Avg
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework with a first order Taylor expansion.
///
/// This class expands the PMNS_Fast class including the use of a first order
/// Taylor expansion to calculate the average on bins faster.
///
/// The model assumes a first order expansion over neutrino energy and angle
/// for both dynamical variables at the same time or for only one.
///
/// This is the first version of this class. A second version will be release
/// with a better implementation with the other classes.
///
/// Reference: https://doi.org/10.48550/arXiv.2308.00037
///
/// \sa PMNS_Fast
///
/// \author jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_Avg_H
#define PMNS_Avg_H

#include "PMNS_Sterile.h"

namespace OscProb {

  class PMNS_Avg : public PMNS_Sterile {
    public:
      PMNS_Avg(int numNus); ///< Constructor
      virtual ~PMNS_Avg();  ///< Destructor
      // Get probability averaged over a bin
      using PMNS_Base::AvgProb;
      using PMNS_Base::AvgProbLoE;
      using PMNS_Base::AvgProbVector;
      using PMNS_Base::AvgProbVectorLoE;

      // Get probability
    protected:
      //virtual double AvgFormula(
      //    int flvi, int flvf, double dbin, vectorD flambda,
       //   matrixC fV); ///< Formula for the average probability over a bin
       //                ///< of width dbin

     // virtual double AvgAlgo(
     //     int flvi, int flvf, double LoE, double dLoE,
     //     double L); ///< Algorithm for the compute of the average
     //               ///< probability over a bin of LoE



      // Attributes

                                 ///< and angle for the entire path



  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
