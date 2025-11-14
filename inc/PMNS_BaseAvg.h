#include <Eigen/Core>
#include "PMNS_Base.h"

#ifndef PMNS_BASEAVG_H
#define PMNS_BASEAVG_H

namespace OscProb {

  class PMNS_BaseAvg : public PMNS_Base {
    public:
      PMNS_BaseAvg(int numNus); ///< Constructor
      virtual void SetPremModel(OscProb::PremModel& prem);
    protected:
      virtual void InitializeTaylorsVectors(); ///< Initialize all member
                                               ///< vectors with zeros

      virtual void SetwidthBin(double dE,
                               double dcosT); ///< Set bin's widths for
                                              ///< energy and angle
      
      virtual void SetCosT(double cosT); ///< Set neutrino angle.

     // Construction of the K matrices
      virtual void BuildKE(
          double L); ///< build K matrix for the inverse of energy in mass basis

      /// Currently missing BuildKCosT due to that Updating the Ham, which is made 
      /// a level up
          
      virtual double LnDerivative(); ///< Compute the derivation of one layer's
                                     ///< length
      //////////
      ///
      ///

      double fcosT;  ///<  Cosine of neutrino angle

      double fdInvE; ///< Bin's width for the inverse of energy in GeV-1
      double fdcosT; ///< Bin's width for angle

      vectorD flambdaInvE;     ///< Eigenvectors of K_invE
      vectorD flambdaCosT; ///< Eigenvectors of K_cosTheta
      matrixC fVInvE;          ///< Eigenvalues of K_invE

      Eigen::MatrixXcd fKInvE; ///< K matrix for the inverse of energy in GeV
                               ///< for the entire path
      Eigen::MatrixXcd fKcosT;      ///< K matrix for neutrino angle for the entire path

      matrixC fVcosT;      ///< Eigenvalues of K_cosTheta
      matrixC fevolutionMatrixS; ///< Evolution matrix S for reference energy
                                 ///< and angle for the entire path
      matrixC fSflavor; ///< S matrix for one layer
      matrixC fKmass;   ///< K matrix in mass basis for one layer
      matrixC fKflavor; ///< K matrix in flavor basis for one layer

      // Variables for the compute of the derivation of one layer's length
      int    fLayer;
      int    fdl;
      double fDetRadius;

      // Copy of the earth model used
      OscProb::PremModel fPrem;

      // Variables for the compute of the derivation of one layer's length
      double fminRsq;


  };

} // namespace OscProb

#endif

