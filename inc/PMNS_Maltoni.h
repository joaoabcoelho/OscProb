#include <Eigen/Core>
#include "PMNS_Base.h"

#ifndef PMNS_BASEAVG_H
#define PMNS_BASEAVG_H

namespace OscProb {

  class PMNS_Maltoni : public PMNS_Base {
    public:
      PMNS_Maltoni(int numNus); ///< Constructor
      virtual void SetPremModel(OscProb::PremModel& prem);

      // Get probability averaged over a bin
      using PMNS_Base::AvgProb;
      using PMNS_Base::AvgProbLoE;
      // Get probability averaged over a bin
      virtual double AvgProb(
          int flvi, int flvf, double E,
          double dE); ///< Compute the average probability over
                      ///< a bin of energy with a Taylor expansion

     virtual double AvgProb(
          int flvi, int flvf, double E, double cosT,
          double dcosT); ///< Compute the average probability over a
                         ///< bin of cosTheta with a Taylor expansion
			 ///
      virtual double AvgProb(
          int flvi, int flvf, double E, double dE,
                             double cosT,
                             double dcosT); ///< Compute the average probability
                                            ///< over a bin of cosTheta and
                                            ///< energy with a Taylor expansion

      virtual double AvgProbLoE(
              int flvi, int flvf, double LoE,
              double dLoE); ///< Compute the average probability over
                        ///< a bin of LoE with a Taylor expansion
     virtual double AvgProbLoE(
		  int flvi, int flvf, double LoE, double dLoE, double cosT,
		  double dcosT); ///< Compute the average probability over a
				 ///< bin of cosTheta and LoE with a Taylor
				 ///< expansion


	      // Get probability vector averaged over a bin
      virtual vectorD AvgProbVector(
		  int flvi, double E,
		  double dE); ///< Compute the average probability vector over a bin
			      ///< of energy using a Taylor expansion
      virtual vectorD AvgProbVectorLoE(
		  int flvi, double LoE,
		  double dLoE); ///< Compute the average probability vector over a
				///< bin of L/E using a Taylor expansion
	       
      virtual matrixD AvgProbMatrix(int nflvi, int nflvf, double E, double dE);
      virtual matrixD AvgProbMatrixLoE(int nflvi, int nflvf, double LoE,
					       double dLoE);

      virtual double ExtrapolationProb(
             int flvi, int flvf, double E,
               double dE); ///< Compute the probability of flvi going to
                      ///< flvf for an energy E+dE

      virtual double ExtrapolationProbLoE(
               int flvi, int flvf, double LoE,
               double dLoE); ///< Compute the probability of flvi going to
                        ///< flvf at LoE+dLoE

      virtual double ExtrapolationProbCosT(
              int flvi, int flvf, double cosT,
              double dcosT); ///< Compute the probability of flvi going to
                         ///< flvf for an angle cosT+dcosT

      virtual void UseOscProbAverage(bool AverageFlag); // Set flag for which averaging to use


	    protected:

          virtual vectorD GetSamplePointsAvgClass(
		  double LoE,
		  double dLoE); ///< Compute the sample points for
				///< a bin of L/E with width dLoE
	
          virtual vectorD GetSamplePointsAvgClass(
                double E, double cosT,
                double dcosT); ///< Compute the sample points for
                         ///< a bin of cosT with width dcosT


              virtual matrixC GetSamplePointsAvgClass(
                  double InvE, double dInvE, double cosT,
                  double dcosT); ///< Compute the sample
                         ///< points for a bin
                         ///< of E and cosT with
                         ///< width dE and dcosT


	      /// Build the full Hamiltonian
	      virtual void UpdateHam() = 0;


	      virtual void InitializeTaylorsVectors(); ///< Initialize all member
						       ///< vectors with zeros

	      virtual void SetwidthBin(double dE,
				       double dcosT); ///< Set bin's widths for
						      ///< energy and angle
	      
	      virtual void SetCosT(double cosT); ///< Set neutrino angle.

	     // Construction of the K matrices
	      virtual void BuildKE(
		  double L); ///< build K matrix for the inverse of energy in mass basis

	      virtual void BuildKcosT(
		  double L); ///< build K matrix for angle in flavor basis
	 
	      virtual double LnDerivative(); ///< Compute the derivation of one layer's
					     ///< length

	      // Propagating
	      virtual void PropagatePathTaylor(
		  NuPath p); ///< Propagate neutrino through a single path
	      virtual void PropagateTaylor(); ///< Propagate neutrino through full path

	      // Rotation from mass to flavor basis
	      virtual void rotateS(); ///< Rotate the S matrix
	      virtual void rotateK(); ///< Rotate one K matrix

	      // Multiplication rule between two pairs of (S,K)
	      virtual void MultiplicationRuleS(); ///< Product between two S matrices
	      virtual void MultiplicationRuleK(          
		  Eigen::MatrixXcd& K); ///< Product between two K matrices

	      /// Solve the K matrix
	      void SolveK(Eigen::MatrixXcd& K, vectorD& lambda,
			  matrixC& V); ///< Solve the K matrix for
				       ///< eigenvectors and eigenvalues
	      template <typename T>
	      void TemplateSolver(Eigen::MatrixXcd& K, vectorD& lambda, matrixC& V);

	    // Propagating
		  virtual double AvgFormula(
		  int flvi, int flvf, double dbin, vectorD flambda,
		  matrixC fV); ///< Formula for the average probability over a bin
			       ///< of width dbin

               virtual double AvgFormulaExtrapolation(
                   int flvi, int flvf, double dbin, vectorD flambda,
                   matrixC fV); ///< Formula for the extrapolation of probability


	      virtual double AvgAlgo(
		  int flvi, int flvf, double LoE, double dLoE,
		  double L); ///< Algorithm for the compute of the average
			     ///< probability over a bin of LoE

               virtual double AvgAlgo(
                int flvi, int flvf, double InvE, double dInvE, double cosT,
                  double dcosT); ///< Algorithm for the compute of the average
                         ///< probability over a bin of 1oE and cosT


	      virtual double AvgAlgoCosT(
          int flvi, int flvf, double E, double cosT,
          double dcosT); ///< Algorithm for the compute of the average
                         ///< probability over a bin of cosT


           virtual double AlgorithmDensityMatrix(
          int flvi, int flvf); ///< Algorithm for the transformations on the
                               ///< density matrix

      virtual void RotateDensityM(
          bool to_basis, matrixC V); ///< Apply rotation to the density matrix

      virtual void HadamardProduct(vectorD lambda,
                                   double  dbin); ///< Apply an Hadamard Product
                                                 ///< to the density matrix


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
      
      Eigen::MatrixXcd fHam; ///< The full Hamiltonian
			     ///
      bool gAverageFlag = false;   /// Flag to call OscProb default or Maltoni average

  };

} // namespace OscProb

#endif

