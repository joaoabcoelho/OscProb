///////////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_Base
///
/// \brief Base class implementing general functions for computing
///        neutrino oscillations.
///
/// This is an abstract class implementing the general functions needed
/// for setting up an oscillation calculator. The method for solving the
/// eigensystem for the Hamiltonian must be defined in the derived classes.
///
/// \sa PMNS_Fast PMNS_NSI PMNS_Sterile PMNS_SNSI
///     PMNS_Decay PMNS_Deco PMNS_LIV PMNS_Iter
///
/// \author Joao Coelho - jcoelho\@apc.in2p3.fr
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_BASE_H
#define PMNS_BASE_H

#include <unordered_set>

#include "EigenPoint.h"

namespace OscProb {

  class PMNS_Base {
    public:
      PMNS_Base(int numNus = 3); ///< Constructor
      virtual ~PMNS_Base();      ///< Destructor

      // Get the oscillation probability
      virtual double Prob(
          vectorC nu_in,
          int     flvf); ///< Compute the probability of nu_in going to flvf
      virtual double Prob(vectorC nu_in, int flvf,
                          double E); ///< Compute the probability of nu_in going
                                     ///< to flvf for energy E
      virtual double Prob(vectorC nu_in, int flvf, double E,
                          double L); ///< Compute the probability of nu_in going
                                     ///< to flvf for energy E and distance L

      // Get the oscillation probability
      virtual double Prob(
          int flvi,
          int flvf); ///< Compute the probability of flvi going to flvf
      virtual double Prob(int flvi, int flvf,
                          double E); ///< Compute the probability of flvi going
                                     ///< to flvf for energy E
      virtual double Prob(int flvi, int flvf, double E,
                          double L); ///< Compute the probability of flvi going
                                     ///< to flvf for energy E and distance L

      // Get the oscillation probability vector
      virtual vectorD ProbVector(
          vectorC nu_in); ///< Compute the probabilities of nu_in going to all
                          ///< flavours
      virtual vectorD ProbVector(
          vectorC nu_in,
          double  E); ///< Compute the probabilities of nu_in going to all
                     ///< flavours for energy E
      virtual vectorD ProbVector(
          vectorC nu_in, double E,
          double L); ///< Compute the probabilities of nu_in going to all
                     ///< flavours for energy E and distance L

      // Get the oscillation probability vector
      virtual vectorD ProbVector(int flvi); ///< Compute the probabilities of
                                            ///< flvi going to all flavours
      virtual vectorD ProbVector(
          int flvi, double E); ///< Compute the probabilities of flvi going to
                               ///< all flavours for energy E
      virtual vectorD ProbVector(
          int flvi, double E,
          double L); ///< Compute the probabilities of flvi going to all
                     ///< flavours for energy E and distance L

      // Get the oscillation probability matrix
      virtual matrixD ProbMatrix(int nflvi,
                                 int nflvf); ///< Compute the probability matrix
      virtual matrixD ProbMatrix(
          int nflvi, int nflvf,
          double E); ///< Compute the probability matrix for energy E
      virtual matrixD ProbMatrix(int nflvi, int nflvf, double E,
                                 double L); ///< Compute the probability matrix
                                            ///< for energy E and distance L

      // Get probability averaged over a bin
      virtual double AvgProb(
          vectorC nu_in, int flvf, double E,
          double dE =
              0); ///< Compute the average probability over a bin of energy
      virtual double AvgProbLoE(
          vectorC nu_in, int flvf, double LoE,
          double dLoE =
              0); ///< Compute the average probability over a bin of L/E

      // Get probability averaged over a bin
      virtual double AvgProb(
          int flvi, int flvf, double E,
          double dE =
              0); ///< Compute the average probability over a bin of energy
      virtual double AvgProbLoE(
          int flvi, int flvf, double LoE,
          double dLoE =
              0); ///< Compute the average probability over a bin of L/E

      // Get probability vector averaged over a bin
      virtual vectorD AvgProbVector(
          vectorC nu_in, double E,
          double dE = 0); ///< Compute the average probability vector over a bin
                          ///< of energy
      virtual vectorD AvgProbVectorLoE(
          vectorC nu_in, double LoE,
          double dLoE =
              0); ///< Compute the average probability vector over a bin of L/E

      // Get probability vector averaged over a bin
      virtual vectorD AvgProbVector(
          int flvi, double E,
          double dE = 0); ///< Compute the average probability vector over a bin
                          ///< of energy
      virtual vectorD AvgProbVectorLoE(
          int flvi, double LoE,
          double dLoE =
              0); ///< Compute the average probability vector over a bin of L/E

      // Get probability vector averaged over a bin
      virtual matrixD AvgProbMatrix(
          int nflvi, int nflvf, double E,
          double dE = 0); ///< Compute the average probability matrix over a bin
                          ///< of energy
      virtual matrixD AvgProbMatrixLoE(
          int nflvi, int nflvf, double LoE,
          double dLoE =
              0); ///< Compute the average probability matrix over a bin of L/E

      virtual vectorC GetMassEigenstate(
          int mi); ///< Get a neutrino mass eigenstate

      // Set the oscillation parameters
      virtual void SetAngle(int i, int j,
                            double th); ///< Set the mixing angle theta_ij
      virtual void SetDelta(int i, int j,
                            double delta); ///< Set the CP phase delta_ij
      virtual void SetDm(int    j,
                         double dm); ///< Set the mass-splitting dm_j1 in eV^2

      // Get the oscillation parameters
      virtual double GetAngle(int i, int j); ///< Get the mixing angle theta_ij
      virtual double GetDelta(int i, int j); ///< Get the CP phase delta_ij
      virtual double GetDm(int j); ///< Get the mass-splitting dm_j1 in eV^2

      // Get the effective oscillation parameters
      virtual double GetDmEff(
          int j); ///< Get the effective mass-splitting dm_j1 in eV^2

      // Set default oscillation parameters
      virtual void SetStdPars(); ///< Set PDG 3-flavor parameters

      // Set energy and anti-neutrino flag
      virtual void SetEnergy(double E);      ///< Set the neutrino energy in GeV
      virtual void SetIsNuBar(bool isNuBar); ///< Set the anti-neutrino flag

      // Get energy and anti-neutrino flag
      virtual double GetEnergy();  ///< Get the neutrino energy in GeV
      virtual bool   GetIsNuBar(); ///< Get the anti-neutrino flag

      // Set the neutrino path
      virtual void SetPath(NuPath p); ///< Set a single path
      virtual void SetPath(double length, double density, double zoa = 0.5,
                           int layer = 0); ///< Set a single path

      virtual void SetPath(std::vector<NuPath> paths); ///< Set a path sequence

      virtual void AddPath(NuPath p); ///< Add a path to the sequence
      virtual void AddPath(double length, double density, double zoa = 0.5,
                           int layer = 0); ///< Add a path to the sequence

      virtual void ClearPath(); ///< Clear the path vector

      virtual void SetLength(double L); ///< Set a single path lentgh in km
      virtual void SetDensity(
          double rho);                 ///< Set single path density in g/cm^3
      virtual void SetZoA(double zoa); ///< Set Z/A value for single path

      virtual void SetLength(vectorD L);    ///< Set multiple path lengths
      virtual void SetDensity(vectorD rho); ///< Set multiple path densities
      virtual void SetZoA(vectorD zoa);     ///< Set multiple path Z/A values
      virtual void SetLayers(
          std::vector<int> lay); ///< Set multiple path layer indices

      // Set a default neutrino path
      virtual void SetStdPath(); ///< Set standard neutrino path

      // Get the neutrino path
      virtual std::vector<NuPath> GetPath(); ///< Get the neutrino path sequence

      /// Compute the sample points for a bin of L/E with width dLoE
      virtual vectorD GetSamplePoints(double LoE, double dLoE);

      // Setup the caching system
      virtual void SetUseCache(bool u = true); ///< Set caching on/off
      virtual void ClearCache();               ///< Clear the cache
      virtual void SetMaxCache(int mc = 1e6);  ///< Set max cache size

      /// Set the AvgProb precision
      virtual void SetAvgProbPrec(double prec);

      /// Use the one-loop correction to Fermi constant
      virtual void SetUseOneLoopGF(bool u = true);

      /// Set Fermi constant
      virtual void SetGF(double Gf = kGmu);

    protected:
      // Some useful complex numbers
      static const complexD zero; ///< zero in complex
      static const complexD one;  ///< one in complex

      // Unit conversion constants
      static const double kKm2eV;  ///< km to eV^-1
      static const double kK2;     ///< mol/GeV^2/cm^3 to eV
      static const double kGeV2eV; ///< GeV to eV
      static const double kNA;     ///< Avogadro constant

      static const double kGmu; ///< G_mu in units of GeV^-2
      double              kGf;  ///< User value of G_F

      virtual void InitializeVectors(); ///< Initialize all member vectors with
                                        ///< zeros

      // Internal caching functions
      virtual bool TryCache();  ///< Try to find a cached eigensystem
      virtual void FillCache(); ///< Cache the current eigensystem

      // Auxiliary path functions
      virtual void SetCurPath(
          NuPath p); ///< Set the path currently in use by the class

      virtual void SetAtt(double att,
                          int    idx); ///< Set one of the path attributes
      virtual void SetAtt(vectorD att,
                          int     idx); ///< Set all values of a path attribute

      // Building and solving
      virtual void RotateH(
          int i, int j,
          matrixC& Ham); ///< Rotate the Hamiltonian by theta_ij and delta_ij
      virtual void RotateState(
          int i, int j); ///< Rotate the neutrino state by theta_ij and delta_ij

      virtual void BuildHms(); ///< Build the matrix of masses squared.

      ///
      /// Solve the full Hamiltonian
      ///
      /// Not implemented in base class.
      ///
      virtual void SolveHam() = 0; ///< Solve the full Hamiltonian for
                                   ///< eigenvectors and eigenvalues

      // Resetting and propagating
      virtual void ResetToFlavour(
          int flv); ///< Reset neutrino state to pure flavour flv
      virtual void SetPureState(
          vectorC nu_in); ///< Set the initial state from a pure state

      virtual void PropagatePath(
          NuPath p);            ///< Propagate neutrino through a single path
      virtual void Propagate(); ///< Propagate neutrino through full path

      virtual double P(
          int flv); ///< Return the probability of final state in flavour flv
      virtual vectorD GetProbVector(); ///< Return vector of probabilities from
                                       ///< final state

      virtual std::vector<int> GetSortedIndices(
          const vectorD x); ///< Get indices that sort a vector

      virtual vectorD ConvertEtoLoE(double E, double dE);

      // Attributes

      int fNumNus; ///< Number of neutrino flavours

      vectorD fDm;    ///< m^2_i - m^2_1 in vacuum
      matrixD fTheta; ///< theta[i][j] mixing angle
      matrixD fDelta; ///< delta[i][j] CP violating phase

      vectorC fNuState; ///< The neutrino current state
      matrixC fHms;     ///< matrix H*2E in eV^2

      vectorC fPhases; ///< Buffer for oscillation phases
      vectorC fBuffer; ///< Buffer for neutrino state tranformations

      vectorD fEval; ///< Eigenvalues of the Hamiltonian
      matrixC fEvec; ///< Eigenvectors of the Hamiltonian

      double fEnergy;  ///< Neutrino energy
      bool   fIsNuBar; ///< Anti-neutrino flag

      std::vector<NuPath> fNuPaths; ///< Vector of neutrino paths
      NuPath              fPath;    ///< Current neutrino path

      bool fBuiltHms; ///< Tag to avoid rebuilding Hms
      bool fGotES;    ///< Tag to avoid recalculating eigensystem

      bool   fUseCache;  ///< Flag for whether to use caching
      double fCachePrec; ///< Precision of cache matching
      int    fMaxCache;  ///< Maximum cache size

      double fAvgProbPrec; ///< AvgProb precision

      std::unordered_set<EigenPoint> fMixCache; ///< Caching set of eigensystems
      EigenPoint                     fProbe;    ///< EigenpPoint to try
  };

  /// An index sorting comparator
  struct IdxCompare {
      /// Take in a target vector
      IdxCompare(const vectorD& target) : target(target) {}

      /// Compare elements a and b of target vector
      bool operator()(int a, int b) const { return target[a] < target[b]; }

    private:
      /// Attribute to store the target vector
      const vectorD target;
  };

} // namespace OscProb

#endif
