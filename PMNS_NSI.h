////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_NSI
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework with NSI. 
///
/// This class expands the PMNS_Fast class including a general
/// matter potential matrix describing Non-Standard Interactions (NSI).
///
/// The matter potential is parametrized by dimensionless quantities
/// epsilon which quantify the intensity of the NSI with respect to the
/// standard matter effects from coherent forward scattering with electrons. 
///
/// \sa PMNS_Fast
///
/// \author coelho\@lal.in2p3.fr
////////////////////////////////////////////////////////////////////////

#ifndef PMNS_NSI_H
#define PMNS_NSI_H

#include "PMNS_Fast.h"

namespace OscProb {
  class PMNS_NSI : public PMNS_Fast {

  public:

    PMNS_NSI();          ///< Constructor
    virtual ~PMNS_NSI(); ///< Destructor

    /// Set any given NSI parameter
    virtual void SetEps(int flvi, int flvj, double val, double phase);

    /// Get any given NSI parameter
    virtual complexD GetEps(int flvi, int flvj);

    /// Set the NSI parameters all at once    
    void SetNSI(double eps_ee,      double eps_emu,      double eps_etau,
                double eps_mumu,    double eps_mutau,    double eps_tautau,
                double delta_emu=0, double delta_etau=0, double delta_mutau=0);

    // Set the diagonal real NSI parameters
    virtual void SetEps_ee    (double a); ///< Set eps_ee parameter
    virtual void SetEps_mumu  (double a); ///< Set eps_mumu parameter
    virtual void SetEps_tautau(double a); ///< Set eps_tautau parameter

    // Set the off-diagonal complex NSI parameters
    virtual void SetEps_emu  (double a, double phi); ///< Set eps_emu parameter
    virtual void SetEps_etau (double a, double phi); ///< Set eps_etau parameter
    virtual void SetEps_mutau(double a, double phi); ///< Set eps_mutau parameter

    // Set relative NSI couplings
    virtual void SetElecCoup(double e); ///< Set electron coupling
    virtual void SetUpCoup(double u);   ///< Set u-quark couling
    virtual void SetDownCoup(double d); ///< Set d-quark couling
    virtual void SetFermCoup(double e, double u, double d); ///< Set all fermion couplings
    virtual void SetCoupByIndex(double c, int i); ///< Set a given fermion coupling

    // Get relative NSI couplings
    virtual double GetElecCoup(); ///< Get electron coupling
    virtual double GetUpCoup();   ///< Get u-quark couling
    virtual double GetDownCoup(); ///< Get d-quark couling
    virtual double GetZoACoup();  ///< Get effective Z/A coupling

  protected:

    /// Build the full Hamiltonian
    virtual void UpdateHam();
    
    complexD fEps[3][3]; ///< Stores each NSI parameter
    
    double fNSIcoup[3]; ///< Relative NSI couplings

  };

}
#endif
////////////////////////////////////////////////////////////////////////
