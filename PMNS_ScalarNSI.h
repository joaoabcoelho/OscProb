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

#ifndef PMNS_ScalarNSI_H
#define PMNS_ScalarNSI_H

#include "PMNS_Fast.h"

namespace OscProb {
  class PMNS_ScalarNSI : public PMNS_Fast {

  public:

    PMNS_ScalarNSI();          ///< Constructor
    virtual ~PMNS_ScalarNSI(); ///< Destructor

    /// Set any given NSI parameter
    virtual void SetEta(int flvi, int flvj, double val, double phase);

    /// Get any given NSI parameter
    virtual complexD GetEta(int flvi, int flvj);

    /// Set the NSI parameters all at once    
    void SetScalarNSI(double eta_ee,      double eta_emu,      double eta_etau,
                double eta_mumu,    double eta_mutau,    double eta_tautau,
                double delta_emu=0, double delta_etau=0, double delta_mutau=0);

    // Set the diagonal real NSI parameters
    virtual void SetEta_ee    (double a); ///< Set eta_ee parameter
    virtual void SetEta_mumu  (double a); ///< Set eta_mumu parameter
    virtual void SetEta_tautau(double a); ///< Set eta_tautau parameter

    // Set the off-diagonal complex NSI parameters
    virtual void SetEta_emu  (double a, double phi=0); ///< Set eta_emu parameter
    virtual void SetEta_etau (double a, double phi=0); ///< Set eta_etau parameter
    virtual void SetEta_mutau(double a, double phi=0); ///< Set eta_mutau parameter

    // Set lowest neutrino mass
    

    virtual void SetM(double m); ///< Set m3
    
   // virtual void SetLowestMass(double m); ///< Set lightest neutrino mass
    
    // Get lowest neutrino mass
    virtual double GetM(); ///< Get m3
    
    
  protected:

    /// Build the full Hamiltonian
    virtual void UpdateHam();
    
    complexD fEta[3][3]; ///< Stores each NSI parameter
    
    double fM; ///< Lowest neutrinom mass

  };

}
#endif
////////////////////////////////////////////////////////////////////////
