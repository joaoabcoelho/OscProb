#ifndef PMNS_NSIs_IO_H
#define PMNS_NSIs_IO_H

#include "PMNS_Fast.h"

namespace OscProb {
  class PMNS_NSIs_IO : public PMNS_Fast {

  public:

    PMNS_NSIs_IO();          ///< Constructor
    virtual ~PMNS_NSIs_IO(); ///< Destructor

    /// Set any given NSI parameter
    virtual void SetEta(int flvi, int flvj, double val, double phase);

    /// Get any given NSI parameter
    virtual complexD GetEta(int flvi, int flvj);

    /// Set the NSI parameters all at once    
    void SetNSIs(double eta_ee,      double eta_emu,      double eta_etau,
                double eta_mumu,    double eta_mutau,    double eta_tautau,
                double delta_emu=0, double delta_etau=0, double delta_mutau=0);

    // Set the diagonal real NSI parameters
    virtual void SetEta_ee    (double a); ///< Set eta_ee parameter
    virtual void SetEta_mumu  (double a); ///< Set eta_mumu parameter
    virtual void SetEta_tautau(double a); ///< Set eta_tautau parameter

    // Set the off-diagonal complex NSI parameters
    virtual void SetEta_emu  (double a, double phi); ///< Set eta_emu parameter
    virtual void SetEta_etau (double a, double phi); ///< Set eta_etau parameter
    virtual void SetEta_mutau(double a, double phi); ///< Set eta_mutau parameter

    // Set lowest neutrino mass
    

    virtual void SetM3(double m3); ///< Set m3
    
    virtual void SetLowestMass(double m3); ///< Set lightest neutrino mass
    
    // Get lowest neutrino mass
    virtual double GetM3(); ///< Get m3
    
    
  protected:

    /// Build the full Hamiltonian
    virtual void UpdateHam();
    
    complexD fEta[3][3]; ///< Stores each NSI parameter
    
    double fM3; ///< Lowest neutrinom mass

  };

}
#endif
////////////////////////////////////////////////////////////////////////

