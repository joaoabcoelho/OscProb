////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_ScalarNSI
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework with scalar NSI. 
///
/// This class reimplements the PMNS_NSI class to a scalar NSI scenario.
///
/// The matter potential is parametrized by epsilon parameters with 
/// dimensions of mass in eV.
///
/// \sa PMNS_NSI
///
/// \author jcoelho\@apc.in2p3.fr and urahaman\@km3net.de
////////////////////////////////////////////////////////////////////////

#ifndef PMNS_ScalarNSI_H
#define PMNS_ScalarNSI_H

#include "PMNS_NSI.h"

namespace OscProb {
  class PMNS_ScalarNSI : public PMNS_NSI {

  public:

    PMNS_ScalarNSI();          ///< Constructor
    virtual ~PMNS_ScalarNSI(); ///< Destructor

    virtual void SetLowestMass(double m); ///< Set lightest neutrino mass
    virtual double GetLowestMass(); ///< Get lightest neutrino mass
    
  protected:

    /// Build the full Hamiltonian
    virtual void UpdateHam();
    virtual void BuildHms();

    double fM; ///< Lightest neutrino mass

  };

}
#endif
////////////////////////////////////////////////////////////////////////
