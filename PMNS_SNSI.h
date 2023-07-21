////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_SNSI
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

#ifndef PMNS_SNSI_H
#define PMNS_SNSI_H

#include "PMNS_NSI.h"

namespace OscProb {
  class PMNS_SNSI : public PMNS_NSI {

  public:

    PMNS_SNSI();          ///< Constructor
    virtual ~PMNS_SNSI(); ///< Destructor

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
