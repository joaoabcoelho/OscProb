///////////////////////////////////////////////////////////////////////////////
/// class OscProb::PMNS_LIV
///
/// Brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework with LIV.
///
///\author Joao Coelho - coelho\@lal.in2p3.fr and
/// \colaborator Nafis R. K. Chowdhury - nrkhanchowdhury\@km3net.de 
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_LIV_H
#define PMNS_LIV_H

#include "PMNS_Fast.h"

namespace OscProb {
  class PMNS_LIV : public PMNS_Fast {

  public:

    PMNS_LIV();          ///< Constructor
    virtual ~PMNS_LIV(); ///< Destructor

    /// Set any given LIV parameter
    virtual void SetaT(int flvi, int flvj, double val, double phase);
    virtual void SetcT(int flvi, int flvj, double val, double phase);

    /// Set the LIV parameters all at once
    void SetLIV(double aT_ee,     double aT_mumu,      double aT_tautau,
            double aT_emu,     double aT_etau,      double aT_mutau,
            double cT_ee,     double CT_mumu,      double CT_tautau,
            double cT_emu,     double cT_etau,      double cT_mutau,
            double delta_aT_emu=0, double delta_aT_etau=0, double delta_aT_mutau=0,
	    double delta_cT_emu=0, double delta_cT_etau=0, double delta_cT_mutau=0);

    //Set diagonal LIV pars
    virtual void SetaT_ee    (double a); ///< Set eps_ee parameter
    virtual void SetaT_mumu  (double a); ///< Set eps_mumu parameter
    virtual void SetaT_tautau(double a); ///< Set eps_tautau parameter
    
    virtual void SetcT_ee    (double a); ///< Set eps_ee parameter
    virtual void SetcT_mumu  (double a); ///< Set eps_mumu parameter
    virtual void SetcT_tautau(double a); ///< Set eps_tautau parameter


    /// Set diagonal LIV pars
    virtual void SetaT_emu  (double a, double phi);
    virtual void SetaT_etau (double a, double phi);
    virtual void SetaT_mutau(double a, double phi);


    virtual void SetcT_emu  (double a, double phi);
    virtual void SetcT_etau (double a, double phi);
    virtual void SetcT_mutau(double a, double phi);

    /// Get any given LIV par
    virtual complexD GetaT(int flvi, int flvj);
    virtual complexD GetcT(int flvi, int flvj);

  protected:

    /// Build the full Hamiltonian
    virtual void UpdateHam();

    complexD faT[3][3]; ///< Stores each LIV1 parameter
    complexD fcT[3][3]; ///< Stores each LIV1 parameter

  };

}
#endif
////////////////////////////////////////////////////////////////////////
