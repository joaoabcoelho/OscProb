///////////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_LIV1
///
/// \brief Implementation of oscillations of neutrinos in matter in a
///        three-neutrino framework with LIV1.
///
/// This class expands the PMNS_Fast class including an additional
/// Hamiltonian part matrix describing Lorentz Invariance Violation (LIV1).
///
/// The additional part is parametrized by six quantities (aT, cTT).
/// see paper arxiv: arXiv:1410.4267v2 [hep-ex] for details of the model.
///
/// \sa PMNS_Fast
///
/// \author M.L. Abdelali
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_LIV1_H
#define PMNS_LIV1_H

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


    /// Set any given LIV1 parameter
    virtual void SetaT_emu  (double a, double phi);
    virtual void SetaT_etau (double a, double phi);
    virtual void SetaT_mutau(double a, double phi);

    virtual void SetcT_ee    (double a); ///< Set eps_ee parameter
    virtual void SetcT_mumu  (double a); ///< Set eps_mumu parameter
    virtual void SetcT_tautau(double a); ///< Set eps_tautau parameter

    virtual void SetcT_emu  (double a, double phi);
    virtual void SetcT_etau (double a, double phi);
    virtual void SetcT_mutau(double a, double phi);

    /// Get any given LIV1 parameter
    virtual complex GetaT(int flvi, int flvj);
    virtual complex GetcT(int flvi, int flvj);

  protected:

    /// Build the full Hamiltonian
    virtual void UpdateHam();

    complex faT[3][3]; ///< Stores each LIV1 parameter
	complex fcT[3][3]; ///< Stores each LIV1 parameter

  };

}
#endif
////////////////////////////////////////////////////////////////////////
