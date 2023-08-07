///////////////////////////////////////////////////////////////////////////////
/// class OscProb::PMNS_LIV
///
/// Implementation of neutrino oscillations in matter in a
/// three-neutrino framework with LIV as modelled by the SME.
/// The SME coefficients are included up to the 8th order,
/// following the approach described in
/// https://doi.org/10.1103/PhysRevD.85.096005.
///
/// This developement is part of the QGRANT project with
/// ID: 101068013,
/// founded by the HORIZON-MSCA-2021-PF-01-01 programme.
///
/// \author Nafis R. K. Chowdhury - nrkhanchowdhury\@km3net.de
/// \author Joao Coelho - jcoelho\@apc.in2p3.fr
/// \author Alba Domi - alba.domi\@fau.de
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_LIV_H
#define PMNS_LIV_H

#include "PMNS_Fast.h"

namespace OscProb {

  class PMNS_LIV : public PMNS_Fast {

    public:

      PMNS_LIV();          ///< Constructor
      virtual ~PMNS_LIV(); ///< Destructor

      /// Set any given LIV parameter of a chosen dimension.
      /// The flavour convention is:
      ///   e   -> 0
      ///   mu  -> 1
      ///   tau -> 2
      virtual void SetaT(int flvi, int flvj, int dim, double val, double phase);
      virtual void SetcT(int flvi, int flvj, int dim, double val, double phase);

      /// Get any given LIV parameter of a chosen dimension.
      /// The flavour convention is:
      ///   e   -> 0
      ///   mu  -> 1
      ///   tau -> 2
      virtual complexD GetaT(int flvi, int flvj, int dim=3);
      virtual complexD GetcT(int flvi, int flvj, int dim=4);

    protected:

      /// Build the full Hamiltonian
      virtual void UpdateHam();

      complexD faT[3][3][3]; ///< Stores each aT LIV parameter of dimension 3,5,7
      complexD fcT[3][3][3]; ///< Stores each cT LIV parameter of dimension 4,6,8


    ///////////////////////////////////////////////////////////////////////////
    //
    // Obsolete functions for backward compatibility...
    //
    ///////////////////////////////////////////////////////////////////////////

    public:

      /// Set any given LIV parameter
      virtual void SetaT(int flvi, int flvj, double val, double phase);
      virtual void SetcT(int flvi, int flvj, double val, double phase);

      /// Set the LIV parameters all at once
      void SetLIV(double aT_ee,          double aT_mumu,         double aT_tautau,
                  double aT_emu,         double aT_etau,         double aT_mutau,
                  double cT_ee,          double CT_mumu,         double CT_tautau,
                  double cT_emu,         double cT_etau,         double cT_mutau,
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

  };

}

#endif

///////////////////////////////////////////////////////////////////////////////
