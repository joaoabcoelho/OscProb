///////////////////////////////////////////////////////////////////////////////
/// class OscProb::PMNS_LIV
///
/// Implementation of neutrino oscillations in matter in a
/// three-neutrino framework with LIV as modelled by the SME.
/// The SME coefficients are included up to the 8th order.
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

  };

}
#endif
