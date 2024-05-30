///////////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_LIVS
///
/// \brief Implements oscillations with LIV as modelled by SME
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
/// \author Alba Domi - alba.domi\@fau.de
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_LIVS_H
#define PMNS_LIVS_H

#include "PMNS_Fast.h"
#include "PMNS_Base.h"

namespace OscProb {

  class PMNS_LIVS : public PMNS_Fast {
  public:
    PMNS_LIVS();          ///< Constructor
    virtual ~PMNS_LIVS(); ///< Destructor

    /// Set any given LIVS parameter of a chosen dimension.
    /// The flavour convention is:
    ///   e   -> 0
    ///   mu  -> 1
    ///   tau -> 2
    virtual void SetA(int flvi, int flvj, int coord, double val);
    virtual void SetC(int flvi, int flvj, int coord1, int coord2, double val);
    virtual void SetNeutrinoDirection(double zenith, double azimuth, double chi);
    virtual void SetTime(double hours);

    /// Get any given LIVS parameter of a chosen dimension.
    /// The flavour convention is:
    ///   e   -> 0
    ///   mu  -> 1
    ///   tau -> 2
    virtual double GetA(int flvi, int flvj, int coord);
    virtual double GetC(int flvi, int flvj, int coord1, int coord2);
    
    
  protected:
    /// Build the full Hamiltonian
    virtual void UpdateHam();
    virtual void SolveHam();
    virtual void FillCache() {} ///< Deactivate cache

    double fa[3][3][3]; ///< Stores each aT LIVS parameter of dimension 3
    double fc[3][3][3][3]; ///< Stores each cT LIVS parameter of dimension 4

    double N[3]; ///< neutrino directional factors NX,NY,NZ

    double zenith;
    double azimuth;
    double chi;
    
    double omega; // 25 h 56 min
    double T; // hours
    
    
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
