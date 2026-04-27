///////////////////////////////////////////////////////////////////////////////
/// \class OscProb::PMNS_SiderealLIV
///
/// \brief Implements oscillations with Sidereal LIV as modelled by SME
///
/// Implementation of neutrino oscillations in matter in a
/// three-neutrino framework with LIV as modelled by the SME.
/// Implements the minimal SME sector (aT, mass dimension 3, and
/// cT, mass dimension 4), following the approach described in
/// https://doi.org/10.1103/PhysRevD.85.096005.
///
/// This developement is part of the QGRANT project with
/// ID: 101068013,
/// founded by the HORIZON-MSCA-2021-PF-01-01 programme.
///
/// \author Alba Domi - alba.domi\@fau.de
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_SiderealLIV_H
#define PMNS_SiderealLIV_H

#include "PMNS_Fast.h"
#include "PMNS_Base.h"

namespace OscProb {

  class PMNS_SiderealLIV : public PMNS_Fast {
  public:
    PMNS_SiderealLIV();          ///< Constructor
    virtual ~PMNS_SiderealLIV(); ///< Destructor

    /// Set any given LIVS parameter of a chosen dimension.
    /// The flavour convention is:
    ///   e   -> 0
    ///   mu  -> 1
    ///   tau -> 2
    virtual void SetA(int flvi, int flvj, int coord, double val);
    virtual void SetC(int flvi, int flvj, int coord1, int coord2, double val);

    /// Set the neutrino arrival direction and detector colatitude.
    /// All angles are in radians.
    virtual void SetNeutrinoDirection(double zenith, double azimuth,
                                      double colatitude);

    /// Set the local sidereal time at which the probability is evaluated.
    /// @param hours - Local sidereal time, in hours.
    virtual void SetTimeHours(double hours);

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

    /// Validate and reorder a flavour pair so that flvi <= flvj.
    /// Emits a warning on cerr if reordering or rejection occurs.
    /// @param flvi  - first flavour index (modified in place)
    /// @param flvj  - second flavour index (modified in place)
    /// @param name  - name of the parameter (used in warning text, e.g. "aT")
    /// @return true if (flvi, flvj) is a valid pair, false otherwise.
    static bool ValidateAndOrder(int& flvi, int& flvj, const char* name);

    /// Validate that a Lorentz spatial index is in [0, 2].
    /// @param coord     - the index to check
    /// @param coordName - human-readable name for the warning ("coord", "coord1", ...)
    /// @param paramName - parameter name for the warning ("aT", "cT", ...)
    /// @return true if in range, false otherwise.
    static bool ValidateCoord(int coord, const char* coordName,
                              const char* paramName);

    // ------------------------------------------------------------------------
    // Sidereal-frequency constants
    // ------------------------------------------------------------------------
    /// Length of one sidereal day, in hours (IERS/IAU conventional value).
    static constexpr double kSiderealDayHours = 23.9344696;

    /// Angular sidereal frequency, in rad/h.
    static constexpr double kOmegaSidereal = 2.0 * M_PI / kSiderealDayHours;

    // ------------------------------------------------------------------------
    // SME coefficients (real; only the upper triangle in flavour indices is used)
    // ------------------------------------------------------------------------
    double fa[3][3][3];    ///< aT LIVS parameters (mass dimension 3)
    double fc[3][3][3][3]; ///< cT LIVS parameters (mass dimension 4)

    // ------------------------------------------------------------------------
    // Geometry / time state
    // ------------------------------------------------------------------------
    double fN[3];   ///< Neutrino directional factors (NX, NY, NZ)
    double fTime;   ///< Local sidereal time, in hours
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
