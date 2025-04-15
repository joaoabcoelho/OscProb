///////////////////////////////////////////////////////////////////////////////
/// mettre info et biblio
///////////////////////////////////////////////////////////////////////////////

#ifndef PMNS_TaylorExp_H
#define PMNS_TaylorExp_H

#include "PMNS_Fast.h"

namespace OscProb {

  class PMNS_TaylorExp : public PMNS_Fast {
    public:
      PMNS_TaylorExp();          ///< Constructor
      virtual ~PMNS_TaylorExp(); ///< Destructor

    protected:
      double fGamma[3]; ///< Stores each decoherence parameter
  };

} // namespace OscProb

#endif

///////////////////////////////////////////////////////////////////////////////
