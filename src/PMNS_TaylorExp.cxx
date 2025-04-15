///////////////////////////////////////////////////////////////////////////////
//info
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "PMNS_TaylorExp.h"

using namespace OscProb;

using namespace std;

//.............................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
/// This class is restricted to 3 neutrino flavours.
///
PMNS_TaylorExp::PMNS_TaylorExp() : PMNS_Fast(), fGamma() {}

//.............................................................................
///
/// Nothing to clean.
///
PMNS_TaylorExp::~PMNS_TaylorExp() {}