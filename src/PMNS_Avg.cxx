///////////////////////////////////////////////////////////////////////////////
//
// Implementation of oscillations of neutrinos in matter in a
// three-neutrino framework with a first order Taylor expansion.
//
// This  class inherits from the PMNS_Fast class
//
///////////////////////////////////////////////////////////////////////////////

#include "PMNS_Avg.h"
#include <algorithm>
#include <iostream>
#include <Eigen/Eigenvalues>

using namespace OscProb;

using namespace std;

//.............................................................................
///
/// Constructor. \sa PMNS_Base::PMNS_Base
///
///
///
/// This bit would copy over to PMNS_MyClass, as the myOsc = ROOT.OscProb.PMNS_MyClass(4)
PMNS_Avg::PMNS_Avg(int numNus) : PMNS_Sterile(numNus)
{}
//.............................................................................
///
/// Nothing to clean.
///
PMNS_Avg::~PMNS_Avg() {}


