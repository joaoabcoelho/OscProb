#ifndef COMPLEXSOLVER_H
#define COMPLEXSOLVER_H

#include "Eigenvalues"

#include "Definitions.h"

void complexsolver(OscProb::matrixC  A,
                   OscProb::matrixC& Q,
                   OscProb::matrixC& Qinv,
                   OscProb::vectorC& w);

#endif
