#ifndef COMPLEXSOLVER_H
#define COMPLEXSOLVER_H

#include <Eigen/Core>

#include "Definitions.h"

void complexsolver(const Eigen::Matrix3cd& A, OscProb::vectorD& w);

#endif
