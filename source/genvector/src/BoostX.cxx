// @(#)root/mathcore:$Id$
// Authors:  M. Fischler  2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 , LCG ROOT FNAL MathLib Team                    *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class BoostX, a 4x4 symmetric matrix representation of
// an axial Lorentz transformation
//
// Created by: Mark Fischler Mon Nov 1  2005
//
#include "Math/GenVector/BoostX.h"
#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/GenVector_exception.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PxPyPzE4D.h"

#include <algorithm>
#include <cmath>

namespace ROOT {

namespace Math {

// ========== I/O =====================

std::ostream &operator<<(std::ostream &os, const BoostX &b) {
  os << " BoostX( beta: " << b.Beta() << ", gamma: " << b.Gamma() << " ) ";
  return os;
}

} // namespace Math
} // namespace ROOT
