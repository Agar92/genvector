// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 , LCG ROOT FNAL MathLib Team                    *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class Rotation in 3 dimensions, represented by 3x3 matrix
//
// Created by: Mark Fischler Tues July 5 2005
//
#include "Math/GenVector/Rotation3D.h"

#include <algorithm>
#include <cmath>

#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"

namespace ROOT
{

namespace Math
{

// ========== Constructors and Assignment =====================

std::ostream &operator<<(std::ostream &os, const Rotation3D &r)
{
  // TODO - this will need changing for machine-readable issues
  //        and even the human readable form needs formatting improvements
  double m[9];
  r.GetComponents(m, m + 9);
  os << "\n" << m[0] << "  " << m[1] << "  " << m[2];
  os << "\n" << m[3] << "  " << m[4] << "  " << m[5];
  os << "\n" << m[6] << "  " << m[7] << "  " << m[8] << "\n";
  return os;
}

} // namespace Math
} // namespace ROOT
