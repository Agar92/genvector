// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 , LCG ROOT FNAL MathLib Team                    *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class LorentzRotation, a 4x4 matrix representation of
// a general Lorentz transformation
//
// Created by: Mark Fischler Mon Aug 8  2005
//

#include "Math/GenVector/GenVectorIO.h"

#include "Math/GenVector/GenVector_exception.h"
#include "Math/GenVector/LorentzRotation.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PxPyPzE4D.h"

#include <algorithm>
#include <cmath>

namespace ROOT
{

namespace Math
{

std::ostream &operator<<(std::ostream &os, const LorentzRotation &r)
{
  // TODO - this will need changing for machine-readable issues
  //        and even the human readable form needs formatiing improvements
  double m[16];
  r.GetComponents(m, m + 16);
  os << "\n" << m[0] << "  " << m[1] << "  " << m[2] << "  " << m[3];
  os << "\n" << m[4] << "  " << m[5] << "  " << m[6] << "  " << m[7];
  os << "\n" << m[8] << "  " << m[9] << "  " << m[10] << "  " << m[11];
  os << "\n"
     << m[12] << "  " << m[13] << "  " << m[14] << "  " << m[15] << "\n";
  return os;
}

} // namespace Math
} // namespace ROOT
