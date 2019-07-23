// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 , LCG ROOT FNAL MathLib Team                    *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Implementation file for rotation in 3 dimensions, represented by quaternion
//
// Created by: Mark Fischler Thurs June 9  2005
//
// Last update: $Id$
//
#include "Math/GenVector/Quaternion.h"

#include <cmath>

#include "Math/GenVector/Quaternion.h"

namespace ROOT
{

namespace Math
{

// ========== I/O =====================

std::ostream &operator<<(std::ostream &os, const Quaternion &q)
{
  // TODO - this will need changing for machine-readable issues
  //        and even the human readable form may need formatiing improvements
  os << "\n{" << q.U() << "   " << q.I() << "   " << q.J() << "   " << q.K()
     << "}\n";
  return os;
}

} // namespace Math
} // namespace ROOT
