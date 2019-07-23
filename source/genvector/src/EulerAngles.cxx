// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 , LCG ROOT FNAL MathLib Team                    *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Implementation file for rotation in 3 dimensions, represented by EulerAngles
//
// Created by: Mark Fischler Thurs June 9  2005
//
// Last update: $Id$
//
#include "Math/GenVector/EulerAngles.h"

#include <cmath>

#include "Math/GenVector/EulerAngles.h"

namespace ROOT
{

namespace Math
{

// ========== I/O =====================

std::ostream &operator<<(std::ostream &os, const EulerAngles &e)
{
  // TODO - this will need changing for machine-readable issues
  //        and even the human readable form may need formatiing improvements
  os << "\n{phi: " << e.Phi() << "   theta: " << e.Theta()
     << "   psi: " << e.Psi() << "}\n";
  return os;
}

} // namespace Math
} // namespace ROOT
