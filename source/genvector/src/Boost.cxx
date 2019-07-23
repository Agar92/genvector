// @(#)root/mathcore:$Id$
// Authors:  M. Fischler  2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 , LCG ROOT FNAL MathLib Team                    *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class Boost, a 4x4 symmetric matrix representation of
// an axial Lorentz transformation
//
// Created by: Mark Fischler Mon Nov 1  2005
//
#include "Math/GenVector/Boost.h"
#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/GenVector_exception.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PxPyPzE4D.h"

#include <algorithm>
#include <cmath>

//#ifdef TEX
/**

   A variable names bgamma appears in several places in this file. A few
   words of elaboration are needed to make its meaning clear.  On page 69
   of Misner, Thorne and Wheeler, (Exercise 2.7) the elements of the matrix
   for a general Lorentz boost are given as

   \f[   \Lambda^{j'}_k = \Lambda^{k'}_j
              = (\gamma - 1) n^j n^k + \delta^{jk}  \f]

   where the n^i are unit vectors in the direction of the three spatial
   axes.  Using the definitions, \f$ n^i = \beta_i/\beta \f$ , then, for
   example,

   \f[   \Lambda_{xy} = (\gamma - 1) n_x n_y
              = (\gamma - 1) \beta_x \beta_y/\beta^2  \f]

   By definition, \f[   \gamma^2 = 1/(1 - \beta^2)  \f]

   so that   \f[   \gamma^2 \beta^2 = \gamma^2 - 1  \f]

   or   \f[   \beta^2 = (\gamma^2 - 1)/\gamma^2  \f]

   If we insert this into the expression for \f$ \Lambda_{xy} \f$, we get

   \f[   \Lambda_{xy} = (\gamma - 1) \gamma^2/(\gamma^2 - 1) \beta_x \beta_y \f]

   or, finally

   \f[   \Lambda_{xy} = \gamma^2/(\gamma+1) \beta_x \beta_y  \f]

   The expression \f$ \gamma^2/(\gamma+1) \f$ is what we call <em>bgamma</em> in
   the code below.

   \class ROOT::Math::Boost
*/
//#endif

namespace ROOT {

namespace Math {

// ========== I/O =====================

std::ostream &operator<<(std::ostream &os, const Boost &b) {
  // TODO - this will need changing for machine-readable issues
  //        and even the human readable form needs formatiing improvements
  double m[16];
  b.GetLorentzRotation(m);
  os << "\n" << m[0] << "  " << m[1] << "  " << m[2] << "  " << m[3];
  os << "\n"
     << "\t"
     << "  " << m[5] << "  " << m[6] << "  " << m[7];
  os << "\n"
     << "\t"
     << "  "
     << "\t"
     << "  " << m[10] << "  " << m[11];
  os << "\n"
     << "\t"
     << "  "
     << "\t"
     << "  "
     << "\t"
     << "  " << m[15] << "\n";
  return os;
}

} // namespace Math
} // namespace ROOT
