// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 ROOT FNAL MathLib Team                          *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for Boost
//
// Created by: Mark Fischler  Mon Nov 1  2005
//
// Last update: $Id$
//
#ifndef ROOT_Math_GenVector_Boost
#define ROOT_Math_GenVector_Boost 1

#include <array>

#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PxPyPzE4D.h"

#include "Math/GenVector/BoostX.h"
#include "Math/GenVector/BoostY.h"
#include "Math/GenVector/BoostZ.h"

namespace ROOT
{

namespace Math
{

//__________________________________________________________________________________________
/**
   Lorentz boost class with the (4D) transformation represented internally
   by a 4x4 orthosymplectic matrix.
   See also BoostX, BoostY and BoostZ for classes representing
   specialized Lorentz boosts.
   Also, the 3-D rotation classes can be considered to be special Lorentz
   transformations which do not mix space and time components.

   @ingroup GenVector

*/

class Boost
{

public:
  typedef double Scalar;

  enum ELorentzRotationMatrixIndex
  {
    kLXX = 0,
    kLXY = 1,
    kLXZ = 2,
    kLXT = 3,
    kLYX = 4,
    kLYY = 5,
    kLYZ = 6,
    kLYT = 7,
    kLZX = 8,
    kLZY = 9,
    kLZZ = 10,
    kLZT = 11,
    kLTX = 12,
    kLTY = 13,
    kLTZ = 14,
    kLTT = 15
  };

  enum EBoostMatrixIndex
  {
    kXX = 0,
    kXY = 1,
    kXZ = 2,
    kXT = 3,
    kYY = 4,
    kYZ = 5,
    kYT = 6,
    kZZ = 7,
    kZT = 8,
    kTT = 9
  };

  // ========== Constructors and Assignment =====================

  /**
      Default constructor (identity transformation)
  */
  Boost() { SetIdentity(); }

  /**
     Construct given a three Scalars beta_x, beta_y, and beta_z
   */
  Boost(Scalar beta_x, Scalar beta_y, Scalar beta_z)
  {
    SetComponents(beta_x, beta_y, beta_z);
  }

  /**
     Construct given a beta vector (which must have methods x(), y(), z())
   */
  template <class Avector> explicit Boost(const Avector &beta)
  {
    SetComponents(beta);
  }

  /**
     Construct given a pair of pointers or iterators defining the
     beginning and end of an array of three Scalars to use as beta_x, _y, and _z
   */
  template <class IT> Boost(IT begin, IT end) { SetComponents(begin, end); }

  /**
     copy constructor
  */
  Boost(Boost const &b) = default;
  /**
     move constructor
  */
  Boost(Boost &&b) = default;

  /**
     Construct from an axial boost
  */

  explicit Boost(BoostX const &bx) { SetComponents(bx.BetaVector()); }
  explicit Boost(BoostY const &by) { SetComponents(by.BetaVector()); }
  explicit Boost(BoostZ const &bz) { SetComponents(bz.BetaVector()); }

  // The compiler-generated copy ctor, copy assignment, and dtor are OK.

  /**
     Assignment operator
   */
  Boost &operator=(Boost const &rhs) = default;
  Boost &operator=(Boost &&rhs) = default;

  /**
     Assign from an axial pure boost
  */
  Boost &operator=(BoostX const &bx) { return operator=(Boost(bx)); }
  Boost &operator=(BoostY const &by) { return operator=(Boost(by)); }
  Boost &operator=(BoostZ const &bz) { return operator=(Boost(bz)); }

  /**
     Re-adjust components to eliminate small deviations from a perfect
     orthosyplectic matrix.
   */
  void Rectify();

  // ======== Components ==============

  /**
     Set components from beta_x, beta_y, and beta_z
  */
  inline void SetComponents(Scalar beta_x, Scalar beta_y, Scalar beta_z);

  /**
     Get components into beta_x, beta_y, and beta_z
  */
  void GetComponents(Scalar &beta_x, Scalar &beta_y, Scalar &beta_z) const;

  /**
     Set components from a beta vector
  */
  template <class Avector> inline void SetComponents(const Avector &beta)
  {
    SetComponents(beta.x(), beta.y(), beta.z());
  }

  /**
     Set given a pair of pointers or iterators defining the beginning and end of
     an array of three Scalars to use as beta_x,beta _y, and beta_z
   */
  template <class IT>
#ifndef NDEBUG
  void SetComponents(IT begin, IT end)
  {
#else
  void SetComponents(IT begin, IT)
  {
#endif
    IT a = begin;
    IT b = ++begin;
    IT c = ++begin;
    assert(++begin == end);
    SetComponents(*a, *b, *c);
  }

  /**
     Get given a pair of pointers or iterators defining the beginning and end of
     an array of three Scalars into which to place beta_x, beta_y, and beta_z
   */
  template <class IT>
#ifndef NDEBUG
  void GetComponents(IT begin, IT end) const
  {
#else
  void GetComponents(IT begin, IT) const
  {
#endif
    IT a = begin;
    IT b = ++begin;
    IT c = ++begin;
    assert(++begin == end);
    GetComponents(*a, *b, *c);
  }

  /**
     Get given a pointer or an iterator defining the beginning of
     an array into which to place beta_x, beta_y, and beta_z
   */
  template <class IT> void GetComponents(IT begin) const
  {
    double bx, by, bz = 0;
    GetComponents(bx, by, bz);
    *begin++ = bx;
    *begin++ = by;
    *begin   = bz;
  }

  /**
     The beta vector for this boost
   */
  typedef DisplacementVector3D<Cartesian3D<double>, DefaultCoordinateSystemTag>
      XYZVector;
  XYZVector BetaVector() const;

  /**
     Get elements of internal 4x4 symmetric representation, into a data
     array suitable for direct use as the components of a LorentzRotation
     Note -- 16 Scalars will be written into the array; if the array is not
     that large, then this will lead to undefined behavior.
  */
  void GetLorentzRotation(Scalar r[]) const;

  // =========== operations ==============

  /**
     Lorentz transformation operation on a Minkowski ('Cartesian')
     LorentzVector
  */
  LorentzVector<ROOT::Math::PxPyPzE4D<double>>
  operator()(const LorentzVector<ROOT::Math::PxPyPzE4D<double>> &v) const;

  /**
     Lorentz transformation operation on a LorentzVector in any
     coordinate system
   */
  template <class CoordSystem>
  LorentzVector<CoordSystem>
  operator()(const LorentzVector<CoordSystem> &v) const
  {
    LorentzVector<PxPyPzE4D<double>> xyzt(v);
    LorentzVector<PxPyPzE4D<double>> r_xyzt = operator()(xyzt);
    return LorentzVector<CoordSystem>(r_xyzt);
  }

  /**
     Lorentz transformation operation on an arbitrary 4-vector v.
     Preconditions:  v must implement methods x(), y(), z(), and t()
     and the arbitrary vector type must have a constructor taking (x,y,z,t)
   */
  template <class Foreign4Vector>
  Foreign4Vector operator()(const Foreign4Vector &v) const
  {
    LorentzVector<PxPyPzE4D<double>> xyzt(v);
    LorentzVector<PxPyPzE4D<double>> r_xyzt = operator()(xyzt);
    return Foreign4Vector(r_xyzt.X(), r_xyzt.Y(), r_xyzt.Z(), r_xyzt.T());
  }

  /**
     Overload operator * for boost on a vector
   */
  template <class A4Vector> inline A4Vector operator*(const A4Vector &v) const
  {
    return operator()(v);
  }

  /**
      Invert a Boost in place
   */
  void Invert();

  /**
      Return inverse of  a boost
   */
  Boost Inverse() const;

  /**
     Equality/inequality operators
   */
  bool operator==(const Boost &rhs) const
  {
    for (unsigned int i = 0; i < 10; ++i)
    {
      if (fM[i] != rhs.fM[i])
        return false;
    }
    return true;
  }
  bool operator!=(const Boost &rhs) const { return !operator==(rhs); }

protected:
  void SetIdentity();

private:
  std::array<Scalar, 10> fM;

}; // Boost

// ============ Class Boost ends here ============

/**
   Stream Output and Input
 */
// TODO - I/O should be put in the manipulator form

std::ostream &operator<<(std::ostream &os, const Boost &b);

// ============ Declarations end here ============

inline void Boost::SetIdentity()
{
  // set identity boost
  fM[kXX] = 1.0;
  fM[kXY] = 0.0;
  fM[kXZ] = 0.0;
  fM[kXT] = 0.0;
  fM[kYY] = 1.0;
  fM[kYZ] = 0.0;
  fM[kYT] = 0.0;
  fM[kZZ] = 1.0;
  fM[kZT] = 0.0;
  fM[kTT] = 1.0;
}

inline void Boost::SetComponents(Scalar bx, Scalar by, Scalar bz)
{
  // set the boost beta as 3 components
  Scalar const bp2 = bx * bx + by * by + bz * bz;
  if (bp2 >= 1)
  {
    GenVector::Throw("Beta Vector supplied to set Boost represents speed >= c");
    // SetIdentity();
    return;
  }
  Scalar const gamma  = 1.0 / sqrt(1.0 - bp2);
  Scalar const bgamma = gamma * gamma / (1.0 + gamma);
  fM[kXX]             = 1.0 + bgamma * bx * bx;
  fM[kYY]             = 1.0 + bgamma * by * by;
  fM[kZZ]             = 1.0 + bgamma * bz * bz;
  fM[kXY]             = bgamma * bx * by;
  fM[kXZ]             = bgamma * bx * bz;
  fM[kYZ]             = bgamma * by * bz;
  fM[kXT]             = gamma * bx;
  fM[kYT]             = gamma * by;
  fM[kZT]             = gamma * bz;
  fM[kTT]             = gamma;
}

inline void Boost::GetComponents(Scalar &bx, Scalar &by, Scalar &bz) const
{
  // get beta of the boots as 3 components
  Scalar gaminv = 1.0 / fM[kTT];
  bx            = fM[kXT] * gaminv;
  by            = fM[kYT] * gaminv;
  bz            = fM[kZT] * gaminv;
}

inline DisplacementVector3D<Cartesian3D<Boost::Scalar>>
Boost::BetaVector() const
{
  // get boost beta vector
  Scalar gaminv = 1.0 / fM[kTT];
  return DisplacementVector3D<Cartesian3D<Scalar>>(
      fM[kXT] * gaminv, fM[kYT] * gaminv, fM[kZT] * gaminv);
}

inline void Boost::GetLorentzRotation(Scalar r[]) const
{
  // get Lorentz rotation corresponding to this boost as an array of 16 values
  r[kLXX] = fM[kXX];
  r[kLXY] = fM[kXY];
  r[kLXZ] = fM[kXZ];
  r[kLXT] = fM[kXT];
  r[kLYX] = fM[kXY];
  r[kLYY] = fM[kYY];
  r[kLYZ] = fM[kYZ];
  r[kLYT] = fM[kYT];
  r[kLZX] = fM[kXZ];
  r[kLZY] = fM[kYZ];
  r[kLZZ] = fM[kZZ];
  r[kLZT] = fM[kZT];
  r[kLTX] = fM[kXT];
  r[kLTY] = fM[kYT];
  r[kLTZ] = fM[kZT];
  r[kLTT] = fM[kTT];
}

inline void Boost::Rectify()
{
  // Assuming the representation of this is close to a true Lorentz Rotation,
  // but may have drifted due to round-off error from many operations,
  // this forms an "exact" orthosymplectic matrix for the Lorentz Rotation
  // again.

  if (fM[kTT] <= 0)
  {
    GenVector::Throw("Attempt to rectify a boost with non-positive gamma");
    return;
  }
  DisplacementVector3D<Cartesian3D<Scalar>> beta(fM[kXT], fM[kYT], fM[kZT]);
  beta /= fM[kTT];
  if (beta.mag2() >= 1)
  {
    beta /= (beta.R() * (1.0 + 1.0e-16));
  }
  SetComponents(beta);
}

inline LorentzVector<PxPyPzE4D<double>> Boost::
operator()(const LorentzVector<PxPyPzE4D<double>> &v) const
{
  // apply bosost to a PxPyPzE LorentzVector
  Scalar x = v.Px();
  Scalar y = v.Py();
  Scalar z = v.Pz();
  Scalar t = v.E();
  return LorentzVector<PxPyPzE4D<double>>(
      fM[kXX] * x + fM[kXY] * y + fM[kXZ] * z + fM[kXT] * t,
      fM[kXY] * x + fM[kYY] * y + fM[kYZ] * z + fM[kYT] * t,
      fM[kXZ] * x + fM[kYZ] * y + fM[kZZ] * z + fM[kZT] * t,
      fM[kXT] * x + fM[kYT] * y + fM[kZT] * z + fM[kTT] * t);
}

inline void Boost::Invert()
{
  // invert in place boost (modifying the object)
  fM[kXT] = -fM[kXT];
  fM[kYT] = -fM[kYT];
  fM[kZT] = -fM[kZT];
}

inline Boost Boost::Inverse() const
{
  // return inverse of boost
  Boost tmp(*this);
  tmp.Invert();
  return tmp;
}

} // namespace Math
} // namespace ROOT

#endif /* ROOT_Math_GenVector_Boost  */
