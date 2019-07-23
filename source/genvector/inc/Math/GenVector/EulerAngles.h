// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class EulerAngles
//
// Created by: Lorenzo Moneta  at Tue May 10 17:55:10 2005
//
// Last update: Tue May 10 17:55:10 2005
//
#ifndef ROOT_Math_GenVector_EulerAngles
#define ROOT_Math_GenVector_EulerAngles 1

#include "Math/GenVector/Allfwd.h"
#include "Math/Math.h"


#include <algorithm>
#include <cassert>
#include <ostream>

namespace ROOT
{
namespace Math
{

//__________________________________________________________________________________________
/**
   EulerAngles class describing rotation as three angles (Euler Angles).
   The Euler angles definition matches that of Classical Mechanics (Goldstein).
   It is also the same convention defined in
   <A HREF="http://mathworld.wolfram.com/EulerAngles.html">mathworld</A>
   and used in Mathematica and CLHEP. Note that the ROOT class TRotation defines
   a slightly different convention.

   @ingroup GenVector
*/
class EulerAngles
{

public:
  typedef double Scalar;

  /**
     Default constructor
  */
  EulerAngles() : fPhi(0.0), fTheta(0.0), fPsi(0.0) {}

  /**
     Constructor from phi, theta and psi
  */
  EulerAngles(Scalar phi, Scalar theta, Scalar psi)
      : fPhi(phi), fTheta(theta), fPsi(psi)
  {
    Rectify();
  } // Added 27 Jan. 06   JMM

  /**
     Construct given a pair of pointers or iterators defining the
     beginning and end of an array of three Scalars, to be treated as
     the angles phi, theta and psi.
  */
  template <class IT> EulerAngles(IT begin, IT end)
  {
    SetComponents(begin, end);
  }

  // The compiler-generated copy ctor, copy assignment, and dtor are OK.

  /**
     Re-adjust components place angles in canonical ranges
  */
  void Rectify();

  // ======== Construction and assignement from any other rotation
  // ==================

  /**
     Create from any other supported rotation (see gv_detail::convert )
   */
  template <class OtherRotation> explicit EulerAngles(const OtherRotation &r);

  /**
     Assign from any other rotation (see gv_detail::convert )
  */
  template <class OtherRotation> EulerAngles &operator=(OtherRotation const &r);

#ifdef OLD
  explicit EulerAngles(const Rotation3D &r) { gv_detail::convert(r, *this); }

  /**
     Construct from a rotation matrix
  */
  explicit EulerAngles(const Rotation3D &r) { gv_detail::convert(r, *this); }

  /**
     Construct from a rotation represented by a Quaternion
  */
  explicit EulerAngles(const Quaternion &q) { gv_detail::convert(q, *this); }

  /**
     Construct from an AxisAngle
  */
  explicit EulerAngles(const AxisAngle &a) { gv_detail::convert(a, *this); }

  /**
     Construct from an axial rotation
  */
  explicit EulerAngles(RotationZ const &r) { gv_detail::convert(r, *this); }
  explicit EulerAngles(RotationY const &r) { gv_detail::convert(r, *this); }
  explicit EulerAngles(RotationX const &r) { gv_detail::convert(r, *this); }

  /**
     Assign from an AxisAngle
  */
  EulerAngles &operator=(AxisAngle const &a)
  {
    return operator=(EulerAngles(a));
  }

  /**
     Assign from a Quaternion
  */
  EulerAngles &operator=(Quaternion const &q)
  {
    return operator=(EulerAngles(q));
  }

  /**
     Assign from an axial rotation
  */
  EulerAngles &operator=(RotationZ const &r)
  {
    return operator=(EulerAngles(r));
  }
  EulerAngles &operator=(RotationY const &r)
  {
    return operator=(EulerAngles(r));
  }
  EulerAngles &operator=(RotationX const &r)
  {
    return operator=(EulerAngles(r));
  }

#endif

  // ======== Components ==============

  /**
     Set the three Euler angles given a pair of pointers or iterators
     defining the beginning and end of an array of three Scalars.
  */
  template <class IT>
#ifndef NDEBUG
  void SetComponents(IT begin, IT end)
  {
#else
  void SetComponents(IT begin, IT)
  {
#endif
    fPhi   = *begin++;
    fTheta = *begin++;
    fPsi   = *begin++;
    assert(begin == end);
    Rectify(); // Added 27 Jan. 06   JMM
  }

  /**
     Get the axis and then the angle into data specified by an iterator begin
     and another to the end of the desired data (4 past start).
  */
  template <class IT>
#ifndef NDEBUG
  void GetComponents(IT begin, IT end) const
  {
#else
  void GetComponents(IT begin, IT) const
  {
#endif
    *begin++ = fPhi;
    *begin++ = fTheta;
    *begin++ = fPsi;
    assert(begin == end);
  }

  /**
     Get the axis and then the angle into data specified by an iterator begin
  */
  template <class IT> void GetComponents(IT begin) const
  {
    *begin++ = fPhi;
    *begin++ = fTheta;
    *begin   = fPsi;
  }

  /**
     Set the components phi, theta, psi based on three Scalars.
  */
  void SetComponents(Scalar phi, Scalar theta, Scalar psi)
  {
    fPhi   = phi;
    fTheta = theta;
    fPsi   = psi;
    Rectify(); // Added 27 Jan. 06   JMM
  }

  /**
     Get the components phi, theta, psi into three Scalars.
  */
  void GetComponents(Scalar &phi, Scalar &theta, Scalar &psi) const
  {
    phi   = fPhi;
    theta = fTheta;
    psi   = fPsi;
  }

  /**
     Set Phi Euler angle // JMM 30 Jan. 2006
  */
  void SetPhi(Scalar phi)
  {
    fPhi = phi;
    Rectify();
  }

  /**
     Return Phi Euler angle
  */
  Scalar Phi() const { return fPhi; }

  /**
     Set Theta Euler angle // JMM 30 Jan. 2006
  */
  void SetTheta(Scalar theta)
  {
    fTheta = theta;
    Rectify();
  }

  /**
     Return Theta Euler angle
  */
  Scalar Theta() const { return fTheta; }

  /**
     Set Psi Euler angle // JMM 30 Jan. 2006
  */
  void SetPsi(Scalar psi)
  {
    fPsi = psi;
    Rectify();
  }

  /**
     Return Psi Euler angle
  */
  Scalar Psi() const { return fPsi; }

  // =========== operations ==============

  /**
     Rotation operation on a displacement vector in any coordinate system and
     tag
  */
  template <class CoordSystem, class U>
  DisplacementVector3D<CoordSystem, U>
  operator()(const DisplacementVector3D<CoordSystem, U> &v) const;

  /**
     Rotation operation on a position vector in any coordinate system
  */
  template <class CoordSystem, class U>
  PositionVector3D<CoordSystem, U>
  operator()(const PositionVector3D<CoordSystem, U> &v) const
  {
    DisplacementVector3D<Cartesian3D<double>, U> xyz(v);
    DisplacementVector3D<Cartesian3D<double>, U> rxyz = operator()(xyz);
    return PositionVector3D<CoordSystem, U>(rxyz);
  }

  /**
     Rotation operation on a Lorentz vector in any 4D coordinate system
  */
  template <class CoordSystem>
  inline LorentzVector<CoordSystem>
  operator()(const LorentzVector<CoordSystem> &v) const;

  /**
     Rotation operation on an arbitrary vector v.
     Preconditions:  v must implement methods x(), y(), and z()
     and the arbitrary vector type must have a constructor taking (x,y,z)
  */
  template <class ForeignVector>
  inline ForeignVector operator()(const ForeignVector &v) const;

  /**
     Overload operator * for rotation on a vector
  */
  template <class AVector> inline AVector operator*(const AVector &v) const
  {
    return operator()(v);
  }

  /**
     Invert a rotation in place
  */
  // theta stays the same and negative rotation in Theta is done via a rotation
  // of + PI in phi and Psi
  void Invert()
  {
    Scalar tmp = -fPhi;
    fPhi       = -fPsi + Pi();
    fPsi       = tmp + Pi();
  }

  /**
     Return inverse of a rotation
  */
  EulerAngles Inverse() const
  {
    return EulerAngles(-fPsi + Pi(), fTheta, -fPhi + Pi());
  }

  // ========= Multi-Rotation Operations ===============

  /**
     Multiply (combine) two rotations
  */
  EulerAngles operator*(const Rotation3D &r) const;
  EulerAngles operator*(const AxisAngle &a) const;
  EulerAngles operator*(const EulerAngles &e) const;
  EulerAngles operator*(const Quaternion &q) const;
  EulerAngles operator*(const RotationX &rx) const;
  EulerAngles operator*(const RotationY &ry) const;
  EulerAngles operator*(const RotationZ &rz) const;

  /**
     Post-Multiply (on right) by another rotation :  T = T*R
  */
  template <class R> EulerAngles &operator*=(const R &r)
  {
    return *this = (*this) * r;
  }

  /**
     Distance between two rotations
  */
  template <class R> Scalar Distance(const R &r) const;

  /**
     Equality/inequality operators
  */
  bool operator==(const EulerAngles &rhs) const
  {
    if (fPhi != rhs.fPhi)
      return false;
    if (fTheta != rhs.fTheta)
      return false;
    if (fPsi != rhs.fPsi)
      return false;
    return true;
  }
  bool operator!=(const EulerAngles &rhs) const { return !operator==(rhs); }

private:
  double fPhi;   // Z rotation angle (first)  defined in [-PI,PI]
  double fTheta; // X rotation angle (second) defined only [0,PI]
  double fPsi;   // Z rotation angle (third)  defined in [-PI,PI]

  static double Pi() { return M_PI; }

}; // EulerAngles

/**
   Distance between two rotations
 */
template <class R>
inline typename EulerAngles::Scalar Distance(const EulerAngles &r1,
                                             const R &r2);

/**
   Multiplication of an axial rotation by an AxisAngle
 */
EulerAngles operator*(RotationX const &r1, EulerAngles const &r2);
EulerAngles operator*(RotationY const &r1, EulerAngles const &r2);
EulerAngles operator*(RotationZ const &r1, EulerAngles const &r2);

/**
   Stream Output and Input
 */
// TODO - I/O should be put in the manipulator form

std::ostream &operator<<(std::ostream &os, const EulerAngles &e);

} // namespace Math

} // namespace ROOT

#include "Math/GenVector/3DDistances.h"
#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/Rotation3D.h"

namespace ROOT
{
namespace Math


{

template <class CoordSystem, class U>
DisplacementVector3D<CoordSystem, U>
EulerAngles::operator()(const DisplacementVector3D<CoordSystem, U> &v) const
{
return Rotation3D(*this)(v);
}
    
    
    
/**
   Distance between two rotations
 */
template <class R>
inline typename EulerAngles::Scalar Distance(const EulerAngles &r1, const R &r2)
{
  return gv_detail::dist(r1, r2);
}

template <class R> EulerAngles::Scalar EulerAngles::Distance(const R &r) const
{
  return gv_detail::dist(*this, r);
}

template <class CoordSystem>
LorentzVector<CoordSystem> EulerAngles::
operator()(const LorentzVector<CoordSystem> &v) const
{
  DisplacementVector3D<Cartesian3D<double>> xyz(v.Vect());
  xyz = operator()(xyz);
  LorentzVector<PxPyPzE4D<double>> xyzt(xyz.X(), xyz.Y(), xyz.Z(), v.E());
  return LorentzVector<CoordSystem>(xyzt);
}

/**
   Rotation operation on an arbitrary vector v.
   Preconditions:  v must implement methods x(), y(), and z()
   and the arbitrary vector type must have a constructor taking (x,y,z)
*/
template <class ForeignVector>
ForeignVector EulerAngles::operator()(const ForeignVector &v) const
{
  DisplacementVector3D<Cartesian3D<double>> xyz(v);
  DisplacementVector3D<Cartesian3D<double>> rxyz = operator()(xyz);
  return ForeignVector(rxyz.X(), rxyz.Y(), rxyz.Z());
}
} // namespace Math
} // namespace ROOT

#include "Math/GenVector/3DConversions.h"

namespace ROOT
{
namespace Math
{

/**
   Create from any other supported rotation (see gv_detail::convert )
 */
template <class OtherRotation>
inline EulerAngles::EulerAngles(const OtherRotation &r)
{
  gv_detail::convert(r, *this);
}

/**
   Assign from any other rotation (see gv_detail::convert )
*/
template <class OtherRotation>
inline EulerAngles &EulerAngles::operator=(OtherRotation const &r)
{
  gv_detail::convert(r, *this);
  return *this;
}

// ========== Constructors and Assignment =====================

inline void EulerAngles::Rectify()
{
  // rectify
  if (fTheta < 0 || fTheta > Pi())
  {
    Scalar t = fTheta - std::floor(fTheta / (2 * Pi())) * 2 * Pi();
    if (t <= Pi())
    {
      fTheta = t;
    }
    else
    {
      fTheta = 2 * Pi() - t;
      fPhi   = fPhi + Pi();
      fPsi   = fPsi + Pi();
    }
  }

  if (fPhi <= -Pi() || fPhi > Pi())
  {
    fPhi = fPhi - std::floor(fPhi / (2 * Pi()) + .5) * 2 * Pi();
  }

  if (fPsi <= -Pi() || fPsi > Pi())
  {
    fPsi = fPsi - std::floor(fPsi / (2 * Pi()) + .5) * 2 * Pi();
  }

} // Rectify()

} // namespace Math
} // namespace ROOT

#include "Math/GenVector/3DConversions.h"
#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PositionVector3D.h"
#include "Math/GenVector/Quaternion.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/RotationX.h"
#include "Math/GenVector/RotationY.h"
#include "Math/GenVector/RotationZ.h"

namespace ROOT
{
namespace Math
{

// ========== Operations =====================

// DisplacementVector3D< Cartesian3D<double> >
// EulerAngles::
// operator() (const DisplacementVector3D< Cartesian3D<double> > & v) const
// {
//   return Rotation3D(*this)(v);
// }

inline EulerAngles EulerAngles::operator*(const Rotation3D &r) const
{
  // combine with a Rotation3D
  return EulerAngles(Rotation3D(*this) * r);
}

inline EulerAngles EulerAngles::operator*(const AxisAngle &a) const
{
  // combine with a AxisAngle
  return EulerAngles(Quaternion(*this) * Quaternion(a));
}

inline EulerAngles EulerAngles::operator*(const EulerAngles &e) const
{
  // combine with a EulerAngles
  return EulerAngles(Quaternion(*this) * Quaternion(e));
}

inline EulerAngles EulerAngles::operator*(const Quaternion &q) const
{
  // combination with a Quaternion
  return EulerAngles(Quaternion(*this) * q);
}

inline EulerAngles EulerAngles::operator*(const RotationX &r) const
{
  // combine with a RotationX
  return EulerAngles(Quaternion(*this) * r);
}

inline EulerAngles EulerAngles::operator*(const RotationY &r) const
{
  // combine with a RotationY
  return EulerAngles(Quaternion(*this) * r);
}

inline EulerAngles EulerAngles::operator*(const RotationZ &r) const
{
  // combine with a RotationZ
  // TODO -- this can be made much faster because it merely adds
  //         the r.Angle() to phi.
  Scalar newPhi = fPhi + r.Angle();
  if (newPhi <= -Pi() || newPhi > Pi())
  {
    newPhi = newPhi - std::floor(newPhi / (2 * Pi()) + .5) * 2 * Pi();
  }
  return EulerAngles(newPhi, fTheta, fPsi);
}

inline EulerAngles operator*(RotationX const &r, EulerAngles const &e)
{
  return EulerAngles(r) * e; // TODO: improve performance
}

inline EulerAngles operator*(RotationY const &r, EulerAngles const &e)
{
  return EulerAngles(r) * e; // TODO: improve performance
}

inline EulerAngles operator*(RotationZ const &r, EulerAngles const &e)
{
  return EulerAngles(r) * e; // TODO: improve performance
}

} // namespace Math
} // namespace ROOT

#endif /* ROOT_Math_GenVector_EulerAngles  */
