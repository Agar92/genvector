// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 ROOT FNAL MathLib Team                          *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for BoostZ
//
// Created by: Mark Fischler  Mon Nov 1  2005
//
// Last update: $Id$
//
#ifndef ROOT_Math_GenVector_BoostZ
#define ROOT_Math_GenVector_BoostZ 1

#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PxPyPzE4D.h"

namespace ROOT {

namespace Math {

//__________________________________________________________________________________________
/**
   Class representing a Lorentz Boost along the Z axis, by beta.
   For efficiency, gamma is held as well.

   @ingroup GenVector
*/

class BoostZ {

public:
  typedef double Scalar;

  enum ELorentzRotationMatrixIndex {
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

  enum EBoostMatrixIndex {
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
  BoostZ();

  /**
     Construct given a Scalar beta_z
  */
  explicit BoostZ(Scalar beta_z) { SetComponents(beta_z); }

  // The compiler-generated copy ctor, copy assignment, and dtor are OK.

  /**
     Re-adjust components to eliminate small deviations from a perfect
     orthosyplectic matrix.
  */
  void Rectify();

  // ======== Components ==============

  /**
     Set components from a Scalar beta_z
  */
  void SetComponents(Scalar beta_z);

  /**
     Get components into a Scalar beta_z
  */
  void GetComponents(Scalar &beta_z) const;

  /**
      Retrieve the beta of the Boost
  */
  Scalar Beta() const { return fBeta; }

  /**
      Retrieve the gamma of the Boost
  */
  Scalar Gamma() const { return fGamma; }

  /**
      Set the given beta of the Boost
  */
  void SetBeta(Scalar beta) { SetComponents(beta); }

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
  operator()(const LorentzVector<CoordSystem> &v) const {
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
  Foreign4Vector operator()(const Foreign4Vector &v) const {
    LorentzVector<PxPyPzE4D<double>> xyzt(v);
    LorentzVector<PxPyPzE4D<double>> r_xyzt = operator()(xyzt);
    return Foreign4Vector(r_xyzt.X(), r_xyzt.Y(), r_xyzt.Z(), r_xyzt.T());
  }

  /**
     Overload operator * for boost on a vector
  */
  template <class A4Vector> inline A4Vector operator*(const A4Vector &v) const {
    return operator()(v);
  }

  /**
     Invert a BoostZ in place
  */
  void Invert();

  /**
     Return inverse of  a BoostZ
  */
  BoostZ Inverse() const;

  /**
     Equality/inequality operators
  */
  bool operator==(const BoostZ &rhs) const {
    if (fBeta != rhs.fBeta)
      return false;
    if (fGamma != rhs.fGamma)
      return false;
    return true;
  }
  bool operator!=(const BoostZ &rhs) const { return !operator==(rhs); }

private:
  Scalar fBeta;  // boost beta z
  Scalar fGamma; // boost gamma

}; // BoostZ

// ============ Class BoostZ ends here ============

/**
   Stream Output and Input
 */
// TODO - I/O should be put in the manipulator form

std::ostream &operator<<(std::ostream &os, const BoostZ &b);

// ============ Declarations end here ============

inline BoostZ::BoostZ() : fBeta(0.0), fGamma(1.0) {}

inline void BoostZ::SetComponents(Scalar bz) {
  // set component
  Scalar bp2 = bz * bz;
  if (bp2 >= 1) {
    GenVector::Throw(
        "Beta Vector supplied to set BoostZ represents speed >= c");
    return;
  }
  fBeta = bz;
  fGamma = 1.0 / sqrt(1.0 - bp2);
}

inline void BoostZ::GetComponents(Scalar &bz) const {
  // get component
  bz = fBeta;
}

inline DisplacementVector3D<Cartesian3D<BoostZ::Scalar>>
BoostZ::BetaVector() const {
  // return beta vector
  return DisplacementVector3D<Cartesian3D<Scalar>>(0.0, 0.0, fBeta);
}

inline void BoostZ::GetLorentzRotation(Scalar r[]) const {
  // get corresponding LorentzRotation
  r[kLXX] = 1.0;
  r[kLXY] = 0.0;
  r[kLXZ] = 0.0;
  r[kLXT] = 0.0;
  r[kLYX] = 0.0;
  r[kLYY] = 1.0;
  r[kLYZ] = 0.0;
  r[kLYT] = 0.0;
  r[kLZX] = 0.0;
  r[kLZY] = 0.0;
  r[kLZZ] = fGamma;
  r[kLZT] = fGamma * fBeta;
  r[kLTX] = 0.0;
  r[kLTY] = 0.0;
  r[kLTZ] = fGamma * fBeta;
  r[kLTT] = fGamma;
}

inline void BoostZ::Rectify() {
  // Assuming the representation of this is close to a true Lorentz Rotation,
  // but may have drifted due to round-off error from many operations,
  // this forms an "exact" orthosymplectic matrix for the Lorentz Rotation
  // again.

  if (fGamma <= 0) {
    GenVector::Throw("Attempt to rectify a boost with non-positive gamma");
    return;
  }
  Scalar beta = fBeta;
  if (beta >= 1) {
    beta /= (beta * (1.0 + 1.0e-16));
  }
  SetComponents(beta);
}

inline LorentzVector<PxPyPzE4D<double>> BoostZ::
operator()(const LorentzVector<PxPyPzE4D<double>> &v) const {
  // apply boost to a LV
  Scalar z = v.Pz();
  Scalar t = v.E();
  return LorentzVector<PxPyPzE4D<double>>(v.Px(), v.Py(),
                                          fGamma * z + fGamma * fBeta * t,
                                          fGamma * fBeta * z + fGamma * t);
}

inline void BoostZ::Invert() {
  // invert
  fBeta = -fBeta;
}

inline BoostZ BoostZ::Inverse() const {
  // return an inverse boostZ
  BoostZ tmp(*this);
  tmp.Invert();
  return tmp;
}

} // namespace Math
} // namespace ROOT

#endif /* ROOT_Math_GenVector_BoostZ  */
