#include "Math/LorentzRotation.h"
#include "Math/LorentzVector.h"
#include <cstdint>

using IntVector    = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<uint64_t>>;
using DoubleVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>;

int main()
{
  constexpr IntVector aVector;
  constexpr auto m2 = aVector.M2();
  static_assert((m2 == 0), "Must compute M2 of a zero vector at compile time");
  constexpr DoubleVector aDVector;
  constexpr auto m2d = aDVector.M2();
  static_assert((m2d == 0.),
                "Must compute M2 of a zero floating vector at compile time");
  constexpr auto threeVector = aDVector.Vect();
  static_assert(
      (threeVector.Mag2() == 0.),
      "Must compute length of a zero floating three vector at compile time");
  constexpr DoubleVector v{15., 1., 3., 100.};
  constexpr DoubleVector u{30., 1., 3., 70.};
  constexpr DoubleVector w{3. * v + 4. * u};
  static_assert(!std::isnan(w.M2()),
                "Must compute M2 of a vector at compile time");

  return 0;
}
