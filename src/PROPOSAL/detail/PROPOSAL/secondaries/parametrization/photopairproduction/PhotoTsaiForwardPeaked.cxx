#include "PROPOSAL/Constants.h"
#include "PROPOSAL/secondaries/parametrization/photopairproduction/PhotoTsaiForwardPeaked.h"

#include <cmath>
#include <tuple>

using namespace PROPOSAL;

std::tuple<Cartesian3D, Cartesian3D>
secondaries::PhotoTsaiForwardPeaked::CalculateDirections(const Vector3D& dir,
    double energy, double, const Component&, std::vector<double> rnd)
{
    auto k = energy/ ME;

    auto cosphi0 = std::cos(1. / k);
    auto cosphi1 = std::cos(1. / k);

    auto theta0 = rnd[4] * 2. * PI;
    auto theta1 = std::fmod(theta0 + PI, 2. * PI);
    // TODO: Sometimes the intergration fails and -1 instead of 1 is returned...
    // :(
    if (cosphi0 == -1.)
        cosphi0 *= (-1);
    if (cosphi1 == -1.)
        cosphi1 *= (-1);
    auto dir_0 = Cartesian3D(dir);
    dir_0.deflect(cosphi0, theta0);
    auto dir_1 = Cartesian3D(dir);
    dir_1.deflect(cosphi1, theta1);
    return std::make_tuple(dir_0, dir_1);
}
