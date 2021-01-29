/*! \file   Scattering.cxx
 *   \brief  Source filefor the Scattering bug routines.
 *
 *   This version has a major bug and produces too small scattering angles.
 *
 *   \date   2013.08.19
 *   \author Tomasz Fuchs
 **/

#include <cmath>
#include <tuple>
#include <string>
#include <sstream>
#include <iostream>


#include "PROPOSAL/methods.h"
#include "PROPOSAL/scattering/multiple_scattering/Parametrization.h"

using namespace PROPOSAL::multiple_scattering;

Parametrization::Parametrization(double _mass)
    : mass(_mass)
{
}

bool Parametrization::operator==(const Parametrization& scattering) const
{
    if (mass != scattering.mass)
        return false;
    else
        return this->compare(scattering);
}

namespace PROPOSAL {
std::ostream& operator<<(
    std::ostream& os, multiple_scattering::Parametrization const& scattering)
{
    std::stringstream ss;
    ss << " Multiple scattering (" << &scattering << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    scattering.print(os);

    os << Helper::Centered(60, "");

    return os;
}
} // namespace PROPOSAL
