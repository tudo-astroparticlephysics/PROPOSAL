
#pragma once

#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"

namespace PROPOSAL {

class Utility;
class UtilityDecorator;


/**
 * \brief Class containing the functions to randomize the continuous energy losses
 *
 */
class ContinuousRandomizer
{
    public:
    ContinuousRandomizer(Utility&);
    ContinuousRandomizer(const ContinuousRandomizer&);
    ~ContinuousRandomizer();

    // bool operator==(const ContinuousRandomizer& scattering) const;
    // bool operator!=(const ContinuousRandomizer& scattering) const;
    // void swap(ContinuousRandomizer& scattering);

    double Randomize(double ei, double ef, double rnd);

    private:
    ContinuousRandomizer& operator=(const ContinuousRandomizer&); // Undefined & not allowed

    UtilityDecorator* DE2de;
};

}
