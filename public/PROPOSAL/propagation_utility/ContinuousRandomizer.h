
#pragma once


namespace PROPOSAL {

class Utility;
class UtilityDecorator;
struct InterpolationDef;


/**
 * \brief Class containing the functions to randomize the continuous energy losses
 *
 */
class ContinuousRandomizer
{
    public:
    ContinuousRandomizer(Utility&);
    ContinuousRandomizer(Utility&, InterpolationDef);

    // Copy constructors
    ContinuousRandomizer(const Utility&, const ContinuousRandomizer&);
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
