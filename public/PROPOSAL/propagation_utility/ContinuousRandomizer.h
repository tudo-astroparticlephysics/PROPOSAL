
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
    ContinuousRandomizer(const Utility&);
    ContinuousRandomizer(const Utility&, const InterpolationDef);

    // Copy constructors
    ContinuousRandomizer(const Utility&, const ContinuousRandomizer&);
    ContinuousRandomizer(const ContinuousRandomizer&);

    ~ContinuousRandomizer();

    bool operator==(const ContinuousRandomizer&) const;
    bool operator!=(const ContinuousRandomizer&) const;

    double Randomize(double ei, double ef, double rnd);

private:
    ContinuousRandomizer& operator=(const ContinuousRandomizer&); // Undefined & not allowed

    UtilityDecorator* DE2de;
};

} // namespace PROPOSAL
