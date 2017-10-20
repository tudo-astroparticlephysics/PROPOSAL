
#pragma once

#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/math/Integral.h"

#define UTILITY_INTEGRAL_DEC(cls)                                                                                      \
    class UtilityIntegral##cls : public UtilityIntegral                                                                \
    {                                                                                                                  \
        public:                                                                                                        \
        UtilityIntegral##cls(const Utility&);                                                                          \
                                                                                                                       \
        UtilityIntegral##cls(const Utility&, const UtilityIntegral##cls&);                                             \
        UtilityIntegral##cls(const UtilityIntegral##cls&);                                                             \
        virtual ~UtilityIntegral##cls();                                                                               \
                                                                                                                       \
        virtual UtilityDecorator* clone(const Utility& utility) const                                                  \
        {                                                                                                              \
            return new UtilityIntegral##cls(utility, *this);                                                           \
        }                                                                                                              \
                                                                                                                       \
        double FunctionToIntegral(double energy);                                                                      \
        virtual double Calculate(double ei, double ef, double rnd);                                                    \
                                                                                                                       \
        private:                                                                                                       \
        UtilityDecorator& operator=(const UtilityDecorator&);                                                          \
    };

namespace PROPOSAL {

class UtilityIntegral : public UtilityDecorator
{
    public:
    UtilityIntegral(const Utility&);

    // Copy constructors
    UtilityIntegral(const Utility&, const UtilityIntegral&);
    UtilityIntegral(const UtilityIntegral&);
    virtual UtilityDecorator* clone(const Utility&) const = 0;

    virtual ~UtilityIntegral();

    // Methods
    virtual double Calculate(double ei, double ef, double rnd) = 0;
    virtual double GetUpperLimit(double ei, double rnd);

    protected:
    UtilityIntegral& operator=(const UtilityIntegral&); // Undefined & not allowed

    virtual bool compare(const UtilityDecorator&) const;

    Integral integral_;

};

class UtilityIntegralDisplacement: public UtilityIntegral
{
    public:
    UtilityIntegralDisplacement(const Utility&);

    // Copy constructors
    UtilityIntegralDisplacement(const Utility&, const UtilityIntegralDisplacement&);
    UtilityIntegralDisplacement(const UtilityIntegralDisplacement&);
    virtual UtilityIntegral* clone(const Utility& utility) const { return new UtilityIntegralDisplacement(utility, *this); }

    virtual ~UtilityIntegralDisplacement();

    // Methods
    virtual double Calculate(double ei, double ef, double rnd);

    private:
    UtilityDecorator& operator=(const UtilityDecorator&); // Undefined & not allowed
};

UTILITY_INTEGRAL_DEC(Interaction)
UTILITY_INTEGRAL_DEC(Decay)
UTILITY_INTEGRAL_DEC(Time)
UTILITY_INTEGRAL_DEC(ContRand)
UTILITY_INTEGRAL_DEC(Scattering)

#undef UTILITY_INTEGRAL_DEC

} /* PROPOSAL */
