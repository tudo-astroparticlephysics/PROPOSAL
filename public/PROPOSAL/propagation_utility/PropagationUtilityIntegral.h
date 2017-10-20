
#pragma once

#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/math/Integral.h"

#define UTILITY_INTEGRAL_DEC(cls)                                                                                      \
    class UtilityIntegral##cls : public UtilityIntegral                                                                \
    {                                                                                                                  \
        public:                                                                                                        \
        UtilityIntegral##cls(const Utility&);                                                                          \
        UtilityIntegral##cls(const UtilityIntegral##cls&);                                                             \
        virtual ~UtilityIntegral##cls();                                                                               \
                                                                                                                       \
        virtual UtilityDecorator* clone() const { return new UtilityIntegral##cls(*this); }                            \
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
    UtilityIntegral(const UtilityIntegral&);
    virtual ~UtilityIntegral();

    virtual UtilityDecorator* clone() const = 0;

    virtual double Calculate(double ei, double ef, double rnd) = 0;
    virtual double GetUpperLimit(double ei, double rnd);

    protected:
    UtilityIntegral& operator=(const UtilityIntegral&); // Undefined & not allowed

    bool compare(const UtilityDecorator&) const;

    Integral integral_;

};

class UtilityIntegralDisplacement: public UtilityIntegral
{
    public:
    UtilityIntegralDisplacement(const Utility&);
    UtilityIntegralDisplacement(const UtilityIntegralDisplacement&);
    virtual ~UtilityIntegralDisplacement();

    virtual UtilityIntegral* clone() const { return new UtilityIntegralDisplacement(*this); }

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
