
#pragma once

#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
// #include "PROPOSAL/math/Interpolant.h"

#define UTILITY_INTERPOLANT_DEC(cls)                                                                                   \
    class UtilityInterpolant##cls : public UtilityInterpolant                                                          \
    {                                                                                                                  \
        public:                                                                                                        \
        UtilityInterpolant##cls(const Utility&, InterpolationDef);                                                     \
        UtilityInterpolant##cls(const UtilityInterpolant##cls&);                                                       \
        ~UtilityInterpolant##cls();                                                                                    \
                                                                                                                       \
        UtilityDecorator* clone() const { return new UtilityInterpolant##cls(*this); }                                 \
                                                                                                                       \
        double Calculate(double ei, double ef, double rnd);                                                            \
        double GetUpperLimit(double ei, double rnd);                                                                   \
                                                                                                                       \
        private:                                                                                                       \
        UtilityInterpolant##cls& operator=(const UtilityInterpolant##cls&);                                            \
                                                                                                                       \
        double BuildInterpolant(double, UtilityIntegral&, Integral&);                                                  \
        void InitInterpolation(const std::string&, UtilityIntegral&, int number_of_sampling_points);                   \
    };

namespace PROPOSAL {

class Interpolant;

class UtilityInterpolant : public UtilityDecorator
{
    public:
    UtilityInterpolant(const Utility&, InterpolationDef);
    UtilityInterpolant(const UtilityInterpolant&);
    virtual ~UtilityInterpolant();

    virtual UtilityDecorator* clone() const = 0;

    virtual double Calculate(double ei, double ef, double rnd) = 0;
    virtual double GetUpperLimit(double ei, double rnd);

    protected:
    UtilityInterpolant& operator=(const UtilityInterpolant&); // Undefined & not allowed

    virtual double BuildInterpolant(double, UtilityIntegral&, Integral&) = 0;
    virtual void InitInterpolation(const std::string&, UtilityIntegral&, int number_of_sampling_points) = 0;

    double stored_result_;
    Interpolant* interpolant_;
    Interpolant* interpolant_diff_;

    InterpolationDef interpolation_def_;
};

class UtilityInterpolantInteraction: public UtilityInterpolant
{
    public:
    UtilityInterpolantInteraction(const Utility&, InterpolationDef);
    UtilityInterpolantInteraction(const UtilityInterpolantInteraction&);
    virtual ~UtilityInterpolantInteraction();

    virtual UtilityInterpolant* clone() const { return new UtilityInterpolantInteraction(*this); }

    double Calculate(double ei, double ef, double rnd);
    double GetUpperLimit(double ei, double rnd);

    private:
    UtilityInterpolantInteraction& operator=(const UtilityInterpolantInteraction&); // Undefined & not allowed

    double BuildInterpolant(double, UtilityIntegral&, Integral&);
    void InitInterpolation(const std::string&, UtilityIntegral&, int number_of_sampling_points);

    double big_low_;
    double up_;
};

class UtilityInterpolantDecay: public UtilityInterpolant
{
    public:
    UtilityInterpolantDecay(const Utility&, InterpolationDef);
    UtilityInterpolantDecay(const UtilityInterpolantDecay&);
    virtual ~UtilityInterpolantDecay();

    virtual UtilityInterpolant* clone() const { return new UtilityInterpolantDecay(*this); }

    double Calculate(double ei, double ef, double rnd);
    double GetUpperLimit(double ei, double rnd);

    private:
    UtilityInterpolantDecay& operator=(const UtilityInterpolantDecay&); // Undefined & not allowed

    double BuildInterpolant(double, UtilityIntegral&, Integral&);
    void InitInterpolation(const std::string&, UtilityIntegral&, int number_of_sampling_points);

    double big_low_;
    double up_;
};

UTILITY_INTERPOLANT_DEC(Displacement)
UTILITY_INTERPOLANT_DEC(Time)
UTILITY_INTERPOLANT_DEC(ContRand)
UTILITY_INTERPOLANT_DEC(Scattering)

} /* PROPOSAL */
