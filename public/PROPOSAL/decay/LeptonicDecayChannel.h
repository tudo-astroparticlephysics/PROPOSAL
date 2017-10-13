
#pragma once

// #include <vector>

#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/math/RootFinder.h"

namespace PROPOSAL
{

class Particle;

class LeptonicDecayChannel : public DecayChannel
{
    public:
    LeptonicDecayChannel();
    LeptonicDecayChannel(const LeptonicDecayChannel& mode);
    virtual ~LeptonicDecayChannel();
    // No copy and assignemnt -> done by clone
    DecayChannel* clone() { return new LeptonicDecayChannel(*this); }

    DecayProducts Decay(Particle*);

    private:
    LeptonicDecayChannel& operator=(const LeptonicDecayChannel&); // Undefined & not allowed

    bool compare(const DecayChannel&) const;

    // ----------------------------------------------------------------------------
    /// @brief Function for electron energy calculation - interface to FindRoot
    // ----------------------------------------------------------------------------
    double DecayRate(double, double);

    // ----------------------------------------------------------------------------
    /// @brief Function for electron energy calculation - interface to FindRoot
    // ----------------------------------------------------------------------------
    double DifferentialDecayRate(double);

    std::pair<double, double> function_and_derivative(double x, double right_side);

    double FindRootBoost(double min, double right_side);

    RootFinder root_finder_;
};

} /* PROPOSAL */

