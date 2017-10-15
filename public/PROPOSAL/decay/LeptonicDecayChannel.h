
#pragma once

// #include <vector>

#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/particle/ParticleDef.h"
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

    DecayProducts Decay(Particle&);

    const std::string& GetName() const { return name_; }

    private:
    LeptonicDecayChannel& operator=(const LeptonicDecayChannel&); // Undefined & not allowed

    bool compare(const DecayChannel&) const;
    void print(std::ostream&) const;

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

    // ParticleDef massive_lepton_;
    // ParticleDef first_neutrino_;
    // ParticleDef second_neutrino_;

    RootFinder root_finder_;
    static const std::string name_;
};

} /* PROPOSAL */

