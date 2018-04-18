
#pragma once

// #include <vector>

#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/particle/ParticleDef.h"

namespace PROPOSAL {

class Particle;

class LeptonicDecayChannel : public DecayChannel
{
public:
    LeptonicDecayChannel(const ParticleDef&, const ParticleDef&, const ParticleDef&);
    LeptonicDecayChannel(const LeptonicDecayChannel& mode);
    virtual ~LeptonicDecayChannel();
    // No copy and assignemnt -> done by clone
    DecayChannel* clone() const { return new LeptonicDecayChannel(*this); }

    DecayProducts Decay(const Particle&);

    const std::string& GetName() const { return name_; }

private:
    ParticleDef massive_lepton_;
    ParticleDef neutrino_;
    ParticleDef anti_neutrino_;
    static const std::string name_;

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
};

} // namespace PROPOSAL
