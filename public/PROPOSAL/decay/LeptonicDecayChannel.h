
#pragma once

// #include <vector>

#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/math/RootFinder.h"

namespace PROPOSAL
{

class PROPOSALParticle;

class LeptonicDecayChannel : public DecayChannel
{
    public:
    LeptonicDecayChannel();
    LeptonicDecayChannel(const LeptonicDecayChannel& mode);
    virtual ~LeptonicDecayChannel();
    // No copy and assignemnt -> done by clone
    DecayChannel* clone() { return new LeptonicDecayChannel(*this); }

    DecayProducts Decay(PROPOSALParticle*);

    const std::string& GetName() const { return name_; }

    private:
    LeptonicDecayChannel& operator=(const LeptonicDecayChannel&); // Undefined & not allowed

    bool compare(const DecayChannel&) const;
    void print(std::ostream&) const;

    // ----------------------------------------------------------------------------
    /// @brief Function for electron energy calculation - interface to FindRoot
    // ----------------------------------------------------------------------------
    double DecayRate(double);

    // ----------------------------------------------------------------------------
    /// @brief Function for electron energy calculation - interface to FindRoot
    // ----------------------------------------------------------------------------
    double DifferentialDecayRate(double);

    RootFinder root_finder_;
    static const std::string name_;
};

} /* PROPOSAL */

