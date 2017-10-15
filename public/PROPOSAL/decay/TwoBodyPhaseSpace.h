
#pragma once

// #include <vector>

#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/particle/ParticleDef.h"

namespace PROPOSAL
{

class Particle;

class TwoBodyPhaseSpace : public DecayChannel
{
    public:
    TwoBodyPhaseSpace(const ParticleDef&, const ParticleDef&);
    TwoBodyPhaseSpace(const TwoBodyPhaseSpace& mode);
    virtual ~TwoBodyPhaseSpace();
    // No copy and assignemnt -> done by clone
    DecayChannel* clone() { return new TwoBodyPhaseSpace(*this); }

    DecayProducts Decay(Particle&);

    const std::string& GetName() const { return name_; }

    private:
    TwoBodyPhaseSpace& operator=(const TwoBodyPhaseSpace&); // Undefined & not allowed

    bool compare(const DecayChannel&) const;
    void print(std::ostream&) const;

    ParticleDef first_daughter_;
    ParticleDef second_daughter_;

    static const std::string name_;
};

} /* PROPOSAL */

