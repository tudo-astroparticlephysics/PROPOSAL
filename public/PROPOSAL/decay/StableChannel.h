
#pragma once

// #include <vector>

#include "PROPOSAL/decay/DecayChannel.h"

namespace PROPOSAL
{

class Particle;


class StableChannel : public DecayChannel
{
    public:
    StableChannel();
    StableChannel(const StableChannel& mode);
    virtual ~StableChannel();
    // No copy and assignemnt -> done by clone
    DecayChannel* clone() { return new StableChannel(*this); }

    DecayProducts Decay(Particle*);

    private:
    StableChannel& operator=(const StableChannel&); // Undefined & not allowed

    bool compare(const DecayChannel&) const;
};

} /* PROPOSAL */
