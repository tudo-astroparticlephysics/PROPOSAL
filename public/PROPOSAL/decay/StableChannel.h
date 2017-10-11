
#pragma once

// #include <vector>

#include "PROPOSAL/decay/DecayChannel.h"

namespace PROPOSAL
{

class PROPOSALParticle;


class StableChannel : public DecayChannel
{
    public:
    StableChannel();
    StableChannel(const StableChannel& mode);
    virtual ~StableChannel();
    // No copy and assignemnt -> done by clone
    DecayChannel* clone() { return new StableChannel(*this); }

    DecayProducts Decay(PROPOSALParticle*);

    private:
    StableChannel& operator=(const StableChannel&); // Undefined & not allowed

    bool compare(const DecayChannel&) const;
};

} /* PROPOSAL */
