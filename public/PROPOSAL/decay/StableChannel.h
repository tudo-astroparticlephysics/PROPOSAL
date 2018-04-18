
#pragma once

// #include <vector>

#include "PROPOSAL/decay/DecayChannel.h"

namespace PROPOSAL {

class Particle;

class StableChannel : public DecayChannel
{
public:
    StableChannel();
    StableChannel(const StableChannel& mode);
    virtual ~StableChannel();
    // No copy and assignemnt -> done by clone
    DecayChannel* clone() const { return new StableChannel(*this); }

    DecayProducts Decay(const Particle&);

    const std::string& GetName() const { return name_; }

private:
    StableChannel& operator=(const StableChannel&); // Undefined & not allowed

    bool compare(const DecayChannel&) const;

    static const std::string name_;
};

} // namespace PROPOSAL
