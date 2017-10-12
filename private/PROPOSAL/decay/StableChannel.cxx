
#include "PROPOSAL/decay/StableChannel.h"

using namespace PROPOSAL;

StableChannel::StableChannel()
    : DecayChannel()
{
}

StableChannel::~StableChannel()
{
}

StableChannel::StableChannel(const StableChannel& mode)
    : DecayChannel(mode)
{
}

bool StableChannel::compare(const DecayChannel& channel) const
{
    const StableChannel* stable = dynamic_cast<const StableChannel*>(&channel);

    if (!stable)
        return false;
    else
        return true;
}

DecayChannel::DecayProducts StableChannel::Decay(Particle*)
{
    // return empty vector;
    DecayProducts products;
    return products;
}
