
#include "PROPOSAL/decay/StableChannel.h"
#include <ostream>

using namespace PROPOSAL;

const std::string StableChannel::name_ = "StableChannel";

StableChannel::StableChannel()
    : DecayChannel()
{
}

StableChannel::~StableChannel() {}

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

DecayChannel::DecayProducts StableChannel::Decay(const Particle&)
{
    // return empty vector;
    return DecayProducts();
}
