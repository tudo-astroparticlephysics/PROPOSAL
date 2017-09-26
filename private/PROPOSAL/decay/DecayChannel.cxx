
#include "PROPOSAL/decay/DecayChannel.h"
// #include "PROPOSAL/math/RootFinder.h"

using namespace PROPOSAL;

bool DecayChannel::operator==(const DecayChannel& table) const
{
    return this->compare(table);
}

bool DecayChannel::operator!=(const DecayChannel& def) const
{
    return !(*this == def);
}
