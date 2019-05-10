
#include <sstream>

#include "PROPOSAL/decay/DecayTable.h"
#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/decay/LeptonicDecayChannel.h"
#include "PROPOSAL/decay/StableChannel.h"
#include "PROPOSAL/decay/TwoBodyPhaseSpace.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructors & destructor
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
DecayTable::DecayTable()
    : channels_()
{
}

// ------------------------------------------------------------------------- //
DecayTable::DecayTable(const DecayTable& table)
{
    clearTable();

    for (DecayMap::const_iterator iter = table.channels_.begin(); iter != table.channels_.end(); ++iter)
    {
        channels_[iter->first] = iter->second->clone();
    }
}

// ------------------------------------------------------------------------- //
DecayTable::~DecayTable()
{
    clearTable();
}

DecayTable& DecayTable::operator=(const DecayTable& table)
{
    if (this != &table)
    {
        DecayTable tmp(table);
        swap(*this, tmp);
    }
    return *this;
}

void PROPOSAL::swap(DecayTable& first, DecayTable& second)
{
    using std::swap;
    swap(first.channels_, second.channels_);
}

bool DecayTable::operator==(const DecayTable& table) const
{
    if (channels_.size() != table.channels_.size())
    {
        return false;
    } else
    {
        DecayMap::const_iterator i, j;
        for (i = channels_.begin(), j = table.channels_.begin(); i != channels_.end(); ++i, ++j)
        {
            if (i->first != j->first)
            {
                return false;
            } else if (*i->second != *j->second)
            {
                return false;
            }
        }
        return true;
    }
}

bool DecayTable::operator!=(const DecayTable& def) const
{
    return !(*this == def);
}

// ------------------------------------------------------------------------- //
std::ostream& PROPOSAL::operator<<(std::ostream& os, DecayTable const& table)
{
    std::stringstream ss;
    ss << " DecayTable (" << &table << ") ";

    os << Helper::Centered(60, ss.str()) << '\n';

    os << "Branching Ratio"
       << "\t\t"
       << "Channel" << '\n';

    for (DecayTable::DecayMap::const_iterator iter = table.channels_.begin(); iter != table.channels_.end(); ++iter)
    {
        os << iter->first << "\t\t" << iter->second->GetName() << '\n';
    }

    os << Helper::Centered(60, "");
    return os;
}
// ------------------------------------------------------------------------- //
// Methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
DecayChannel& DecayTable::SelectChannel() const
{
    for (int i = 0; i < 1000; ++i)
    {
        double sumBranchingRatio = 0.0;
        double random            = RandomGenerator::Get().RandomDouble();

        for (DecayMap::const_iterator iter = channels_.begin(); iter != channels_.end(); ++iter)
        {
            sumBranchingRatio += iter->first;

            if (random < sumBranchingRatio)
            {
                return *iter->second;
            }
        }
    }

    log_fatal("No decay channel found. If your particle is stable, call \"SetStable\"!");
    return *channels_.begin()->second; // return first channel just to prevent warnings
}

// ------------------------------------------------------------------------- //
void DecayTable::SetStable()
{
    clearTable();

    // TODO(mario): Find better way Wed 2017/08/23
    // A stable channel which alwas will be selected
    channels_[1.1] = new StableChannel();
}

// ------------------------------------------------------------------------- //
DecayTable& DecayTable::addChannel(double Br, const DecayChannel& dc)
{
    channels_[Br] = dc.clone();
    return *this;
}

void DecayTable::SetUniformSampling(bool uniform) const
{
    for (DecayMap::const_iterator iter = channels_.begin(); iter != channels_.end(); ++iter)
    {
        iter->second->SetUniformSampling(uniform);
    }
}

// ------------------------------------------------------------------------- //
// private methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
void DecayTable::clearTable()
{
    if (!channels_.empty())
    {
        for (DecayMap::const_iterator iter = channels_.begin(); iter != channels_.end(); ++iter)
        {
            delete iter->second;
        }
    }

    channels_.clear();
}
