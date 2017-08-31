
#include "PROPOSAL/decay/DecayTable.h"
#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/decay/StableChannel.h"
#include "PROPOSAL/math/MathModel.h"
#include "PROPOSAL/Output.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructors & destructor
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
DecayTable::~DecayTable()
{
    clearTable();
}

// ------------------------------------------------------------------------- //
DecayTable::DecayTable(const DecayTable& table)
{
    clearTable();

    for (DecayMap::const_iterator iter = table.channels.begin(); iter != table.channels.end(); ++iter)
    {
        channels[iter->first] = iter->second->clone();
    }
}

DecayTable& DecayTable::operator=(const DecayTable& table)
{
    if (this != &table)
    {
      DecayTable tmp(table);
      swap(tmp);
    }
    return *this;
}

void DecayTable::swap(DecayTable& table)
{
    using std::swap;
    swap(channels, table.channels);
}

bool DecayTable::operator==(const DecayTable& table) const
{
    if (channels.size() != table.channels.size())
    {
        return false;
    }
    else
    {
        DecayMap::const_iterator i, j;
        for (i = channels.begin(), j = table.channels.begin(); i != channels.end(); ++i, ++j)
        {
            if (i->first != j->first)
            {
                return false;
            }
            else if (*i->second != *j->second)
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
// Methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
DecayChannel& DecayTable::SelectChannel() const
{
    for (int i = 0; i < 1000; ++i)
    {
        double sumBranchingRatio = 0.0;
        for (DecayMap::const_iterator iter = channels.begin(); iter != channels.end(); ++iter)
        {
            sumBranchingRatio += iter->first;

            if (RandomGenerator::Get().RandomDouble() < sumBranchingRatio)
            {
                return *iter->second;
            }
        }
    }

    log_fatal("No decay channel found. If your particle is stable, call \"SetStable\"!");
}

void DecayTable::SetStable()
{
    clearTable();

    //TODO(mario): Find better way Wed 2017/08/23
    // A stable channel which alwas will be selected
    channels[1.1] = new StableChannel();
}

// ------------------------------------------------------------------------- //
void DecayTable::addChannel(double Br, DecayChannel& dc)
{
    channels[Br] = dc.clone();
}

void DecayTable::clearTable()
{
    if (!channels.empty())
    {
        for (DecayMap::const_iterator iter = channels.begin(); iter != channels.end(); ++iter)
        {
            delete iter->second;
        }
    }

    channels.clear();
}
