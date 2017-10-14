
#include "PROPOSAL/decay/DecayTable.h"
#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/decay/StableChannel.h"
#include "PROPOSAL/decay/LeptonicDecayChannel.h"
#include "PROPOSAL/decay/TwoBodyPhaseSpace.h"
#include "PROPOSAL/math/MathModel.h"

#include "PROPOSAL/Output.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructors & destructor
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
DecayTable::DecayTable()
    : registered_mode_()
    , registered_str_()
    , channels_()
{
    Register("leptonic_chanel", LeptonicDecay);
    Register("two_body_phase_space", TwoBodyDecay);
    Register("stable", Stable);
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

    registered_mode_.clear();
    registered_str_.clear();
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
    }
    else
    {
        DecayMap::const_iterator i, j;
        for (i = channels_.begin(), j = table.channels_.begin(); i != channels_.end(); ++i, ++j)
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
std::ostream& PROPOSAL::operator<<(std::ostream& os, DecayTable const& table)
{
    std::stringstream ss;
    ss << " DecayTable (" << &table << ") ";

    os << Helper::Centered(60, ss.str()) << '\n';

    os << "Branching Ratio" << "\t\t"  << "Channel" << '\n';
    for (DecayTable::DecayMap::const_iterator iter = table.channels_.begin(); iter != table.channels_.end(); ++iter) {
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
        for (DecayMap::const_iterator iter = channels_.begin(); iter != channels_.end(); ++iter)
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
    channels_[1.1] = new StableChannel();
}

// ------------------------------------------------------------------------- //
void DecayTable::addChannel(double Br, DecayChannel& dc)
{
    channels_[Br] = dc.clone();
}

// ------------------------------------------------------------------------- //
// void DecayTable::addChannel(double Br, Mode mode, double parent_mass, const std::vector<double>& masses)
// {
//     for(unsigned int i = 0; i < masses.size(); ++i)
//     {
//         std::cout << masses[i] << std::endl;
//     }
//     std::vector<Mode>::const_iterator iter;
//     iter = std::find(registered_mode_.begin(), registered_mode_.end(), mode);
//
//     if (iter != registered_mode_.end())
//     {
//         if (*iter == LeptonicDecay)
//         {
//             channels_[Br] = new LeptonicDecayChannel();
//         }
//         else if (*iter == TwoBodyDecay)
//         {
//             try
//             {
//                 channels_[Br] = new TwoBodyPhaseSpace(parent_mass, masses.at(0));
//             }
//             catch (const std::out_of_range& ex)
//             {
//                 log_fatal("No secondary masses are given!");
//             }
//         }
//         else if (*iter == Stable)
//         {
//             channels_[Br] = new StableChannel();
//         }
//         else
//         {
//             log_fatal("DecayChannel %s not registerd!", typeid(mode).name());
//         }
//     }
//     else
//     {
//         log_fatal("DecayChannel %s not registerd!", typeid(mode).name());
//     }
// }

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

// ------------------------------------------------------------------------- //
void DecayTable::Register(const std::string& name, Mode model)
{
    registered_str_.push_back(name);
    registered_mode_.push_back(model);
}
