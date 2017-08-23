
#pragma once

#include <iostream>
#include <map>

#include "PROPOSAL/Decay.h"


namespace PROPOSAL
{

class DecayTable
{
    typedef std::map<double, DecayChannel*> DecayMap;

    public:
        DecayTable() {}
        DecayTable(const DecayTable& table);
        virtual ~DecayTable();

        DecayTable& operator=(const DecayTable&);
        void swap(DecayTable&);

        bool operator==(const DecayTable&) const;
        bool operator!=(const DecayTable&) const;

        // ----------------------------------------------------------------------------
        /// @brief Get a decay channel
        ///
        /// The Decay channels will be sampled from the previous given branching ratios
        ///
        /// @return Sampled Decay channel
        // ----------------------------------------------------------------------------
        DecayChannel& SelectChannel() const;

        // ----------------------------------------------------------------------------
        /// @brief Add decay channels to the decay table
        ///
        /// The given decaychannel will be copied and stored in a map along with
        /// the branching ration.
        ///
        /// @param Br Branching ration of the channel
        /// @param dc the decay channel
        // ----------------------------------------------------------------------------
        void addChannel(double Br, DecayChannel& dc);

        // ----------------------------------------------------------------------------
        /// @brief Provide a decay table for stable particles
        ///
        /// Clears the decay table and adds one stable decay channel with
        /// branching ration 1.1 therefor it always will be selected.
        /// The decay channel will return an empty decay product vector.
        // ----------------------------------------------------------------------------
        void SetStable();

    private:
        void clearTable();
        DecayMap channels;
};

} /* PROPOSAL */
