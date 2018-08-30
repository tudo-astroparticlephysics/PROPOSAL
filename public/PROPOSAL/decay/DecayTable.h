
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/


#pragma once

#include <map>
#include <ostream>
#include <vector>

namespace PROPOSAL {

class DecayChannel;
class DecayTable;

void swap(DecayTable&, DecayTable&);

class DecayTable
{
    typedef std::map<double, DecayChannel*> DecayMap;

public:
    DecayTable();
    DecayTable(const DecayTable& table);
    virtual ~DecayTable();

    DecayTable& operator=(const DecayTable&);
    friend void swap(DecayTable&, DecayTable&);

    bool operator==(const DecayTable&) const;
    bool operator!=(const DecayTable&) const;

    friend std::ostream& operator<<(std::ostream&, DecayTable const&);

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
    DecayTable& addChannel(double Br, const DecayChannel& dc);

    // ----------------------------------------------------------------------------
    /// @brief Provide a decay table for stable particles
    ///
    /// Clears the decay table and adds one stable decay channel with
    /// branching ration 1.1 therefor it always will be selected.
    /// The decay channel will return an empty decay product vector.
    // ----------------------------------------------------------------------------
    void SetStable();

    // ----------------------------------------------------------------------------
    /// @brief Sets the uniform flag in the ManyBodyPhaseSpace channels
    ///
    /// If uniform is true, the momenta will be sampled uniform in the phase space.
    /// This is done by rejection, since the pure raubold lynch algorithm does not
    /// create a uniform phase space. So enabling uniform sampling comes in with
    /// a cost of performance.
    ///
    /// @param uniform
    // ----------------------------------------------------------------------------
    void SetUniformSampling(bool uniform) const;

private:
    void clearTable();

    DecayMap channels_;
};

std::ostream& operator<<(std::ostream&, PROPOSAL::DecayTable const&);

} // namespace PROPOSAL
