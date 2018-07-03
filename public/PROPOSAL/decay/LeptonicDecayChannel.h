
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

#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/particle/ParticleDef.h"

namespace PROPOSAL {

class Particle;

class LeptonicDecayChannelApprox : public DecayChannel
{
public:
    LeptonicDecayChannelApprox(const ParticleDef&, const ParticleDef&, const ParticleDef&);
    LeptonicDecayChannelApprox(const LeptonicDecayChannelApprox& mode);
    virtual ~LeptonicDecayChannelApprox();
    // No copy and assignemnt -> done by clone
    DecayChannel* clone() const { return new LeptonicDecayChannelApprox(*this); }

    virtual DecayProducts Decay(const Particle&);

    const std::string& GetName() const { return name_; }

protected:
    ParticleDef massive_lepton_;
    ParticleDef neutrino_;
    ParticleDef anti_neutrino_;
    static const std::string name_;

    LeptonicDecayChannelApprox& operator=(const LeptonicDecayChannelApprox&); // Undefined & not allowed

    bool compare(const DecayChannel&) const;
    void print(std::ostream&) const;

    // ----------------------------------------------------------------------------
    /// @brief Function for electron energy calculation - interface to FindRoot
    // ----------------------------------------------------------------------------
    virtual double DecayRate(double x, double parent_mass, double E_max, double right_side);

    // ----------------------------------------------------------------------------
    /// @brief Function for electron energy calculation - interface to FindRoot
    // ----------------------------------------------------------------------------
    virtual double DifferentialDecayRate(double x, double parent_mass, double E_max);

    std::pair<double, double> function_and_derivative(double x, double parent_mass, double E_max, double right_side);

    double FindRootBoost(double min, double parent_mass, double E_max, double right_side);
};

class LeptonicDecayChannel : public LeptonicDecayChannelApprox
{
public:
    LeptonicDecayChannel(const ParticleDef&, const ParticleDef&, const ParticleDef&);
    LeptonicDecayChannel(const LeptonicDecayChannel& mode);
    virtual ~LeptonicDecayChannel();

    // No copy and assignemnt -> done by clone
    DecayChannel* clone() const { return new LeptonicDecayChannel(*this); }

    const std::string& GetName() const { return name_; }

private:
    double DecayRate(double x, double parent_mass, double E_max, double right_side);
    double DifferentialDecayRate(double x, double parent_mass, double E_max);

    static const std::string name_;
};

} // namespace PROPOSAL
