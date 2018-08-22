
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

#include <string>
#include <vector>

#include "PROPOSAL/decay/DecayTable.h"
#include "PROPOSAL/methods.h"

#define PARTICLE_DEF(cls)                                                                                              \
    class cls##Def : public ParticleDef                                                                                \
    {                                                                                                                  \
    public:                                                                                                            \
        static const cls##Def& Get()                                                                                   \
        {                                                                                                              \
            static const cls##Def instance;                                                                            \
            return instance;                                                                                           \
        }                                                                                                              \
                                                                                                                       \
    private:                                                                                                           \
        cls##Def();                                                                                                    \
        ~cls##Def();                                                                                                   \
    };

namespace PROPOSAL {

namespace HardComponentTables {

typedef std::vector<std::vector<double> > VecType;

// ----------------------------------------------------------------------------
/// @brief Table used for photonuclear interaction
///
/// Used in parametrizations with real photon assumption.
/// This table is needed when propagation muons.
///
// ----------------------------------------------------------------------------
extern const VecType MuonTable;

// ----------------------------------------------------------------------------
/// @brief Table used for photonuclear interaction
///
/// Used in parametrizations with real photon assumption.
/// This table is needed when propagation taus.
///
// ----------------------------------------------------------------------------
extern const VecType TauTable;

// ----------------------------------------------------------------------------
/// @brief Table used for photonuclear interaction
///
/// Used in parametrizations with real photon assumption.
/// This table is empty and indicates to do no hard component calculations.
///
// ----------------------------------------------------------------------------
extern const VecType EmptyTable;

} // namespace HardComponentTables

// ----------------------------------------------------------------------------
/// @brief Struct to define Basic Particle Properties
///
/// Used to construct Particles
// ----------------------------------------------------------------------------
struct ParticleDef
{
    class Builder;

    const std::string name;
    const double mass;
    const double low; //!< energy below which the particle is lost [MeV]
    const double lifetime;
    const double charge;
    const HardComponentTables::VecType& hard_component_table;
    const DecayTable decay_table;

    ParticleDef();
    ParticleDef(std::string name,
                double mass,
                double low,
                double lifetime,
                double charge,
                const HardComponentTables::VecType& table,
                const DecayTable&);

    ParticleDef(const ParticleDef&);
    virtual ~ParticleDef();

    ParticleDef* clone() const { return new ParticleDef(*this); }

    bool operator==(const ParticleDef&) const;
    bool operator!=(const ParticleDef&) const;

    friend std::ostream& operator<<(std::ostream&, ParticleDef const&);

private:
    ParticleDef& operator=(const ParticleDef&); // Undefined & not allowed
};


std::ostream& operator<<(std::ostream&, PROPOSAL::ParticleDef const&);

class ParticleDef::Builder
{
public:
    Builder();

    // --------------------------------------------------------------------- //
    // Setter
    // --------------------------------------------------------------------- //

    Builder& SetName(const std::string var)
    {
        name = var;
        return *this;
    }
    Builder& SetMass(const double var)
    {
        mass = var;
        return *this;
    }
    Builder& SetLow(const double var)
    {
        low = var;
        return *this;
    }
    Builder& SetLifetime(const double var)
    {
        lifetime = var;
        return *this;
    }
    Builder& SetCharge(const double var)
    {
        charge = var;
        return *this;
    }
    Builder& SetHardComponentTable(const HardComponentTables::VecType& var)
    {
        hard_component_table = &var;
        return *this;
    }
    Builder& SetDecayTable(const DecayTable& var)
    {
        decay_table = var;
        return *this;
    }
    Builder& SetParticleDef(const ParticleDef& var)
    {
        name                 = var.name;
        mass                 = var.mass;
        low                  = var.low;
        lifetime             = var.lifetime;
        charge               = var.charge;
        hard_component_table = &var.hard_component_table;
        decay_table          = var.decay_table;
        return *this;
    }

    ParticleDef build()
    {

        if (low < mass)
        {
            low = mass;
        }

        return ParticleDef(name, mass, low, lifetime, charge, *hard_component_table, decay_table);
    }

private:
    std::string name;
    double mass;
    double low;
    double lifetime;
    double charge;
    const HardComponentTables::VecType* hard_component_table;
    DecayTable decay_table;
};

// ------------------------------------------------------------------------- //
// Predefined particle definitions
// ------------------------------------------------------------------------- //

PARTICLE_DEF(MuMinus)
PARTICLE_DEF(MuPlus)

PARTICLE_DEF(EMinus)
PARTICLE_DEF(EPlus)

PARTICLE_DEF(TauMinus)
PARTICLE_DEF(TauPlus)

PARTICLE_DEF(StauMinus)
PARTICLE_DEF(StauPlus)

PARTICLE_DEF(Pi0)
PARTICLE_DEF(PiMinus)
PARTICLE_DEF(PiPlus)

// ----------------------------------------------------------------------------
/// @brief Particle definition of K0
///
/// K0 is defined with a liftime of -1, meaning stable.
/// The reason is, that PROPOSAL uses this particle soly as
/// return value from decays and does therefore not decide between
/// K short and K long.
///
/// Used to construct Particles
// ----------------------------------------------------------------------------
PARTICLE_DEF(K0)

PARTICLE_DEF(KMinus)
PARTICLE_DEF(KPlus)

PARTICLE_DEF(PMinus)
PARTICLE_DEF(PPlus)

PARTICLE_DEF(NuE)
PARTICLE_DEF(NuEBar)

PARTICLE_DEF(NuMu)
PARTICLE_DEF(NuMuBar)

PARTICLE_DEF(NuTau)
PARTICLE_DEF(NuTauBar)

PARTICLE_DEF(Monopole)
PARTICLE_DEF(Gamma)

PARTICLE_DEF(SMPMinus)
PARTICLE_DEF(SMPPlus)

} // namespace PROPOSAL


PROPOSAL_MAKE_HASHABLE(PROPOSAL::ParticleDef, t.mass, t.lifetime, t.charge)

#undef PARTICLE_DEF
