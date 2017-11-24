/*! \file   ParticleDef.h
*   \brief  Header file for definition of the ParticleDef class object.
*
*   For more details see the class documentation.
*
*   \date   Sat Aug  5 14:47:16 CEST 2017
*   \author Mario Dunsch
*/

#pragma once

#include <boost/functional/hash.hpp>

#include <string>
#include <vector>

#include "PROPOSAL/decay/DecayTable.h"

#define PARTICLE_DEF(cls)                                                                                              \
    class cls##Def : public ParticleDef                                                                                \
    {                                                                                                                  \
        public:                                                                                                        \
        static const cls##Def& Get()                                                                                   \
        {                                                                                                              \
            static const cls##Def instance;                                                                            \
            return instance;                                                                                           \
        }                                                                                                              \
                                                                                                                       \
        private:                                                                                                       \
        cls##Def();                                                                                                    \
        ~cls##Def();                                                                                                   \
};

namespace PROPOSAL {

// class DecayTable;

namespace HardBBTables {

typedef std::vector<std::vector<double> > VecType;

extern const double muon[8][7];
extern const double tau[8][7];

VecType getBBVector(const double table[8][7]);

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
/// This table is empty and indicates to do no hardbb calculations.
///
// ----------------------------------------------------------------------------
extern const VecType EmptyTable;

} // namespace HardBBTable

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
    const double low;       //!< energy below which the particle is lost [MeV]
    const double lifetime;
    const double charge;
    const HardBBTables::VecType& hardbb_table;
    const DecayTable decay_table;

    ParticleDef();
    ParticleDef(std::string name,
                double mass,
                double low,
                double lifetime,
                double charge,
                const HardBBTables::VecType& table,
                const DecayTable&);

    ParticleDef(const ParticleDef&);
    virtual ~ParticleDef();

    ParticleDef* clone() const { return new ParticleDef(*this); }

    // ParticleDef& operator=(const ParticleDef&);
    // void swap(ParticleDef&);

    bool operator==(const ParticleDef&) const;
    bool operator!=(const ParticleDef&) const;

    friend std::ostream& operator<<(std::ostream&, ParticleDef const&);

    private:
    ParticleDef& operator=(const ParticleDef&); // Undefined & not allowed
};

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
    Builder& SetHardBBTable(const HardBBTables::VecType& var)
    {
        hardbb_table = &var;
        return *this;
    }
    Builder& SetDecayTable(const DecayTable& var)
    {
        decay_table = var;
        return *this;
    }
    Builder& SetParticleDef(const ParticleDef& var)
    {
        name = var.name;
        mass = var.mass;
        low = var.low;
        lifetime = var.lifetime;
        charge = var.charge;
        hardbb_table = &var.hardbb_table;
        decay_table = var.decay_table;
        return *this;
    }

    ParticleDef build()
    {

        if (low < mass)
        {
            low = mass;
        }

        return ParticleDef(name, mass, low, lifetime, charge, *hardbb_table, decay_table);
    }

    private:
    std::string name;
    double mass;
    double low;
    double lifetime;
    double charge;
    const HardBBTables::VecType* hardbb_table;
    DecayTable decay_table;
};

// ----------------------------------------------------------------------------
/// @brief Used to hash ParticleDefs in PropagatorService hash table
//
/// Hash will currently only done by mass, lifetime and charge
// ----------------------------------------------------------------------------
std::size_t hash_value(ParticleDef const& particle_def);

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
PARTICLE_DEF(StableMassiveParticle)

} // namespace PROPOSAL

#undef PARTICLE_DEF
