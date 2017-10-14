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

    // ParticleDef& operator=(const ParticleDef&);
    // void swap(ParticleDef&);

    bool operator==(const ParticleDef&) const;
    bool operator!=(const ParticleDef&) const;

    friend std::ostream& operator<<(std::ostream&, ParticleDef const&);
};

class ParticleDef::Builder
{
    public:
    // default values for variables
    static const std::string default_name;
    static const double default_mass;
    static const double default_low;
    static const double default_lifetime;
    static const double default_charge;
    static const HardBBTables::VecType* default_hardbb_table;
    // static const DecayTable default_decay_table;

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

    // --------------------------------------------------------------------- //
    // Predefined particle
    // --------------------------------------------------------------------- //

    Builder& SetMuMinus();
    Builder& SetMuPlus();

    Builder& SetEMinus();
    Builder& SetEPlus();

    Builder& SetTauMinus();
    Builder& SetTauPlus();

    Builder& SetStauMinus();
    Builder& SetStauPlus();

    Builder& SetP0();
    Builder& SetPiMinus();
    Builder& SetPiPlus();

    Builder& SetKMinus();
    Builder& SetKPlus();

    Builder& SetPMinus();
    Builder& SetPPlus();

    Builder& SetNuE();
    Builder& SetNuEBar();

    Builder& SetNuMu();
    Builder& SetNuMuBar();

    Builder& SetNuTau();
    Builder& SetNuTauBar();

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

} // namespace PROPOSAL
