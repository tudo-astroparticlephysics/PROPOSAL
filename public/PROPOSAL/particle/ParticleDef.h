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

} // namespace HardBBTable

// ----------------------------------------------------------------------------
/// @brief Struct to define Basic Particle Properties
///
/// Used to construct Particles
// ----------------------------------------------------------------------------
struct ParticleDef
{
    std::string name;
    double mass;
    double low;       //!< energy below which the particle is lost [MeV]
    double lifetime;
    double charge;
    const HardBBTables::VecType* hardbb_table;
    DecayTable decay_table;

    ParticleDef();
    ParticleDef(std::string name,
                double mass,
                double lifetime,
                double charge,
                const HardBBTables::VecType* table = NULL);

    ParticleDef(const ParticleDef&);
    virtual ~ParticleDef();

    ParticleDef& operator=(const ParticleDef&);
    void swap(ParticleDef&);

    bool operator==(const ParticleDef&) const;
    bool operator!=(const ParticleDef&) const;

    friend std::ostream& operator<<(std::ostream&, ParticleDef const&);
};


// ----------------------------------------------------------------------------
/// @brief Used to hash ParticleDefs in PropagatorService hash table
//
/// Hash will currently only done by mass, lifetime and charge
// ----------------------------------------------------------------------------
std::size_t hash_value(ParticleDef const& particle_def);

// ------------------------------------------------------------------------- //
// Create default particle definitions
// ------------------------------------------------------------------------- //


PARTICLE_DEF(MuMinus)
PARTICLE_DEF(MuPlus)

PARTICLE_DEF(EMinus)
PARTICLE_DEF(EPlus)

PARTICLE_DEF(TauMinus)
PARTICLE_DEF(TauPlus)

PARTICLE_DEF(StauMinus)
PARTICLE_DEF(StauPlus)

PARTICLE_DEF(P0)
PARTICLE_DEF(PiMinus)
PARTICLE_DEF(PiPlus)

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

PARTICLE_DEF(DeltaE)
PARTICLE_DEF(Brems)
PARTICLE_DEF(NuclInt)
PARTICLE_DEF(Hadrons)
PARTICLE_DEF(Epair)
PARTICLE_DEF(Mupair)

PARTICLE_DEF(ContinuousEnergyLoss)
PARTICLE_DEF(Monopole)
PARTICLE_DEF(Gamma)
PARTICLE_DEF(StableMassiveParticle)

#undef PARTICLE_DEF
} // namespace PROPOSAL
