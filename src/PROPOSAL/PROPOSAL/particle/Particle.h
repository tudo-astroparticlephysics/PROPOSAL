
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
#include <memory>
#include <string>

#include "PROPOSAL/math/Cartesian3D.h"
#include "PROPOSAL/particle/ParticleDef.h"

namespace PROPOSAL {
enum class InteractionType : int {
    Undefined = 0,
    Particle = 1000000001,
    Brems = 1000000002,
    Ioniz = 1000000003,
    Epair = 1000000004,
    Photonuclear = 1000000005,
    MuPair = 1000000006,
    Hadrons = 1000000007,
    ContinuousEnergyLoss = 1000000008,
    WeakInt = 1000000009,
    Compton = 1000000010,
    Decay = 1000000011,
    Annihilation = 1000000012,
    Photopair = 1000000013,
};
struct InteractionType_hash {
    template <class T> std::size_t operator()(const T& type) const
    {
        return static_cast<int>(type);
    }
};
} // namespace PROPOSAL


namespace PROPOSAL {
static const std::unordered_map<InteractionType, std::string,
    InteractionType_hash>
    Type_Interaction_Name_Map{
        { InteractionType::Particle, "Particle" },
        { InteractionType::Brems, "Brems" },
        { InteractionType::Ioniz, "Ioniz" },
        { InteractionType::Epair, "Epair" },
        { InteractionType::Photonuclear, "Photonuclear" },
        { InteractionType::MuPair, "MuPair" },
        { InteractionType::Hadrons, "Hadrons" },
        { InteractionType::ContinuousEnergyLoss, "ContinuousEnergyLoss" },
        { InteractionType::WeakInt, "WeakInt" },
        { InteractionType::Compton, "Compton" },
        { InteractionType::Decay, "Decay" },
        { InteractionType::Annihilation, "Annihilation" },
        { InteractionType::Photopair, "Photopair" },
    };
} // namespace PROPOSAL

namespace PROPOSAL {

struct ParticleState {
public:
    ParticleState();
    ParticleState(const Vector3D&, const Vector3D&, const double&,
                  const double&, const double&);
    ParticleState(const ParticleType&, const Vector3D&, const Vector3D&, const double&,
                  const double&, const double&);
    ~ParticleState() = default;
    bool operator==(const ParticleState&) const;
    bool operator!=(const ParticleState&) const;
    friend std::ostream& operator<<(std::ostream&, ParticleState const&);

    int type;
    Cartesian3D position;  //!< position coordinates [cm]
    Cartesian3D direction; //!< direction vector, angles in [rad]
    double energy;                 //!< energy [MeV]
    double time;                   //!< age [sec]
    double propagated_distance;    //!< propagation distance [cm]

    void SetType(ParticleType particle_type) { type = static_cast<int>(particle_type); }
    void SetMomentum(double momentum);
    double GetMomentum() const;
    ParticleDef GetParticleDef() const;

private:
    void print(std::ostream&) const {}

};

struct Loss {
    Loss(int type, double energy, double parent_particle_energy)
        : type(type), energy(energy), parent_particle_energy(parent_particle_energy) {};
    int type;
    double energy;
    double parent_particle_energy;
};

struct StochasticLoss : public Loss {
    StochasticLoss(int, double, const Vector3D&, const Vector3D&, double, double, double);
    Cartesian3D position;
    Cartesian3D direction;
    double time;
    double propagated_distance;
};

struct ContinuousLoss : public Loss {
    ContinuousLoss(double, double, const Vector3D&, double, const Vector3D&,
            const Vector3D&, double, double);
    Cartesian3D start_position;
    double length;
    Cartesian3D direction_initial;
    Cartesian3D direction_final;
    double time_initial;
    double time_final;
};
} // namespace PROPOSAL
