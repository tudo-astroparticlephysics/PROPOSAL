
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

#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/particle/ParticleDef.h"

namespace PROPOSAL {
enum class InteractionType : int {
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
class UtilityDecorator;

class DynamicData {
public:
    DynamicData();
    DynamicData(const int&);
    DynamicData(const int&, const Vector3D&, const Vector3D&, const double&,
        const double&, const double&, const double&);

    friend std::ostream& operator<<(std::ostream&, DynamicData const&);
    DynamicData& operator=(const DynamicData&);
    bool operator==(const DynamicData&) const;
    bool operator!=(const DynamicData&) const;

    // --------------------------------------------------------------------- //
    // Getter & Setter
    // --------------------------------------------------------------------- //

    // Setter
    void SetPosition(const Vector3D& position) { position_ = position; }
    void SetDirection(const Vector3D& direction) { direction_ = direction; }

    void SetParentParticleEnergy(double parent_particle_energy)
    {
        parent_particle_energy_ = parent_particle_energy;
    }
    void SetTime(double time) { time_ = time; }
    void SetPropagatedDistance(double prop_dist)
    {
        propagated_distance_ = prop_dist;
    }
    void SetType(InteractionType type) { type_ = static_cast<int>(type); }
    void SetType(int type) { type_ = type; }

    // Getter
    int GetType() const { return type_; }

    Vector3D GetPosition() const { return position_; }
    Vector3D GetDirection() const { return direction_; }

    void SetEnergy(double energy);
    void SetMomentum(double momentum);

    double GetEnergy() const { return energy_; }
    double GetMomentum() const;
    double GetParentParticleEnergy() const { return parent_particle_energy_; }
    double GetTime() const { return time_; }
    double GetPropagatedDistance() const { return propagated_distance_; }
    std::string GetName() const;

    // Deflect the direction of a particle by cosphi_deflect with an azimuth of
    // theta_deflect
    void DeflectDirection(double cosphi_deflect, double theta_deflect);

protected:
    virtual void print(std::ostream&) const {}

    int type_;

    Vector3D position_;  //!< position coordinates [cm]
    Vector3D direction_; //!< direction vector, angles in [rad]

    double energy_;                 //!< energy [MeV]
    double parent_particle_energy_; //!< energy of the parent particle
    double time_;                   //!< age [sec]
    double propagated_distance_;    //!< propagation distance [cm]

    std::shared_ptr<UtilityDecorator> utility_decorator;
};

struct Loss {
    enum { TYPE, POSITION, DIRECTION, ENERGY, TIME };
    using secondary_t = std::tuple<int, Vector3D, Vector3D, double, double>;
};
} // namespace PROPOSAL
