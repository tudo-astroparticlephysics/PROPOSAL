
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

#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/particle/ParticleDef.h"

namespace PROPOSAL {

class DynamicData
{
public:
    enum Type
    {
        None = 0,
        Particle,
        Brems,
        DeltaE,
        Epair,
        NuclInt,
        MuPair,
        Hadrons,
        ContinuousEnergyLoss,
        WeakInt
    };

public:
    DynamicData(DynamicData::Type);
    DynamicData(const DynamicData&);
    virtual ~DynamicData();

    friend std::ostream& operator<<(std::ostream&, DynamicData const&);

    // --------------------------------------------------------------------- //
    // Getter & Setter
    // --------------------------------------------------------------------- //

    // Setter
    void SetPosition(const Vector3D& position) { position_ = position; }
    void SetDirection(const Vector3D& direction) { direction_ = direction; }

    virtual void SetEnergy(double energy) { energy_ = energy; }
    void SetParentParticleEnergy(double parent_particle_energy) { parent_particle_energy_ = parent_particle_energy; }
    void SetTime(double time) { time_ = time; }
    void SetPropagatedDistance(double prop_dist) { propagated_distance_ = prop_dist; }

    // Getter
    Type GetTypeId() const { return type_id_; }
    static std::string GetNameFromType(Type);

    Vector3D GetPosition() const { return position_; }
    Vector3D GetDirection() const { return direction_; }

    double GetEnergy() const { return energy_; }
    double GetParentParticleEnergy() const { return parent_particle_energy_; }
    double GetTime() const { return time_; }
    double GetPropagatedDistance() const { return propagated_distance_; }

protected:
    virtual void print(std::ostream&) const {}

    const Type type_id_;

    Vector3D position_;  //!< position coordinates [cm]
    Vector3D direction_; //!< direction vector, angles in [rad]

    double energy_;                 //!< energy [MeV]
    double parent_particle_energy_; //!< energy of the parent particle
    double time_;                   //!< age [sec]
    double propagated_distance_;    //!< propagation distance [cm]
};

// ----------------------------------------------------------------------------
/// @brief This class provides the main particle properties and functions.
///
/// All coordinates, angles and physical values are stored in this class.
// ----------------------------------------------------------------------------
class Particle : public DynamicData
{
public:
    Particle();
    Particle(const Particle&);
    Particle(const ParticleDef&);

    // destructors
    virtual ~Particle() {}

    // Operators
    bool operator==(const Particle& particle) const;
    bool operator!=(const Particle& particle) const;

    // --------------------------------------------------------------------- //
    // Methods
    // --------------------------------------------------------------------- //

    // ----------------------------------------------------------------------------
    /// @brief Copies state data of the given particle
    ///
    /// Copies:
    ///  - energy
    ///  - parent particle energy
    ///  - time
    ///  - propagated distance
    ///  - position, direction
    ///  - position, energy, time at detector entry, exit and closest approach
    ///    points
    ///  - energy lost in the detector
    ///
    /// @param Particle
    // ----------------------------------------------------------------------------
    void InjectState(const Particle&);

    // --------------------------------------------------------------------- //
    // Getter & Setter
    // --------------------------------------------------------------------- //

    // Setter
    void SetEnergy(double energy);
    void SetMomentum(double momentum);

    void SetParentParticleId(int parent_particle_id) { parent_particle_id_ = parent_particle_id; }
    void SetParticleId(int particle_id) { particle_id_ = particle_id; }

    void SetEntryPoint(const Vector3D& entry_point) { entry_point_ = entry_point; }
    void SetEntryTime(double entry_time) { entry_time_ = entry_time; }
    void SetEntryEnergy(double entry_energy) { entry_energy_ = entry_energy; }

    void SetExitPoint(const Vector3D& exit_point) { exit_point_ = exit_point; }
    void SetExitTime(const double exit_time) { exit_time_ = exit_time; }
    void SetExitEnergy(const double exit_energy) { exit_energy_ = exit_energy; }

    void SetClosestApproachPoint(const Vector3D& closest_approach_point)
    {
        closest_approach_point_ = closest_approach_point;
    }
    void SetClosestApproachTime(const double closest_approach_time) { closest_approach_time_ = closest_approach_time; }
    void SetClosestApproachEnergy(const double closest_approach_energy)
    {
        closest_approach_energy_ = closest_approach_energy;
    }

    void SetElost(const double elost) { elost_ = elost; }

    // Getter
    const ParticleDef& GetParticleDef() const { return particle_def_; }
    const DecayTable& GetDecayTable() const { return particle_def_.decay_table; }

    double GetMomentum() const { return momentum_; }
    double GetLow() const { return particle_def_.low; }

    double GetMass() const { return particle_def_.mass; }
    double GetLifetime() const { return particle_def_.lifetime; }
    double GetCharge() const { return particle_def_.charge; }
    std::string GetName() const { return particle_def_.name; }

    // ----------------------------------------------------------------------------
    /// @brief Return interpolation tables used for Photonuclear Crossection
    ///
    /// Only for muons and taus a non empty vector will be returned.
    ///
    /// @return 2d double vector
    // ----------------------------------------------------------------------------
    const HardComponentTables::VecType& getHardComponent() const { return particle_def_.hard_component_table; }

    int GetParentParticleId() const { return parent_particle_id_; }
    int GetParticleId() const { return particle_id_; }

    Vector3D GetEntryPoint() const { return entry_point_; }
    double GetEntryTime() const { return entry_time_; }
    double GetEntryEnergy() const { return entry_energy_; }

    Vector3D GetExitPoint() const { return exit_point_; }
    double GetExitTime() const { return exit_time_; }
    double GetExitEnergy() const { return exit_energy_; }

    Vector3D GetClosestApproachPoint() const { return closest_approach_point_; }
    double GetClosestApproachTime() const { return closest_approach_time_; }
    double GetClosestApproachEnergy() const { return closest_approach_energy_; }

    double GetElost() const { return elost_; }

    // --------------------------------------------------------------------- //

private:
    Particle& operator=(const Particle&);
    virtual void print(std::ostream&) const;

    const ParticleDef particle_def_; //!< static defenitions of the particle

    double momentum_;        //!< momentum [MeV]

    int parent_particle_id_; //!< parent particle id
    int particle_id_;        //!< particle id

    Vector3D entry_point_; //!< entry point coordinates [cm]
    double entry_time_;    //!< time-coordinate entry Point [sec]
    double entry_energy_;  //!< energy at entry point [MeV]

    Vector3D exit_point_; //!< exit point coordinates [cm]
    double exit_time_;    //!< time-coordinate exit Point [sec]
    double exit_energy_;  //!< energy at exit point [MeV]

    Vector3D closest_approach_point_; // point of closest approach (to geometry center) [cm]
    double closest_approach_time_;    //!< time-coordinate at point of closest approach [sec]
    double closest_approach_energy_;  //!< energy at at point of closest approach [MeV]

    double elost_; //!< energy lost in the detector volume [MeV]
};

std::ostream& operator<<(std::ostream&, PROPOSAL::DynamicData const&);

} // namespace PROPOSAL
