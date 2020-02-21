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

#include <memory>
#include <string>
#include <vector>

#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/math/Vector3D.h"

namespace PROPOSAL {

class Geometry;

class Secondaries {

public:
    Secondaries();
    Secondaries(std::shared_ptr<ParticleDef>);

    void reserve(size_t number_secondaries);
    void clear() { secondaries_.clear(); };

    DynamicData& operator[](std::size_t idx) { return secondaries_[idx]; };

    void push_back(const DynamicData& continuous_loss);
    void emplace_back(const int& type);
    void emplace_back(const int& type, const Vector3D& position,
        const Vector3D& direction, const double& energy,
        const double& parent_particle_energy, const double& time,
        const double& distance);

    void append(Secondaries& secondaries);

    Secondaries Query(const int&) const;
    Secondaries Query(const std::string&) const;
    Secondaries Query(const Geometry& geometry) const;

    void DoDecay();

    std::vector<Vector3D> GetPosition() const;
    std::vector<Vector3D> GetDirection() const;
    std::vector<double> GetEnergy() const;
    std::vector<double> GetParentParticleEnergy() const;
    std::vector<double> GetTime() const;
    std::vector<double> GetPropagatedDistance() const;
    std::vector<DynamicData> GetSecondaries() const { return secondaries_; };
    std::vector<DynamicData>& GetModifyableSecondaries()
    {
        return secondaries_;
    };
    unsigned int GetNumberOfParticles() const { return secondaries_.size(); };
    Secondaries GetOnlyLostInsideDetector() const;

    // TODO: Prelimary, see note below
    double GetELost() const;
    DynamicData GetEntryPoint() const;
    DynamicData GetExitPoint() const;
    DynamicData GetClosestApproachPoint() const;
    void SetEntryPoint(const DynamicData& entry_point);
    void SetExitPoint(const DynamicData& exit_point);
    void SetClosestApproachPoint(const DynamicData& closest_approach_point);

private:
    std::vector<DynamicData> secondaries_;
    std::shared_ptr<ParticleDef> primary_def_;

    // TODO: Entry and Exit point must not necessary be saved.
    // It can be calculated by a given structure
    std::unique_ptr<DynamicData> entry_point_;
    std::unique_ptr<DynamicData> exit_point_;
    std::unique_ptr<DynamicData> closest_approach_point_;
};

} // namespace PROPOSAL
