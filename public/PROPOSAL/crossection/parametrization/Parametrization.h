
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

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/json.hpp"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include <functional>

namespace PROPOSAL {

namespace Components {
    class Component;
}

class Parametrization {
public:
    Parametrization(
        const ParticleDef&, std::shared_ptr<const Medium>, double multiplier);
    Parametrization(const Parametrization&);
    virtual ~Parametrization();

    bool operator==(const Parametrization& parametrization) const;
    bool operator!=(const Parametrization& parametrization) const;

    virtual Parametrization* clone() const {};

    friend std::ostream& operator<<(std::ostream&, Parametrization const&);

    // bounds of integration
    struct KinematicLimits {
        double vMin; //!< lower physical, kinematic limit
        double vMax; //!< upper physical, kinematic limit
    };

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    virtual double DifferentialCrossSection(double energy, double v) = 0;

    double FunctionToDNdxIntegral(double energy, double v);
    virtual double FunctionToDEdxIntegral(double energy, double v);
    double FunctionToDE2dxIntegral(double energy, double v);

    virtual double Calculaterho(
        double energy, double v, double rnd1, double rnd2)
    {
        (void)energy;
        (void)v;
        (void)rnd1;
        (void)rnd2;
        return 0;
    }

    virtual KinematicLimits GetKinematicLimits(double energy) = 0;

    // ----------------------------------------------------------------- //
    // Getter
    // ----------------------------------------------------------------- //

    virtual const std::string& GetName() const = 0; //{ return name_; }
    std::shared_ptr<const Medium> GetMedium() const { return medium_; }
    double GetParticleMass() const { return particle_mass_; }
    double GetParticleCharge() const { return particle_charge_; }
    double GetParticleLow() const { return particle_low_; }

    virtual InteractionType GetInteractionType() const = 0;
    double GetMultiplier() const { return multiplier_; }

    virtual size_t GetHash() const;

    // ----------------------------------------------------------------- //
    // Setter
    // ----------------------------------------------------------------- //

    // void SetCurrentComponent(Components::Component* component)
    // {current_component_ = component;}
    void SetCurrentComponent(int index) { component_index_ = index; }

protected:
    virtual bool compare(const Parametrization&) const;
    virtual void print(std::ostream&) const {};

    // const std::string name_;

    double particle_mass_;
    double particle_charge_;
    double particle_low_;
    std::shared_ptr<const Medium> medium_;

    // const Components::Component* current_component_;
    const std::vector<Components::Component>& components_;
    int component_index_;

    double multiplier_;
};

class Parametrization_builder : public Parametrization {
    std::function<double(const Parametrization&, double, double)> diff_cross;
    std::function<KinematicLimits(const Parametrization&, double)> kinematic_limits;

public:
    Parametrization_builder(const ParticleDef& p_def,
        std::shared_ptr<const Medium> medium, double multiplier)
        : Parametrization(p_def, medium, multiplier)
        , diff_cross(nullptr)
        , kinematic_limits(nullptr){};

    double DifferentialCrossSection(double energy, double v) override;
    KinematicLimits GetKinematicLimits(double energy) override;
};

std::ostream& operator<<(std::ostream&, PROPOSAL::Parametrization const&);

size_t GetHash(const std::vector<Parametrization*>& params);

} // namespace PROPOSAL
