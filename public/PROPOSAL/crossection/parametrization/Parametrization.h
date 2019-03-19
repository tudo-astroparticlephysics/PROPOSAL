
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
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/medium/Medium.h"

namespace PROPOSAL {

namespace Components {
class Component;
}

class Medium;

class Parametrization
{
public:
    Parametrization(const ParticleDef&, const Medium&, const EnergyCutSettings&, double multiplier);
    Parametrization(const Parametrization&);
    virtual ~Parametrization();

    bool operator==(const Parametrization& parametrization) const;
    bool operator!=(const Parametrization& parametrization) const;

    virtual Parametrization* clone() const = 0;

    friend std::ostream& operator<<(std::ostream&, Parametrization const&);

    // bounds of integration
    struct IntegralLimits
    {
        double vMax; //!< upper bound of integration
        double vUp;  //!< lower bound of integration
        double vMin; //!< lowest physical possible bound of integration
    };

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    virtual double DifferentialCrossSection(double energy, double v) = 0;

    double FunctionToDNdxIntegral(double energy, double v);
    virtual double FunctionToDEdxIntegral(double energy, double v);
    double FunctionToDE2dxIntegral(double energy, double v);

    virtual double Calculaterho(double energy, double v, double rnd1, double rnd2){
        (void)energy; (void)v; (void)rnd1; (void)rnd2; return 0;}

    virtual IntegralLimits GetIntegralLimits(double energy) = 0;

    // ----------------------------------------------------------------- //
    // Getter
    // ----------------------------------------------------------------- //

    virtual const std::string& GetName() const = 0; //{ return name_; }

    const ParticleDef& GetParticleDef() const { return particle_def_; }
    const Medium& GetMedium() const { return *medium_; }
    const EnergyCutSettings& GetEnergyCuts() const { return cut_settings_; }
    double GetMultiplier() const { return multiplier_; }

    virtual size_t GetHash() const;

    // ----------------------------------------------------------------- //
    // Setter
    // ----------------------------------------------------------------- //

    // void SetCurrentComponent(Components::Component* component) {current_component_ = component;}
    void SetCurrentComponent(int index) { component_index_ = index; }

protected:
    typedef std::vector<Components::Component*> ComponentVec;

    virtual bool compare(const Parametrization&) const;
    virtual void print(std::ostream&) const {};

    // const std::string name_;

    const ParticleDef particle_def_;
    const Medium* medium_;
    const EnergyCutSettings cut_settings_;

    // const Components::Component* current_component_;
    const ComponentVec& components_;
    int component_index_;

    double multiplier_;
};

std::ostream& operator<<(std::ostream&, PROPOSAL::Parametrization const&);

} // namespace PROPOSAL
