
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

#include "PROPOSAL/medium/Components.h"

using std::string;

namespace PROPOSAL {
class ParticleDef;

class Parametrization {
protected:
    const string param_name_;
    const double particle_mass_;
    const double particle_charge_;
    const component_list components_;
    Components::Component current_component_;
    const double lower_energy_lim_;

public:
    Parametrization(const string&, const ParticleDef&, const component_list&, double);
    virtual ~Parametrization() = default;

    struct KinematicLimits {
        double vMin;
        double vMax;

        KinematicLimits(double, double);
    };

    virtual double DifferentialCrossSection(double, double) = 0;
    virtual double FunctionToDNdxIntegral(double, double);
    virtual double FunctionToDEdxIntegral(double, double);
    virtual double FunctionToDE2dxIntegral(double, double);

    virtual KinematicLimits GetKinematicLimits(double energy) = 0;
    double GetLowerEnergyLim() const { return lower_energy_lim_; }
    string GetName() const { return param_name_; }
    component_list GetComponents() const { return components_; }
    virtual size_t GetHash() const;

    void SetCurrentComponent(Components::Component&);
};


} // namespace PROPOSAL
