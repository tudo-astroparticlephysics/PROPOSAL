
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

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/scattering/multiple_scattering/Highland.h"

namespace PROPOSAL {
namespace multiple_scattering {

    template <class T, class Cross> class HighlandIntegral : public Highland {
        T highland_integral;

        T BuildHighlandIntegral(Cross&&);
        double CalculateTheta0(double, double, double) override;

    public:
        HighlandIntegral(const ParticleDef&, Medium const&, Cross&&);

        std::unique_ptr<Parametrization> clone() const override
        {
            return std::unique_ptr<Parametrization>(
                std::make_unique<HighlandIntegral<T, Cross>>(*this));
        }

        static Interpolant1DBuilder::Definition interpol_def;

        double Integral(Displacement&, Cross, double);
    };

    template <class T, class Cross>
    HighlandIntegral<T, Cross>::HighlandIntegral(
        const ParticleDef& p_def, Medium const& medium, Cross&& cross)
        : Highland(p_def, medium)
        , highland_integral(BuildHighlandIntegral(cross))
    {
    }

    template <class T> struct HighlandIntegral<T, std::nullptr_t> {
        HighlandIntegral(const ParticleDef&, Medium const&, std::nullptr_t)
        {
            throw std::invalid_argument("CrossSectionVector needs to be passed "
                                        "to use scattering_highland.");
        }
        std::unique_ptr<Parametrization> clone() const { return nullptr; }
    };

    template <class T, class Cross>
    T HighlandIntegral<T, Cross>::BuildHighlandIntegral(Cross&& cross)
    {
        auto disp
            = std::shared_ptr<Displacement>(make_displacement(cross, false));
        auto higland_integral_func = [this, disp, &cross](double energy) {
            return Integral(*disp, cross, energy);
        };
        T decay_integral(
            higland_integral_func, CrossSectionVector::GetLowerLim(cross));
        if (typeid(T) == typeid(UtilityInterpolant)) {
            auto hash = CrossSectionVector::GetHash(cross);
            decay_integral.BuildTables("highland", hash, interpol_def, true);
        };
        return decay_integral;
    }

    template <class T, class Cross>
    double HighlandIntegral<T, Cross>::Integral(
        Displacement& disp, Cross, double energy)
    {
        auto square_momentum = (energy - mass) * (energy + mass);
        auto aux = energy / square_momentum;
        return disp.FunctionToIntegral(energy) * aux * aux;
    }

    template <class T, class Cross>
    double HighlandIntegral<T, Cross>::CalculateTheta0(
        double grammage, double ei, double ef)
    {
        auto aux = 13.6
            * std::sqrt(highland_integral.Calculate(ei, ef) / radiation_length)
            * std::abs(charge);
        aux *= std::max(1 + 0.038 * std::log(grammage / radiation_length), 0.0);
        return std::min(aux, 1.0);
    }

    template <class T, class Cross>
    Interpolant1DBuilder::Definition HighlandIntegral<T, Cross>::interpol_def
        = { 200 };
} // namespace multiple_scattering
} // namespace PROPOSAL
