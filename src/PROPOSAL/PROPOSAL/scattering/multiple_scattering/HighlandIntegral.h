
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
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include "PROPOSAL/scattering/multiple_scattering/Highland.h"

namespace PROPOSAL {
namespace multiple_scattering {
    template <class T> class HighlandIntegral : public Highland {
        T highland_integral;

        inline double CalculateTheta0(
            double grammage, double ei, double ef) final
        {
            auto aux = 13.6
                * std::sqrt(
                    highland_integral.Calculate(ei, ef) / radiation_length)
                * std::abs(charge);
            aux *= std::max(
                1. + 0.038 * std::log(grammage / radiation_length), 0.0);
            return std::min(aux, 1.0);
        }

        inline double Integral(Displacement& disp, double energy)
        {
            auto square_momentum = (energy - mass) * (energy + mass);
            auto aux = energy / square_momentum;
            return disp.FunctionToIntegral(energy) * aux * aux;
        }

        inline std::unique_ptr<Parametrization> clone() const override
        {
            return std::unique_ptr<Parametrization>(
                std::make_unique<HighlandIntegral<T>>(*this));
        }

    public:
        HighlandIntegral(const ParticleDef& p, Medium const& m,
            std::shared_ptr<Displacement> disp)
            : Highland(p, m)
            , highland_integral(
                  [this, disp](double E) { return Integral(*disp, E); },
                  disp->GetLowerLim(), disp->GetHash()) {};
    };
} // namespace multiple_scattering

inline auto make_highland_integral(ParticleDef const& p, Medium const& m,
    std::shared_ptr<Displacement> disp, bool interpol = false)
{
    auto scatter = std::unique_ptr<multiple_scattering::Parametrization>();
    if (interpol)
        scatter.reset(
            new multiple_scattering::HighlandIntegral<UtilityInterpolant>(
                p, m, disp));
    else
        scatter.reset(
            new multiple_scattering::HighlandIntegral<UtilityIntegral>(
                p, m, disp));
    return scatter;
}

template <typename Cross>
inline auto make_highland_integral(
    ParticleDef const& p, Medium const& m, Cross&& c, bool interpol = false)
{
    auto disp = std::shared_ptr<Displacement>(make_displacement(c, false));
    return make_highland_integral(p, m, disp, interpol);
}
} // namespace PROPOSAL
