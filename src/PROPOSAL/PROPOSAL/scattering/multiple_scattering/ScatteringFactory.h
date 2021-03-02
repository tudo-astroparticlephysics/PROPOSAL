
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

#include <algorithm>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include "PROPOSAL/scattering/multiple_scattering/Highland.h"
#include "PROPOSAL/scattering/multiple_scattering/HighlandIntegral.h"
#include "PROPOSAL/scattering/multiple_scattering/Moliere.h"

namespace PROPOSAL {

enum class ScatteringType : int {
    NoScattering,
    Moliere,
    Highland,
    HighlandIntegral
};

static const std::unordered_map<std::string, ScatteringType>
    MultipleScatteringTable = { { "moliere", ScatteringType::Moliere },
        { "highland", ScatteringType::Highland },
        { "highlandintegral", ScatteringType::HighlandIntegral },
        { "noscattering", ScatteringType::NoScattering } };

namespace detail {

    constexpr static char multiple_scattering_msg[]
        = "This constructor is not provided.";

    inline std::unique_ptr<multiple_scattering::Parametrization> make_multiple_scattering(ScatteringType t)
    {
        switch (t) {
            case ScatteringType::NoScattering:
                return std::unique_ptr<multiple_scattering::Parametrization>(
                        nullptr);
            default:
                throw std::out_of_range(multiple_scattering_msg);
        }
    }

    inline std::unique_ptr<multiple_scattering::Parametrization> make_multiple_scattering(
            ScatteringType t, ParticleDef const& p, Medium const& m)
    {
        switch (t) {
            case ScatteringType::Highland:
                return make_highland(p, m);
            case ScatteringType::Moliere:
                return make_moliere(p, m);
            default:
                return make_multiple_scattering(t);
        }
    }

    template <typename... Args>
    inline std::unique_ptr<multiple_scattering::Parametrization> make_multiple_scattering(
        ScatteringType t, ParticleDef const& p, Medium const& m, Args&&... args)
    {
        switch (t) {
        case ScatteringType::HighlandIntegral:
            return make_highland_integral(p, m, std::forward<Args>(args)...);
        default:
            return make_multiple_scattering(t, p, m);
        }
    }
}

template <typename... Args>
auto make_multiple_scattering(std::string const& name, Args... args)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    auto it = MultipleScatteringTable.find(name_lower);
    if (it != MultipleScatteringTable.end()) {
        return detail::make_multiple_scattering(it->second, args...);
    }
    throw std::out_of_range("This scattering model is not provided.");
}

} // namespace PROPOSAL