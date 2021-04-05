
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

#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/scattering/stochastic_deflection/Parametrization.h"
#include "PROPOSAL/scattering/stochastic_deflection/bremsstrahlung/TsaiApproximationBremsstrahlung.h"
#include "PROPOSAL/scattering/stochastic_deflection/ionization/NaivIonization.h"
#include "PROPOSAL/scattering/stochastic_deflection/nuclearInteraction/BorogPetrukhinNuclearInteraction.h"
#include "PROPOSAL/scattering/stochastic_deflection/pairProd/KelnerPairProduction.h"

namespace PROPOSAL {

/**
 * @brief Creates a default deflection of specific type.
 *
 * @tparam Args InteractionType and stochastic_deflection::Parametrization
 * cstr. arguments
 * @param args interaction type and the arguments of the corresponding cstr.
 */
template <typename... Args>
inline auto make_default_stochastic_deflection(Args... args)
{
    return DefaultFactory<stochastic_deflection::Parametrization>::Create(
        args...);
}

/**
 * @brief Creates a vector of default deflection for a vector of types.
 *
 * @tparam Args std::vector<InteractionType> and
 * stochastic_deflection::Parametrization cstr. arguments
 * @param types list of interaction types
 * @param args the corresponding cstr. argmuents
 */
template <typename... Args>
inline auto make_default_stochastic_deflection(
    std::vector<InteractionType> const& types, Args... args)
{
    using param_ptr = std::unique_ptr<stochastic_deflection::Parametrization>;
    auto v = std::vector<param_ptr>();
    for (auto t : types)
        v.emplace_back(make_default_stochastic_deflection(t, args...));
    return v;
}

using func_ptr = std::unique_ptr<stochastic_deflection::Parametrization> (*)(
    const ParticleDef&, const Medium&);

template <typename Param>
std::unique_ptr<stochastic_deflection::Parametrization> create_deflection(
    const ParticleDef& p_def, const Medium& medium)
{
    return make_unique<Param>(p_def, medium);
}

static const std::map<std::string, func_ptr> DeflectionTable
    = { { "tsaiapproximationbremsstrahlung",
            create_deflection<stochastic_deflection::TsaiApproximationBremsstrahlung> },
        { "naivionization",
              create_deflection<stochastic_deflection::NaivIonization> },
        {"borogpetrukhinnuclearinteraction",
                create_deflection<stochastic_deflection::BorogPetrukhinNuclearInteraction>},
        {"kelnerpairproduction",
                create_deflection<stochastic_deflection::KelnerPairProduction>}

      };

template <typename Cross = std::nullptr_t>
std::unique_ptr<stochastic_deflection::Parametrization>
make_stochastic_deflection(
    std::string const& name, ParticleDef const& p_def, Medium const& medium)
{
    std::string name_lower = name;
    std::transform(name.begin(), name.end(), name_lower.begin(), ::tolower);

    auto it = DeflectionTable.find(name_lower);
    if (it != DeflectionTable.end()) {
        return it->second(p_def, medium);
    }
    throw std::out_of_range(
        "This stochastic deflection model is not provided.");
}

} // namespace PROPOSAL
