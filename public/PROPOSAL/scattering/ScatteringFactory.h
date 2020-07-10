
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
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include "PROPOSAL/scattering/Scattering.h"
#include "PROPOSAL/scattering/ScatteringHighland.h"
#include "PROPOSAL/scattering/ScatteringHighlandIntegral.h"
#include "PROPOSAL/scattering/ScatteringMoliere.h"

using std::make_shared;
namespace PROPOSAL {

enum class ScatteringType : int { Moliere, Highland, HighlandIntegral };

static const std::unordered_map<std::string, ScatteringType> ScatteringTable
    = { { "moliere", ScatteringType::Moliere },
          { "highland", ScatteringType::Highland },
          { "highland_integral", ScatteringType::HighlandIntegral } };

template <typename Cross = std::nullptr_t>
unique_ptr<Scattering> make_scattering(std::string const& name,
    ParticleDef const& p_def, Medium const& medium, Cross&& cross = nullptr)
{
    auto m = make_shared<const Medium>(medium);
    auto it = ScatteringTable.find(name);
    if (it != ScatteringTable.end()) {
        switch (it->second) {
        case ScatteringType::HighlandIntegral:
            return unique_ptr<Scattering>(
                new ScatteringHighlandIntegral<UtilityInterpolant, Cross>(
                    p_def, m, std::forward<Cross>(cross)));
        case ScatteringType::Highland:
            return unique_ptr<Scattering>(new ScatteringHighland(p_def, m));
        case ScatteringType::Moliere:
            return unique_ptr<Scattering>(new ScatteringMoliere(p_def, m));
        default:
            throw std::out_of_range("This constructor is not provided.");
        }
    }
    throw std::out_of_range("This scattering model is not provided.");
}

} // namespace PROPOSAL
