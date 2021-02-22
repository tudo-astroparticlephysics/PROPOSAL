
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

#include "CubicInterpolation/CubicSplines.h"
#include "CubicInterpolation/Interpolant.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"

#include <functional>
#include <memory>

namespace PROPOSAL {
class Integral;
class Interpolant;
struct InterpolationDef;

class UtilityInterpolant : public UtilityIntegral {
    using interpolant_t
        = cubic_splines::Interpolant<cubic_splines::CubicSplines<double>>;
    using interpolant_ptr = std::shared_ptr<interpolant_t>;

    double lower_lim;
    interpolant_ptr interpolant_;

    const double bias = 1.; // bias to be able to use log substitution for f. should cancel out in every calculation.

    // maybe interpolate function to integral will give a performance boost.
    // in general this function should be underfrequently called
    // Interpolant1DBuilder builder_diff;
    // std::unique_ptr<Interpolant> interpolant_diff_;

    std::string gen_path() const;
    std::string gen_name(std::string prefix) const;

public:
    UtilityInterpolant(std::function<double(double)>, double, size_t);
    virtual ~UtilityInterpolant() = default;

    void BuildTables(const std::string prefix, bool reverse = false)
    {
        auto def = cubic_splines::CubicSplines<double>::Definition();
        auto reference_x = lower_lim;
        if (reverse)
            reference_x = 1.e14;

        hash_combine(this->hash, reverse);

        def.f = [&](double energy) {
            return UtilityIntegral::Calculate(energy, reference_x) + bias;
        };

        def.f_trafo = std::make_unique<cubic_splines::ExpAxis<double>>(1., 0.);
        def.axis = std::make_unique<cubic_splines::ExpAxis<double>>(
            lower_lim, 1e14, (size_t)100);

        interpolant_ = std::make_shared<interpolant_t>(
            std::move(def), gen_path(), gen_name(prefix));
    }

    double Calculate(double, double) final;
    double GetUpperLimit(double, double) final;
};
} // namespace PROPOSAL
