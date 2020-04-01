
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

#include "PROPOSAL/crossection/CrossSection.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/propagation_utility/Displacement.h"
#include "PROPOSAL/scattering/ScatteringHighland.h"

namespace PROPOSAL {

template <class T>
class ScatteringHighlandIntegral : public ScatteringHighland {
public:
    ScatteringHighlandIntegral<T>(const ParticleDef& p_def,
        std::shared_ptr<const Medium> medium, CrossSectionList cross)
        : ScatteringHighland(p_def, medium)
        , displacement(new DisplacementBuilder<UtilityIntegral>(cross))
        , integral(std::bind(&ScatteringHighlandIntegral::HighlandIntegral,
              this, std::placeholders::_1))
    {
        if (typeid(T) == typeid(UtilityInterpolant)) {
            size_t hash_digest = 0;
            for (const auto& crosssection : cross)
                hash_combine(hash_digest, crosssection->GetHash());
            interpol_def.function1d = [this](double energy) {
                return reinterpret_cast<UtilityIntegral*>(&integral)->Calculate(interpol_def.xmax, energy, 0);
            };

            integral.BuildTables(name, hash_digest, interpol_def);
        }
    }

    double HighlandIntegral(double energy)
    {
        double square_momentum = (energy - mass) * (energy + mass);
        double aux = energy / square_momentum;

        return displacement->FunctionToIntegral(energy) * aux * aux;
    }

    static Interpolant1DBuilder::Definition interpol_def;

private:
    ScatteringHighlandIntegral& operator=(const ScatteringHighlandIntegral&)
        = delete;

    bool compare(const Scattering& scattering) const override
    {
        // TODO: Add ScatteringHighlandIntegral comparison operator
        throw std::logic_error(
            "This comparison function has not been implemented yet. "
            "The developers need to put on the thinking caps to fix it.");
    }

    void print(std::ostream&) const override {}

    double CalculateTheta0(
        double dr, double ei, double ef, const Vector3D& pos) override
    {
        double aux = integral.Calculate(ei, ef, 0.0)
            * medium_->GetDensityDistribution().Evaluate(pos);
        double cutoff = 1;
        double radiation_length = medium_->GetRadiationLength(pos);

        aux = 13.6 * std::sqrt(std::max(aux, 0.0) / radiation_length)
            * std::abs(charge);
        aux *= std::max(1 + 0.038 * std::log(dr / radiation_length), 0.0);

        return std::min(aux, cutoff);
    }

    std::string name = "scattering";
    T integral;
    std::unique_ptr<Displacement> displacement;
};

Interpolant1DBuilder::Definition interpol_def(nullptr, 200, 0., 1e14, 5, false, false, true, 5, false, false, false);

} // namespace PROPOSAL
