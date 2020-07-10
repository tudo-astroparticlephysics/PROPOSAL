
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
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/scattering/ScatteringHighland.h"

using std::shared_ptr;

namespace PROPOSAL {

template <class T, class Cross, class Enable = void>
class ScatteringHighlandIntegral : public ScatteringHighland {
    T highland_integral;

    T BuildHighlandIntegral(Cross&&);
    double CalculateTheta0(double, double, double) override;

public:
    ScatteringHighlandIntegral(
        const ParticleDef&, shared_ptr<const Medium>, Cross&&);

    static Interpolant1DBuilder::Definition interpol_def;

    template <typename Disp> double HighlandIntegral(Disp&&, double);
};

template <class T, class Cross>
class ScatteringHighlandIntegral<T, Cross,
    typename std::enable_if<
        std::is_null_pointer<typename decay<Cross>::type>::value>::type>
    : public Scattering {
public:
    ScatteringHighlandIntegral(
        const ParticleDef&, shared_ptr<const Medium>, Cross&&)
    {
        throw std::invalid_argument("FUUU");
    };

    bool compare(const Scattering&) const {};
    void print(std::ostream&) const {};

    RandomAngles CalculateRandomAngle(double grammage, double ei, double ef,
        const std::array<double, 4>& rnd){};
};

template <class T, class Cross, class Enable>
ScatteringHighlandIntegral<T, Cross, Enable>::ScatteringHighlandIntegral(
    const ParticleDef& p_def, shared_ptr<const Medium> medium, Cross&& cross)
    : ScatteringHighland(p_def, medium)
    , highland_integral(BuildHighlandIntegral(cross))
{
}

template <class T, class Cross, class Enable>
T ScatteringHighlandIntegral<T, Cross, Enable>::BuildHighlandIntegral(
    Cross&& cross)
{
    auto disp = DisplacementBuilder<UtilityIntegral, Cross>(cross);
    auto higland_integral_func = [this, &disp](double energy) {
        return FunctionToIntegral(disp, energy);
    };
    T decay_integral(higland_integral_func, disp.GetLowerLim(cross));
    if (typeid(T) == typeid(UtilityInterpolant)) {
        auto hash = disp.GetHash(cross);
        decay_integral.BuildTables("highland", hash, interpol_def);
    };
    return decay_integral;
}

template <class T, class Cross, class Enable>
template <typename Disp>
double ScatteringHighlandIntegral<T, Cross, Enable>::HighlandIntegral(
    Disp&& disp, double energy)
{
    auto square_momentum = (energy - mass) * (energy + mass);
    auto aux = energy / square_momentum;
    return disp.FunctionToIntegral(energy) * aux * aux;
}

template <class T, class Cross, class Enable>
double ScatteringHighlandIntegral<T, Cross, Enable>::CalculateTheta0(
    double grammage, double ei, double ef)
{
    auto radiation_length = medium_->GetRadiationLength();
    auto aux = 13.6
        * std::sqrt(highland_integral.Calculate(ei, ef) / radiation_length)
        * std::abs(charge);
    aux *= std::max(1 + 0.038 * std::log(grammage / radiation_length), 0.0);
    return std::min(aux, 1.f);
}

/* template <class T> */
/* Interpolant1DBuilder::Definition ScatteringHighlandIntegral<T>::interpol_def(
 */
/*     nullptr, 200, 0., 1e14, 5, false, false, true, 5, false, false, false);
 */

} // namespace PROPOSAL
