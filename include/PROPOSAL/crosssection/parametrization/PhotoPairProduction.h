
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

#include "PROPOSAL/crosssection/parametrization/Parametrization.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/crosssection/CrossSection.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/crosssection/CrossSectionMultiplier.h"

namespace PROPOSAL {
namespace crosssection {
class PhotoPairProduction : public Parametrization {
public:
    PhotoPairProduction();
    virtual ~PhotoPairProduction() = default;

    using only_stochastic = std::true_type;
    using component_wise = std::true_type;

    double GetLowerEnergyLim(const ParticleDef&) const noexcept override;
    std::tuple<double, double> GetKinematicLimits(
        const ParticleDef&, const Component&, double) const noexcept override;
};

struct PhotoPairTsai : public PhotoPairProduction {
    PhotoPairTsai() = default;
    using base_param_t = PhotoPairProduction;

    virtual double DifferentialCrossSection(
        const ParticleDef&, const Component&, double, double) const;
};

/* class PhotoAngleDistribution { */
/* public: */
/*     PhotoAngleDistribution() = default; */
/*     virtual ~PhotoAngleDistribution() = default; */

/*     struct DeflectionAngles { */
/*         double cosphi0, theta0, cosphi1, theta1; */
/*     }; */

/*     virtual DeflectionAngles SampleAngles( */
/*         const Component& comp, double energy, double rho) */
/*         = 0; */
/* }; */

/* class PhotoAngleTsaiIntegral : public PhotoAngleDistribution { */
/*     Integral integral_; */
/*     double FunctionToIntegral(const Component&, double, double, double); */

/* public: */
/*     PhotoAngleTsaiIntegral(); */

/*     DeflectionAngles SampleAngles(const Component&, double, double) override; */
/* }; */

/* class PhotoAngleNoDeflection : public PhotoAngleDistribution { */
/* public: */
/*     PhotoAngleNoDeflection(); */

/*     DeflectionAngles SampleAngles(const Component&, double, double) override; */
/* }; */

/* class PhotoAngleEGS : public PhotoAngleDistribution { */
/* public: */
/*     PhotoAngleEGS(); */

/*     DeflectionAngles SampleAngles(const Component&, double, double) override; */
/* }; */

// Factory pattern functions

template <typename P, typename M>
using photopair_func_ptr = cross_t_ptr<P, M>(*)(P, M, bool, double);

template <typename Param, typename P, typename M>
cross_t_ptr<P, M> create_photopair(P p_def, M medium, bool interpol,
                                   double multiplier = 1.0) {
    auto param = Param();
    auto cross =  make_crosssection(param, p_def, medium, nullptr, interpol);
    if (multiplier == 1.0)
        return cross;
    return make_crosssection_multiplier(std::shared_ptr<crosssection_t<P, M>>(
            std::move(cross)), multiplier);
}

template<typename P, typename M>
static std::map<std::string, photopair_func_ptr<P, M>> photopair_map = {
        {"tsai", create_photopair<PhotoPairTsai, P, M>}
};

template<typename P, typename M>
cross_t_ptr<P, M> make_photopairproduction(P p_def, M medium, bool interpol,
                                           const std::string& param_name,
                                           double multiplier = 1.0){
    std::string name = param_name;
    std::transform(param_name.begin(), param_name.end(), name.begin(), ::tolower);
    auto it = photopair_map<P, M>.find(name);
    if (it == photopair_map<P, M>.end())
        throw std::logic_error("Unknown parametrization for photopairproduction");

    return it->second(p_def, medium, interpol, multiplier);
}

template<typename P, typename M>
cross_t_ptr<P, M> make_photopairproduction(P p_def, M medium, bool interpol,
                                    const nlohmann::json& config){
    if (!config.contains("parametrization"))
        throw std::logic_error("No parametrization passed for photopairproduction");
    std::string param_name = config["parametrization"];
    double multiplier = config.value("multiplier", 1.0);
    return make_photopairproduction(p_def, medium, interpol, param_name,
                                    multiplier);
}

} //namespace crosssection
} // namespace PROPOSAL
