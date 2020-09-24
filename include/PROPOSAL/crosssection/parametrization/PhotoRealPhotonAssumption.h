
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

#include "PROPOSAL/crosssection/parametrization/Photonuclear.h"
#include <memory>
#include <unordered_map>
using std::unordered_map;

#define PHOTO_PARAM_REAL_DEC(param, parent)                                    \
    class Photo##param : public Photo##parent {                                \
    public:                                                                    \
        Photo##param(bool hard_component);                                     \
        using base_param_t = Photonuclear;                                     \
                                                                               \
        virtual double CalculateParametrization(                               \
            const Component&, double nu) const;                                \
    };

namespace PROPOSAL {
namespace crosssection {
class PhotoRealPhotonAssumption : public Photonuclear {
protected:
    bool hard_component_;
    unordered_map<size_t, std::shared_ptr<RealPhoton>> hard_component_map;

public:
    PhotoRealPhotonAssumption(bool hard_component);
    virtual ~PhotoRealPhotonAssumption() = default;

    virtual double DifferentialCrossSection(
        const ParticleDef&, const Component&, double energy, double v) const;
    virtual double CalculateParametrization(
        const Component&, double nu) const = 0;
    double NucleusCrossSectionCaldwell(double nu) const;
};

PHOTO_PARAM_REAL_DEC(Zeus, RealPhotonAssumption)
PHOTO_PARAM_REAL_DEC(BezrukovBugaev, RealPhotonAssumption)
PHOTO_PARAM_REAL_DEC(Kokoulin, BezrukovBugaev)

class PhotoRhode : public PhotoRealPhotonAssumption {
    std::shared_ptr<Interpolant> interpolant_;

    double MeasuredSgN(double e) const;

public:
    PhotoRhode(bool hard_component);
    using base_param_t = Photonuclear;
    double CalculateParametrization(const Component&, double nu) const override;
};

#undef PHOTO_PARAM_REAL_DEC

// Factory pattern functions

template <typename P, typename M>
using photoreal_func_ptr = cross_t_ptr<P, M>(*)(P, M, std::shared_ptr<const
        EnergyCutSettings>, bool, bool);

template <typename Param, typename P, typename M>
cross_t_ptr<P, M> create_photoreal(P p_def, M medium,std::shared_ptr<const
        EnergyCutSettings> cuts, bool hard_component, bool interpol) {
    auto param = Param(hard_component);
    return make_crosssection(param, p_def, medium, cuts, interpol);
}

template<typename P, typename M>
static std::map<std::string, photoreal_func_ptr<P, M>> photoreal_map = {
        {"Zeus", create_photoreal<PhotoZeus, P, M>},
        {"BezrukovBugaev", create_photoreal<PhotoBezrukovBugaev, P, M>},
        {"Kokoulin", create_photoreal<PhotoKokoulin, P, M>},
        {"Rhode", create_photoreal<PhotoRhode, P, M>},
};

template<typename P, typename M>
cross_t_ptr<P, M> make_photonuclearreal(P p_def, M medium, std::shared_ptr<const
        EnergyCutSettings> cuts, bool interpol, const nlohmann::json& config){
    if (!config.contains("parametrization"))
        throw std::logic_error("No parametrization passed for photonuclear");
    std::string param_name = config["parametrization"];

    auto it = photoreal_map<P, M>.find(param_name);
    if (it != photoreal_map<P, M>.end()) {
        bool hard_component = config.value("hard_component", true);
        return it->second(p_def, medium, cuts, hard_component, interpol);
    }
    throw std::invalid_argument("Unknown parametrization for photonuclear");
}

} // namespace crosssection
} // namespace PROPOSAL