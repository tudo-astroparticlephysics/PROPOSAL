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
#include <memory>
#include "PROPOSAL/crosssection/CrossSection.h"

#define BREMSSTRAHLUNG_DEF(param)                                              \
    struct Brems##param : public Bremsstrahlung {                              \
        Brems##param(bool);                                                    \
        using base_param_t = Bremsstrahlung;                                   \
                                                                               \
        double CalculateParametrization(const ParticleDef&, const Component&,  \
            double energy, double v) const override;                                 \
    };

using std::unique_ptr;

namespace PROPOSAL {
class Interpolant;
} // namespace PROPOSAL

namespace PROPOSAL {
namespace crosssection {
class Bremsstrahlung : public Parametrization {

protected:
    bool lorenz_;       // enable lorenz cut
    double lorenz_cut_; // in [MeV]
    bool init_lpm_effect_;
    bool lpm_;
    double eLpm_;

    double lpm(const ParticleDef&, const Component&, double, double, double);

public:
    Bremsstrahlung(bool);
    virtual ~Bremsstrahlung() = default;

    using only_stochastic = std::false_type;
    using component_wise = std::true_type;
    using base_param_t = Bremsstrahlung;

    double DifferentialCrossSection(
        const ParticleDef&, const Component&, double, double) const override;
    virtual double CalculateParametrization(
        const ParticleDef&, const Component&, double, double) const
        = 0;

    double GetLowerEnergyLim(const ParticleDef&) const noexcept override;
    tuple<double, double> GetKinematicLimits(
        const ParticleDef&, const Component&, double) const noexcept override;
    size_t GetHash() const noexcept override;
};

BREMSSTRAHLUNG_DEF(PetrukhinShestakov)
BREMSSTRAHLUNG_DEF(KelnerKokoulinPetrukhin)
BREMSSTRAHLUNG_DEF(CompleteScreening)
BREMSSTRAHLUNG_DEF(AndreevBezrukovBugaev)
BREMSSTRAHLUNG_DEF(SandrockSoedingreksoRhode)

class BremsElectronScreening : public Bremsstrahlung {
    std::shared_ptr<Interpolant> interpolant_;

public:
    BremsElectronScreening(bool);
    using base_param_t = Bremsstrahlung;

    double CalculateParametrization(
        const ParticleDef&, const Component&, double energy, double v) const override;
    double DifferentialCrossSection(
        const ParticleDef&, const Component&, double energy, double v) const override;
};

// Factory pattern functions

template <typename P, typename M>
using brems_func_ptr = cross_t_ptr<P, M>(*)(P, M, std::shared_ptr<const
        EnergyCutSettings>, bool, bool);

template <typename Param, typename P, typename M>
cross_t_ptr<P, M> create_brems(P p_def, M medium,std::shared_ptr<const
        EnergyCutSettings> cuts, bool lpm, bool interpol) {
    auto param = Param(lpm);
    return make_crosssection(param, p_def, medium, cuts, interpol);
}

template<typename P, typename M>
static std::map<std::string, brems_func_ptr<P, M>> brems_map = {
        {"KelnerKokoulinPetrukhin", create_brems<BremsKelnerKokoulinPetrukhin, P, M>},
        {"PetrukhinShestakov", create_brems<BremsPetrukhinShestakov, P, M>},
        {"CompleteScreening", create_brems<BremsCompleteScreening, P, M>},
        {"AndreevBezrukovBugaev", create_brems<BremsAndreevBezrukovBugaev, P, M>},
        {"SandrockSoedingreksoRhode", create_brems<BremsSandrockSoedingreksoRhode, P, M>},
        {"ElectronScreening", create_brems<BremsElectronScreening, P, M>}
};

template<typename P, typename M>
cross_t_ptr<P, M> make_bremsstrahlung(P p_def, M medium, std::shared_ptr<const
        EnergyCutSettings> cuts, bool interpol, const nlohmann::json& config){
    if (!config.contains("parametrization"))
        throw std::logic_error("No parametrization passed for bremsstrahlung");

    std::string param_name = config["parametrization"];
    auto it = brems_map<P, M>.find(param_name);
    if (it == brems_map<P, M>.end())
        throw std::logic_error("Unknown parametrization for bremsstrahlung");

    bool lpm = config.value("lpm", true);
    return it->second(p_def, medium, cuts, lpm, interpol);
}

} // namespace crosssection
} // namespace PROPOSAL

#undef BREMSSTRAHLUNG_DEF
