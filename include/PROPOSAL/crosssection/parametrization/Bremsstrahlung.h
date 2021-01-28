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

#define BREMSSTRAHLUNG_DEF(param)                                              \
    struct Brems##param : public Bremsstrahlung {                              \
        Brems##param(bool lpm = false);                                        \
        Brems##param(bool lpm, const ParticleDef&, const Medium&,              \
            double density_correction = 1.0);                                  \
                                                                               \
        double CalculateParametrization(const ParticleDef&, const Component&,  \
            double energy, double v) const final;                              \
    };

namespace PROPOSAL {
class ParticleDef;
class Medium;
class Component;
class Interpolant;
} // namespace PROPOSAL

namespace PROPOSAL {
namespace crosssection {
    class BremsLPM;

    class Bremsstrahlung : public Parametrization<Component> {

    protected:
        bool lorenz_;       // enable lorenz cut
        double lorenz_cut_; // in [MeV]
        std::shared_ptr<BremsLPM> lpm_;

    public:
        Bremsstrahlung() = default;

        virtual ~Bremsstrahlung() = default;

        double DifferentialCrossSection(const ParticleDef&, const Component&,
            double, double) const override;
        virtual double CalculateParametrization(
            const ParticleDef&, const Component&, double, double) const = 0;

        double GetLowerEnergyLim(const ParticleDef&) const final;
        KinematicLimits GetKinematicLimits(
            const ParticleDef&, const Component&, double) const final;
    };

    BREMSSTRAHLUNG_DEF(PetrukhinShestakov)
    BREMSSTRAHLUNG_DEF(KelnerKokoulinPetrukhin)
    BREMSSTRAHLUNG_DEF(CompleteScreening)
    BREMSSTRAHLUNG_DEF(AndreevBezrukovBugaev)
    BREMSSTRAHLUNG_DEF(SandrockSoedingreksoRhode)

    class BremsElectronScreening : public Bremsstrahlung {
        std::shared_ptr<Interpolant> interpolant_;

    public:
        BremsElectronScreening(bool lpm = false);
        BremsElectronScreening(bool lpm, const ParticleDef&, const Medium&,
            double density_correction = 1.0);

        double CalculateParametrization(const ParticleDef&, const Component&,
            double energy, double v) const final;
        double DifferentialCrossSection(const ParticleDef&, const Component&,
            double energy, double v) const final;
    };

    // LPM effect object
    class BremsLPM {
        size_t hash;
        double mass_;
        double mol_density_;
        double mass_density_;
        double sum_charge_;
        double density_correction_;
        double eLpm_;

    public:
        BremsLPM(const ParticleDef&, const Medium&, const Bremsstrahlung&,
            double density_correction = 1.0);
        double suppression_factor(
            double energy, double v, const Component&) const;
        size_t GetHash() const noexcept { return hash; }
    };

    /* // Factory pattern functions */

    /* template <typename P, typename M> */
    /* using brems_func_ptr = cross_t_ptr<P, M>(*)(P, M, std::shared_ptr<const
     */
    /*         EnergyCutSettings>, bool, bool, double); */

    /* template <typename Param, typename P, typename M> */
    /* cross_t_ptr<P, M> create_brems(P p_def, M medium, */
    /*                                std::shared_ptr<const EnergyCutSettings>
     * cuts, */
    /*                                bool lpm, bool interpol, */
    /*                                double density_correction = 1.0) { */
    /*         auto param = Param(lpm, p_def, medium, density_correction); */
    /*         return make_crosssection(param, p_def, medium, cuts, interpol);
     */
    /* } */

    /* template<typename P, typename M> */
    /* static std::map<std::string, brems_func_ptr<P, M>> brems_map = { */
    /*         {"kelnerkokoulinpetrukhin",
     * create_brems<BremsKelnerKokoulinPetrukhin, P, M>}, */
    /*         {"petrukhinshestakov", create_brems<BremsPetrukhinShestakov, P,
     * M>}, */
    /*         {"completescreening", create_brems<BremsCompleteScreening, P,
     * M>}, */
    /*         {"andreevbezrukovbugaev",
     * create_brems<BremsAndreevBezrukovBugaev, P, M>}, */
    /*         {"sandrocksoedingreksorhode",
     * create_brems<BremsSandrockSoedingreksoRhode, P, M>}, */
    /*         {"electronscreening", create_brems<BremsElectronScreening, P, M>}
     */
    /* }; */

    /* template<typename P, typename M> */
    /* cross_t_ptr<P, M> make_bremsstrahlung(P p_def, M medium,
     * std::shared_ptr<const */
    /*         EnergyCutSettings> cuts, bool interpol, bool lpm, std::string
     * param_name, */
    /*         double density_correction = 1.0){ */
    /*     std::string name = param_name; */
    /*     std::transform(param_name.begin(), param_name.end(), name.begin(),
     * ::tolower); */
    /*     auto it = brems_map<P, M>.find(name); */
    /*     if (it == brems_map<P, M>.end()) */
    /*         throw std::logic_error("Unknown parametrization for
     * bremsstrahlung"); */

    /*     return it->second(p_def, medium, cuts, lpm, interpol,
     * density_correction); */
    /* } */

    /* template<typename P, typename M> */
    /* cross_t_ptr<P, M> make_bremsstrahlung(P p_def, M medium,
     * std::shared_ptr<const */
    /*         EnergyCutSettings> cuts, bool interpol, const nlohmann::json&
     * config, */
    /*         double density_correction = 1.0){ */
    /*     if (!config.contains("parametrization")) */
    /*         throw std::logic_error("No parametrization passed for
     * bremsstrahlung"); */

    /*     std::string param_name = config["parametrization"]; */
    /*     bool lpm = config.value("lpm", true); */

    /*     return make_bremsstrahlung(p_def, medium, cuts, interpol, lpm,
     * param_name, */
    /*                                density_correction); */
    /* } */

} // namespace crosssection
} // namespace PROPOSAL

#undef BREMSSTRAHLUNG_DEF
