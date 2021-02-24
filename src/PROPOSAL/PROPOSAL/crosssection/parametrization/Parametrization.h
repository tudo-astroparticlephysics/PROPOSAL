
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

#include <memory>
#include <type_traits>

namespace PROPOSAL {
struct ParticleDef;
class Medium;
class Component;
enum class InteractionType;
} // namespace PROPOSAL

namespace PROPOSAL {
namespace crosssection {

    struct KinematicLimits {
        double v_min;
        double v_max;
    };

    template <typename T> struct ParametrizationName {
        static constexpr auto value = "name_not_available";
    };

    template <typename T> struct ParametrizationId {
        static constexpr size_t value = 0;
    };


    template <typename T> struct is_component_wise : std::true_type {
    };

    template <typename T> struct is_only_stochastic : std::false_type {
    };

    template <typename Target> class Parametrization {
    protected:
        size_t hash;

    public:
        Parametrization()
            : hash(0)
        {
        }

        virtual ~Parametrization() = default;

        virtual std::unique_ptr<Parametrization<Target>> clone() const = 0;

        virtual double DifferentialCrossSection(
            ParticleDef const&, Target const&, double, double) const = 0;

        virtual KinematicLimits GetKinematicLimits(
            ParticleDef const&, Target const&, double) const = 0;

        virtual double FunctionToDEdxIntegral(
            ParticleDef const& p, const Target& t, double E, double v) const
        {
            return v * DifferentialCrossSection(p, t, E, v);
        }

        virtual double FunctionToDE2dxIntegral(
            ParticleDef const& p, Target const& t, double E, double v) const
        {
            return v * v * DifferentialCrossSection(p, t, E, v);
        }

        virtual double GetLowerEnergyLim(ParticleDef const&) const noexcept = 0;

        inline size_t GetHash() const noexcept { return hash; };
    };

    template class Parametrization<Medium>;
    template class Parametrization<Component>;

} // namespace crosssection
} // namespace PROPOSAL
