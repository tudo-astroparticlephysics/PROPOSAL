
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

#include "PROPOSAL/propagation_utility/Interaction.h"
#include <vector>
#include <memory>
#include <tuple>
#include <functional>

namespace PROPOSAL {
class Component;
class Displacement;
class Time;
class Scattering;
class Decay;
class Cartesian3D;
struct ContRand;
class Vector3D;
enum class InteractionType;
}

namespace PROPOSAL {
class PropagationUtility {
public:
    PropagationUtility() = default;

    virtual Interaction::Loss EnergyStochasticloss(double, double) const = 0;
    virtual double EnergyDecay(double, std::function<double()>, double) const = 0;
    virtual std::pair<double, double> EnergyDistanceStochasticInteraction(double, std::function<double()>) const = 0;
    virtual double EnergyRandomize(double, double, std::function<double()>) const = 0;
    virtual double EnergyDistance(double, double) const = 0;
    virtual double LengthContinuous(double, double) const = 0;
    virtual double TimeElapsed(double, double, double, double) const = 0;

    virtual std::tuple<Cartesian3D, Cartesian3D> DirectionsScatter(
            double, double, double, const Vector3D&, std::function<double()>) const = 0;
    virtual Cartesian3D DirectionDeflect(InteractionType, double, double,
                                 const Vector3D&, std::function<double()>) const = 0;

    virtual double GetLowerPropagationLim() const = 0;
};
}

namespace PROPOSAL {
class PropagationUtilityContinuous : public PropagationUtility {
public:
    struct Collection {

        bool operator==(const Collection& lhs);

        // obligatory pointers
        std::shared_ptr<Interaction> interaction_calc;
        std::shared_ptr<Displacement> displacement_calc;
        std::shared_ptr<Time> time_calc;

        // optional pointers
        std::shared_ptr<Scattering> scattering;
        std::shared_ptr<Decay> decay_calc;
        std::shared_ptr<ContRand> cont_rand;
    };

    PropagationUtilityContinuous(Collection const& collection);

    Interaction::Loss EnergyStochasticloss(double, double) const override;
    double EnergyDecay(double, std::function<double()>, double) const override;
    std::pair<double, double> EnergyDistanceStochasticInteraction(double, std::function<double()>) const override;
    double EnergyRandomize(double, double, std::function<double()>) const override;
    double EnergyDistance(double, double) const override;
    double LengthContinuous(double, double) const override;
    double TimeElapsed(double, double, double, double) const override;

    // TODO: return value doesn't tell what it include. Maybe it would be better
    // to give a tuple of two directions back. One is the mean over the
    // displacement and the other is the actual direction. With a get method
    // there could be a possible access with the position of the object stored
    // in an enum.

    std::tuple<Cartesian3D, Cartesian3D> DirectionsScatter(
        double, double, double, const Vector3D&, std::function<double()>) const override;
    Cartesian3D DirectionDeflect(InteractionType, double, double,
                                 const Vector3D&, std::function<double()>) const override;

    double GetLowerPropagationLim() const override;


    Collection collection;
};
} // namespace PROPOSAL
