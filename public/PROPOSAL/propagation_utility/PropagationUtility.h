
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

#include <vector>

#include "PROPOSAL/EnergyCutSettings.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/scattering/Scattering.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

using std::tuple;
namespace PROPOSAL {

class CrossSection;
class UtilityDecorator;
struct InterpolationDef;

typedef std::vector<std::shared_ptr<CrossSection>> CrossSectionList;
//typedef const double random_number;
//typedef const double energy;
//typedef const double distance;

class Utility {
public:
    struct Definition {
        std::shared_ptr<InterpolationDef> inter_def;  // integration used
        std::shared_ptr<Scattering> scattering;       // no scattering
        CrossSectionList crosssections; // copy is expensive

        Definition(CrossSectionList,
            std::shared_ptr<Scattering>,
            std::shared_ptr<InterpolationDef>);
        ~Definition();
    };

    Utility(std::unique_ptr<Definition> utility_def)
        : utility_def(std::move(utility_def)){};

    std::shared_ptr<CrossSection> TypeInteraction(
        double, const std::array<double, 2>&);
    double EnergyStochasticloss(
        const CrossSection&, double, const std::array<double, 2>&);
    double EnergyDecay(double, double);
    double EnergyInteraction(double, double);
    double EnergyRandomize(double, double, double);
    double LengthContinuous(double, double, double);
    double ElapsedTime(double, double, double);

    // TODO: return value doesn't tell what it include. Maybe it would be better
    // to give a tuple of two directions back. One is the mean over the
    // displacement and the other is the actual direction. With a get method
    // there could be a possible access with the position of the object stored
    // in an enum.
    tuple<Vector3D, Vector3D> DirectionsScatter(double, double, double, const Vector3D&, const Vector3D&, const std::array<double, 4>&);
    Vector3D DirectionDeflect(double, double, double, const Vector3D&, const Vector3D&);

private:
    Definition utility_def;

    std::unique_ptr<UtilityDecorator> displacement_calc = nullptr;
    std::unique_ptr<UtilityDecorator> interaction_calc = nullptr;
    std::unique_ptr<UtilityDecorator> decay_calc = nullptr;
    std::unique_ptr<UtilityDecorator> cont_rand = nullptr;
    std::unique_ptr<UtilityDecorator> exact_time = nullptr;

    double mass;
};
} // namespace PROPOSAL

namespace PROPOSAL {
class UtilityDecorator {

public:
    UtilityDecorator(CrossSectionList cross)
        : crosssections(cross){};

    virtual double FunctionToIntegral(double energy) = 0;
    virtual double Calculate(double ei, double ef, double rnd) = 0;
    virtual double GetUpperLimit(double ei, double rnd) = 0;

protected:
    CrossSectionList crosssections;
};
} // namespace PROPOSAL
