
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
#include "PROPOSAL/propagation_utility/Time.h"
#include "PROPOSAL/propagation_utility/Decay.h"
#include "PROPOSAL/propagation_utility/ContRand.h"
#include "PROPOSAL/propagation_utility/Displacement.h"
#include "PROPOSAL/propagation_utility/Interaction.h"

using std::tuple;
namespace PROPOSAL {

class Utility {
public:
    struct Definition {
        CrossSectionList cross;              // copy is expensive
        std::shared_ptr<Scattering> scattering;      // no scattering
        std::shared_ptr<InterpolationDef> inter_def; // integration used

        Definition(CrossSectionList, const ParticleDef&,
            std::shared_ptr<Scattering>, std::shared_ptr<InterpolationDef>);
        ~Definition();

        std::unique_ptr<Displacement> displacement_calc = nullptr;
        std::unique_ptr<Interaction> interaction_calc = nullptr;
        std::unique_ptr<Decay> decay_calc = nullptr;
        std::unique_ptr<ContRand> cont_rand = nullptr;
        std::unique_ptr<Time> time_calc = nullptr;
    };

    double EnergyStochasticloss(
         CrossSection&, double, const std::array<double, 2>&);
    double EnergyDecay(double, double);
    double EnergyInteraction(double, double);
    double EnergyRandomize(double, double, double);
    double LengthContinuous(double, double);
    double TimeElapsed(double, double);

    /* // TODO: return value doesn't tell what it include. Maybe it would be better */
    /* // to give a tuple of two directions back. One is the mean over the */
    /* // displacement and the other is the actual direction. With a get method */
    /* // there could be a possible access with the position of the object stored */
    /* // in an enum. */
    tuple<Vector3D, Vector3D> DirectionsScatter(double, double, double,
        const Vector3D&, const Vector3D&, const std::array<double, 4>&);
    std::pair<double, double> DirectionDeflect(CrossSection& , double particle_energy, double loss_energy);

private:
    std::unique_ptr<Definition> utility_def;
};
} // namespace PROPOSAL
