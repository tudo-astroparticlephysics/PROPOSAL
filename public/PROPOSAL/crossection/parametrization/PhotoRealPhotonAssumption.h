
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

#include "PROPOSAL/crossection/parametrization/Photonuclear.h"
#include <memory>
#include <unordered_map>
using std::unordered_map;

#define PHOTO_PARAM_REAL_DEC(param, parent)                                    \
    class Photo##param : public Photo##parent {                                \
    public:                                                                    \
        Photo##param(bool hard_component);                                     \
                                                                               \
        virtual double CalculateParametrization(const Component&, double nu);  \
    };

namespace PROPOSAL {

class PhotoRealPhotonAssumption : public Photonuclear {
protected:
    bool hard_component_;
    unordered_map<size_t, std::unique_ptr<RealPhoton>> hard_component_map;

public:
    PhotoRealPhotonAssumption(bool hard_component);
    virtual ~PhotoRealPhotonAssumption() = default;

    virtual double DifferentialCrossSection(
        const ParticleDef&, const Component&, double energy, double v);
    virtual double CalculateParametrization(const Component&, double nu) = 0;
    double NucleusCrossSectionCaldwell(double nu);
};

PHOTO_PARAM_REAL_DEC(Zeus, RealPhotonAssumption)
PHOTO_PARAM_REAL_DEC(BezrukovBugaev, RealPhotonAssumption)
PHOTO_PARAM_REAL_DEC(Kokoulin, BezrukovBugaev)

class PhotoRhode : public PhotoRealPhotonAssumption {
    std::unique_ptr<Interpolant> interpolant_;

    double MeasuredSgN(double e);

public:
    PhotoRhode(bool hard_component);

    double CalculateParametrization(const Component&, double nu) override;
};

#undef Q2_PHOTO_PARAM_INTEGRAL_DEC

} // namespace PROPOSAL
