
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

#include <functional>

#include "PROPOSAL/crosssection/parametrization/Photonuclear.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/methods.h"

#define Q2_PHOTO_PARAM_INTEGRAL_DEC(param)                                     \
    struct Photo##param : public PhotoQ2Integral {                             \
        Photo##param(shared_ptr<ShadowEffect>);                                \
        using base_param_t = Photonuclear;                                     \
        double FunctionToQ2Integral(const ParticleDef&, const Component&,      \
            double energy, double v, double Q2) const;                         \
    };

namespace PROPOSAL {
namespace crosssection {
class PhotoQ2Integral : public Photonuclear {
public:
    PhotoQ2Integral(shared_ptr<ShadowEffect>);
    virtual ~PhotoQ2Integral() = default;

    virtual double DifferentialCrossSection(
        const ParticleDef&, const Component&, double energy, double v) const;
    virtual double FunctionToQ2Integral(const ParticleDef&, const Component&,
        double energy, double v, double Q2) const = 0;

    std::shared_ptr<ShadowEffect> shadow_effect_;
};

Q2_PHOTO_PARAM_INTEGRAL_DEC(AbramowiczLevinLevyMaor91)
Q2_PHOTO_PARAM_INTEGRAL_DEC(AbramowiczLevinLevyMaor97)
Q2_PHOTO_PARAM_INTEGRAL_DEC(ButkevichMikhailov)
Q2_PHOTO_PARAM_INTEGRAL_DEC(RenoSarcevicSu)
} // namespace crosssection
} // namespace PROPOSAL

#undef Q2_PHOTO_PARAM_INTEGRAL_DEC
