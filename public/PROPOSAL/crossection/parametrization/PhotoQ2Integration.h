
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

#include "PROPOSAL/crossection/parametrization/Photonuclear.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/methods.h"
#include "PROPOSAL/Logging.h"

#define Q2_PHOTO_PARAM_INTEGRAL_DEC(param)                                                                             \
    class Photo##param : public PhotoQ2Integral                                                                        \
    {                                                                                                                  \
    public:                                                                                                            \
        Photo##param(const ParticleDef&,                                                                               \
                     std::shared_ptr<const Medium>,                                                                                    \
                     double multiplier,                                                                                \
                     const ShadowEffect& shadow_effect);                                                               \
        Photo##param(const Photo##param&);                                                                             \
        virtual ~Photo##param();                                                                                       \
                                                                                                                       \
        virtual Parametrization* clone() const { return new Photo##param(*this); }                                     \
        static Photonuclear* create(const ParticleDef& particle_def,                                                   \
                                    std::shared_ptr<const Medium> medium,                                                              \
                                    double multiplier,                                                                 \
                                    const ShadowEffect& shadow_effect)                                                 \
        {                                                                                                              \
            return new Photo##param(particle_def, medium, multiplier, shadow_effect);                                  \
        }                                                                                                              \
                                                                                                                       \
        double FunctionToQ2Integral(double energy, double v, double Q2);                                               \
                                                                                                                       \
        const std::string& GetName() const { return name_; }                                                           \
                                                                                                                       \
    protected:                                                                                                         \
        static const std::string name_;                                                                                \
    };


namespace PROPOSAL {

// class Interpolant;

/******************************************************************************
 *                            Photo Q2 Integration                            *
 ******************************************************************************/

class PhotoQ2Integral : public Photonuclear
{
public:
    PhotoQ2Integral(const ParticleDef&,
                    std::shared_ptr<const Medium>,
                    double multiplier,
                    const ShadowEffect&);
    PhotoQ2Integral(const PhotoQ2Integral&);
    virtual ~PhotoQ2Integral();

    virtual Parametrization* clone() const = 0;

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    virtual double DifferentialCrossSection(double energy, double v);

    virtual double FunctionToQ2Integral(double energy, double v, double Q2) = 0;

    // --------------------------------------------------------------------- //
    // Getter
    // --------------------------------------------------------------------- //

    virtual size_t GetHash() const;

protected:
    virtual bool compare(const Parametrization&) const;
    virtual void print(std::ostream&) const;

    ShadowEffect* shadow_effect_;
    Integral integral_;
};

/******************************************************************************
 *                     Declare Integral Parametrizations                      *
 ******************************************************************************/

Q2_PHOTO_PARAM_INTEGRAL_DEC(AbramowiczLevinLevyMaor91)
Q2_PHOTO_PARAM_INTEGRAL_DEC(AbramowiczLevinLevyMaor97)
Q2_PHOTO_PARAM_INTEGRAL_DEC(ButkevichMikhailov)
Q2_PHOTO_PARAM_INTEGRAL_DEC(RenoSarcevicSu)

#undef Q2_PHOTO_PARAM_INTEGRAL_DEC

} // namespace PROPOSAL
