
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

#define PHOTO_PARAM_REAL_DEC(param, parent)                                                                            \
    class Photo##param : public Photo##parent                                                                          \
    {                                                                                                                  \
    public:                                                                                                            \
        Photo##param(const ParticleDef&,                                                                               \
                     const Medium&,                                                                                    \
                     const EnergyCutSettings&,                                                                         \
                     double multiplier,                                                                                \
                     bool hard_component);                                                                             \
        Photo##param(const Photo##param&);                                                                             \
        virtual ~Photo##param();                                                                                       \
                                                                                                                       \
        Parametrization* clone() const { return new Photo##param(*this); }                                             \
        static Photonuclear* create(const ParticleDef& particle_def,                                                   \
                                    const Medium& medium,                                                              \
                                    const EnergyCutSettings& cuts,                                                     \
                                    double multiplier,                                                                 \
                                    bool hard_component)                                                               \
        {                                                                                                              \
            return new Photo##param(particle_def, medium, cuts, multiplier, hard_component);                           \
        }                                                                                                              \
                                                                                                                       \
        virtual double CalculateParametrization(double nu);                                                            \
                                                                                                                       \
        const std::string& GetName() const { return name_; }                                                           \
                                                                                                                       \
    private:                                                                                                           \
        static const std::string name_;                                                                                \
    };

namespace PROPOSAL {

/******************************************************************************
 *                         PhotoRealPhotonAssumption                           *
 ******************************************************************************/

class PhotoRealPhotonAssumption : public Photonuclear
{
public:
    PhotoRealPhotonAssumption(const ParticleDef&,
                              const Medium&,
                              const EnergyCutSettings&,
                              double multiplier,
                              bool hard_component);
    PhotoRealPhotonAssumption(const PhotoRealPhotonAssumption&);
    virtual ~PhotoRealPhotonAssumption();

    virtual Parametrization* clone() const = 0;

    // ----------------------------------------------------------------- //
    // Public methods
    // ----------------------------------------------------------------- //

    virtual double DifferentialCrossSection(double energy, double v);

    virtual double CalculateParametrization(double nu) = 0;
    double NucleusCrossSectionCaldwell(double nu);

    // --------------------------------------------------------------------- //
    // Getter
    // --------------------------------------------------------------------- //

    virtual size_t GetHash() const;

protected:
    virtual bool compare(const Parametrization&) const;
    virtual void print(std::ostream&) const;

    RealPhoton* hard_component_;
};

/******************************************************************************
 *                       Zeus, BezrukovBugaev, Kokoulin                       *
 ******************************************************************************/

// Signature: (new class, parent class)
PHOTO_PARAM_REAL_DEC(Zeus, RealPhotonAssumption)
PHOTO_PARAM_REAL_DEC(BezrukovBugaev, RealPhotonAssumption)
PHOTO_PARAM_REAL_DEC(Kokoulin, BezrukovBugaev) // Kokoulin derives from BezrukovBugaev

/******************************************************************************
 *                           Rhode Parametrization                            *
 ******************************************************************************/

class PhotoRhode : public PhotoRealPhotonAssumption
{
public:
    PhotoRhode(const ParticleDef&, const Medium&, const EnergyCutSettings&, double multiplier, bool hard_component);
    PhotoRhode(const PhotoRhode&);
    virtual ~PhotoRhode();

    Parametrization* clone() const { return new PhotoRhode(*this); }
    static Photonuclear* create(const ParticleDef&,
                                const Medium&,
                                const EnergyCutSettings&,
                                double multiplier,
                                bool hard_component);

    double CalculateParametrization(double nu);

    const std::string& GetName() const { return name_; }

private:
    virtual bool compare(const Parametrization&) const;

    double MeasuredSgN(double e);

    static const std::string name_;
    Interpolant* interpolant_;
};

#undef Q2_PHOTO_PARAM_INTEGRAL_DEC

} // namespace PROPOSAL
