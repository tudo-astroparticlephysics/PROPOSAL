
#pragma once

#include <boost/bind.hpp>

#include "PROPOSAL/crossection/parametrization/Photonuclear.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/math/Interpolant.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"

#define Q2_PHOTO_PARAM_INTEGRAL_DEC(param)                                                                             \
    class Photo##param : public PhotoQ2Integral                                                                        \
    {                                                                                                                  \
        public:                                                                                                        \
        Photo##param(const ParticleDef&,                                                                               \
                     const Medium&,                                                                                    \
                     const EnergyCutSettings&,                                                                         \
                     const ShadowEffect& shadow_effect,                                                                \
                     double multiplier);                                                                               \
        Photo##param(const Photo##param&);                                                                             \
        virtual ~Photo##param();                                                                                       \
                                                                                                                       \
        virtual Parametrization* clone() const { return new Photo##param(*this); }                                     \
        static Photonuclear* create(const ParticleDef& particle_def,                                                   \
                                    const Medium& medium,                                                              \
                                    const EnergyCutSettings& cuts,                                                     \
                                    const ShadowEffect& shadow_effect,                                                 \
                                    double multiplier)                                                                 \
        {                                                                                                              \
            return new Photo##param(particle_def, medium, cuts, shadow_effect, multiplier);                            \
        }                                                                                                              \
                                                                                                                       \
        double FunctionToQ2Integral(double energy, double v, double Q2);                                               \
                                                                                                                       \
        const std::string& GetName() const { return name_; }                                                           \
                                                                                                                       \
        protected:                                                                                                     \
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
                    const Medium&,
                    const EnergyCutSettings&,
                    const ShadowEffect&,
                    double multiplier);
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

/******************************************************************************
*                    Declare Interpolant Parametrizations                    *
******************************************************************************/

template <class Param = PhotoAbramowiczLevinLevyMaor97>
class PhotoQ2Interpolant : public Param
{
    public:
    typedef std::vector<Interpolant*> InterpolantVec;

    public:
    PhotoQ2Interpolant(const ParticleDef&,
                              const Medium&,
                              const EnergyCutSettings&,
                              const ShadowEffect&,
                              double multiplier,
                              InterpolationDef def = InterpolationDef());
    PhotoQ2Interpolant(const PhotoQ2Interpolant&);
    virtual ~PhotoQ2Interpolant();

    Parametrization* clone() const { return new PhotoQ2Interpolant<Param>(*this); }
    static Photonuclear* create(const ParticleDef& particle_def,
                                const Medium& medium,
                                const EnergyCutSettings& cuts,
                                const ShadowEffect& shadow_effect,
                                double multiplier,
                                InterpolationDef def = InterpolationDef())
    {
        return new PhotoQ2Interpolant<Param>(particle_def, medium, cuts, shadow_effect, multiplier, def);
    }

    double DifferentialCrossSection(double energy, double v);

    protected:
    double FunctionToBuildPhotoInterpolant(double energy, double v, int component);

    InterpolantVec interpolant_;
};

template <class Param>
PhotoQ2Interpolant<Param>::PhotoQ2Interpolant(const ParticleDef& particle_def,
                                                     const Medium& medium,
                                                     const EnergyCutSettings& cuts,
                                                     const ShadowEffect& shadow_effect,
                                                     double multiplier,
                                                     InterpolationDef def)
    : Param(particle_def, medium, cuts, shadow_effect, multiplier)
    , interpolant_(this->medium_->GetNumComponents(), NULL)
{
    std::vector<Interpolant2DBuilder> builder2d(this->components_.size());
    Helper::InterpolantBuilderContainer builder_container2d(this->components_.size());

    for (unsigned int i = 0; i < this->components_.size(); ++i)
    {
        builder2d[i]
            .SetMax1(NUM1)
            .SetX1Min(this->particle_def_.low)
            .SetX1Max(BIGENERGY)
            .SetMax2(NUM1)
            .SetX2Min(0.0)
            .SetX2Max(1.0)
            .SetRomberg1(def.order_of_interpolation)
            .SetRational1(false)
            .SetRelative1(false)
            .SetIsLog1(true)
            .SetRomberg2(def.order_of_interpolation)
            .SetRational2(false)
            .SetRelative2(false)
            .SetIsLog2(false)
            .SetRombergY(def.order_of_interpolation)
            .SetRationalY(false)
            .SetRelativeY(false)
            .SetLogSubst(false)
            .SetFunction2D(
                boost::bind(&PhotoQ2Interpolant::FunctionToBuildPhotoInterpolant, this, _1, _2, i));

        builder_container2d[i].first  = &builder2d[i];
        builder_container2d[i].second = &interpolant_[i];
    }

    Helper::InitializeInterpolation("Photo", builder_container2d, std::vector<Parametrization*>(1, this), def);
}

template <class Param>
PhotoQ2Interpolant<Param>::PhotoQ2Interpolant(const PhotoQ2Interpolant& photo)
    : Param(photo)
    , interpolant_()
{
    interpolant_.resize(photo.interpolant_.size());

    for(unsigned int i = 0; i < photo.interpolant_.size(); ++i)
    {
        interpolant_[i] = new Interpolant(*photo.interpolant_[i]);
    }
}

template <class Param>
PhotoQ2Interpolant<Param>::~PhotoQ2Interpolant()
{
}

template <class Param>
double PhotoQ2Interpolant<Param>::DifferentialCrossSection(double energy, double v)
{
    Parametrization::IntegralLimits limits = this->GetIntegralLimits(energy);

    if (v >= limits.vUp)
    {
        return std::max(interpolant_[this->component_index_]->Interpolate(
                            energy, log(v / limits.vUp) / log(limits.vMax / limits.vUp)),
                        0.0);
    }

    return Param::DifferentialCrossSection(energy, v);
}

template <class Param>
double PhotoQ2Interpolant<Param>::FunctionToBuildPhotoInterpolant(double energy, double v, int component)
{
    this->component_index_      = component;
    Parametrization::IntegralLimits limits = this->GetIntegralLimits(energy);

    if (limits.vUp == limits.vMax)
    {
        return 0;
    }

    v = limits.vUp * exp(v * log(limits.vMax / limits.vUp));

    return Param::DifferentialCrossSection(energy, v);
}

#undef Q2_PHOTO_PARAM_INTEGRAL_DEC

} /* PROPOSAL */
