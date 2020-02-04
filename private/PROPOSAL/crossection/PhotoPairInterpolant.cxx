
#include <functional>
#include <cmath>

#include "PROPOSAL/crossection/PhotoPairIntegral.h"
#include "PROPOSAL/crossection/PhotoPairInterpolant.h"
#include "PROPOSAL/crossection/parametrization/PhotoPairProduction.h"

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

PhotoPairInterpolant::PhotoPairInterpolant(const PhotoPairProduction& param, const PhotoAngleDistribution& photoangle, InterpolationDef def)
        : CrossSectionInterpolant(InteractionType::Particle, param)
        , photoangle_(photoangle.clone()), rndc_(-1.){
    // Use own initialization
    PhotoPairInterpolant::InitdNdxInterpolation(def);
    eminus_def_ = &EMinusDef::Get();
    eplus_def_ = &EPlusDef::Get();
}


PhotoPairInterpolant::PhotoPairInterpolant(const PhotoPairInterpolant& param)
        : CrossSectionInterpolant(param), photoangle_(param.GetPhotoAngleDistribution().clone()), rndc_(param.rndc_)
        , eminus_def_(param.eminus_def_), eplus_def_(param.eplus_def_)
{
}

PhotoPairInterpolant::~PhotoPairInterpolant() {
        delete photoangle_;
}

bool PhotoPairInterpolant::compare(const CrossSection& cross_section) const
{
    // to avoid errors in CrossSectionInterpolant (due to non-defined dEdx interpolant objects) comparison is done here
    const PhotoPairInterpolant* cross_section_interpolant =
            static_cast<const PhotoPairInterpolant*>(&cross_section);

    if (dndx_interpolant_1d_.size() != cross_section_interpolant->dndx_interpolant_1d_.size())
        return false;
    else if (dndx_interpolant_2d_.size() != cross_section_interpolant->dndx_interpolant_2d_.size())
        return false;

    for (unsigned int i = 0; i < dndx_interpolant_1d_.size(); ++i)
    {
        if (*dndx_interpolant_1d_[i] != *cross_section_interpolant->dndx_interpolant_1d_[i])
            return false;
    }
    for (unsigned int i = 0; i < dndx_interpolant_2d_.size(); ++i)
    {
        if (*dndx_interpolant_2d_[i] != *cross_section_interpolant->dndx_interpolant_2d_[i])
            return false;
    }
    if(*photoangle_!=*cross_section_interpolant->photoangle_)
        return false;

    return true;
}

// ------------------------------------------------------------------------- //
double PhotoPairInterpolant::CalculatedNdx(double energy) {
    if(energy < 2. * ME){
        return 0;
    } else
        return CrossSectionInterpolant::CalculatedNdx(energy);
}

// ------------------------------------------------------------------------- //
double PhotoPairInterpolant::CalculatedNdx(double energy, double rnd) {
    if(energy < 2. * ME){
        return 0;
    } else
        return CrossSectionInterpolant::CalculatedNdx(energy, rnd);
}

// ------------------------------------------------------------------------- //
double PhotoPairInterpolant::CalculateStochasticLoss(double energy, double rnd1, double rnd2)
{
    if(rnd1 != rnd_){
        CalculatedNdx(energy, rnd1);
    }
    rndc_ = rnd2; // save random number for component sampling
    return energy; // losses are always catastrophic

}


// ------------------------------------------------------------------------- //
std::pair<std::vector<DynamicData>, bool> PhotoPairInterpolant::CalculateProducedParticles(double energy, double energy_loss, const Vector3D& initial_direction){
    (void)energy_loss;
    double rnd;
    double rsum;
    double rho;

    std::vector<DynamicData> particle_list{};

    if(rndc_<0){
        //CalculateStochasticLoss has never been called before, return empty list
        //TODO: find a better way of checking the random numbers
        log_warn("CalculateProducedParticles has been called with no call of CalculateStochasticLoss for PhotoPairProduction");
        return std::make_pair(particle_list, true);
    }

    particle_list.push_back(DynamicData(eplus_def_->particle_type));
    particle_list.push_back(DynamicData(eminus_def_->particle_type));

    rnd  = rndc_ * sum_of_rates_;
    rsum = 0;

    for (size_t i = 0; i < components_.size(); ++i)
    {
        rsum += prob_for_component_[i];

        if (rsum > rnd)
        {
            parametrization_->SetCurrentComponent(i);
            Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);
            rho = (limits.vUp * std::exp(dndx_interpolant_2d_.at(i)->FindLimit(energy, rnd_ * prob_for_component_[i]) *
                                         std::log(limits.vMax / limits.vUp)));

            particle_list[0].SetEnergy(energy * (1-rho));
            particle_list[1].SetEnergy(energy * rho);

            PhotoAngleDistribution::DeflectionAngles angles;
            angles = photoangle_->SampleAngles(energy, rho, i);
            particle_list[0].SetDirection(initial_direction);
            particle_list[1].SetDirection(initial_direction);
            particle_list[0].DeflectDirection(angles.cosphi0, angles.theta0);
            particle_list[1].DeflectDirection(angles.cosphi1, angles.theta1);

            return std::make_pair(particle_list, true);
        }
    }

    log_fatal("could not sample ProducedParticles for PhotoPairProduction!");
    return std::make_pair(particle_list, true);
}

// ------------------------------------------------------------------------- //
void PhotoPairInterpolant::InitdNdxInterpolation(const InterpolationDef& def)
{
    // --------------------------------------------------------------------- //
    // Builder for dNdx with no logarithmic energy
    // --------------------------------------------------------------------- //

    std::vector<Interpolant1DBuilder> builder1d(components_.size());
    std::vector<Interpolant2DBuilder> builder2d(components_.size());

    Helper::InterpolantBuilderContainer builder_container1d(components_.size());
    Helper::InterpolantBuilderContainer builder_container2d(components_.size());
    Helper::InterpolantBuilderContainer builder_return;

    Integral integral(IROMB, IMAXS, IPREC);

    for (unsigned int i = 0; i < components_.size(); ++i)
    {
        // !!! IMPORTANT !!!
        // Order of builder matter because the functions needed for 1d interpolation
        // needs the already intitialized 2d interpolants.
        builder2d[i]
                .SetMax1(def.nodes_cross_section)
                .SetX1Min(ME)
                .SetX1Max(def.max_node_energy)
                .SetMax2(def.nodes_cross_section)
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
                .SetRationalY(true)
                .SetRelativeY(false)
                .SetLogSubst(false)
                .SetFunction2D(std::bind(
                        &CrossSectionInterpolant::FunctionToBuildDNdxInterpolant2D,
                        this,
                        std::placeholders::_1,
                        std::placeholders::_2,
                        std::ref(integral),
                        i));

        builder_container2d[i].first  = &builder2d[i];
        builder_container2d[i].second = &dndx_interpolant_2d_[i];

        builder1d[i]
                .SetMax(def.nodes_cross_section)
                .SetXMin(ME)
                .SetXMax(def.max_node_energy)
                .SetRomberg(def.order_of_interpolation)
                .SetRational(false)
                .SetRelative(false)
                .SetIsLog(true)
                .SetRombergY(def.order_of_interpolation)
                .SetRationalY(true)
                .SetRelativeY(false)
                .SetLogSubst(false)
                .SetFunction1D(std::bind(&CrossSectionInterpolant::FunctionToBuildDNdxInterpolant, this, std::placeholders::_1, i));

        builder_container1d[i].first  = &builder1d[i];
        builder_container1d[i].second = &dndx_interpolant_1d_[i];
    }

    builder_return.insert(builder_return.end(), builder_container2d.begin(), builder_container2d.end());
    builder_return.insert(builder_return.end(), builder_container1d.begin(), builder_container1d.end());
    // builder2d.insert(builder2d.end(), builder1d.begin(), builder1d.end());

    Helper::InitializeInterpolation("dNdx", builder_return, std::vector<Parametrization*>(1, parametrization_), def);
}
