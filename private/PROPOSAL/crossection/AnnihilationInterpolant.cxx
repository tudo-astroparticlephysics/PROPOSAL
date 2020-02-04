
#include <functional>

#include "PROPOSAL/crossection/AnnihilationIntegral.h"
#include "PROPOSAL/crossection/AnnihilationInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Annihilation.h"

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/RandomGenerator.h"


using namespace PROPOSAL;

AnnihilationInterpolant::AnnihilationInterpolant(const Annihilation& param, InterpolationDef def)
        : CrossSectionInterpolant(InteractionType::Particle, param), rndc_(-1.) {
    // Use parent CrossSecition dNdx interpolation
    InitdNdxInterpolation(def);
    gamma_def_ = &GammaDef::Get();
}

AnnihilationInterpolant::AnnihilationInterpolant(const AnnihilationInterpolant& annihilation)
        : CrossSectionInterpolant(annihilation), rndc_(annihilation.rndc_), gamma_def_(annihilation.gamma_def_)
{
}

AnnihilationInterpolant::~AnnihilationInterpolant() {}

bool AnnihilationInterpolant::compare(const CrossSection& cross_section) const
{
    const AnnihilationInterpolant* cross_section_interpolant =
            static_cast<const AnnihilationInterpolant*>(&cross_section);

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

    return true;
}

double AnnihilationInterpolant::CalculateStochasticLoss(double energy, double rnd1, double rnd2) {
    if(rnd1 != rnd_){
        CalculatedNdx(energy, rnd1);
    }
    rndc_ = rnd2; // save random number for component sampling
    return energy; // losses are always catastrophic
}

std::pair<std::vector<DynamicData>, bool> AnnihilationInterpolant::CalculateProducedParticles(double energy,
                                                                                            double energy_loss,
                                                                                            const Vector3D& initial_direction) {
    (void)energy_loss;
    double rnd, rsum, rho;

    std::vector<DynamicData> particle_list{};

    if(rndc_<0){
        //CalculateStochasticLoss has never been called before, return empty list
        //TODO: find a better way of checking the random numbers
        log_warn("CalculateProducedParticles has been called with no call of CalculateStochasticLoss for PhotoPairProduction");
        return std::make_pair(particle_list, true);
    }

    particle_list.push_back(DynamicData(gamma_def_->particle_type));
    particle_list.push_back(DynamicData(gamma_def_->particle_type));

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

            // The available energy is the positron energy plus the mass of the electron
            particle_list[0].SetEnergy((energy + ME) * (1-rho));
            particle_list[1].SetEnergy((energy + ME) * rho);

            particle_list[0].SetDirection(initial_direction);
            particle_list[1].SetDirection(initial_direction);

            double cosphi0 = ((energy + ME) * (1. - rho) - ME)/( (1. - rho) * std::sqrt( (energy + ME) * (energy - ME) ) );
            double cosphi1 = ((energy + ME) * rho - ME)/( rho * std::sqrt((energy + ME) * (energy - ME)));

            double rndtheta = RandomGenerator::Get().RandomDouble();

            particle_list[0].DeflectDirection(cosphi0, rndtheta * 2. * PI);
            particle_list[1].DeflectDirection(cosphi1, std::fmod(rndtheta * 2. * PI + PI, 2. * PI));

            return std::make_pair(particle_list, true);
        }
    }


    log_fatal("could not sample ProducedParticles for PhotoPairProduction!");
    return std::make_pair(particle_list, true);

}
