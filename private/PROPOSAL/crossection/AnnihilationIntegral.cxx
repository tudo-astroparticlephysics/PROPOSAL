
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/AnnihilationIntegral.h"
#include "PROPOSAL/crossection/parametrization/Annihilation.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/RandomGenerator.h"

using namespace PROPOSAL;

AnnihilationIntegral::AnnihilationIntegral(const Annihilation& param)
        : CrossSectionIntegral(InteractionType::Particle, param), rndc_(-1.)
{
    gamma_def_ = &GammaDef::Get();
}

AnnihilationIntegral::AnnihilationIntegral(const AnnihilationIntegral& annihilation)
        : CrossSectionIntegral(annihilation), rndc_(annihilation.rndc_), gamma_def_(annihilation.gamma_def_)
{
}

AnnihilationIntegral::~AnnihilationIntegral() {}

std::pair<std::vector<DynamicData>, bool> AnnihilationIntegral::CalculateProducedParticles(double energy,
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
            rho = dndx_integral_[i].GetUpperLimit();

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

double AnnihilationIntegral::CalculateStochasticLoss(double energy, double rnd1, double rnd2) {
    if(rnd1 != rnd_){
        CalculatedNdx(energy, rnd1);
    }
    rndc_ = rnd2; // save random number for component sampling
    return energy; // losses are always catastrophic
}
