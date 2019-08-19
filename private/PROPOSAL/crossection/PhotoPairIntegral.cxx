
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/PhotoPairIntegral.h"
#include "PROPOSAL/crossection/parametrization/PhotoPairProduction.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Logging.h"


using namespace PROPOSAL;

PhotoPairIntegral::PhotoPairIntegral(const PhotoPairProduction& param)
        : CrossSectionIntegral(DynamicData::Epair, param), rndc_(-1)
{
}

PhotoPairIntegral::PhotoPairIntegral(const PhotoPairIntegral& photo)
        : CrossSectionIntegral(photo), rndc_(photo.rndc_)
{
}

PhotoPairIntegral::~PhotoPairIntegral() {}



// ------------------------------------------------------------------------- //
double PhotoPairIntegral::CalculateStochasticLoss(double energy, double rnd1, double rnd2)
{
    if(rnd1 != rnd_){
        CalculatedNdx(energy, rnd1);
    }
    rndc_ = rnd2; // save random number for component sampling
    return energy; // losses are always catastrophic

}

// ------------------------------------------------------------------------- //
std::pair<std::vector<Particle*>, bool> PhotoPairIntegral::CalculateProducedParticles(double energy, double energy_loss){
    (void)energy_loss;
    double rnd;
    double rsum;
    double rho;

    std::vector<Particle*> particle_list{};

    if(rndc_<0){
        //CalculateStochasticLoss has never been called before, return empty list
        //TODO: find a better way of checking the random numbers
        log_debug("CalculateProducedParticles has been called with no call of CalculateStochasticLoss for PhotoPairProduction");
        return std::make_pair(particle_list, true);
    }

    particle_list.push_back(new Particle(EPlusDef::Get()));
    particle_list.push_back(new Particle(EMinusDef::Get()));

    rnd  = rndc_ * sum_of_rates_;
    rsum = 0;

    for (size_t i = 0; i < components_.size(); ++i)
    {
        rsum += prob_for_component_[i];

        if (rsum > rnd)
        {
            parametrization_->SetCurrentComponent(i);
            rho = dndx_integral_[i].GetUpperLimit();

            particle_list[0]->SetEnergy(energy * (1-rho));
            particle_list[1]->SetEnergy(energy * rho);

            return std::make_pair(particle_list, true);
        }
    }


    log_fatal("could not sample ProducedParticles for PhotoPairProduction!");
    return std::make_pair(particle_list, true);
}