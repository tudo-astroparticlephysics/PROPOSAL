
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/PhotoPairIntegral.h"
#include "PROPOSAL/crossection/parametrization/PhotoPairProduction.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Logging.h"


using namespace PROPOSAL;

PhotoPairIntegral::PhotoPairIntegral(const PhotoPairProduction& param, const PhotoAngleDistribution& photoangle)
        : CrossSectionIntegral(InteractionType::Particle, param), photoangle_(photoangle.clone()), rndc_(-1)
{
    eminus_def_ = &EMinusDef::Get();
    eplus_def_ = &EPlusDef::Get();
}

PhotoPairIntegral::PhotoPairIntegral(const PhotoPairIntegral& photo)
        : CrossSectionIntegral(photo), photoangle_(photo.GetPhotoAngleDistribution().clone()), rndc_(photo.rndc_)
        , eminus_def_(photo.eminus_def_), eplus_def_(photo.eplus_def_)
{
}

PhotoPairIntegral::~PhotoPairIntegral() {
    delete photoangle_;
}

bool PhotoPairIntegral::compare(const CrossSection& cross_section) const {
    const PhotoPairIntegral* cross_section_integral =
            static_cast<const PhotoPairIntegral*>(&cross_section);


    if(*photoangle_!=*cross_section_integral->photoangle_)
        return false;

    return CrossSectionIntegral::compare(cross_section);
}

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
std::pair<std::vector<DynamicData>, bool> PhotoPairIntegral::CalculateProducedParticles(double energy, double energy_loss, const Vector3D& initial_direction){
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
            rho = dndx_integral_[i].GetUpperLimit();

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
