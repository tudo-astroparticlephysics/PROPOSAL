
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/MupairIntegral.h"
#include "PROPOSAL/crossection/parametrization/MupairProduction.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/math/RandomGenerator.h"


using namespace PROPOSAL;

MupairIntegral::MupairIntegral(const MupairProduction& param)
    : CrossSectionIntegral(GetType(param), param)
{
    muminus_def_ = &MuMinusDef::Get();
    muplus_def_ = &MuPlusDef::Get();
}

MupairIntegral::MupairIntegral(const MupairIntegral& mupair)
    : CrossSectionIntegral(mupair)
    , muminus_def_(mupair.muminus_def_)
    , muplus_def_(mupair.muplus_def_)
{
}

MupairIntegral::~MupairIntegral() {}

// ----------------------------------------------------------------- //
// Public methods
// ----------------------------------------------------------------- //

double MupairIntegral::CalculatedEdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return parametrization_->GetMultiplier() * MupairIntegral::CalculatedEdxWithoutMultiplier(energy);
}

double MupairIntegral::CalculatedEdxWithoutMultiplier(double energy)
{
    double sum = 0;

    for (int i = 0; i < parametrization_->GetMedium()->GetNumComponents(); i++)
    {
        parametrization_->SetCurrentComponent(i);
        Parametrization::IntegralLimits limits = parametrization_->GetIntegralLimits(energy);

        sum += dedx_integral_.Integrate(
            limits.vMin,
            limits.vUp,
            std::bind(&Parametrization::FunctionToDEdxIntegral, parametrization_, energy, std::placeholders::_1),
            4);
    }

    return energy * sum;
}

std::pair<std::vector<DynamicData>, bool> MupairIntegral::CalculateProducedParticles(double energy, double energy_loss, const Vector3D& initial_direction){
    std::vector<DynamicData> mupair;

    if(parametrization_->IsParticleOutputEnabled() == false){
        return std::make_pair(mupair, false);
    }

    //Create MuPair particles
    mupair.push_back(DynamicData(muminus_def_->particle_type));
    mupair.push_back(DynamicData(muplus_def_->particle_type));

    //Sample random numbers
    double rnd1 = RandomGenerator::Get().RandomDouble();
    double rnd2 = RandomGenerator::Get().RandomDouble();

    //Sample and assign energies
    double rho = parametrization_->Calculaterho(energy, energy_loss/energy, rnd1, rnd2);

    mupair[0].SetEnergy(0.5*energy_loss*(1 + rho));
    mupair[1].SetEnergy(0.5*energy_loss*(1 - rho));
    mupair[0].SetDirection(initial_direction);
    mupair[1].SetDirection(initial_direction);
    return std::make_pair(mupair, false);

}

InteractionType MupairIntegral::GetType(const MupairProduction& param){
    if(param.IsParticleOutputEnabled()){
        return InteractionType::Particle;
    }
    else{
        return InteractionType::MuPair;
    }
}
