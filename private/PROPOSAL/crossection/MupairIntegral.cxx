
#include <functional>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/MupairIntegral.h"
#include "PROPOSAL/crossection/parametrization/MupairProduction.h"
#include "PROPOSAL/medium/Medium.h"

using namespace PROPOSAL;

MupairIntegral::MupairIntegral(const MupairProduction& param)
    : CrossSectionIntegral(DynamicData::MuPair, param)
{
}

MupairIntegral::MupairIntegral(const MupairIntegral& mupair)
    : CrossSectionIntegral(mupair)
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

    for (int i = 0; i < parametrization_->GetMedium().GetNumComponents(); i++)
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

std::vector<Particle*> MupairIntegral::CalculateProducedParticles(double energy, double energy_loss, double rnd1, double rnd2){

    //Create MuPair particles
    std::vector<Particle*> mupair;
    mupair.push_back(new Particle(MuMinusDef::Get()));
    mupair.push_back(new Particle(MuPlusDef::Get()));

    //Sample and assign energies
    double rho = parametrization_->Calculaterho(energy, energy_loss/energy, rnd1, rnd2);

    mupair[0]->SetEnergy(0.5*energy_loss*(1 + rho));
    mupair[1]->SetEnergy(0.5*energy_loss*(1 - rho));
    return mupair;

}

