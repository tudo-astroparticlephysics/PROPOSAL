
#include <functional>

#include "PROPOSAL/crossection/MupairIntegral.h"
#include "PROPOSAL/crossection/MupairInterpolant.h"
#include "PROPOSAL/crossection/parametrization/MupairProduction.h"

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"

#include "PROPOSAL/math/RandomGenerator.h"


using namespace PROPOSAL;

MupairInterpolant::MupairInterpolant(const MupairProduction& param, InterpolationDef def)
    : CrossSectionInterpolant(GetType(param), param)
{
    // Use parent CrossSecition dNdx interpolation
    InitdNdxInterpolation(def);

    // --------------------------------------------------------------------- //
    // Builder for DEdx
    // --------------------------------------------------------------------- //

    Interpolant1DBuilder builder1d;
    Helper::InterpolantBuilderContainer builder_container;

    // Needed for CalculatedEdx integration
    MupairIntegral mupair(param);

    builder1d.SetMax(def.nodes_cross_section)
        .SetXMin(param.GetParticleDef().mass)
        .SetXMax(def.max_node_energy)
        .SetRomberg(def.order_of_interpolation)
        .SetRational(true)
        .SetRelative(false)
        .SetIsLog(true)
        .SetRombergY(def.order_of_interpolation)
        .SetRationalY(false)
        .SetRelativeY(false)
        .SetLogSubst(true)
        .SetFunction1D(std::bind(&CrossSectionIntegral::CalculatedEdxWithoutMultiplier, &mupair, std::placeholders::_1));

    builder_container.push_back(std::make_pair(&builder1d, &dedx_interpolant_));

    // --------------------------------------------------------------------- //
    // Builder for DE2dx
    // --------------------------------------------------------------------- //

    Interpolant1DBuilder builder_de2dx;
    Helper::InterpolantBuilderContainer builder_container_de2dx;

    builder_de2dx.SetMax(def.nodes_continous_randomization)
        .SetXMin(param.GetParticleDef().mass)
        .SetXMax(def.max_node_energy)
        .SetRomberg(def.order_of_interpolation)
        .SetRational(false)
        .SetRelative(false)
        .SetIsLog(true)
        .SetRombergY(def.order_of_interpolation)
        .SetRationalY(false)
        .SetRelativeY(false)
        .SetLogSubst(false)
        .SetFunction1D(std::bind(&MupairIntegral::CalculatedE2dxWithoutMultiplier, &mupair, std::placeholders::_1));

    builder_container_de2dx.push_back(std::make_pair(&builder_de2dx, &de2dx_interpolant_));

    Helper::InitializeInterpolation("dEdx", builder_container, std::vector<Parametrization*>(1, parametrization_), def);
    Helper::InitializeInterpolation(
        "dE2dx", builder_container_de2dx, std::vector<Parametrization*>(1, parametrization_), def);

    muminus_def_ = &MuMinusDef::Get();
    muplus_def_ = &MuPlusDef::Get();
}

MupairInterpolant::MupairInterpolant(const MupairInterpolant& mupair)
    : CrossSectionInterpolant(mupair)
    , muminus_def_(mupair.muminus_def_)
    , muplus_def_(mupair.muplus_def_)
{
}

MupairInterpolant::~MupairInterpolant() {}

// ----------------------------------------------------------------- //
// Public methods
// ----------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double MupairInterpolant::CalculatedEdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return parametrization_->GetMultiplier() * std::max(dedx_interpolant_->Interpolate(energy), 0.0);
}

std::pair<std::vector<DynamicData>, bool> MupairInterpolant::CalculateProducedParticles(double energy, double energy_loss, const Vector3D& initial_direction){
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

InteractionType MupairInterpolant::GetType(const MupairProduction& param){
    if(param.IsParticleOutputEnabled()){
        return InteractionType::Particle;
    }
    else{
        return InteractionType::MuPair;
    }
}
