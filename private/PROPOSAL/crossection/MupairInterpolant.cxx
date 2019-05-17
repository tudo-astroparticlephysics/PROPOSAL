
#include <functional>

#include "PROPOSAL/crossection/MupairIntegral.h"
#include "PROPOSAL/crossection/MupairInterpolant.h"
#include "PROPOSAL/crossection/parametrization/MupairProduction.h"

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

MupairInterpolant::MupairInterpolant(const MupairProduction& param, InterpolationDef def)
    : CrossSectionInterpolant(DynamicData::MuPair, param)
{
    // Use parent CrossSecition dNdx interpolation
    InitdNdxInerpolation(def);

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
}

MupairInterpolant::MupairInterpolant(const MupairInterpolant& mupair)
    : CrossSectionInterpolant(mupair)
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

std::vector<Particle*> MupairInterpolant::CalculateProducedParticles(double energy, double energy_loss, double rnd1, double rnd2){

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
