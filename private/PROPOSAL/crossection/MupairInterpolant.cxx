
#include <functional>
#include <PROPOSAL/crossection/EpairInterpolant.h>

#include "PROPOSAL/crossection/MupairIntegral.h"
#include "PROPOSAL/crossection/MupairInterpolant.h"
#include "PROPOSAL/crossection/parametrization/MupairProduction.h"

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"

#include "PROPOSAL/math/RandomGenerator.h"


using namespace PROPOSAL;

MupairInterpolant::MupairInterpolant(const MupairProduction& param, std::shared_ptr<const EnergyCutSettings> cuts, const InterpolationDef& def)
    : CrossSectionInterpolant(param, cuts)
    , interpolant_(param.GetMedium()->GetNumComponents())
{
    // Interpolate differential cross section for dNdx interpolation (if necessary)
    if(interpolateDifferentialCrossSection[param.GetName()]){
        InitializeDifferentialCrossSectionInterpolation(def);
        functiondNdxIntegral = std::bind(&MupairInterpolant::InterpolatedCrossSection, this,
                                         std::placeholders::_1, std::placeholders::_2);
    }
    else{
        functiondNdxIntegral = std::bind(&Parametrization::FunctionToDNdxIntegral, &this->GetParametrization(),
                                         std::placeholders::_1, std::placeholders::_2);
    }

    // Use parent CrossSecition dNdx interpolation
    InitdNdxInterpolation(def);

    // --------------------------------------------------------------------- //
    // Builder for DEdx
    // --------------------------------------------------------------------- //

    Interpolant1DBuilder builder1d;

    // Needed for CalculatedEdx integration
    MupairIntegral mupair(param, cuts);

    builder1d.SetMax(def.nodes_cross_section)
        .SetXMin(param.GetParticleMass())
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

    // --------------------------------------------------------------------- //
    // Builder for DE2dx
    // --------------------------------------------------------------------- //

    Interpolant1DBuilder builder_de2dx;

    builder_de2dx.SetMax(def.nodes_continous_randomization)
        .SetXMin(param.GetParticleMass())
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

    dedx_interpolant_ = Helper::InitializeInterpolation("dEdx", builder1d, parametrization_->GetHash(), def);
    de2dx_interpolant_ = Helper::InitializeInterpolation("dE2dx", builder_de2dx, parametrization_->GetHash(), def);

    muminus_def_ = &MuMinusDef::Get();
    muplus_def_ = &MuPlusDef::Get();
}

/*MupairInterpolant::MupairInterpolant(const MupairInterpolant& mupair)
    : CrossSectionInterpolant(mupair)
    , muminus_def_(mupair.muminus_def_)
    , muplus_def_(mupair.muplus_def_)
{
}*/

MupairInterpolant::~MupairInterpolant() {}

bool MupairInterpolant::compare(const CrossSection& cross_section) const{
    const MupairInterpolant* mupair = static_cast<const MupairInterpolant*>(&cross_section);

    if (interpolant_.size() != mupair->interpolant_.size())
        return false;

    for (unsigned int i = 0; i < interpolant_.size(); ++i)
    {
        if (*interpolant_[i] != *mupair->interpolant_[i])
            return false;
    }

    return CrossSectionInterpolant::compare(cross_section);
}

// ----------------------------------------------------------------- //
// Public methods
// ----------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double MupairInterpolant::CalculatedEdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0) {
        return 0;
    }

    return parametrization_->GetMultiplier() * std::max(dedx_interpolant_->Interpolate(energy), 0.0);
}

std::pair<std::vector<DynamicData>, bool> MupairInterpolant::CalculateProducedParticles(double energy, double energy_loss, const Vector3D& initial_direction){
    std::vector<DynamicData> mupair;

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

//----------------------------------------------------------------------------//
double MupairInterpolant::FunctionToBuildDNdxInterpolant2D(double energy,
                                                          double v,
                                                          Integral& integral,
                                                          int component)
{
    parametrization_->SetCurrentComponent(component);
    Parametrization::KinematicLimits limits = parametrization_->GetKinematicLimits(energy);

    double vUp = GetEnergyCut(energy);

    if (vUp == limits.vMax)
    {
        return 0;
    }

    v = vUp * std::exp(v * std::log(limits.vMax / vUp));

    return integral.Integrate(
            vUp, v, std::bind(functiondNdxIntegral, energy, std::placeholders::_1), 4);
}

//----------------------------------------------------------------------------//
double MupairInterpolant::FunctionToBuildDiffCrossSectionInterpolant(double energy, double v, int component) {
    parametrization_->SetCurrentComponent(component);
    Parametrization::KinematicLimits limits = parametrization_->GetKinematicLimits(energy);
    double cut = GetEnergyCut(energy);

    if (cut == limits.vMax)
    {
        return 0;
    }

    v = cut * std::exp(v * std::log(limits.vMax / cut));
    return parametrization_->DifferentialCrossSection(energy, v);
}

//----------------------------------------------------------------------------//
double MupairInterpolant::InterpolatedCrossSection(double energy, double v) {
    Parametrization::KinematicLimits limits = parametrization_->GetKinematicLimits(energy);
    double cut = GetEnergyCut(energy);
    if (v >= cut)
    {
        return std::max(
                interpolant_.at(parametrization_->GetCurrentComponent())->Interpolate(energy, std::log(v / cut) / std::log(limits.vMax / cut)),
                0.0);
    } else
    {
        return parametrization_->DifferentialCrossSection(energy, v);
    }
}

//----------------------------------------------------------------------------//
void MupairInterpolant::InitializeDifferentialCrossSectionInterpolation(const InterpolationDef &def) {
    std::vector<Interpolant2DBuilder> builder2d(this->components_.size());
    Helper::InterpolantBuilderContainer builder_container2d(this->components_.size());

    for (unsigned int i = 0; i < this->components_.size(); ++i)
    {
        builder2d[i]
                .SetMax1(def.nodes_cross_section)
                .SetX1Min(parametrization_->GetParticleMass())
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
                .SetRationalY(false)
                .SetRelativeY(false)
                .SetLogSubst(false)
                .SetFunction2D(std::bind(&MupairInterpolant::FunctionToBuildDiffCrossSectionInterpolant, this, std::placeholders::_1, std::placeholders::_2, i));

        builder_container2d[i]  = &builder2d[i];
    }

    interpolant_ = Helper::InitializeInterpolation("Mupair", builder_container2d, this->GetHash(), def);
}