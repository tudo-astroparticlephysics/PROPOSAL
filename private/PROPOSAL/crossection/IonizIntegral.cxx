
#include <functional>

#include <cmath>

#include "PROPOSAL/crossection/IonizIntegral.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"

#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

IonizIntegral::IonizIntegral(const Ionization& param, std::shared_ptr<EnergyCutSettings> cuts)
    : CrossSectionIntegral(param, cuts)
{
}

IonizIntegral::IonizIntegral(const IonizIntegral& brems)
    : CrossSectionIntegral(brems)
{
}

IonizIntegral::~IonizIntegral() {}

// ----------------------------------------------------------------- //
// Public methods
// ----------------------------------------------------------------- //

double IonizIntegral::CalculatedEdxWithoutMultiplier(double energy)
{
    Parametrization::KinematicLimits limits = parametrization_->GetKinematicLimits(energy);

    if(parametrization_->GetName() != "IonizBergerSeltzerBhabha" && parametrization_->GetName() != "IonizBergerSeltzerMoller"){
        return energy * dedx_integral_.Integrate(
                limits.vMin,
                cuts_.GetCut(energy),
                std::bind(&Parametrization::FunctionToDEdxIntegral, parametrization_, energy, std::placeholders::_1),
                4);
    }
    else{
        return parametrization_->FunctionToDEdxIntegral(energy, 0);
    }

}

double IonizIntegral::CalculatedEdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return parametrization_->GetMultiplier() * IonizIntegral::CalculatedEdxWithoutMultiplier(energy);
}

// ------------------------------------------------------------------------- //
double IonizIntegral::CalculatedE2dx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    return parametrization_->GetMultiplier() * IonizIntegral::CalculatedE2dxWithoutMultiplier(energy);
}

double IonizIntegral::CalculatedE2dxWithoutMultiplier(double energy)
{
    Parametrization::KinematicLimits limits = parametrization_->GetKinematicLimits(energy);

    return de2dx_integral_.Integrate(
        limits.vMin,
        cuts_.GetCut(energy),
        std::bind(&Parametrization::FunctionToDE2dxIntegral, parametrization_, energy, std::placeholders::_1),
        2);
}

// ------------------------------------------------------------------------- //
double IonizIntegral::CalculatedNdx(double energy)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    Parametrization::KinematicLimits limits = parametrization_->GetKinematicLimits(energy);
    ;
    sum_of_rates_ =
        dndx_integral_[0].Integrate(cuts_.GetCut(energy),
                                    limits.vMax,
                                    std::bind(&Parametrization::FunctionToDNdxIntegral, parametrization_, energy, std::placeholders::_1),
                                    3,
                                    1);

    return parametrization_->GetMultiplier() * sum_of_rates_;
}

// ------------------------------------------------------------------------- //
double IonizIntegral::CalculatedNdx(double energy, double rnd)
{
    if (parametrization_->GetMultiplier() <= 0)
    {
        return 0;
    }

    rnd_ = rnd;

    Parametrization::KinematicLimits limits = parametrization_->GetKinematicLimits(energy);
    ;
    sum_of_rates_ = dndx_integral_[0].IntegrateWithRandomRatio(
        cuts_.GetCut(energy),
        limits.vMax,
        std::bind(&Parametrization::FunctionToDNdxIntegral, parametrization_, energy, std::placeholders::_1),
        3,
        rnd,
        1);

    return parametrization_->GetMultiplier() * sum_of_rates_;
}

// ------------------------------------------------------------------------- //
double IonizIntegral::CalculateStochasticLoss(double energy, double rnd1)
{
    double rnd, rsum;
    const Medium& medium = parametrization_->GetMedium();

    rnd  = medium.GetSumCharge() * rnd1;
    rsum = 0;

    for (unsigned int i = 0; i < components_.size(); i++)
    {
        rsum += components_[i]->GetAtomInMolecule() * components_[i]->GetNucCharge();

        if (rsum > rnd)
        {
            return energy * dndx_integral_[0].GetUpperLimit();
        }
    }

    log_fatal("SumCharge of medium was not initialized correctly");

    return 0;
}
