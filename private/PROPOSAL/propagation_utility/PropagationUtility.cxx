#include <boost/bind.hpp>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/math/MathModel.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/propagation_utility/PropagationUtility.h"

#include "PROPOSAL/crossection/CrossSection.h"

#include "PROPOSAL/crossection/IonizIntegral.h"
#include "PROPOSAL/crossection/EpairIntegral.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"
#include "PROPOSAL/crossection/parametrization/EpairProduction.h"

#include "PROPOSAL/crossection/IonizInterpolant.h"
#include "PROPOSAL/crossection/EpairInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"
#include "PROPOSAL/crossection/parametrization/EpairProduction.h"

#include "PROPOSAL/methods.h"

using namespace std;
using namespace PROPOSAL;

/******************************************************************************
*                            Propagation utility                              *
******************************************************************************/

Utility::Definition::Definition()
    : do_interpolation(true)
    , brems_multiplier(1.0)
    , photo_multiplier(1.0)
    , ioniz_multiplier(1.0)
    , epair_multiplier(1.0)
    , photo_parametrization(PhotonuclearFactory::AbramowiczLevinLevyMaor97)
    , photo_shadow(PhotonuclearFactory::ShadowButkevichMikhailov)
    , hardbb_enabled(true)
    , brems_parametrization(BremsstrahlungFactory::KelnerKokoulinPetrukhin)
    , lpm_effect_enabled(false)
    , order_of_interpolation(5)
    , raw(true)
    , path_to_tables("")
{
}

Utility::Definition::~Definition()
{
}

// ------------------------------------------------------------------------- //
// Constructors
// ------------------------------------------------------------------------- //

// Standard constructor
Utility::Utility(const ParticleDef& particle_def)
      //TODO(mario): init different Fri 2017/09/01
    : utility_def_()
    , particle_def_(particle_def)
    , medium_(new Water())
    , cut_settings_()
    , crosssections_()
{
}

Utility::Utility(const ParticleDef& particle_def, const Medium& medium,
                       const EnergyCutSettings& cut_settings,
                       const Definition& def)
    : utility_def_(def)
    , particle_def_(particle_def)
    , medium_(medium.clone())
    , cut_settings_(cut_settings)
    , crosssections_()
{
    Parametrization::Definition param_def;

    param_def.multiplier         = utility_def_.brems_multiplier;
    param_def.lpm_effect_enabled = utility_def_.lpm_effect_enabled;

    if (utility_def_.do_interpolation)
    {
        param_def.path_to_tables         = utility_def_.path_to_tables;
        param_def.raw                    = def.raw;
        param_def.order_of_interpolation = utility_def_.order_of_interpolation;

        crosssections_.push_back(BremsstrahlungFactory::Get().CreateBremsstrahlung(
            utility_def_.brems_parametrization, particle_def_, *medium_, cut_settings_, param_def, true));

        param_def.multiplier = utility_def_.photo_multiplier;

        crosssections_.push_back(PhotonuclearFactory::Get().CreatePhotonuclear(
            utility_def_.photo_parametrization,
            particle_def_,
            *medium_,
            cut_settings_,
            utility_def_.photo_shadow,
            utility_def_.hardbb_enabled,
            param_def,
            true));

        param_def.multiplier = utility_def_.ioniz_multiplier;

        crosssections_.push_back(
            new IonizInterpolant(Ionization(particle_def_, *medium_, cut_settings_, param_def)));

        param_def.multiplier = utility_def_.epair_multiplier;

        crosssections_.push_back(new EpairInterpolant(
            EpairProductionRhoInterpolant(particle_def_, *medium_, cut_settings_, param_def)));
    }
    else
    {
        crosssections_.push_back(BremsstrahlungFactory::Get().CreateBremsstrahlung(
            utility_def_.brems_parametrization, particle_def_, *medium_, cut_settings_, param_def, false));

        param_def.multiplier = def.photo_multiplier;

        crosssections_.push_back(PhotonuclearFactory::Get().CreatePhotonuclear(utility_def_.photo_parametrization,
                                                                               particle_def_,
                                                                               *medium_,
                                                                               cut_settings_,
                                                                               utility_def_.photo_shadow,
                                                                               utility_def_.hardbb_enabled,
                                                                               param_def,
                                                                               false));

        param_def.multiplier = utility_def_.ioniz_multiplier;

        crosssections_.push_back(new IonizIntegral(Ionization(particle_def_, *medium_, cut_settings_, param_def)));

        param_def.multiplier = utility_def_.epair_multiplier;

        crosssections_.push_back(
            new EpairIntegral(EpairProductionRhoIntegral(particle_def_, *medium_, cut_settings_, param_def)));
    }
}

Utility::Utility(const Utility& collection)
    :utility_def_(collection.utility_def_)
    ,particle_def_(collection.particle_def_)
    ,medium_(collection.medium_->clone())
    ,cut_settings_(collection.cut_settings_)
{
    crosssections_.resize(collection.crosssections_.size());

    for (unsigned int i = 0; i < crosssections_.size(); ++i)
    {
        crosssections_[i] = collection.crosssections_[i]->clone();
    }
}

Utility::~Utility()
{
    delete medium_;

    for(std::vector<CrossSection*>::const_iterator iter = crosssections_.begin(); iter != crosssections_.end(); ++iter)
    {
        delete *iter;
    }

    crosssections_.clear();
}

/******************************************************************************
*                            Utility Decorator                            *
******************************************************************************/

UtilityDecorator::UtilityDecorator(const Utility& utility)
    : utility_(utility)
{
}

UtilityDecorator::UtilityDecorator(const UtilityDecorator& decorator)
    : utility_(decorator.utility_)
{

}

UtilityDecorator::~UtilityDecorator()
{
}

// ------------------------------------------------------------------------- //
double UtilityDecorator::FunctionToIntegral( double energy)
{
    double result = 0.0;

    const std::vector<CrossSection*> crosssections = utility_.GetCrosssections();

    for(std::vector<CrossSection*>::const_iterator iter = crosssections.begin(); iter != crosssections.end(); ++iter)
    {
        result += (*iter)->CalculatedEdx(energy);
    }

    return -1.0 / result;
}
