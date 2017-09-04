
#include "PROPOSAL/crossection/parametrization/Parametrization.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// CrossSection Definition
// ------------------------------------------------------------------------- //


Parametrization::Definition::Definition()
    : multiplier(1.0)
    , lpm_effect_enabled(false)
    , order_of_interpolation(5)
    , path_to_tables("")
    , raw(true)
{
}

Parametrization::Definition::~Definition()
{
}

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

Parametrization::Parametrization(const ParticleDef& particle_def,
                                 const Medium& medium,
                                 const EnergyCutSettings& cuts,
                                 Definition param_def)
    : particle_def_(particle_def)
    , medium_(medium.clone())
    , cut_settings_(cuts)
    , current_component_(medium_->GetComponents().at(0))
    , param_def_(param_def)
    // , lpm_effect_enabled_(false)
    , init_lpm_effect_(false)
{
}

Parametrization::Parametrization(const Parametrization& param)
    : particle_def_(param.particle_def_)
    , medium_(param.medium_->clone())
    , cut_settings_(param.cut_settings_)
    , current_component_(medium_->GetComponents().at(0)) // //TODO(mario): Check better way Mon 2017/09/04
    , param_def_(param.param_def_)
    // , multiplier_(param.multiplier_)
    // , lpm_effect_enabled_(param.lpm_effect_enabled_)
    , init_lpm_effect_(param.init_lpm_effect_)
{
}

Parametrization::~Parametrization()
{
    delete medium_;
}
