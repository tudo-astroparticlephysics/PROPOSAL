#include <boost/bind.hpp>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/math/MathModel.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/propagation_utility/PropagationUtility.h"

#include "PROPOSAL/crossection/CrossSection.h"

#include "PROPOSAL/methods.h"

using namespace std;
using namespace PROPOSAL;

/******************************************************************************
*                                 Sector                                 *
******************************************************************************/

PropagationUtility::Definition::Definition()
    : brems_multiplier(1.0)
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

PropagationUtility::Definition::~Definition()
{
}

// ------------------------------------------------------------------------- //
// Constructors
// ------------------------------------------------------------------------- //

// Standard constructor
PropagationUtility::PropagationUtility(const ParticleDef& particle_def)
    : ini_(0)
      //TODO(mario): init different Fri 2017/09/01
    , utility_def_()
    , particle_def_(particle_def)
    , medium_(new Water())
    , cut_settings_()
    , crosssections_()
{
    // crosssections_.push_back(new Bremsstrahlung(particle_, medium_, &cut_settings_));
    // crosssections_.push_back(new Epairproduction(particle_, medium_, &cut_settings_));
    // crosssections_.push_back(new Photonuclear(particle_, medium_, &cut_settings_));
    // crosssections_.push_back(new Ionization(particle_, medium_, &cut_settings_));

    //TODO(mario): Polymorphic initilaization in collections childs  Sun 2017/08/27
    // if (utility_def_.do_continuous_randomization)
    // {
    //     randomizer_ = new ContinuousRandomization(particle_);
    // }

}

PropagationUtility::PropagationUtility(const ParticleDef& particle_def, const Medium& medium,
                       const EnergyCutSettings& cut_settings,
                       const Definition& def)
    : ini_(0)
    , utility_def_(def)
    , particle_def_(particle_def)
    , medium_(medium.clone())
    , cut_settings_(cut_settings)
    , crosssections_()
{
}

PropagationUtility::PropagationUtility(const PropagationUtility& collection)
    :ini_(collection.ini_)
    ,utility_def_(collection.utility_def_)
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

PropagationUtility::~PropagationUtility()
{
    delete medium_;
}


// ------------------------------------------------------------------------- //
// Integral functions
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double PropagationUtility::FunctionToIntegral( double energy)
{
    double result = 0.0;

    for(std::vector<CrossSection*>::iterator iter = crosssections_.begin(); iter != crosssections_.end(); ++iter)
    {
        result += (*iter)->CalculatedEdx(energy);

        //TODO(mario): Redo Fri 2017/09/15
        // log_debug("energy %f , dE/dx = %f", energy, aux);
    }

    return -1.0 / result;
}

// ------------------------------------------------------------------------- //
double PropagationUtility::FunctionToPropIntegralDecay( double energy)
{
    log_debug(" + %f", energy);

    return FunctionToIntegral(energy) * MakeDecay(energy);
}

// ------------------------------------------------------------------------- //
double PropagationUtility::FunctionToPropIntegralInteraction( double energy)
{
    double total_rate = 0.0;

    for(std::vector<CrossSection*>::iterator iter = crosssections_.begin(); iter != crosssections_.end(); ++iter)
    {
        //TODO(mario): name Wed 2017/09/06
        // log_debug("Rate for %s = %f", crosssections_.at(i)->GetName().c_str(), rate);

        total_rate += (*iter)->CalculatedNdx(energy);
    }

    return FunctionToIntegral(energy) * total_rate;
}
