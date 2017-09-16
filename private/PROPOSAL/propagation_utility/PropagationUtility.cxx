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
*                            Propagation utility                              *
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
    , do_exact_time_calculation(false)
    , do_continuous_randomization(false)
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

    for(std::vector<CrossSection*>::const_iterator iter = crosssections_.begin(); iter != crosssections_.end(); ++iter)
    {
        delete *iter;
    }

    crosssections_.clear();
}

// ------------------------------------------------------------------------- //
// Public already defined member functions
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double PropagationUtility::MakeDecay(double energy)
{
    // TODO(mario): multiplier? Was not used before Fri 2017/08/25
    // if(multiplier_ <= 0 || particle_.GetLifetime() < 0)
    // {
    //     return 0;
    // }
    //
    // return multiplier_/max((particle_.GetMomentum()/particle_.GetMass())*particle_.GetLifetime()*SPEED, XRES);
    if (particle_def_.lifetime < 0)
    {
        return 0;
    }

    //TODO(mario): Better way? Sat 2017/09/02
    double square_momentum = energy * energy - particle_def_.mass * particle_def_.mass;
    double particle_momentum = sqrt(max(square_momentum, 0.0));

    // return multiplier / max((particle_momentum / particle_.GetMass()) * particle_.GetLifetime() * SPEED, XRES);
    return 1.0 / max((particle_momentum / particle_def_.mass) * particle_def_.lifetime * SPEED, XRES);
}

double PropagationUtility::Randomize(double initial_energy, double final_energy, double rnd)
{
    double sigma, xhi, xlo, rndtmp;

    // this happens if small distances are propagated and the
    // energy loss is so small that it is smaller than the precision
    // which is checked for in during the calculation.
    if (initial_energy == final_energy)
    {
        return final_energy;
    }

    sigma = sqrt(CalculateDE2de(initial_energy, final_energy));

    // It is not drawn from the real gaus distribution but rather from the
    // area which is possible due to the limits of the initial energy and the
    // particle mass. Another possibility would be to draw again but that would be
    // more expensive.
    //
    // calculate the allowed region
    xhi = 0.5 + boost::math::erf((initial_energy - final_energy) / (SQRT2 * sigma)) / 2;
    xlo = 0.5 + boost::math::erf((particle_def_.low - final_energy) / (SQRT2 * sigma)) / 2;

    // draw random number from the allowed region.
    rndtmp = xlo + (xhi - xlo) * rnd;

    // Calculate and return the needed value.
    return SQRT2 * sigma * erfInv(2 * (rndtmp - 0.5)) + final_energy;
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

// ------------------------------------------------------------------------- //
double PropagationUtility::FunctionToTimeIntegral(double energy)
{
    //TODO(mario): Better way? Sat 2017/09/02
    double square_momentum = energy * energy - particle_def_.mass * particle_def_.mass;
    double particle_momentum = sqrt(max(square_momentum, 0.0));

    return energy * SPEED / particle_momentum * FunctionToIntegral(energy);
}

// ------------------------------------------------------------------------- //
double PropagationUtility::FunctionToDE2deIntegral(double energy)
{
    double sum = 0.0;

    for(std::vector<CrossSection*>::iterator iter = crosssections_.begin(); iter != crosssections_.end(); ++iter)
    {
        sum += (*iter)->CalculatedE2dx(energy);
    }

    return FunctionToIntegral(energy) * sum;
}
