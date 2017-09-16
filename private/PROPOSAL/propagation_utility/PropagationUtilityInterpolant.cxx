
#include <boost/bind.hpp>

#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/crossection/IonizInterpolant.h"
#include "PROPOSAL/crossection/EpairInterpolant.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"
#include "PROPOSAL/crossection/parametrization/EpairProduction.h"

#include "PROPOSAL/math/InterpolantBuilder.h"

using namespace PROPOSAL;
using namespace std;

PropagationUtilityInterpolant::PropagationUtilityInterpolant(const ParticleDef& particle_def)
    : PropagationUtility(particle_def)
    , up_(false)
    , big_low_interatction_(0)
    , big_low_decay_(0)
    , store_dif_interaction_(0)
    , store_dif_decay_(0)
    , interpolant_(NULL)
    , interpolant_diff_(NULL)
    , interpol_time_particle_(NULL)
    , interpol_time_particle_diff_(NULL)
    , interpol_prop_decay_(NULL)
    , interpol_prop_decay_diff_(NULL)
    , interpol_prop_interaction_(NULL)
    , interpol_prop_interaction_diff_(NULL)
    , dE2de_interpolant_(NULL)
    , dE2de_interpolant_diff_(NULL)
{
    InitInterpolation();
}

PropagationUtilityInterpolant::PropagationUtilityInterpolant(const ParticleDef& particle,  const Medium& medium,
                                             const EnergyCutSettings& cut_settings,
                                             const Definition& def)
    :PropagationUtility(particle, medium, cut_settings, def)
    , up_(false)
    , big_low_interatction_(0)
    , big_low_decay_(0)
    , store_dif_interaction_(0)
    , store_dif_decay_(0)
    , interpolant_(NULL)
    , interpolant_diff_(NULL)
    , interpol_time_particle_(NULL)
    , interpol_time_particle_diff_(NULL)
    , interpol_prop_decay_(NULL)
    , interpol_prop_decay_diff_(NULL)
    , interpol_prop_interaction_(NULL)
    , interpol_prop_interaction_diff_(NULL)
    , dE2de_interpolant_(NULL)
    , dE2de_interpolant_diff_(NULL)
{
    Parametrization::Definition param_def;

    param_def.multiplier             = utility_def_.brems_multiplier;
    param_def.path_to_tables         = utility_def_.path_to_tables;
    param_def.raw                    = def.raw;
    param_def.lpm_effect_enabled     = utility_def_.lpm_effect_enabled;
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

    InitInterpolation();
}

PropagationUtilityInterpolant::~PropagationUtilityInterpolant()
{
    delete interpolant_;
    delete interpolant_diff_;
    delete interpol_time_particle_;
    delete interpol_time_particle_diff_;
    delete interpol_prop_decay_;
    delete interpol_prop_decay_diff_;
    delete interpol_prop_interaction_;
    delete interpol_prop_interaction_diff_;
    delete dE2de_interpolant_;
    delete dE2de_interpolant_diff_;
}

PropagationUtilityInterpolant::PropagationUtilityInterpolant(const PropagationUtilityInterpolant& collection)
    : PropagationUtility(collection)
    , up_(collection.up_)
    , big_low_interatction_(collection.big_low_interatction_)
    , big_low_decay_(collection.big_low_decay_)
    , store_dif_interaction_(collection.store_dif_interaction_)
    , store_dif_decay_(collection.store_dif_decay_)
    , interpolant_(NULL)
    , interpolant_diff_(NULL)
    , interpol_time_particle_(NULL)
    , interpol_time_particle_diff_(NULL)
    , interpol_prop_decay_(NULL)
    , interpol_prop_decay_diff_(NULL)
    , interpol_prop_interaction_(NULL)
    , interpol_prop_interaction_diff_(NULL)
    , dE2de_interpolant_(NULL)
    , dE2de_interpolant_diff_(NULL)
{
    if (collection.interpolant_ != NULL)
        interpolant_ = new Interpolant(*collection.interpolant_);
    if (collection.interpolant_diff_ != NULL)
        interpolant_diff_ = new Interpolant(*collection.interpolant_diff_);
    if (collection.interpol_time_particle_ != NULL)
        interpol_time_particle_ = new Interpolant(*collection.interpol_time_particle_);
    if (collection.interpol_time_particle_diff_ != NULL)
        interpol_time_particle_diff_ = new Interpolant(*collection.interpol_time_particle_diff_);
    if (collection.interpol_prop_decay_ != NULL)
        interpol_prop_decay_ = new Interpolant(*collection.interpol_prop_decay_);
    if (collection.interpol_prop_decay_diff_ != NULL)
        interpol_prop_decay_diff_ = new Interpolant(*collection.interpol_prop_decay_diff_);
    if (collection.interpol_prop_interaction_ != NULL)
        interpol_prop_interaction_ = new Interpolant(*collection.interpol_prop_interaction_);
    if (collection.interpol_prop_interaction_diff_ != NULL)
        interpol_prop_interaction_diff_ = new Interpolant(*collection.interpol_prop_interaction_diff_);
    if (collection.dE2de_interpolant_ != NULL)
        dE2de_interpolant_ = new Interpolant(*collection.dE2de_interpolant_);
    if (collection.dE2de_interpolant_diff_ != NULL)
        dE2de_interpolant_diff_ = new Interpolant(*collection.dE2de_interpolant_diff_);
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double PropagationUtilityInterpolant::CalculateDisplacement( double ei, double ef, double dist)
{
    (void)dist;

    if (fabs(ei - ef) > fabs(ei) * HALF_PRECISION)
    {
        double aux;

        ini_ = interpolant_->Interpolate(ei);
        aux  = ini_ - interpolant_->Interpolate(ef);

        if (fabs(aux) > fabs(ini_) * HALF_PRECISION)
        {
            return std::max(aux, 0.0);
        }
    }

    ini_ = 0;

    return std::max((interpolant_diff_->Interpolate((ei + ef) / 2)) * (ef - ei), 0.0);
}

// ------------------------------------------------------------------------- //
double PropagationUtilityInterpolant::CalculateFinalEnergy( double ei, double dist)
{
    if (ini_ != 0)
    {
        double aux;
        aux = interpolant_->FindLimit(ini_ - dist);

        if (fabs(aux) > fabs(ei) * HALF_PRECISION)
        {
            return std::min(std::max(aux, particle_def_.low), ei);
        }
    }

    return std::min(
        std::max(ei + dist / interpolant_diff_->Interpolate(ei + dist / (2 * interpolant_diff_->Interpolate(ei))),
                 particle_def_.low),
        ei);
}

// ------------------------------------------------------------------------- //
double PropagationUtilityInterpolant::CalculateFinalEnergy(
                                                double ei,
                                                double rnd,
                                                bool particle_interaction)
{
    if( particle_interaction )
    {
        if( abs(rnd) > abs(store_dif_interaction_)*HALF_PRECISION)
        {
            double aux;

            if(up_&&particle_interaction)
            {
                if(particle_interaction)
                {
                    aux =   interpol_prop_interaction_->FindLimit(store_dif_interaction_-rnd);
                }
                else
                {
                    aux =   interpol_prop_decay_->FindLimit(store_dif_decay_-rnd);
                }
            }
            else
            {
                if(particle_interaction)
                {
                    aux =   interpol_prop_interaction_->FindLimit(store_dif_interaction_+rnd);
                }
                else
                {
                    aux =   interpol_prop_decay_->FindLimit(store_dif_decay_+rnd);
                }
            }
            if(abs(ei-aux) > abs(ei)*HALF_PRECISION)
            {
                return std::min(std::max(aux, particle_def_.low), ei);
            }
        }
    }
    else
    {
        if(abs(rnd) > abs(store_dif_decay_)*HALF_PRECISION)
        {
            double aux;

            if(up_&&particle_interaction)
            {
                if(particle_interaction)
                {
                    aux =   interpol_prop_interaction_->FindLimit(store_dif_interaction_-rnd);
                }
                else
                {
                    aux =   interpol_prop_decay_->FindLimit(store_dif_decay_-rnd);
                }
            }
            else
            {
                if(particle_interaction)
                {
                    aux =   interpol_prop_interaction_->FindLimit(store_dif_interaction_+rnd);
                }
                else
                {
                    aux =   interpol_prop_decay_->FindLimit(store_dif_decay_+rnd);
                }
            }

            if(abs(ei-aux) > abs(ei)*HALF_PRECISION)
            {
                return std::min(std::max(aux, particle_def_.low), ei);
            }
        }
    }

    if(particle_interaction)
    {
        return std::min(std::max(ei + rnd/interpol_prop_interaction_diff_->Interpolate(ei + rnd/(2*interpol_prop_interaction_diff_->Interpolate(ei))), particle_def_.low), ei);
    }
    else
    {
        return std::min(std::max(ei + rnd/interpol_prop_decay_diff_->Interpolate(ei + rnd/(2*interpol_prop_decay_diff_->Interpolate(ei))), particle_def_.low), ei);
    }
}

// ------------------------------------------------------------------------- //
double PropagationUtilityInterpolant::CalculateTrackingIntegal(
                                                    double initial_energy,
                                                    double rnd,
                                                    bool particle_interaction)
{
    (void) rnd;

    if(particle_interaction)
    {
        store_dif_interaction_ =   interpol_prop_interaction_->Interpolate(initial_energy);
    }
    else
    {
        store_dif_decay_ =   interpol_prop_decay_->Interpolate(initial_energy);
    }

    if(up_&&particle_interaction)
    {
        if(particle_interaction)
        {
            return std::max(store_dif_interaction_, 0.0);
        }
        else
        {
            return std::max(store_dif_decay_, 0.0);
        }
    }
    else
    {
        if(particle_interaction)
        {
            return std::max(big_low_interatction_-store_dif_interaction_, 0.0);
        }
        else
        {
            return std::max(big_low_decay_-store_dif_decay_, 0.0);
        }
    }
}

// ------------------------------------------------------------------------- //
double PropagationUtilityInterpolant::CalculateParticleTime( double ei, double ef)
{
    if (abs(ei - ef) > abs(ei) * HALF_PRECISION)
    {
        double aux  = interpol_time_particle_->Interpolate(ei);
        double aux2 = aux - interpol_time_particle_->Interpolate(ef);

        if (abs(aux2) > abs(aux) * HALF_PRECISION)
        {
            return aux2;
        }
    }

    return interpol_time_particle_diff_->Interpolate((ei + ef) / 2) * (ef - ei);
}

// ------------------------------------------------------------------------- //
double PropagationUtilityInterpolant::CalculateDE2de(double ei, double ef)
{
    if (abs(ei - ef) > abs(ei) * HALF_PRECISION)
    {
        double aux  = dE2de_interpolant_->Interpolate(ei);
        double aux2 = aux - dE2de_interpolant_->Interpolate(ef);

        if (abs(aux2) > abs(aux) * HALF_PRECISION)
        {
            return max(aux2, 0.0);
        }
    }

    else
    {
        return max(dE2de_interpolant_diff_->Interpolate((ei + ef) / 2) * (ef - ei), 0.0);
    }
}

// ------------------------------------------------------------------------- //
// Helper functions to init interpolation
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double PropagationUtilityInterpolant::FunctionToBuildInterpolant( double energy)
{
    Integral integral(IROMB, IMAXS, IPREC2);

    return integral.Integrate(energy, particle_def_.low, boost::bind(&PropagationUtilityInterpolant::FunctionToIntegral, this,  _1),4);
}

// ------------------------------------------------------------------------- //
double PropagationUtilityInterpolant::InterpolPropDecay( double energy)
{
    Integral integral(IROMB, IMAXS, IPREC2);

    return -integral.Integrate(energy, BIGENERGY, boost::bind(&PropagationUtilityInterpolant::FunctionToPropIntegralDecay, this,  _1),4);
}


// ------------------------------------------------------------------------- //
double PropagationUtilityInterpolant::InterpolPropInteraction( double energy)
{
    Integral integral(IROMB, IMAXS, IPREC2);

    if(up_)
    {
        return integral.Integrate(energy, particle_def_.low, boost::bind(&PropagationUtilityInterpolant::FunctionToPropIntegralInteraction, this,  _1),4);
    }
    else
    {
        return -integral.Integrate(energy, BIGENERGY, boost::bind(&PropagationUtilityInterpolant::FunctionToPropIntegralInteraction, this,  _1),4);
    }
}

// ------------------------------------------------------------------------- //
double PropagationUtilityInterpolant::InterpolTimeParticleDiff( double energy)
{
    Integral integral(IROMB, IMAXS, IPREC2);

    return integral.Integrate(energy, particle_def_.low, boost::bind(&PropagationUtilityInterpolant::FunctionToTimeIntegral, this,  _1),4);
}

// ------------------------------------------------------------------------- //
double PropagationUtilityInterpolant::FunctionToBuildDE2deInterpolant(double energy, Integral& integral)
{
    return integral.Integrate(energy, particle_def_.low, boost::bind(&PropagationUtilityInterpolant::FunctionToDE2deIntegral, this, _1) , 4);
}

// ------------------------------------------------------------------------- //
// Init methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
void PropagationUtilityInterpolant::InitInterpolation()
{
    Integral integral(IROMB, IMAXS, IPREC2);

    double a = abs(-integral.Integrate(particle_def_.low, particle_def_.low*10, boost::bind(&PropagationUtilityInterpolant::FunctionToPropIntegralInteraction, this,  _1),4));
    double b = abs(-integral.Integrate(BIGENERGY, BIGENERGY/10, boost::bind(&PropagationUtilityInterpolant::FunctionToPropIntegralInteraction, this,  _1),4));

    if (a < b)
    {
        up_ = true;
    }
    else
    {
        up_ = false;
    }

    std::vector<std::pair<Interpolant**, boost::function<double(double)> > > interpolants;

    interpolants.push_back(std::make_pair(&interpolant_, boost::bind(&PropagationUtilityInterpolant::FunctionToBuildInterpolant, this, _1)));
    interpolants.push_back(std::make_pair(&interpolant_diff_, boost::bind(&PropagationUtilityInterpolant::FunctionToIntegral, this, _1)));
    interpolants.push_back(std::make_pair(&interpol_prop_decay_, boost::bind(&PropagationUtilityInterpolant::InterpolPropDecay, this, _1)));
    interpolants.push_back(std::make_pair(&interpol_prop_decay_diff_, boost::bind(&PropagationUtilityInterpolant::FunctionToPropIntegralDecay, this, _1)));
    interpolants.push_back(std::make_pair(&interpol_prop_interaction_, boost::bind(&PropagationUtilityInterpolant::InterpolPropInteraction, this, _1)));
    interpolants.push_back(std::make_pair(&interpol_prop_interaction_diff_, boost::bind(&PropagationUtilityInterpolant::FunctionToPropIntegralInteraction, this, _1)));

    unsigned int number_of_interpolants = interpolants.size();

    std::vector<Interpolant1DBuilder> builder_vec(number_of_interpolants);

    Helper::InterpolantBuilderContainer builder_container(number_of_interpolants);

    for (unsigned int i = 0; i < number_of_interpolants; ++i)
    {
        builder_vec[i]
            .SetMax(NUM3)
            .SetXMin(particle_def_.low)
            .SetXMax(BIGENERGY)
            .SetRomberg(utility_def_.order_of_interpolation)
            .SetRational(false)
            .SetRelative(false)
            .SetIsLog(true)
            .SetRombergY(utility_def_.order_of_interpolation)
            .SetRationalY(false)
            .SetRelativeY(false)
            .SetLogSubst(false)
            .SetFunction1D(interpolants[i].second);

        builder_container[i] = std::make_pair(&builder_vec[i], interpolants[i].first);
    }

    std::vector<Parametrization*> params(crosssections_.size(), NULL);
    for (unsigned int i = 0; i < crosssections_.size(); ++i)
    {
        params[i] = &crosssections_[i]->GetParametrization();
    }

    Helper::InitializeInterpolation("utility",
                                    builder_container,
                                    params);

    if(utility_def_.do_exact_time_calculation)
    {
        InitTimeInterpolation();
    }

    if(utility_def_.do_continuous_randomization)
    {
        InitRandomizeInterpolation();
    }

    big_low_decay_        = interpol_prop_decay_->Interpolate(particle_def_.low);
    big_low_interatction_ = interpol_prop_interaction_->Interpolate(particle_def_.low);
}

// ------------------------------------------------------------------------- //
void PropagationUtilityInterpolant::InitTimeInterpolation()
{
    std::vector<std::pair<Interpolant**, boost::function<double(double)> > > interpolants;

    interpolants.push_back(std::make_pair(&interpol_time_particle_, boost::bind(&PropagationUtilityInterpolant::FunctionToTimeIntegral, this, _1)));
    interpolants.push_back(std::make_pair(&interpol_time_particle_diff_, boost::bind(&PropagationUtilityInterpolant::InterpolTimeParticleDiff, this, _1)));

    unsigned int number_of_interpolants = interpolants.size();

    std::vector<Interpolant1DBuilder> builder_vec(number_of_interpolants);

    Helper::InterpolantBuilderContainer builder_container(number_of_interpolants);

    for (unsigned int i = 0; i < number_of_interpolants; ++i)
    {
        builder_vec[i].SetMax(NUM3)
            .SetXMin(particle_def_.low)
            .SetXMax(BIGENERGY)
            .SetRomberg(utility_def_.order_of_interpolation)
            .SetRational(false)
            .SetRelative(false)
            .SetIsLog(true)
            .SetRombergY(utility_def_.order_of_interpolation)
            .SetRationalY(false)
            .SetRelativeY(false)
            .SetLogSubst(false)
            .SetFunction1D(interpolants[i].second);

        builder_container[i] = std::make_pair(&builder_vec[i], interpolants[i].first);
    }

    std::vector<Parametrization*> params(crosssections_.size(), NULL);
    for (unsigned int i = 0; i < crosssections_.size(); ++i)
    {
        params[i] = &crosssections_[i]->GetParametrization();
    }

    Helper::InitializeInterpolation("time",
                                    builder_container,
                                    params);
}

// ------------------------------------------------------------------------- //
void PropagationUtilityInterpolant::InitRandomizeInterpolation()
{
    std::vector<std::pair<Interpolant**, boost::function<double(double)> > > interpolants;

    Integral integral(IROMB, IMAXS, IPREC2);

    interpolants.push_back(
        std::make_pair(&dE2de_interpolant_, boost::bind(&PropagationUtilityInterpolant::FunctionToBuildDE2deInterpolant, this, _1, boost::ref(integral))));

    interpolants.push_back(std::make_pair(
        &dE2de_interpolant_diff_,
        boost::bind(&PropagationUtilityInterpolant::FunctionToDE2deIntegral, this, _1)));

    unsigned int number_of_interpolants = interpolants.size();

    std::vector<Interpolant1DBuilder> builder_vec(number_of_interpolants);
    Helper::InterpolantBuilderContainer builder_container(number_of_interpolants);


    for (unsigned int i = 0; i < number_of_interpolants; ++i)
    {
        builder_vec[i].SetMax(NUM2)
            .SetXMin(particle_def_.low)
            .SetXMax(BIGENERGY)
            .SetRomberg(utility_def_.order_of_interpolation)
            .SetRational(false)
            .SetRelative(false)
            .SetIsLog(true)
            .SetRombergY(utility_def_.order_of_interpolation)
            .SetRationalY(false)
            .SetRelativeY(false)
            .SetLogSubst(false)
            .SetFunction1D(interpolants[i].second);

        builder_container[i] = std::make_pair(&builder_vec[i], interpolants[i].first);
    }

    std::vector<Parametrization*> params(crosssections_.size(), NULL);
    for (unsigned int i = 0; i < crosssections_.size(); ++i)
    {
        params[i] = &crosssections_[i]->GetParametrization();
    }

    Helper::InitializeInterpolation("cont",
                                    builder_container,
                                    params);
}
