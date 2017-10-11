
#include <boost/bind.hpp>

#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"
#include "PROPOSAL/math/Interpolant.h"
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

/******************************************************************************
*                              Utility Integral                              *
******************************************************************************/

UtilityInterpolant::UtilityInterpolant(const Utility& utility, InterpolationDef def)
    : UtilityDecorator(utility)
    , stored_result_(0)
    , interpolant_(NULL)
    , interpolant_diff_(NULL)
    , interpolation_def_(def)
{
}

UtilityInterpolant::UtilityInterpolant(const UtilityInterpolant& collection)
    : UtilityDecorator(collection)
    , stored_result_(collection.stored_result_)
    , interpolant_(new Interpolant(*collection.interpolant_))
    , interpolant_diff_(new Interpolant(*collection.interpolant_diff_))
    , interpolation_def_(collection.interpolation_def_)
{
}

UtilityInterpolant::~UtilityInterpolant()
{
    delete interpolant_;
    delete interpolant_diff_;
}

// ------------------------------------------------------------------------- //
double UtilityInterpolant::GetUpperLimit(double ei, double rnd)
{
    return std::min( std::max(ei + rnd / interpolant_diff_->Interpolate(ei + rnd / (2 * interpolant_diff_->Interpolate(ei))), utility_.GetParticleDef().low), ei);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolant::InitInterpolation(const std::string& name, UtilityIntegral& utility, int number_of_sampling_points)
{
    Integral integral(IROMB, IMAXS, IPREC2);
    const ParticleDef& particle_def = utility_.GetParticleDef();

    std::vector<std::pair<Interpolant**, boost::function<double(double)> > > interpolants;

    interpolants.push_back(std::make_pair(&interpolant_, boost::bind(&UtilityInterpolant::BuildInterpolant, this, _1, boost::ref(utility), boost::ref(integral))));
    interpolants.push_back(std::make_pair(&interpolant_diff_, boost::bind(&UtilityIntegral::FunctionToIntegral, &utility, _1)));

    unsigned int number_of_interpolants = interpolants.size();

    std::vector<Interpolant1DBuilder> builder_vec(number_of_interpolants);

    Helper::InterpolantBuilderContainer builder_container(number_of_interpolants);

    for (unsigned int i = 0; i < number_of_interpolants; ++i)
    {
        builder_vec[i]
            .SetMax(number_of_sampling_points)
            .SetXMin(particle_def.low)
            .SetXMax(BIGENERGY)
            .SetRomberg(interpolation_def_.order_of_interpolation)
            .SetRational(false)
            .SetRelative(false)
            .SetIsLog(true)
            .SetRombergY(interpolation_def_.order_of_interpolation)
            .SetRationalY(false)
            .SetRelativeY(false)
            .SetLogSubst(false)
            .SetFunction1D(interpolants[i].second);

        builder_container[i] = std::make_pair(&builder_vec[i], interpolants[i].first);
    }

    std::vector<CrossSection*> crosssections = utility_.GetCrosssections();
    std::vector<Parametrization*> params(crosssections.size(), NULL);

    for (unsigned int i = 0; i < crosssections.size(); ++i)
    {
        params[i] = &crosssections[i]->GetParametrization();
    }

    Helper::InitializeInterpolation(name, builder_container, params, interpolation_def_);
}

/******************************************************************************
*                            Utility Displacement                            *
******************************************************************************/

UtilityInterpolantDisplacement::UtilityInterpolantDisplacement(const Utility& utility, InterpolationDef def)
    : UtilityInterpolant(utility, def)
{
    UtilityIntegralDisplacement utility_disp(utility_);
    InitInterpolation("displacement", utility_disp, NUM3);
}

UtilityInterpolantDisplacement::UtilityInterpolantDisplacement(const UtilityInterpolantDisplacement& collection)
    : UtilityInterpolant(collection)
{
}

UtilityInterpolantDisplacement::~UtilityInterpolantDisplacement()
{
}

double UtilityInterpolantDisplacement::Calculate(double ei, double ef, double rnd)
{
    (void) rnd;

    if (fabs(ei - ef) > fabs(ei) * HALF_PRECISION)
    {
        double aux;

        stored_result_ = interpolant_->Interpolate(ei);
        aux  = stored_result_ - interpolant_->Interpolate(ef);

        if (fabs(aux) > fabs(stored_result_) * HALF_PRECISION)
        {
            return std::max(aux, 0.0);
        }
    }

    stored_result_ = 0;

    return std::max((interpolant_diff_->Interpolate((ei + ef) / 2)) * (ef - ei), 0.0);
}

double UtilityInterpolantDisplacement::GetUpperLimit(double ei, double rnd)
{
    if (stored_result_ != 0)
    {
        double aux;
        aux = interpolant_->FindLimit(stored_result_ - rnd);

        if (fabs(aux) > fabs(ei) * HALF_PRECISION)
        {
            return std::min(std::max(aux, utility_.GetParticleDef().low), ei);
        }
    }

    return UtilityInterpolant::GetUpperLimit(ei, rnd);
}

double UtilityInterpolantDisplacement::BuildInterpolant(double energy, UtilityIntegral& utility, Integral& integral)
{
    return integral.Integrate(energy, utility_.GetParticleDef().low, boost::bind(&UtilityIntegral::FunctionToIntegral, &utility,  _1),4);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantDisplacement::InitInterpolation(const std::string& name, UtilityIntegral& utility, int number_of_sampling_points)
{
    UtilityInterpolant::InitInterpolation(name, utility, number_of_sampling_points);
}

/******************************************************************************
*                            Utility Ineraction                            *
******************************************************************************/

UtilityInterpolantInteraction::UtilityInterpolantInteraction(const Utility& utility, InterpolationDef def)
    : UtilityInterpolant(utility, def)
    , big_low_(0)
    , up_(0)
{
    UtilityIntegralInteraction utility_int(utility_);
    InitInterpolation("interaction", utility_int, NUM3);
}

UtilityInterpolantInteraction::UtilityInterpolantInteraction(const UtilityInterpolantInteraction& collection)
    : UtilityInterpolant(collection)
    , big_low_(collection.big_low_)
    , up_(collection.up_)
{
}

UtilityInterpolantInteraction::~UtilityInterpolantInteraction()
{
}

double UtilityInterpolantInteraction::Calculate(double ei, double ef, double rnd)
{
    (void) rnd;
    (void) ef;

    stored_result_ =   interpolant_->Interpolate(ei);

    if (up_)
    {
        return std::max(stored_result_, 0.0);
    }
    else
    {
        return std::max(big_low_ - stored_result_, 0.0);
    }
}

double UtilityInterpolantInteraction::GetUpperLimit(double ei, double rnd)
{
    if( abs(rnd) > abs(stored_result_)*HALF_PRECISION)
    {
        double aux;

        if(up_)
        {
            aux =   interpolant_->FindLimit(stored_result_ - rnd);
        }
        else
        {
            aux = interpolant_->FindLimit(stored_result_ + rnd);
        }

        if(abs(ei-aux) > abs(ei)*HALF_PRECISION)
        {
            return std::min(std::max(aux, utility_.GetParticleDef().low), ei);
        }
    }

    return UtilityInterpolant::GetUpperLimit(ei, rnd);
}

double UtilityInterpolantInteraction::BuildInterpolant(double energy, UtilityIntegral& utility, Integral& integral)
{
    if(up_)
    {
        return integral.Integrate(energy, utility_.GetParticleDef().low, boost::bind(&UtilityIntegral::FunctionToIntegral, &utility,  _1),4);
    }
    else
    {
        return -integral.Integrate(energy, BIGENERGY, boost::bind(&UtilityIntegral::FunctionToIntegral, &utility,  _1),4);
    }
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantInteraction::InitInterpolation(const std::string& name, UtilityIntegral& utility, int number_of_sampling_points)
{
    Integral integral(IROMB, IMAXS, IPREC2);
    const ParticleDef& particle_def = utility_.GetParticleDef();

    double a = abs(-integral.Integrate(particle_def.low, particle_def.low*10, boost::bind(&UtilityIntegral::FunctionToIntegral, &utility,  _1),4));
    double b = abs(-integral.Integrate(BIGENERGY, BIGENERGY/10, boost::bind(&UtilityIntegral::FunctionToIntegral, &utility,  _1),4));

    if (a < b)
    {
        up_ = true;
    }
    else
    {
        up_ = false;
    }

    UtilityInterpolant::InitInterpolation(name, utility, number_of_sampling_points);

    big_low_= interpolant_->Interpolate(particle_def.low);
}

/******************************************************************************
*                            Utility Decay                            *
******************************************************************************/

UtilityInterpolantDecay::UtilityInterpolantDecay(const Utility& utility, InterpolationDef def)
    : UtilityInterpolant(utility, def)
    , big_low_(0)
    , up_(0)
{
    UtilityIntegralDecay utility_decay(utility_);
    InitInterpolation("decay", utility_decay, NUM3);
}

UtilityInterpolantDecay::UtilityInterpolantDecay(const UtilityInterpolantDecay& collection)
    : UtilityInterpolant(collection)
    , big_low_(collection.big_low_)
    , up_(collection.up_)
{
}

UtilityInterpolantDecay::~UtilityInterpolantDecay()
{
}

double UtilityInterpolantDecay::Calculate(double ei, double ef, double rnd)
{
    (void) rnd;
    (void) ef;

    stored_result_ = interpolant_->Interpolate(ei);

    if (up_)
    {
        return std::max(stored_result_, 0.0);
    }
    else
    {
        return std::max(big_low_ - stored_result_, 0.0);
    }
}

double UtilityInterpolantDecay::GetUpperLimit(double ei, double rnd)
{
    if(abs(rnd) > abs(stored_result_)*HALF_PRECISION)
    {
        double aux;

        aux = interpolant_->FindLimit(stored_result_ + rnd);

        if(abs(ei-aux) > abs(ei)*HALF_PRECISION)
        {
            return std::min(std::max(aux, utility_.GetParticleDef().low), ei);
        }
    }

    return UtilityInterpolant::GetUpperLimit(ei, rnd);
}

double UtilityInterpolantDecay::BuildInterpolant(double energy, UtilityIntegral& utility, Integral& integral)
{
    return -integral.Integrate(energy, BIGENERGY, boost::bind(&UtilityIntegral::FunctionToIntegral, &utility,  _1),4);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantDecay::InitInterpolation(const std::string& name, UtilityIntegral& utility, int number_of_sampling_points)
{
    Integral integral(IROMB, IMAXS, IPREC2);
    const ParticleDef& particle_def = utility_.GetParticleDef();

    double a = abs(-integral.Integrate(particle_def.low, particle_def.low*10, boost::bind(&UtilityIntegral::FunctionToIntegral, &utility,  _1),4));
    double b = abs(-integral.Integrate(BIGENERGY, BIGENERGY/10, boost::bind(&UtilityIntegral::FunctionToIntegral, &utility,  _1),4));

    if (a < b)
    {
        up_ = true;
    }
    else
    {
        up_ = false;
    }

    UtilityInterpolant::InitInterpolation(name, utility, number_of_sampling_points);

    big_low_= interpolant_->Interpolate(particle_def.low);
}


/******************************************************************************
*                            Utility Time                            *
******************************************************************************/

UtilityInterpolantTime::UtilityInterpolantTime(const Utility& utility, InterpolationDef def)
    : UtilityInterpolant(utility, def)
{
    UtilityIntegralTime utility_time(utility_);
    InitInterpolation("time", utility_time, NUM3);
}

UtilityInterpolantTime::UtilityInterpolantTime(const UtilityInterpolantTime& collection)
    : UtilityInterpolant(collection)
{
}

UtilityInterpolantTime::~UtilityInterpolantTime()
{
}

double UtilityInterpolantTime::Calculate(double ei, double ef, double rnd)
{
    (void) rnd;

    if (abs(ei - ef) > abs(ei) * HALF_PRECISION)
    {
        double aux  = interpolant_->Interpolate(ei);
        double aux2 = aux - interpolant_->Interpolate(ef);

        if (abs(aux2) > abs(aux) * HALF_PRECISION)
        {
            return aux2;
        }
    }

    return interpolant_diff_->Interpolate((ei + ef) / 2) * (ef - ei);
}

// ------------------------------------------------------------------------- //
double UtilityInterpolantTime::GetUpperLimit(double ei, double rnd)
{
    (void) ei;
    (void) rnd;

    return 0;
}

double UtilityInterpolantTime::BuildInterpolant(double energy, UtilityIntegral& utility, Integral& integral)
{
    return integral.Integrate(energy, utility_.GetParticleDef().low, boost::bind(&UtilityIntegral::FunctionToIntegral, &utility,  _1),4);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantTime::InitInterpolation(const std::string& name, UtilityIntegral& utility, int number_of_sampling_points)
{
    UtilityInterpolant::InitInterpolation(name, utility, number_of_sampling_points);
}

/******************************************************************************
*                            Utility ContRand                            *
******************************************************************************/

UtilityInterpolantContRand::UtilityInterpolantContRand(const Utility& utility, InterpolationDef def)
    : UtilityInterpolant(utility, def)
{
    UtilityIntegralContRand utility_contrand(utility_);
    InitInterpolation("contrand", utility_contrand, NUM2);
}

UtilityInterpolantContRand::UtilityInterpolantContRand(const UtilityInterpolantContRand& collection)
    : UtilityInterpolant(collection)
{
}

UtilityInterpolantContRand::~UtilityInterpolantContRand()
{
}

double UtilityInterpolantContRand::Calculate(double ei, double ef, double rnd)
{
    (void) rnd;

    if (abs(ei - ef) > abs(ei) * HALF_PRECISION)
    {
        double aux  = interpolant_->Interpolate(ei);
        double aux2 = aux - interpolant_->Interpolate(ef);

        if (abs(aux2) > abs(aux) * HALF_PRECISION)
        {
            return max(aux2, 0.0);
        }
    }
    else
    {
        return max(interpolant_diff_->Interpolate((ei + ef) / 2) * (ef - ei), 0.0);
    }
    //TODO(mario): return Thu 2017/09/21
}

// ------------------------------------------------------------------------- //
double UtilityInterpolantContRand::GetUpperLimit(double ei, double rnd)
{
    (void) ei;
    (void) rnd;

    return 0;
}

double UtilityInterpolantContRand::BuildInterpolant(double energy, UtilityIntegral& utility, Integral& integral)
{
    return integral.Integrate(energy, utility_.GetParticleDef().low, boost::bind(&UtilityIntegral::FunctionToIntegral, &utility,  _1),4);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantContRand::InitInterpolation(const std::string& name, UtilityIntegral& utility, int number_of_sampling_points)
{
    UtilityInterpolant::InitInterpolation(name, utility, number_of_sampling_points);
}


/******************************************************************************
*                            Utility Scattering                            *
******************************************************************************/

UtilityInterpolantScattering::UtilityInterpolantScattering(const Utility& utility, InterpolationDef def)
    : UtilityInterpolant(utility, def)
{
    UtilityIntegralScattering utility_scattering(utility_);
    InitInterpolation("scattering", utility_scattering, NUM2);
}

UtilityInterpolantScattering::UtilityInterpolantScattering(const UtilityInterpolantScattering& collection)
    : UtilityInterpolant(collection)
{
}

UtilityInterpolantScattering::~UtilityInterpolantScattering()
{
}

double UtilityInterpolantScattering::Calculate(double ei, double ef, double rnd)
{
    (void) rnd;

    if (fabs(ei - ef) > fabs(ei) * HALF_PRECISION)
    {
        double aux  = interpolant_->Interpolate(ei);
        double aux2 = aux - interpolant_->Interpolate(ef);

        if (fabs(aux2) > fabs(aux) * HALF_PRECISION)
        {
            return aux2;
        }
        else
        {
            return interpolant_diff_->Interpolate((ei + ef) / 2) * (ef - ei);
        }
    }
    else
    {
        return interpolant_diff_->Interpolate((ei + ef) / 2) * (ef - ei);
    }
}

// ------------------------------------------------------------------------- //
double UtilityInterpolantScattering::GetUpperLimit(double ei, double rnd)
{
    (void) ei;
    (void) rnd;

    return 0;
}

double UtilityInterpolantScattering::BuildInterpolant(double energy, UtilityIntegral& utility, Integral& integral)
{
    return integral.Integrate(energy, BIGENERGY, boost::bind(&UtilityIntegral::FunctionToIntegral, &utility, _1), 4);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantScattering::InitInterpolation(const std::string& name, UtilityIntegral& utility, int number_of_sampling_points)
{
    UtilityInterpolant::InitInterpolation(name, utility, number_of_sampling_points);
}
