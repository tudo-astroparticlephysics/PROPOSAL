
#include <functional>

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"

#include "PROPOSAL/crossection/CrossSection.h"

#include "PROPOSAL/math/InterpolantBuilder.h"

using namespace PROPOSAL;

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

UtilityInterpolant::UtilityInterpolant(const Utility& utility, const UtilityInterpolant& collection)
    : UtilityDecorator(utility)
    , stored_result_(collection.stored_result_)
    , interpolant_(new Interpolant(*collection.interpolant_))
    , interpolant_diff_(new Interpolant(*collection.interpolant_diff_))
    , interpolation_def_(collection.interpolation_def_)
{
    if (utility != collection.GetUtility())
    {
        log_fatal("Utilities of the decorators should have same values!");
    }
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

bool UtilityInterpolant::compare(const UtilityDecorator& utility_decorator) const
{
    const UtilityInterpolant* utility_interpolant = static_cast<const UtilityInterpolant*>(&utility_decorator);

    if (stored_result_ != utility_interpolant->stored_result_)
        return false;
    else if (*interpolant_ != *utility_interpolant->interpolant_)
        return false;
    else if (*interpolant_diff_ != *utility_interpolant->interpolant_diff_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
double UtilityInterpolant::GetUpperLimit(double ei, double rnd)
{
    return std::min(
        std::max(ei + rnd / interpolant_diff_->Interpolate(ei + rnd / (2 * interpolant_diff_->Interpolate(ei))),
                 utility_.GetParticleDef().low),
        ei);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolant::InitInterpolation(const std::string& name,
                                           UtilityIntegral& utility,
                                           int number_of_sampling_points)
{
    Integral integral(IROMB, IMAXS, IPREC2);
    const ParticleDef& particle_def = utility_.GetParticleDef();

    std::vector<std::pair<Interpolant**, std::function<double(double)> > > interpolants;

    interpolants.push_back(std::make_pair(
        &interpolant_,
        std::bind(&UtilityInterpolant::BuildInterpolant, this, std::placeholders::_1, std::ref(utility), std::ref(integral))));
    interpolants.push_back(
        std::make_pair(&interpolant_diff_, std::bind(&UtilityIntegral::FunctionToIntegral, &utility, std::placeholders::_1)));

    unsigned int number_of_interpolants = interpolants.size();

    std::vector<Interpolant1DBuilder> builder_vec(number_of_interpolants);

    Helper::InterpolantBuilderContainer builder_container(number_of_interpolants);

    for (unsigned int i = 0; i < number_of_interpolants; ++i)
    {
        builder_vec[i]
            .SetMax(number_of_sampling_points)
            .SetXMin(particle_def.low)
            .SetXMax(interpolation_def_.max_node_energy)
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
    InitInterpolation("displacement", utility_disp, def.nodes_propagate);
}

UtilityInterpolantDisplacement::UtilityInterpolantDisplacement(const Utility& utility,
                                                               const UtilityInterpolantDisplacement& collection)
    : UtilityInterpolant(utility, collection)
{
}

UtilityInterpolantDisplacement::UtilityInterpolantDisplacement(const UtilityInterpolantDisplacement& collection)
    : UtilityInterpolant(collection)
{
}

UtilityInterpolantDisplacement::~UtilityInterpolantDisplacement() {}

double UtilityInterpolantDisplacement::Calculate(double ei, double ef, double rnd)
{
    (void)rnd;

    if (std::abs(ei - ef) > std::abs(ei) * HALF_PRECISION)
    {
        double aux;

        stored_result_ = interpolant_->Interpolate(ei);
        aux            = stored_result_ - interpolant_->Interpolate(ef);

        if (std::abs(aux) > std::abs(stored_result_) * HALF_PRECISION)
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

        if (std::abs(aux) > std::abs(ei) * HALF_PRECISION)
        {
            return std::min(std::max(aux, utility_.GetParticleDef().low), ei);
        }
    }

    return UtilityInterpolant::GetUpperLimit(ei, rnd);
}

double UtilityInterpolantDisplacement::BuildInterpolant(double energy, UtilityIntegral& utility, Integral& integral)
{
    return integral.Integrate(
        energy, utility_.GetParticleDef().low, std::bind(&UtilityIntegral::FunctionToIntegral, &utility, std::placeholders::_1), 4);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantDisplacement::InitInterpolation(const std::string& name,
                                                       UtilityIntegral& utility,
                                                       int number_of_sampling_points)
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
    InitInterpolation("interaction", utility_int, def.nodes_propagate);
}

UtilityInterpolantInteraction::UtilityInterpolantInteraction(const Utility& utility,
                                                             const UtilityInterpolantInteraction& collection)
    : UtilityInterpolant(utility, collection)
    , big_low_(collection.big_low_)
    , up_(collection.up_)
{
}

UtilityInterpolantInteraction::UtilityInterpolantInteraction(const UtilityInterpolantInteraction& collection)
    : UtilityInterpolant(collection)
    , big_low_(collection.big_low_)
    , up_(collection.up_)
{
}

UtilityInterpolantInteraction::~UtilityInterpolantInteraction() {}

bool UtilityInterpolantInteraction::compare(const UtilityDecorator& utility_decorator) const
{
    const UtilityInterpolantInteraction* utility_interpolant =
        static_cast<const UtilityInterpolantInteraction*>(&utility_decorator);

    if (big_low_ != utility_interpolant->big_low_)
        return false;
    if (up_ != utility_interpolant->up_)
        return false;
    else
        return UtilityInterpolant::compare(utility_decorator);
}

double UtilityInterpolantInteraction::Calculate(double ei, double ef, double rnd)
{
    (void)rnd;
    (void)ef;

    stored_result_ = interpolant_->Interpolate(ei);

    if (up_)
    {
        return std::max(stored_result_, 0.0);
    } else
    {
        return std::max(big_low_ - stored_result_, 0.0);
    }
}

double UtilityInterpolantInteraction::GetUpperLimit(double ei, double rnd)
{
    if (std::abs(rnd) > std::abs(stored_result_) * HALF_PRECISION)
    {
        double aux;

        if (up_)
        {
            aux = interpolant_->FindLimit(stored_result_ - rnd);
        } else
        {
            aux = interpolant_->FindLimit(stored_result_ + rnd);
        }

        if (std::abs(ei - aux) > std::abs(ei) * HALF_PRECISION)
        {
            return std::min(std::max(aux, utility_.GetParticleDef().low), ei);
        }
    }

    return UtilityInterpolant::GetUpperLimit(ei, rnd);
}

double UtilityInterpolantInteraction::BuildInterpolant(double energy, UtilityIntegral& utility, Integral& integral)
{
    if (up_)
    {
        return integral.Integrate(
            energy, utility_.GetParticleDef().low, std::bind(&UtilityIntegral::FunctionToIntegral, &utility, std::placeholders::_1), 4);
    } else
    {
        return -integral.Integrate(
            energy, interpolation_def_.max_node_energy, std::bind(&UtilityIntegral::FunctionToIntegral, &utility, std::placeholders::_1), 4);
    }
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantInteraction::InitInterpolation(const std::string& name,
                                                      UtilityIntegral& utility,
                                                      int number_of_sampling_points)
{
    Integral integral(IROMB, IMAXS, IPREC2);
    const ParticleDef& particle_def = utility_.GetParticleDef();

    double a = std::abs(-integral.Integrate(
        particle_def.low, particle_def.low * 10, std::bind(&UtilityIntegral::FunctionToIntegral, &utility, std::placeholders::_1), 4));
    double b = std::abs(-integral.Integrate(
        interpolation_def_.max_node_energy, interpolation_def_.max_node_energy/ 10, std::bind(&UtilityIntegral::FunctionToIntegral, &utility, std::placeholders::_1), 4));

    if (a < b)
    {
        up_ = true;
    } else
    {
        up_ = false;
    }

    UtilityInterpolant::InitInterpolation(name, utility, number_of_sampling_points);

    big_low_ = interpolant_->Interpolate(particle_def.low);
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
    InitInterpolation("decay", utility_decay, def.nodes_propagate);
}

UtilityInterpolantDecay::UtilityInterpolantDecay(const Utility& utility, const UtilityInterpolantDecay& collection)
    : UtilityInterpolant(utility, collection)
    , big_low_(collection.big_low_)
    , up_(collection.up_)
{
}

UtilityInterpolantDecay::UtilityInterpolantDecay(const UtilityInterpolantDecay& collection)
    : UtilityInterpolant(collection)
    , big_low_(collection.big_low_)
    , up_(collection.up_)
{
}

UtilityInterpolantDecay::~UtilityInterpolantDecay() {}

bool UtilityInterpolantDecay::compare(const UtilityDecorator& utility_decorator) const
{
    const UtilityInterpolantDecay* utility_interpolant =
        static_cast<const UtilityInterpolantDecay*>(&utility_decorator);

    if (big_low_ != utility_interpolant->big_low_)
        return false;
    if (up_ != utility_interpolant->up_)
        return false;
    else
        return UtilityInterpolant::compare(utility_decorator);
}

double UtilityInterpolantDecay::Calculate(double ei, double ef, double rnd)
{
    (void)rnd;
    (void)ef;

    stored_result_ = interpolant_->Interpolate(ei);

    if (up_)
    {
        return std::max(stored_result_, 0.0);
    } else
    {
        return std::max(big_low_ - stored_result_, 0.0);
    }
}

double UtilityInterpolantDecay::GetUpperLimit(double ei, double rnd)
{
    if (std::abs(rnd) > std::abs(stored_result_) * HALF_PRECISION)
    {
        double aux;

        aux = interpolant_->FindLimit(stored_result_ + rnd);

        if (std::abs(ei - aux) > std::abs(ei) * HALF_PRECISION)
        {
            return std::min(std::max(aux, utility_.GetParticleDef().low), ei);
        }
    }

    return UtilityInterpolant::GetUpperLimit(ei, rnd);
}

double UtilityInterpolantDecay::BuildInterpolant(double energy, UtilityIntegral& utility, Integral& integral)
{
    return -integral.Integrate(energy, interpolation_def_.max_node_energy, std::bind(&UtilityIntegral::FunctionToIntegral, &utility, std::placeholders::_1), 4);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantDecay::InitInterpolation(const std::string& name,
                                                UtilityIntegral& utility,
                                                int number_of_sampling_points)
{
    Integral integral(IROMB, IMAXS, IPREC2);
    const ParticleDef& particle_def = utility_.GetParticleDef();

    double a = std::abs(-integral.Integrate(
        particle_def.low, particle_def.low * 10, std::bind(&UtilityIntegral::FunctionToIntegral, &utility, std::placeholders::_1), 4));
    double b = std::abs(-integral.Integrate(
        interpolation_def_.max_node_energy, interpolation_def_.max_node_energy/ 10, std::bind(&UtilityIntegral::FunctionToIntegral, &utility, std::placeholders::_1), 4));

    if (a < b)
    {
        up_ = true;
    } else
    {
        up_ = false;
    }

    UtilityInterpolant::InitInterpolation(name, utility, number_of_sampling_points);

    big_low_ = interpolant_->Interpolate(particle_def.low);
}

/******************************************************************************
 *                            Utility Time                            *
 ******************************************************************************/

UtilityInterpolantTime::UtilityInterpolantTime(const Utility& utility, InterpolationDef def)
    : UtilityInterpolant(utility, def)
{
    UtilityIntegralTime utility_time(utility_);
    InitInterpolation("time", utility_time, def.nodes_propagate);
}

UtilityInterpolantTime::UtilityInterpolantTime(const Utility& utility, const UtilityInterpolantTime& collection)
    : UtilityInterpolant(utility, collection)
{
}

UtilityInterpolantTime::UtilityInterpolantTime(const UtilityInterpolantTime& collection)
    : UtilityInterpolant(collection)
{
}

UtilityInterpolantTime::~UtilityInterpolantTime() {}

double UtilityInterpolantTime::Calculate(double ei, double ef, double rnd)
{
    (void)rnd;

    if (std::abs(ei - ef) > std::abs(ei) * HALF_PRECISION)
    {
        double aux  = interpolant_->Interpolate(ei);
        double aux2 = aux - interpolant_->Interpolate(ef);

        if (std::abs(aux2) > std::abs(aux) * HALF_PRECISION)
        {
            return aux2;
        }
    }

    return interpolant_diff_->Interpolate((ei + ef) / 2) * (ef - ei);
}

// ------------------------------------------------------------------------- //
double UtilityInterpolantTime::GetUpperLimit(double ei, double rnd)
{
    (void)ei;
    (void)rnd;

    return 0;
}

double UtilityInterpolantTime::BuildInterpolant(double energy, UtilityIntegral& utility, Integral& integral)
{
    return integral.Integrate(
        energy, utility_.GetParticleDef().low, std::bind(&UtilityIntegral::FunctionToIntegral, &utility, std::placeholders::_1), 4);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantTime::InitInterpolation(const std::string& name,
                                               UtilityIntegral& utility,
                                               int number_of_sampling_points)
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
    InitInterpolation("contrand", utility_contrand, def.nodes_continous_randomization);
}

UtilityInterpolantContRand::UtilityInterpolantContRand(const Utility& utility,
                                                       const UtilityInterpolantContRand& collection)
    : UtilityInterpolant(utility, collection)
{
}

UtilityInterpolantContRand::UtilityInterpolantContRand(const UtilityInterpolantContRand& collection)
    : UtilityInterpolant(collection)
{
}

UtilityInterpolantContRand::~UtilityInterpolantContRand() {}

double UtilityInterpolantContRand::Calculate(double ei, double ef, double rnd)
{
    if (std::abs(ei - ef) > std::abs(ei) * HALF_PRECISION)
    {
        double aux  = interpolant_->Interpolate(ei);
        double aux2 = aux - interpolant_->Interpolate(ef);

        if (std::abs(aux2) > std::abs(aux) * HALF_PRECISION)
            return std::max(aux2, 0.0);
    } else
    {
        return std::max(interpolant_diff_->Interpolate((ei + ef) / 2) * (ef - ei), 0.0);
    }

    // If the previous conditions do not hold, create a temporary integral.
    // Its ok, because you almost never end up here
    return UtilityIntegralContRand(utility_).Calculate(ei, ef, rnd);
}

// ------------------------------------------------------------------------- //
double UtilityInterpolantContRand::GetUpperLimit(double ei, double rnd)
{
    (void)ei;
    (void)rnd;

    return 0;
}

double UtilityInterpolantContRand::BuildInterpolant(double energy, UtilityIntegral& utility, Integral& integral)
{
    return integral.Integrate(
        energy, utility_.GetParticleDef().low, std::bind(&UtilityIntegral::FunctionToIntegral, &utility, std::placeholders::_1), 4);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantContRand::InitInterpolation(const std::string& name,
                                                   UtilityIntegral& utility,
                                                   int number_of_sampling_points)
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
    InitInterpolation("scattering", utility_scattering, def.nodes_continous_randomization);
}

UtilityInterpolantScattering::UtilityInterpolantScattering(const Utility& utility,
                                                           const UtilityInterpolantScattering& collection)
    : UtilityInterpolant(utility, collection)
{
}

UtilityInterpolantScattering::UtilityInterpolantScattering(const UtilityInterpolantScattering& collection)
    : UtilityInterpolant(collection)
{
}

UtilityInterpolantScattering::~UtilityInterpolantScattering() {}

double UtilityInterpolantScattering::Calculate(double ei, double ef, double rnd)
{
    (void)rnd;

    if (std::abs(ei - ef) > std::abs(ei) * HALF_PRECISION)
    {
        double aux  = interpolant_->Interpolate(ei);
        double aux2 = aux - interpolant_->Interpolate(ef);

        if (std::abs(aux2) > std::abs(aux) * HALF_PRECISION)
        {
            return aux2;
        } else
        {
            return interpolant_diff_->Interpolate((ei + ef) / 2) * (ef - ei);
        }
    } else
    {
        return interpolant_diff_->Interpolate((ei + ef) / 2) * (ef - ei);
    }
}

// ------------------------------------------------------------------------- //
double UtilityInterpolantScattering::GetUpperLimit(double ei, double rnd)
{
    (void)ei;
    (void)rnd;

    return 0;
}

double UtilityInterpolantScattering::BuildInterpolant(double energy, UtilityIntegral& utility, Integral& integral)
{
    return integral.Integrate(energy, interpolation_def_.max_node_energy, std::bind(&UtilityIntegral::FunctionToIntegral, &utility, std::placeholders::_1), 4);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantScattering::InitInterpolation(const std::string& name,
                                                     UtilityIntegral& utility,
                                                     int number_of_sampling_points)
{
    UtilityInterpolant::InitInterpolation(name, utility, number_of_sampling_points);
}
