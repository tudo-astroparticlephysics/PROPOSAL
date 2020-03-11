
#include <cmath>
#include <functional>

#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"

#include "PROPOSAL/crossection/CrossSection.h"

#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/math/MathMethods.h"

using namespace PROPOSAL;

/******************************************************************************
 *                              Utility Integral                              *
 ******************************************************************************/

UtilityInterpolant::UtilityInterpolant(
        CrossSectionList cross, InterpolationDef def)
    : UtilityDecorator(cross)
    , stored_result_(0)
    , interpolant_(NULL)
    , interpolant_diff_(NULL)
    , interpolation_def_(def)
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
    return std::min(
        std::max(ei + rnd / interpolant_diff_->Interpolate(
                          ei + rnd / (2 * interpolant_diff_->Interpolate(ei))),
            utility_.GetParticleDef().low),
        ei);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolant::InitInterpolation(const std::string& name,
    UtilityIntegral& utility, int number_of_sampling_points)
{
    Integral integral(IROMB, IMAXS, IPREC2);
    const ParticleDef& particle_def = utility_.GetParticleDef();

    std::vector<std::pair<Interpolant**, std::function<double(double)>>>
        interpolants;

    interpolants.push_back(std::make_pair(&interpolant_,
        std::bind(&UtilityInterpolant::BuildInterpolant, this,
            std::placeholders::_1, std::ref(utility), std::ref(integral))));
    interpolants.push_back(std::make_pair(&interpolant_diff_,
        std::bind(&UtilityIntegralDisplacement::FunctionToIntegral, &utility,
            std::placeholders::_1)));

    unsigned int number_of_interpolants = interpolants.size();

    std::vector<Interpolant1DBuilder> builder_vec(number_of_interpolants);

    Helper::InterpolantBuilderContainer builder_container(
        number_of_interpolants);

    for (unsigned int i = 0; i < number_of_interpolants; ++i) {
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

        builder_container[i]
            = std::make_pair(&builder_vec[i], interpolants[i].first);
    }

    std::vector<CrossSection*> crosssections = utility_.GetCrosssections();
    std::vector<Parametrization*> params(crosssections.size(), NULL);

    for (unsigned int i = 0; i < crosssections.size(); ++i) {
        params[i] = &crosssections[i]->GetParametrization();
    }

    Helper::InitializeInterpolation(
        name, builder_container, params, interpolation_def_);
}



/******************************************************************************
 *                            Utility Displacement                            *
 ******************************************************************************/

UtilityInterpolantDisplacement::UtilityInterpolantDisplacement(
        CrossSectionList cross, const ParticleDef& p_def, InterpolationDef def)
    : UtilityInterpolant(cross, def)
{
    UtilityIntegralDisplacement utility_disp(cross, p_def);
    InitInterpolation("displacement", utility_disp, def.nodes_propagate);
}

UtilityInterpolantDisplacement::~UtilityInterpolantDisplacement() {}

double UtilityInterpolantDisplacement::Calculate(
    double ei, double ef, double rnd)
{
    (void)rnd;

    if (std::abs(ei - ef) > std::abs(ei) * HALF_PRECISION) {
        double aux;

        stored_result_ = interpolant_->Interpolate(ei);
        aux = stored_result_ - interpolant_->Interpolate(ef);

        if (std::abs(aux) > std::abs(stored_result_) * HALF_PRECISION
            && aux >= 0) {
            return std::max(aux, 0.0);
        }
    }

    stored_result_ = 0;

    return std::max(
        (interpolant_diff_->Interpolate((ei + ef) / 2)) * (ef - ei), 0.0);
}

double UtilityInterpolantDisplacement::GetUpperLimit(double ei, double rnd)
{
    <double(double)>f = [&](double ef) { return Calculate(ei, ef, rnd) - rnd; };
    std::function<double(double)> df = [&](double ef) { return interpolant_diff_->Interpolate(ef); };

    int MaxSteps = 200;
    try {
        return std::max(NewtonRaphson(f, df, 0, ei, ei, MaxSteps,
                            PARTICLE_POSITION_RESOLUTION),
            utility_.GetParticleDef().mass);
    } catch (MathException& e) {
        return utility_.GetParticleDef().mass;
    }
}

double UtilityInterpolantDisplacement::BuildInterpolant(
    double energy, UtilityIntegral& utility, Integral& integral)
{
    return integral.Integrate(energy, utility_.GetParticleDef().low,
        std::bind(&UtilityIntegralDisplacement::FunctionToIntegral, &utility,
            std::placeholders::_1),
        4);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantDisplacement::InitInterpolation(const std::string& name,
    UtilityIntegral& utility, int number_of_sampling_points)
{
    UtilityInterpolant::InitInterpolation(
        name, utility, number_of_sampling_points);
}

/******************************************************************************
 *                            Utility Ineraction                            *
 ******************************************************************************/

UtilityInterpolantInteraction::UtilityInterpolantInteraction(
        CrossSectionList cross, const ParticleDef& p_def, InterpolationDef def)
    : UtilityInterpolant(cross, def)
    , big_low_(0)
    , up_(0)
{
    UtilityIntegralInteraction utility_int(cross, p_def);
    InitInterpolation("interaction", utility_int, def.nodes_propagate);
}

UtilityInterpolantInteraction::~UtilityInterpolantInteraction() {}

double UtilityInterpolantInteraction::Calculate(
    double ei, double ef, double rnd)
{
    (void)rnd;
    (void)ef;

    stored_result_ = interpolant_->Interpolate(ei);

    if (up_) {
        return std::max(stored_result_, 0.0);
    } else {
        return std::max(big_low_ - stored_result_, 0.0);
    }
}

double UtilityInterpolantInteraction::GetUpperLimit(double ei, double rnd)
{
    if (std::abs(rnd) > std::abs(stored_result_) * HALF_PRECISION) {
        double aux;

        if (up_) {
            aux = interpolant_->FindLimit(stored_result_ - rnd);
        } else {
            aux = interpolant_->FindLimit(stored_result_ + rnd);
        }

        if (std::abs(ei - aux) > std::abs(ei) * HALF_PRECISION) {
            return std::min(std::max(aux, utility_.GetParticleDef().low), ei);
        }
    }

    return UtilityInterpolant::GetUpperLimit(ei, rnd);
}

double UtilityInterpolantInteraction::BuildInterpolant(
    double energy, UtilityIntegral& utility, Integral& integral)
{
    if (up_) {
        return integral.Integrate(energy, utility_.GetParticleDef().low,
            std::bind(&UtilityIntegralDisplacement::FunctionToIntegral, &utility,
                std::placeholders::_1),
            4);
    } else {
        return -integral.Integrate(energy, interpolation_def_.max_node_energy,
            std::bind(&UtilityIntegralDisplacement::FunctionToIntegral, &utility,
                std::placeholders::_1),
            4);
    }
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantInteraction::InitInterpolation(const std::string& name,
    UtilityIntegral& utility, int number_of_sampling_points)
{
    Integral integral(IROMB, IMAXS, IPREC2);
    const ParticleDef& particle_def = utility_.GetParticleDef();

    double a
        = std::abs(-integral.Integrate(particle_def.low, particle_def.low * 10,
            std::bind(&UtilityIntegralDisplacement::FunctionToIntegral, &utility,
                std::placeholders::_1),
            4));
    double b = std::abs(-integral.Integrate(interpolation_def_.max_node_energy,
        interpolation_def_.max_node_energy / 10,
        std::bind(&UtilityIntegralDisplacement::FunctionToIntegral, &utility,
            std::placeholders::_1),
        4));

    if (a < b) {
        up_ = true;
    } else {
        up_ = false;
    }

    UtilityInterpolant::InitInterpolation(
        name, utility, number_of_sampling_points);

    big_low_ = interpolant_->Interpolate(particle_def.low);
}

/******************************************************************************
 *                            Utility Decay                            *
 ******************************************************************************/

UtilityInterpolantDecay::UtilityInterpolantDecay(
        CrossSectionList cross, const ParticleDef& p_def, InterpolationDef def)
    : UtilityInterpolant(cross, def)
    , big_low_(0)
    , up_(0)
{
    UtilityIntegralDecay utility_decay(cross, p_def);
    InitInterpolation("decay", utility_decay, def.nodes_propagate);
}

UtilityInterpolantDecay::~UtilityInterpolantDecay() {}

double UtilityInterpolantDecay::Calculate(double ei, double ef, double rnd)
{
    (void)rnd;
    (void)ef;

    stored_result_ = interpolant_->Interpolate(ei);

    if (up_) {
        return std::max(stored_result_, 0.0);
    } else {
        return std::max(big_low_ - stored_result_, 0.0);
    }
}

double UtilityInterpolantDecay::GetUpperLimit(double ei, double rnd)
{
    if (std::abs(rnd) > std::abs(stored_result_) * HALF_PRECISION) {
        double aux;

        aux = interpolant_->FindLimit(stored_result_ + rnd);

        if (std::abs(ei - aux) > std::abs(ei) * HALF_PRECISION) {
            return std::min(std::max(aux, utility_.GetParticleDef().low), ei);
        }
    }

    return UtilityInterpolant::GetUpperLimit(ei, rnd);
}

double UtilityInterpolantDecay::BuildInterpolant(
    double energy, UtilityIntegral& utility, Integral& integral)
{
    return -integral.Integrate(energy, interpolation_def_.max_node_energy,
        std::bind(&UtilityIntegralDisplacement::FunctionToIntegral, &utility,
            std::placeholders::_1),
        4);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantDecay::InitInterpolation(const std::string& name,
    UtilityIntegral& utility, int number_of_sampling_points)
{
    Integral integral(IROMB, IMAXS, IPREC2);
    const ParticleDef& particle_def = utility_.GetParticleDef();

    double a
        = std::abs(-integral.Integrate(particle_def.low, particle_def.low * 10,
            std::bind(&UtilityIntegralDisplacement::FunctionToIntegral, &utility,
                std::placeholders::_1),
            4));
    double b = std::abs(-integral.Integrate(interpolation_def_.max_node_energy,
        interpolation_def_.max_node_energy / 10,
        std::bind(&UtilityIntegralDisplacement::FunctionToIntegral, &utility,
            std::placeholders::_1),
        4));

    if (a < b) {
        up_ = true;
    } else {
        up_ = false;
    }

    UtilityInterpolant::InitInterpolation(
        name, utility, number_of_sampling_points);

    big_low_ = interpolant_->Interpolate(particle_def.low);
}

/******************************************************************************
 *                            Utility Time                            *
 ******************************************************************************/

UtilityInterpolantTime::UtilityInterpolantTime(
        CrossSectionList cross, const ParticleDef& p_def, InterpolationDef def)
    : UtilityInterpolant(cross, def)
{
    UtilityIntegralTime utility_time(cross, p_def);
    InitInterpolation("time", utility_time, def.nodes_propagate);
}

UtilityInterpolantTime::~UtilityInterpolantTime() {}

double UtilityInterpolantTime::Calculate(double ei, double ef, double rnd)
{
    (void)rnd;

    if (std::abs(ei - ef) > std::abs(ei) * HALF_PRECISION) {
        double aux = interpolant_->Interpolate(ei);
        double aux2 = aux - interpolant_->Interpolate(ef);

        if (std::abs(aux2) > std::abs(aux) * HALF_PRECISION) {
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

double UtilityInterpolantTime::BuildInterpolant(
    double energy, UtilityIntegral& utility, Integral& integral)
{
    return integral.Integrate(energy, utility_.GetParticleDef().low,
        std::bind(&UtilityIntegralDisplacement::FunctionToIntegral, &utility,
            std::placeholders::_1),
        4);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantTime::InitInterpolation(const std::string& name,
    UtilityIntegral& utility, int number_of_sampling_points)
{
    UtilityInterpolant::InitInterpolation(
        name, utility, number_of_sampling_points);
}

/******************************************************************************
 *                            Utility ContRand                            *
 ******************************************************************************/

UtilityInterpolantContRand::UtilityInterpolantContRand(
        CrossSectionList cross, const ParticleDef& p_def, InterpolationDef def)
    : UtilityInterpolant(cross, def)
{
    UtilityIntegralContRand utility_contrand(cross, p_def);
    InitInterpolation(
        "contrand", utility_contrand, def.nodes_continous_randomization);
}

UtilityInterpolantContRand::~UtilityInterpolantContRand() {}

double UtilityInterpolantContRand::Calculate(double ei, double ef, double rnd)
{
    if (std::abs(ei - ef) > std::abs(ei) * HALF_PRECISION) {
        double aux = interpolant_->Interpolate(ei);
        double aux2 = aux - interpolant_->Interpolate(ef);

        if (std::abs(aux2) > std::abs(aux) * HALF_PRECISION)
            return std::max(aux2, 0.0);
    } else {
        return std::max(
            interpolant_diff_->Interpolate((ei + ef) / 2) * (ef - ei), 0.0);
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

double UtilityInterpolantContRand::BuildInterpolant(
    double energy, UtilityIntegral& utility, Integral& integral)
{
    return integral.Integrate(energy, utility_.GetParticleDef().low,
        std::bind(&UtilityIntegralDisplacement::FunctionToIntegral, &utility,
            std::placeholders::_1),
        4);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantContRand::InitInterpolation(const std::string& name,
    UtilityIntegral& utility, int number_of_sampling_points)
{
    UtilityInterpolant::InitInterpolation(
        name, utility, number_of_sampling_points);
}

/******************************************************************************
 *                            Utility Scattering                            *
 ******************************************************************************/

UtilityInterpolantScattering::UtilityInterpolantScattering(
        CrossSectionList cross, const ParticleDef& p_def, InterpolationDef def)
    : UtilityInterpolant(cross, def)
{
    UtilityIntegralScattering utility_scattering(cross, p_def);
    InitInterpolation(
        "scattering", utility_scattering, def.nodes_continous_randomization);
}

UtilityInterpolantScattering::~UtilityInterpolantScattering() {}

double UtilityInterpolantScattering::Calculate(double ei, double ef, double rnd)
{
    (void)rnd;

    if (std::abs(ei - ef) > std::abs(ei) * HALF_PRECISION) {
        double aux = interpolant_->Interpolate(ei);
        double aux2 = aux - interpolant_->Interpolate(ef);

        if (std::abs(aux2) > std::abs(aux) * HALF_PRECISION) {
            return aux2;
        } else {
            return interpolant_diff_->Interpolate((ei + ef) / 2) * (ef - ei);
        }
    } else {
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

double UtilityInterpolantScattering::BuildInterpolant(
    double energy, UtilityIntegral& utility, Integral& integral)
{
    return integral.Integrate(energy, interpolation_def_.max_node_energy,
        std::bind(&UtilityIntegralDisplacement::FunctionToIntegral, &utility,
            std::placeholders::_1),
        4);
}

// ------------------------------------------------------------------------- //
void UtilityInterpolantScattering::InitInterpolation(const std::string& name,
    UtilityIntegral& utility, int number_of_sampling_points)
{
    UtilityInterpolant::InitInterpolation(
        name, utility, number_of_sampling_points);
}
