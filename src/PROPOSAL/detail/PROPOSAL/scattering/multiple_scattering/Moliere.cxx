
#include <cmath>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/scattering/multiple_scattering/Coefficients.h"
#include "PROPOSAL/scattering/multiple_scattering/Moliere.h"
#include <CubicInterpolation/Axis.h>
#include "PROPOSAL/Logging.h"

using namespace PROPOSAL::multiple_scattering;

double Moliere::GetPrefactor(double ei, double grammage) {
    double momentum_Sq = (ei - mass) * (ei + mass);
    double beta_Sq = 1. / (1. + mass * mass / momentum_Sq); // beta^2 = (v/c)^2

    // beta^2 p^2 with beta^2 = (v/c)^2 = 1/(1+m^2/p^2)
    double beta_p_Sq = momentum_Sq / ei;
    beta_p_Sq *= beta_p_Sq;

    double chi_0 = 0.;

    std::vector<double> chi_A_Sq; // screening angle^2 in rad^2
    chi_A_Sq.resize(numComp_);

    for (int i = 0; i < numComp_; i++) {
        // Calculate Chi_0 * p
        chi_0 = ME * ALPHA * std::pow(Zi_[i] * 128. / (9. * PI * PI), 1. / 3.);
        // Calculate Chi_a^2
        chi_A_Sq[i] = chi_0 * chi_0 / momentum_Sq
                      * (1.13 + 3.76 * ALPHA * ALPHA * Zi_[i] * Zi_[i] / beta_Sq);
    }

    // Calculate Chi_c^2
    chiCSq_ = ((4. * PI * NA * ALPHA * ALPHA * HBAR * HBAR * SPEED * SPEED)
               * (grammage) / beta_p_Sq)
              * ZSq_A_average_;

    // Calculate B

    for (int i = 0; i < numComp_; i++) {
        // calculate B-ln(B) = ln(chi_c^2/chi_a^2)+1-2*EULER_MASCHERONI via
        // Newton-Raphson method
        double xn = 15.;

        for (int n = 0; n < 6; n++) {
            if (xn < 0)
                return 0; // xn would become nan for further iterations
            xn = xn
                 * ((1. - std::log(xn) - std::log(chiCSq_ / chi_A_Sq[i]) - 1.
                     + 2. * EULER_MASCHERONI)
                    / (1. - xn));
        }

        //  Check for inappropriate values of B. If B < 4.5 it is practical to
        //  assume no deviation.
        if ((xn < 4.5) || xn != xn) {
            return 0;
        }

        B_[i] = xn;
    }

    double pre_factor = std::sqrt(chiCSq_ * B_[max_weight_index_]);
    return pre_factor;
}

ScatteringOffset Moliere::CalculateRandomAngle(
    double grammage, double ei, double ef, const std::array<double, 4>& rnd)
{
    (void)ef;
    ScatteringOffset offsets;

    auto pre_factor = GetPrefactor(ei, grammage);

    if (pre_factor == 0)
        return offsets;

    auto rnd1 = GetRandom(pre_factor, rnd[0]);
    auto rnd2 = GetRandom(pre_factor, rnd[1]);

    offsets.sx = 0.5 * (rnd1 / SQRT3 + rnd2);
    offsets.tx = rnd2;

    rnd1 = GetRandom(pre_factor, rnd[2]);
    rnd2 = GetRandom(pre_factor, rnd[3]);

    offsets.sy = 0.5 * (rnd1 / SQRT3 + rnd2);
    offsets.ty = rnd2;

    return offsets;
}

double Moliere::CalculateScatteringAngle(double grammage, double ei, double ef, double rnd) {
    (void)ef;

    auto pre_factor = GetPrefactor(ei, grammage);

    if (pre_factor == 0)
        return 0;

    return GetRandom(pre_factor, rnd);
};

double Moliere::CalculateScatteringAngle2D(double grammage, double ei, double ef, double rnd1, double rnd2) {
    (void)ef;

    auto pre_factor = GetPrefactor(ei, grammage);

    if (pre_factor == 0)
        return 0;

    auto angle1 = GetRandom(pre_factor, rnd1);
    auto angle2 = GetRandom(pre_factor, rnd2);

    return std::sqrt(angle1 * angle1 + angle2 * angle2);
};

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Moliere::Moliere(const ParticleDef& particle_def, Medium const& medium)
    : Parametrization(particle_def.mass)
    , numComp_(medium.GetNumComponents())
    , ZSq_A_average_(0.0)
    , Zi_(numComp_)
    , weight_ZZ_(numComp_)
    , weight_ZZ_sum_(0.)
    , max_weight_index_(0)
    , chiCSq_(0.0)
    , B_(numComp_)
{
    std::vector<double> Ai(numComp_,
        0); // atomic number of different components
    std::vector<double> ki(
        numComp_, 0); // number of atoms in molecule of different components
    std::vector<double> weight(
        numComp_, 0); // number of atoms in molecule of different components
    double A_sum = 0.;

    for (int i = 0; i < numComp_; i++) {
        Component component = medium.GetComponents().at(i);
        Zi_[i] = component.GetNucCharge();
        ki[i] = component.GetAtomInMolecule();
        Ai[i] = component.GetAtomicNum();

        A_sum += ki[i] * Ai[i];
    }

    double A_average = 0.;
    double ZSq_average = 0.;

    for (int i = 0; i < numComp_; i++) {
        weight[i] = ki[i] * Ai[i] / A_sum;
        A_average += weight[i] * Ai[i];
        weight_ZZ_[i] = weight[i] * Zi_[i] * Zi_[i];
        weight_ZZ_sum_ += weight_ZZ_[i];

        // Calculate Z^2_average for Chi_c^2
        // in case of an electron, replace Z^2 by Z(Z+1) to take into account
        // scatterings on atomic electrons in the medium
        if (mass == ME)
            ZSq_average += weight[i] * Zi_[i] * (Zi_[i] + 1.);
        else
            ZSq_average += weight_ZZ_[i];
    }
    weight_ZZ_sum_ = 1. / weight_ZZ_sum_;
    ZSq_A_average_ = ZSq_average / A_average;

    // guessing an initial value by assuming a gaussian distribution
    // only regarding the component with maximum weight
    // needed in the method GetRandom()
    for (int i = 0; i + 1 < numComp_; i++) {
        if (weight[i + 1] > weight[i])
            max_weight_index_ = i + 1;
    }
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

bool Moliere::compare(const Parametrization& scattering) const
{
    auto sc = dynamic_cast<const Moliere*>(&scattering);

    if (!sc)
        return false;
    else if (numComp_ != sc->numComp_)
        return false;
    else if (ZSq_A_average_ != sc->ZSq_A_average_)
        return false;
    else if (Zi_ != sc->Zi_)
        return false;
    else if (weight_ZZ_ != sc->weight_ZZ_)
        return false;
    else if (weight_ZZ_sum_ != sc->weight_ZZ_sum_)
        return false;
    else if (max_weight_index_ != sc->max_weight_index_)
        return false;
    else if (chiCSq_ != sc->chiCSq_)
        return false;
    else if (B_ != sc->B_)
        return false;
    else
        return true;
}

void Moliere::print(std::ostream& os) const { os << "work in progress\n"; }

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------private member functions---------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//--------------------------calculate distribution----------------------------//
//----------------------------------------------------------------------------//

double Moliere::f1M(double x)
{
    // approximation for large numbers to avoid numerical errors
    if (x > 12.)
        return 0.5 * std::sqrt(PI)
            / (std::pow(x, 1.5) * std::pow(1. - 4.5 / x, 2. / 3.));

    double sum = c1[69];

    // Horner's method
    for (int p = 68; p >= 0; p--)
        sum = sum * x + c1[p];

    return sum;
}

//----------------------------------------------------------------------------//

double f2Mlarge(double x)
{
    double a = 0.00013567765224589459194769192063035;
    double b = -0.0022635502525409950842771866774683;
    double c = 0.0098037758070269476889935233998585;

    // the junction of both parametrizations is smoothed by an interpolation
    // parabola
    if ((x >= 4.25 * 4.25) && (x <= 6.5 * 6.5))
        return a * x + b * std::sqrt(x) + c;

    double sum = 0;

    for (int p = 2; p < 13; p++) {
        sum += PROPOSAL::c2large[p] * (0.5 * std::log(x) + PROPOSAL::s2large[p])
            * std::pow(x, -(p + 0.5));
    }

    return sum;
}

double Moliere::f2M(double x)
{
    // approximation for larger x to avoid numerical errors
    if (x > 4.25 * 4.25)
        return f2Mlarge(x);

    double sum = c2[69];

    for (int p = 68; p >= 0; p--)
        sum = sum * x + c2[p];

    return sum;
}

//----------------------------------------------------------------------------//

double Moliere::f(double theta)
{
    double y1 = 0;

    for (int i = 0; i < numComp_; i++) {
        double x = theta * theta / (chiCSq_ * B_[i]);

        y1 += weight_ZZ_[i] / std::sqrt(chiCSq_ * B_[i] * PI)
            * (std::exp(-x) + f1M(x) / B_[i] + f2M(x) / (B_[i] * B_[i]));
    }

    return y1 * weight_ZZ_sum_;
}

//----------------------------------------------------------------------------//
//----------------------calculate indefinite integral-------------------------//
//----------------------------------------------------------------------------//

double F1Mlarge(double x)
{
    double sum = PROPOSAL::C1large[14];

    // Horner's method
    for (int p = 13; p >= 0; p--)
        sum = PROPOSAL::C1large[p] + sum / x;

    return sum;
}

double Moliere::F1M(double x)
{
    if (x > 12.)
        return F1Mlarge(x);

    double sum = c1[69] / (2. * 69 + 1.);

    // Horner's method
    for (int p = 68; p >= 0; p--)
        sum = sum * x + c1[p] / (2. * p + 1.);

    return sum * std::sqrt(std::max(x, 0.));
}

//----------------------------------------------------------------------------//

double F2Mlarge(double x)
{
    double a = -0.00026360133958801203364619158975302;
    double b = 0.0039965027465457608410459577896745;
    double c = -0.016305842044996649714549974419242;

    // the junction of both parametrizations is smoothed by an interpolation
    // parabola
    if ((x >= 4.25 * 4.25) && (x <= 6.5 * 6.5))
        return a * x + b * std::sqrt(x) + c;

    double sum = 0;

    for (int p = 2; p < 13; p++) {
        sum += -0.5 * PROPOSAL::c2large[p] / p
            * (0.5 / p + 0.5 * std::log(x) + PROPOSAL::s2large[p])
            * std::pow(x, -p);
    }

    return sum;
}

double Moliere::F2M(double x)
{
    if (x > 4.25 * 4.25)
        return F2Mlarge(x);

    double sum = c2[69] / (2. * 69 + 1.);

    for (int p = 68; p >= 0; p--)
        sum = sum * x + c2[p] / (2. * p + 1.);

    return sum * std::sqrt(std::max(x, 0.));
}

//----------------------------------------------------------------------------//

double Moliere::F(double theta)
{
    double y1 = 0;

    for (int i = 0; i < numComp_; i++) {
        double x = theta * theta / (chiCSq_ * B_[i]);

        y1 += weight_ZZ_[i]
            * (0.5 * std::erf(std::sqrt(std::max(x, 0.)))
                + std::sqrt(1. / PI)
                    * (F1M(x) / B_[i] + F2M(x) / (B_[i] * B_[i])));
    }

    return (theta < 0.) ? (-1.) * y1 * weight_ZZ_sum_ : y1 * weight_ZZ_sum_;
}

//----------------------------------------------------------------------------//
//-------------------------generate random angle------------------------------//
//----------------------------------------------------------------------------//

double Moliere::GetRandom(double pre_factor, double rnd)
{
    //  Generate random angles following Moliere's distribution by comparing a
    //  uniformly distributed random number with the integral of the
    //  distribution. Therefore, determine the angle where the integral is equal
    //  to the random number.

    // Newton-Raphson method:
    double theta_n;

    // guessing an initial value by assuming a gaussian distribution
    // only regarding the component with maximum weight
    double theta_np1 = pre_factor * normalppf(rnd) / SQRT2;

    // rnd element of ]-0.5,0.5]
    rnd = rnd - 0.5;

    // iterating until the number of correct digits is greater than 4
    unsigned int i = 0;
    do {
        i++;
        theta_n = theta_np1;
        theta_np1 = theta_n - (F(theta_n) - rnd) / f(theta_n);
        if (i == 100) {
            Logging::Get("proposal.scattering")->warn(
                    "Iteration in Moliere::GetRandom did not converge after 100 iterations. "
                    "Return current value theta = {} (previous iteration: theta = {}", theta_np1, theta_n);
            return theta_np1;
        }
    } while (std::abs((theta_n - theta_np1) / theta_np1) > 1e-4);

    return theta_np1;
}

MoliereInterpol::MoliereInterpol(const ParticleDef& p_def, const Medium& medium) : Moliere(p_def, medium) {
    // initialize interpolation tables
    Logging::Get("proposal.scattering")->debug("Initialize interpooation tables for MoliereInterpol.");

    auto def_f1M = cubic_splines::CubicSplines<double>::Definition();
    def_f1M.f = [&](double x) {return Moliere::f1M(x);};
    def_f1M.axis = std::make_unique<cubic_splines::LinAxis<double>>(0.f, 20, size_t(100));
    f1M_interpolant_ = std::make_shared<interpolant_t>(std::move(def_f1M));

    auto def_f2M = cubic_splines::CubicSplines<double>::Definition();
    def_f2M.f = [&](double x) {return Moliere::f2M(x);};
    def_f2M.axis = std::make_unique<cubic_splines::LinAxis<double>>(0.f, 20, size_t(100));
    f2M_interpolant_ = std::make_shared<interpolant_t>(std::move(def_f2M));

    auto def_F1M = cubic_splines::CubicSplines<double>::Definition();
    def_F1M.f = [&](double x) {return Moliere::F1M(x);};
    def_F1M.axis = std::make_unique<cubic_splines::LinAxis<double>>(0.f, 20, size_t(100));
    F1M_interpolant_ = std::make_shared<interpolant_t>(std::move(def_F1M));

    auto def_F2M = cubic_splines::CubicSplines<double>::Definition();
    def_F2M.f = [&](double x) {return Moliere::F2M(x);};
    def_F2M.axis = std::make_unique<cubic_splines::LinAxis<double>>(0.f, 20, size_t(100));
    F2M_interpolant_ = std::make_shared<interpolant_t>(std::move(def_F2M));
}

double MoliereInterpol::f1M(double x)
{
    if (x > 20.)
        return Moliere::f1M(x); // use analytical evaluation outside range of interpolation tables
    return f1M_interpolant_->evaluate(x);
}

double MoliereInterpol::f2M(double x)
{
    if (x > 20)
        return Moliere::f2M(x);  // use analytical evaluation outside range of interpolation tables
    return f2M_interpolant_->evaluate(x);
}

double MoliereInterpol::F1M(double x)
{
    if (x > 20.)
        return Moliere::F1M(x); // use analytical evaluation outside range of interpolation tables
    return F1M_interpolant_->evaluate(x);
}

double MoliereInterpol::F2M(double x)
{
    if (x > 20)
        return Moliere::F2M(x); // use analytical evaluation outside range of interpolation tables
    return F2M_interpolant_->evaluate(x);
}