
#include <cmath>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/scattering/multiple_scattering/Coefficients.h"
#include "PROPOSAL/scattering/multiple_scattering/Moliere.h"

using namespace PROPOSAL::multiple_scattering;

ScatteringOffset Moliere::CalculateRandomAngle(
    double grammage, double ei, double ef, const std::array<double, 4>& rnd)
{
    (void)ef;
    ScatteringOffset offsets;

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
            xn = xn
                * ((1. - std::log(xn) - std::log(chiCSq_ / chi_A_Sq[i]) - 1.
                       + 2. * EULER_MASCHERONI)
                    / (1. - xn));
        }

        //  Check for inappropriate values of B. If B < 4.5 it is practical to
        //  assume no deviation.
        if ((xn < 4.5) || xn != xn) {
            return offsets;
        }

        B_[i] = xn;
    }

    double pre_factor = std::sqrt(chiCSq_ * B_[max_weight_index_]);

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

    return sum * std::sqrt(x);
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

    return sum * std::sqrt(x);
}

//----------------------------------------------------------------------------//

double Moliere::F(double theta)
{
    double y1 = 0;

    for (int i = 0; i < numComp_; i++) {
        double x = theta * theta / (chiCSq_ * B_[i]);

        y1 += weight_ZZ_[i]
            * (0.5 * std::erf(std::sqrt(x))
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
    do {
        theta_n = theta_np1;
        theta_np1 = theta_n - (F(theta_n) - rnd) / f(theta_n);

    } while (std::abs((theta_n - theta_np1) / theta_np1) > 1e-4);

    return theta_np1;
}
