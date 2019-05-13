
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/scattering/Coefficients.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/scattering/ScatteringMoliere.h"

using namespace PROPOSAL;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------public member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

Scattering::RandomAngles ScatteringMoliere::CalculateRandomAngle(double dr, double ei, double ef)
{
    (void)ei;
    (void)ef;

    double momentum = particle_.GetMomentum(); // momentum in MeV/c
    double mass     = particle_.GetMass();     // mass in MeV/c^2

    double beta_Sq = 1. / (1. + mass * mass / (momentum * momentum)); // beta^2 = (v/c)^2
    double chi_0   = 0.;

    std::vector<double> chi_A_Sq; // screening angle^2 in rad^2
    chi_A_Sq.resize(numComp_);

    for (int i = 0; i < numComp_; i++)
    {
        // Calculate Chi_0
        chi_0 = (ME * ALPHA * std::pow(Zi_[i] * 128. / (9. * PI * PI), 1. / 3.)) / momentum;
        // Calculate Chi_a^2
        chi_A_Sq[i] = chi_0 * chi_0 * (1.13 + 3.76 * ALPHA * ALPHA * Zi_[i] * Zi_[i] / beta_Sq);
    }

    // Calculate Chi_c^2
    chiCSq_ = ((4. * PI * NA * ALPHA * ALPHA * HBAR * HBAR * SPEED * SPEED) *
               (medium_->GetMassDensity() * medium_->GetDensityCorrection() * dr) / (momentum * momentum * beta_Sq)) *
              ZSq_A_average_;

    // Calculate B
    Scattering::RandomAngles random_angles;

    for (int i = 0; i < numComp_; i++)
    {
        // calculate B-ln(B) = ln(chi_c^2/chi_a^2)+1-2*EULER_MASCHERONI via Newton-Raphson method
        double xn = 15.;

        for (int n = 0; n < 6; n++)
        {
            xn = xn * ((1. - std::log(xn) - std::log(chiCSq_ / chi_A_Sq[i]) - 1. + 2. * EULER_MASCHERONI) / (1. - xn));
        }

        //  Check for inappropriate values of B. If B < 4.5 it is practical to assume no deviation.
        if ((xn < 4.5) || xn != xn)
        {
            random_angles.sx = 0;
            random_angles.sy = 0;
            random_angles.tx = 0;
            random_angles.ty = 0;
            return random_angles;
        }

        B_[i] = xn;
    }

    double pre_factor = std::sqrt(chiCSq_ * B_[max_weight_index_]);

    double rnd1, rnd2;

    rnd1 = GetRandom(pre_factor);
    rnd2 = GetRandom(pre_factor);

    random_angles.sx = 0.5 * (rnd1 / SQRT3 + rnd2);
    random_angles.tx = rnd2;

    rnd1 = GetRandom(pre_factor);
    rnd2 = GetRandom(pre_factor);

    random_angles.sy = 0.5 * (rnd1 / SQRT3 + rnd2);
    random_angles.ty = rnd2;

    return random_angles;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

ScatteringMoliere::ScatteringMoliere(Particle& particle, const Medium& medium)
    : Scattering(particle)
    , medium_(medium.clone())
    , numComp_(medium_->GetNumComponents())
    , ZSq_A_average_(0.0)
    , Zi_(numComp_)
    , weight_ZZ_(numComp_)
    , weight_ZZ_sum_(0.)
    , max_weight_index_(0)
    , chiCSq_(0.0)
    , B_(numComp_)
{
    std::vector<double> Ai(numComp_, 0);     // atomic number of different components
    std::vector<double> ki(numComp_, 0);     // number of atoms in molecule of different components
    std::vector<double> weight(numComp_, 0); // number of atoms in molecule of different components
    double A_sum = 0.;

    for (int i = 0; i < numComp_; i++)
    {
        Components::Component* component = medium_->GetComponents().at(i);
        Zi_[i]                           = component->GetNucCharge();
        ki[i]                            = component->GetAtomInMolecule();
        Ai[i]                            = component->GetAtomicNum();

        A_sum += ki[i] * Ai[i];
    }

    double A_average   = 0.;
    double ZSq_average = 0.;

    for (int i = 0; i < numComp_; i++)
    {
        weight[i] = ki[i] * Ai[i] / A_sum;
        A_average += weight[i] * Ai[i];
        weight_ZZ_[i] = weight[i] * Zi_[i] * Zi_[i];
        weight_ZZ_sum_ += weight_ZZ_[i];

        // Calculate Z^2_average for Chi_c^2
        // if case of an electron, replace Z^2 by Z(Z+1) to into account scatterings
        // on atomic electrons in the medium
        if (particle_.GetMass() == ME)
            ZSq_average += weight[i] * Zi_[i] * (Zi_[i] + 1.);
        else
            ZSq_average += weight_ZZ_[i];
    }
    weight_ZZ_sum_ = 1. / weight_ZZ_sum_;
    ZSq_A_average_ = ZSq_average / A_average;

    // guessing an initial value by assuming a gaussian distribution
    // only regarding the component with maximum weight
    // needed in the method GetRandom()
    for (int i = 0; i + 1 < numComp_; i++)
    {
        if (weight[i + 1] > weight[i])
            max_weight_index_ = i + 1;
    }
}

ScatteringMoliere::ScatteringMoliere(const ScatteringMoliere& scattering)
    : Scattering(scattering)
    , medium_(scattering.medium_->clone())
    , numComp_(scattering.numComp_)
    , ZSq_A_average_(scattering.ZSq_A_average_)
    , Zi_(scattering.Zi_)
    , weight_ZZ_(scattering.weight_ZZ_)
    , weight_ZZ_sum_(scattering.weight_ZZ_sum_)
    , max_weight_index_(scattering.max_weight_index_)
    , chiCSq_(scattering.chiCSq_)
    , B_(scattering.B_)
{
}

ScatteringMoliere::ScatteringMoliere(Particle& particle, const ScatteringMoliere& scattering)
    : Scattering(particle)
    , medium_(scattering.medium_->clone())
    , numComp_(scattering.numComp_)
    , ZSq_A_average_(scattering.ZSq_A_average_)
    , Zi_(scattering.Zi_)
    , weight_ZZ_(scattering.weight_ZZ_)
    , weight_ZZ_sum_(scattering.weight_ZZ_sum_)
    , max_weight_index_(scattering.max_weight_index_)
    , chiCSq_(scattering.chiCSq_)
    , B_(scattering.B_)
{
}

ScatteringMoliere::~ScatteringMoliere()
{
    delete medium_;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

bool ScatteringMoliere::compare(const Scattering& scattering) const
{
    const ScatteringMoliere* scatteringMoliere = dynamic_cast<const ScatteringMoliere*>(&scattering);

    if (!scatteringMoliere)
        return false;
    else if (*medium_ != *scatteringMoliere->medium_)
        return false;
    else if (numComp_ != scatteringMoliere->numComp_)
        return false;
    else if (ZSq_A_average_ != scatteringMoliere->ZSq_A_average_)
        return false;
    else if (Zi_ != scatteringMoliere->Zi_)
        return false;
    else if (weight_ZZ_ != scatteringMoliere->weight_ZZ_)
        return false;
    else if (weight_ZZ_sum_ != scatteringMoliere->weight_ZZ_sum_)
        return false;
    else if (max_weight_index_ != scatteringMoliere->max_weight_index_)
        return false;
    else if (chiCSq_ != scatteringMoliere->chiCSq_)
        return false;
    else if (B_ != scatteringMoliere->B_)
        return false;
    else
        return true;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------private member functions---------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//--------------------------calculate distribution----------------------------//
//----------------------------------------------------------------------------//

double ScatteringMoliere::f1M(double x)
{
    // approximation for large numbers to avoid numerical errors
    if (x > 12.)
        return 0.5 * std::sqrt(PI) / (std::pow(x, 1.5) * std::pow(1. - 4.5 / x, 2. / 3.));

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

    // the junction of both parametrizations is smoothed by an interpolation parabola
    if ((x >= 4.25 * 4.25) && (x <= 6.5 * 6.5))
        return a * x + b * std::sqrt(x) + c;

    double sum = 0;

    for (int p = 2; p < 13; p++)
    {
        sum += c2large[p] * (0.5 * std::log(x) + s2large[p]) * std::pow(x, -(p + 0.5));
    }

    return sum;
}

double ScatteringMoliere::f2M(double x)
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

double ScatteringMoliere::f(double theta)
{
    double y1 = 0;

    for (int i = 0; i < numComp_; i++)
    {
        double x = theta * theta / (chiCSq_ * B_[i]);

        y1 += weight_ZZ_[i] / std::sqrt(chiCSq_ * B_[i] * PI) * (std::exp(-x) + f1M(x) / B_[i] + f2M(x) / (B_[i] * B_[i]));
    }

    return y1 * weight_ZZ_sum_;
}

//----------------------------------------------------------------------------//
//----------------------calculate indefinite integral-------------------------//
//----------------------------------------------------------------------------//

double F1Mlarge(double x)
{
    double sum = C1large[14];

    // Horner's method
    for (int p = 13; p >= 0; p--)
        sum = C1large[p] + sum / x;

    return sum;
}

double ScatteringMoliere::F1M(double x)
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

    // the junction of both parametrizations is smoothed by an interpolation parabola
    if ((x >= 4.25 * 4.25) && (x <= 6.5 * 6.5))
        return a * x + b * std::sqrt(x) + c;

    double sum = 0;

    for (int p = 2; p < 13; p++)
    {
        sum += -0.5 * c2large[p] / p * (0.5 / p + 0.5 * std::log(x) + s2large[p]) * std::pow(x, -p);
    }

    return sum;
}

double ScatteringMoliere::F2M(double x)
{
    if (x > 4.25 * 4.25)
        return F2Mlarge(x);

    double sum = c2[69] / (2. * 69 + 1.);

    for (int p = 68; p >= 0; p--)
        sum = sum * x + c2[p] / (2. * p + 1.);

    return sum * std::sqrt(x);
}

//----------------------------------------------------------------------------//

double ScatteringMoliere::F(double theta)
{
    double y1 = 0;

    for (int i = 0; i < numComp_; i++)
    {
        double x = theta * theta / (chiCSq_ * B_[i]);

        y1 += weight_ZZ_[i] *
              (0.5 * std::erf(std::sqrt(x)) + std::sqrt(1. / PI) * (F1M(x) / B_[i] + F2M(x) / (B_[i] * B_[i])));
    }

    return (theta < 0.) ? (-1.) * y1 * weight_ZZ_sum_ : y1 * weight_ZZ_sum_;
}

//----------------------------------------------------------------------------//
//-------------------------generate random angle------------------------------//
//----------------------------------------------------------------------------//

double ScatteringMoliere::GetRandom(double pre_factor)
{
    //  Generate random angles following Moliere's distribution by comparing a
    //  uniformly distributed random number with the integral of the distribution.
    //  Therefore, determine the angle where the integral is equal to the random number.

    double rnd = RandomGenerator::Get().RandomDouble();

    // Newton-Raphson method:
    double theta_n;

    // guessing an initial value by assuming a gaussian distribution
    // only regarding the component with maximum weight
    double theta_np1 = pre_factor * inverseErrorFunction(rnd) / SQRT2;

    // rnd element of ]-0.5,0.5]
    rnd = rnd - 0.5;

    // iterating until the number of correct digits is greater than 4
    do
    {
        theta_n   = theta_np1;
        theta_np1 = theta_n - (F(theta_n) - rnd) / f(theta_n);

    } while (std::abs((theta_n - theta_np1) / theta_np1) > 1e-4);

    return theta_np1;
}
