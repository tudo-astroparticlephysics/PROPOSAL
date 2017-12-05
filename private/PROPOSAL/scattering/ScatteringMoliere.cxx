
#include <boost/math/special_functions/erf.hpp>

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/scattering/Coefficients.h"
#include "PROPOSAL/Constants.h"
// #include "PROPOSAL/Output.h"
#include "PROPOSAL/scattering/ScatteringMoliere.h"
#include "PROPOSAL/methods.h"

using namespace std;
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

    double momentum = particle_.GetMomentum();  // momentum in MeV/c
    double mass = particle_.GetMass();          // mass in MeV/c²

    double beta_Sq = 1./( 1.+mass*mass/(momentum*momentum) ); //beta² = v²/c²
    double chi_0 = 0.;
    double ZSq_average = 0.;

    vector<double> chi_A_Sq; // screening angle² in rad²
    chi_A_Sq.resize(numComp_);

    for(int i = 0; i < numComp_; i++)
    {
        // Calculate Chi_0
        chi_0 = ( ME * ALPHA * pow(Zi_[i] * 128. / (9. * PI * PI), 1./3.) )/momentum ;
        // Calculate Chi_a^2
        chi_A_Sq[i] = chi_0 * chi_0 * ( 1.13 + 3.76 * ALPHA * ALPHA * Zi_[i] * Zi_[i] / beta_Sq );

        // Calculate Z^2_average for Chi_c^2
        //if case of an electron, replace Z² by Z(Z+1) to into account scatterings
        //on atomic electrons in the medium
        if(mass == ME) ZSq_average += weight_[i] * Zi_[i] * (Zi_[i] + 1.);
        else ZSq_average += weight_[i] * Zi_[i] * Zi_[i];
    }
    // Calculate Chi_c^2
    chiCSq_ = ( (4.*PI*NA*ALPHA*ALPHA*HBAR*HBAR*SPEED*SPEED)
                * (medium_->GetMassDensity()*medium_->GetDensityCorrection()*dr)
                / (momentum*momentum*beta_Sq) )
            * ( ZSq_average/A_average_ );


    // Calculate B
    Scattering::RandomAngles random_angles;

    for(int i = 0; i < numComp_; i++)
    {
        //calculate B-ln(B) = ln(chi_c^2/chi_a^2)+1-2*EULER_MASCHERONI via Newton-Raphson method
        double xn = 15.;

        for(int n = 0; n < 6; n++)
        {
            xn = xn * ( (1. - log(xn) - log(chiCSq_ / chi_A_Sq[i]) - 1. + 2. * EULER_MASCHERONI) / (1. - xn) );
        }

        //  Check for inappropriate values of B. If B < 4.5 it is practical to assume no deviation.
        if( (xn < 4.5) || xn != xn )
        {
            random_angles.sx = 0;
            random_angles.sy = 0;
            random_angles.tx = 0;
            random_angles.ty = 0;
            return random_angles;
        }

        B_[i] = xn;
    }

    double rnd1,rnd2;

    rnd1 = GetRandom();
    rnd2 = GetRandom();

    random_angles.sx = (rnd1/SQRT3+rnd2)/2;
    random_angles.tx = rnd2;

    rnd1 = GetRandom();
    rnd2 = GetRandom();

    random_angles.sy = (rnd1/SQRT3+rnd2)/2;
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
    , A_average_(0.0)
    , Zi_(numComp_)
    , weight_(numComp_)
    , chiCSq_(0.0)
    , B_(numComp_)
{
    std::vector<double> Ai(numComp_, 0); // atomic number of different components
    std::vector<double> ki(numComp_, 0); // number of atoms in molecule of different components
    double A = 0.;

    for (int i = 0; i < numComp_; i++)
    {
        Components::Component* component = medium_->GetComponents().at(i);
        Zi_[i] = component->GetNucCharge();
        ki[i] = component->GetAtomInMolecule();
        Ai[i] = component->GetAtomicNum();

        A += ki[i] * Ai[i];
    }

    for(int i = 0; i < numComp_; i++)
    {
        weight_[i] = ki[i] * Ai[i] / A;
        A_average_ += weight_[i] * Ai[i];
    }
}

ScatteringMoliere::ScatteringMoliere(const ScatteringMoliere& scattering)
    : Scattering(scattering)
    , medium_(scattering.medium_->clone())
    , numComp_(scattering.numComp_)
    , A_average_(scattering.A_average_)
    , Zi_(scattering.Zi_)
    , weight_(scattering.weight_)
    , chiCSq_(scattering.chiCSq_)
    , B_(scattering.B_)
{
}

ScatteringMoliere::ScatteringMoliere(Particle& particle, const ScatteringMoliere& scattering)
    : Scattering(particle)
    , medium_(scattering.medium_->clone())
    , numComp_(scattering.numComp_)
    , A_average_(scattering.A_average_)
    , Zi_(scattering.Zi_)
    , weight_(scattering.weight_)
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
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// ScatteringMoliere& ScatteringMoliere::operator=(const ScatteringMoliere &scattering){
//     if (this != &scattering)
//     {
//       ScatteringMoliere tmp(scattering);
//       swap(tmp);
//     }
//     return *this;
// }

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
    else if (A_average_ != scatteringMoliere->A_average_)
        return false;
    else if (Zi_ != scatteringMoliere->Zi_)
        return false;
    else if (weight_ != scatteringMoliere->weight_)
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

// void ScatteringMoliere::swap(ScatteringMoliere &scattering)
// {
//     using std::swap;
//
//     swap(this->Zi_ , scattering.Zi_);
//     swap(this->ki_ , scattering.ki_);
//     swap(this->Ai_ , scattering.Ai_);
//     swap(this->weight_ , scattering.weight_);
//     swap(this->chi0_ , scattering.chi0_);
//     swap(this->chiASq_ , scattering.chiASq_);
//     swap(this->B_ , scattering.B_);
//
//     swap(this->dx_ , scattering.dx_);
//     swap(this->betaSq_ , scattering.betaSq_);
//     swap(this->p_ , scattering.p_);
//     swap(this->m_ , scattering.m_);
//     swap(this->numComp_ , scattering.numComp_);
//     swap(this->chiCSq_ , scattering.chiCSq_);
//
//     if(scattering.medium_ != NULL)
//     {
//         medium_->swap(*scattering.medium_) ;
//     }
//     else
//     {
//         medium_ = NULL;
//     }
// }

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
        return 0.5 * sqrt(PI) / (pow(x, 1.5) * pow(1. - 4.5 / x, 2. / 3.));

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
        return a * x + b * sqrt(x) + c;

    double sum = 0;

    for (int p = 2; p < 13; p++)
    {
        sum += c2large[p] * (0.5 * log(x) + s2large[p]) * pow(x, -(p + 0.5));
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
    double y2 = 0;

    for (int i = 0; i < numComp_; i++)
    {
        double x = theta * theta / (chiCSq_ * B_.at(i));

        y1 += weight_.at(i) * Zi_.at(i) * Zi_.at(i) / sqrt(chiCSq_ * B_.at(i) * PI) *
              (exp(-x) + f1M(x) / B_.at(i) + f2M(x) / (B_.at(i) * B_.at(i)));
        y2 += weight_.at(i) * Zi_.at(i) * Zi_.at(i);
    }

    return y1 / y2;
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

    return sum * sqrt(x);
}

//----------------------------------------------------------------------------//

double F2Mlarge(double x)
{
    double a = -0.00026360133958801203364619158975302;
    double b = 0.0039965027465457608410459577896745;
    double c = -0.016305842044996649714549974419242;

    // the junction of both parametrizations is smoothed by an interpolation parabola
    if ((x >= 4.25 * 4.25) && (x <= 6.5 * 6.5))
        return a * x + b * sqrt(x) + c;

    double sum = 0;

    for (int p = 2; p < 13; p++)
    {
        sum += -0.5 * c2large[p] / p * (0.5 / p + 0.5 * log(x) + s2large[p]) * pow(x, -p);
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

    return sum * sqrt(x);
}

//----------------------------------------------------------------------------//

double ScatteringMoliere::F(double theta)
{
    double y1 = 0;
    double y2 = 0;

    for (int i = 0; i < numComp_; i++)
    {
        double x = theta * theta / (chiCSq_ * B_.at(i));

        y1 += weight_.at(i) * Zi_.at(i) * Zi_.at(i) *
              (0.5 * boost::math::erf(sqrt(x)) + sqrt(1. / PI) * (F1M(x) / B_.at(i) + F2M(x) / (B_.at(i) * B_.at(i))));
        y2 += weight_.at(i) * Zi_.at(i) * Zi_.at(i);
    }

    return (theta < 0.) ? (-1.) * y1 / y2 : y1 / y2;
}

//----------------------------------------------------------------------------//
//-------------------------generate random angle------------------------------//
//----------------------------------------------------------------------------//

double ScatteringMoliere::GetRandom()
{
    //  Generate random angles following Moliere's distribution by comparing a
    //  uniformly distributed random number with the integral of the distribution.
    //  Therefore, determine the angle where the integral is equal to the random number.

    // rndm element of ]-0.5,0.5]
    double rndm = (RandomGenerator::Get().RandomDouble() - 0.5);

    // Newton-Raphson method:
    double theta_n;

    // guessing an initial value by assuming a gaussian distribution
    // only regarding the component j with maximum weight

    int j = 0;
    for (int i = 0; i + 1 < numComp_; i++)
    {
        if (weight_.at(i + 1) > weight_.at(i))
            j = i + 1;
    }
    double theta_np1 = sqrt(chiCSq_ * B_.at(j)) * erfInv(2. * rndm);

    // iterating until the number of correct digits is greater than 4

    do
    {
        theta_n   = theta_np1;
        theta_np1 = theta_n - (F(theta_n) - rndm) / f(theta_n);

    } while (abs((theta_n - theta_np1) / theta_np1) > 1e-4);

    return theta_np1;
}


