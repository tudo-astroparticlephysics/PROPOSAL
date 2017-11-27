
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


//----------------------------------------------------------------------------//
//-----------------------------Coefficients-----------------------------------//
//---for calculating the power series approximation of the moliere function---//
//----------------------------------------------------------------------------//


// const double ScatteringMoliere::c1[100] =
// {
//     0.01824498698928826,
//     -1.054734960967865,
//     1.378945800806554,
//     -0.8101747070430586,
//     0.3020799653590784,
//     -0.08217510264333025,
//     0.01757489395499939,
//     -0.003095373240445565,
//     0.0004633127963647091,
//     -6.029130794154548e-05,
//     6.939349333147516e-06,
//     -7.159829943698265e-07,
//     6.694120779750297e-08,
//     -5.721860009237694e-09,
//     4.504494235551189e-10,
//     -3.286570977596057e-11,
//     2.234424657611573e-12,
//     -1.422140651267253e-13,
//     8.50844668823957e-15,
//     -4.80239726059346e-16,
//     2.56544019782712e-17,
//     -1.3008032373069e-18,
//     6.276721156927733e-20,
//     -2.888980198088948e-21,
//     1.271082178072295e-22,
//     -5.356321836042788e-24,
//     2.165708913690074e-25,
//     -8.415665707369898e-27,
//     3.14768814769592e-28,
//     -1.134804220079789e-29,
//     3.948607075433194e-31,
//     -1.327667573801528e-32,
//     4.318678128634822e-34,
//     -1.360474072770111e-35,
//     4.154710513947732e-37,
//     -1.231145280289021e-38,
//     3.543063949967192e-40,
//     -9.910855129282805e-42,
//     2.69678926347067e-43,
//     -7.143475307428035e-45,
//     1.843336870282648e-46,
//     -4.636847647068766e-48,
//     1.137731434272442e-49,
//     -2.724695331159021e-51,
//     6.372463895718094e-53,
//     -1.4562852800865e-54,
//     3.253589567012378e-56,
//     -7.110068913054863e-58,
//     1.520504345733737e-59,
//     -3.183490667084008e-61,
//     6.528486713873094e-63,
//     -1.311890818170701e-64,
//     2.584252665251881e-66,
//     -4.99221608462446e-68,
//     9.460964789499932e-70,
//     -1.759614514252362e-71,
//     3.212849154606866e-73,
//     -5.761014854313268e-75,
//     1.014807133464239e-76,
//     -1.756624689658179e-78,
//     2.98893079441324e-80,
//     -5.000577728087579e-82,
//     8.228369765510281e-84,
//     -1.332031746609897e-85,
//     2.121957008008998e-87,
//     -3.327287074593805e-89,
//     5.136681841644386e-91,
//     -7.809396838908658e-93,
//     1.169486848348948e-94,
//     -1.72549543381347e-96,
//     2.508809249655915e-98,
//     -3.595413291827729e-100,
//     5.079801219388072e-102,
//     -7.076983722567338e-104,
//     9.723837840810993e-106,
//     -1.317945326279369e-107,
//     1.762410958111165e-109,
//     -2.325652802263962e-111,
//     3.028909076825982e-113,
//     -3.89408163088736e-115,
//     4.942802177984765e-117,
//     -6.195278874139763e-119,
//     7.668956698453442e-121,
//     -9.377048488608183e-123,
//     1.132701624601603e-124,
//     -1.351910188741138e-126,
//     1.594502032991125e-128,
//     -1.858693324663112e-130,
//     2.141681719707223e-132,
//     -2.43963241933955e-134,
//     2.747720630653711e-136,
//     -3.060234006823078e-138,
//     3.370734378882746e-140,
//     -3.672273506135163e-142,
//     3.957653085865992e-144,
//     -4.219715312480803e-146,
//     4.451647266046611e-148,
//     -4.647280664772496e-150,
//     4.80136823945151e-152,
//     -4.909819238915968e-154
// };

// const double ScatteringMoliere::c2[100] =
// {
//     0.3692951315189029,
//     -2.901210618562379,
//     4.763691522462663,
//     -3.668389620520657,
//     1.743233030563621,
//     -0.5857757559172652,
//     0.1507057475725596,
//     -0.03124919421553912,
//     0.005411101880491738,
//     -0.0008029915660482544,
//     0.0001041435915389889,
//     -1.198693445962836e-05,
//     1.239576100587234e-06,
//     -1.163301889847139e-07,
//     9.990755927592503e-09,
//     -7.907851249726336e-10,
//     5.803579436334176e-11,
//     -3.969886777147483e-12,
//     2.542633424164171e-13,
//     -1.530925400639093e-14,
//     8.696260972023862e-16,
//     -4.675164455450479e-17,
//     2.385523336370255e-18,
//     -1.158267999875986e-19,
//     5.363958217246529e-21,
//     -2.374296161332396e-22,
//     1.006470730715194e-23,
//     -4.093159567679296e-25,
//     1.5996354171463e-26,
//     -6.016551229383919e-28,
//     2.180970406636743e-29,
//     -7.629492110705818e-31,
//     2.578781877602269e-32,
//     -8.431429670892704e-34,
//     2.669429048218415e-35,
//     -8.192196123057652e-37,
//     2.439244055889828e-38,
//     -7.052895752459682e-40,
//     1.98198110121125e-41,
//     -5.417437723374357e-43,
//     1.441367960601787e-44,
//     -3.735506587014256e-46,
//     9.436466002585026e-48,
//     -2.325045788414374e-49,
//     5.590871105668466e-51,
//     -1.312817566476688e-52,
//     3.011935532920134e-54,
//     -6.755103130065185e-56,
//     1.481773486168858e-57,
//     -3.180564595352823e-59,
//     6.683429378680969e-61,
//     -1.375495310358957e-62,
//     2.773751380431616e-64,
//     -5.482786499876722e-66,
//     1.062748090934042e-67,
//     -2.020774364918115e-69,
//     3.770691103605001e-71,
//     -6.907025753601122e-73,
//     1.242436024375922e-74,
//     -2.195383992072398e-76,
//     3.811853327900785e-78,
//     -6.505522184103775e-80,
//     1.091628995757091e-81,
//     -1.801512577182476e-83,
//     2.924740086946532e-85,
//     -4.672395489277705e-87,
//     7.346918415137358e-89,
//     -1.137343481536269e-90,
//     1.733816262095051e-92,
//     -2.603397903194606e-94,
//     3.851253561734265e-96,
//     -5.61413579760181e-98,
//     8.066318236042616e-100,
//     -1.142534156297296e-101,
//     1.595701569049596e-103,
//     -2.197898380015261e-105,
//     2.986204180183809e-107,
//     -4.002834291275826e-109,
//     5.294562476317114e-111,
//     -6.911669661645124e-113,
//     8.906363201403711e-115,
//     -1.133064328971227e-116,
//     1.423363576915385e-118,
//     -1.765846157643044e-120,
//     2.163876866238606e-122,
//     -2.619510733347758e-124,
//     3.133137855471063e-126,
//     -3.703162047314654e-128,
//     4.325741534530339e-130,
//     -4.994619251572251e-132,
//     5.701067045150251e-134,
//     -6.433962204999298e-136,
//     7.180006537985007e-138,
//     -7.924088243222738e-140,
//     8.649775953179418e-142,
//     -9.339923463142228e-144,
//     9.977353925282273e-146,
//     -1.054558462375476e-147,
//     1.102954869090941e-149,
//     -1.141626881676896e-151
// };

// const double ScatteringMoliere::c2large[50] =
// {
//     0,
//     0,
//     5.317361552716548,
//     39.88021164537411,
//     279.1614815176188,
//     2093.711111382141,
//     17273.11666890266,
//     157185.3616870142,
//     1571853.616870142,
//     17178114.5272237,
//     203990110.0107814,
//     2617873078.471694,
//     36126648482.90939,
//     533689125315.7068,
//     8405603723722.382,
//     140632216146893.7,
//     2491199257459260,
//     4.658542611448816e+16,
//     9.171505766289856e+17,
//     1.896343692265226e+19,
//     4.108744666574657e+20,
//     9.309550415580998e+21,
//     2.201708673284906e+23,
//     5.425639230594947e+24,
//     1.390936602752523e+26,
//     3.704124648634435e+27,
//     1.023264434185263e+29,
//     2.928582810638222e+30,
//     8.673110631505504e+31,
//     2.65493553219974e+33,
//     8.39149266427418e+34,
//     2.735915970369392e+36,
//     9.192677660441157e+37,
//     3.180369932523594e+39,
//     1.132012922857617e+41,
//     4.142138195001734e+42,
//     1.55695665094477e+44,
//     6.007628448859746e+45,
//     2.378019594340316e+47,
//     9.65026059703239e+48,
//     4.012476774555574e+50,
//     1.708389149781931e+52,
//     7.444305720174763e+53,
//     3.318163098443751e+55,
//     1.512134326290795e+57,
//     7.041974391621668e+58,
//     3.349739182196398e+60,
//     1.626856662820051e+62,
//     8.063550415716772e+63,
//     4.077239907010832e+65
// };

// const double ScatteringMoliere::s2large[50] =
// {
//     0,
//     0,
//     -0.3515783203226215,
//     -0.5515783203226214,
//     -0.6944354631797642,
//     -0.8055465742908754,
//     -0.8964556651999662,
//     -0.9733787421230431,
//     -1.04004540878971,
//     -1.098868938201474,
//     -1.151500517148843,
//     -1.19911956476789,
//     -1.242597825637456,
//     -1.282597825637456,
//     -1.319634862674493,
//     -1.354117621295182,
//     -1.386375685811311,
//     -1.416678716114342,
//     -1.44525014468577,
//     -1.472277171712797,
//     -1.497918197353823,
//     -1.522308441256262,
//     -1.54556425520975,
//     -1.567786477431973,
//     -1.589063073176654,
//     -1.609471236441959,
//     -1.629079079579214,
//     -1.647947004107516,
//     -1.666128822289334,
//     -1.683672681938457,
//     -1.70062183448083,
//     -1.717015277103781,
//     -1.732888292976797,
//     -1.748272908361412,
//     -1.76319828149574,
//     -1.777691035118929,
//     -1.791775542161182,
//     -1.805474172298168,
//     -1.818807505631502,
//     -1.831794518618515,
//     -1.844452746466616,
//     -1.856798425478962,
//     -1.868846618250046,
//     -1.880611324132399,
//     -1.892105577005962,
//     -1.903341532062142,
//     -1.914330543051153,
//     -1.925083231223196,
//     -1.93560954701267,
//     -1.945918825363185
// };

// const double ScatteringMoliere::C1large[50] =
// {
//     0,
//     -0.443113462726379,
//     -0.6646701940895685,
//     -1.661675485223921,
//     -5.815864198283724,
//     -26.17138889227676,
//     -143.9426389075222,
//     -935.6271528988941,
//     -7017.203646741706,
//     -59646.2309973045,
//     -566639.1944743927,
//     -5949711.541981123,
//     -68421682.73278293,
//     -855271034.1597866,
//     -11546158961.15712,
//     -167419304936.7782,
//     -2594999226520.062,
//     -42817487237581.03,
//     -749306026657668,
//     -1.386216149316686e+16,
//     -2.703121491167537e+17,
//     -5.541399056893451e+18,
//     -1.191400797232092e+20,
//     -2.680651793772207e+21,
//     -6.299531715364686e+22,
//     -1.543385270264348e+24,
//     -3.935632439174088e+25,
//     -1.042942596381133e+27,
//     -2.868092140048117e+28,
//     -8.174062599137132e+29,
//     -2.411348466745454e+31,
//     -7.354612823573634e+32,
//     -2.316703039425695e+34,
//     -7.529284878133508e+35,
//     -2.522310434174725e+37,
//     -8.701970997902803e+38,
//     -3.089199704255495e+40,
//     -1.127557892053256e+42,
//     -4.228342095199709e+43,
//     -1.627911706651888e+45,
//     -6.430251241274957e+46,
//     -2.604251752716358e+48,
//     -1.080764477377288e+50,
//     -4.593249028853476e+51,
//     -1.998063327551262e+53,
//     -8.891381807603116e+54,
//     -4.045578722459418e+56,
//     -1.881194105943629e+58,
//     -8.935672003232238e+59,
//     -4.333800921567636e+61
// };

