

#include <boost/functional/hash.hpp>
#include <cmath>

#include "PROPOSAL/crossection/parametrization/PhotoRealPhotonAssumption.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/Components.h"

#define PHOTO_PARAM_REAL_IMPL(param, parent)                                                                           \
    Photo##param::Photo##param(const ParticleDef& particle_def,                                                        \
                               const Medium& medium,                                                                   \
                               const EnergyCutSettings& cuts,                                                          \
                               const RealPhoton& hardBB,                                                               \
                               double multiplier)                                                                      \
        : Photo##parent(particle_def, medium, cuts, hardBB, multiplier)                                                \
    {                                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    Photo##param::Photo##param(const Photo##param& photo)                                                              \
        : Photo##parent(photo)                                                                                         \
    {                                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    Photo##param::~Photo##param() {}                                                                                   \
                                                                                                                       \
    const std::string Photo##param::name_ = "Photo" #param;

using namespace PROPOSAL;

/******************************************************************************
*                         PhotoRealPhotonAssumption                           *
******************************************************************************/

PhotoRealPhotonAssumption::PhotoRealPhotonAssumption(const ParticleDef& particle_def,
                                                     const Medium& medium,
                                                     const EnergyCutSettings& cuts,
                                                     const RealPhoton& hardBB,
                                                     double multiplier)
    : Photonuclear(particle_def, medium, cuts, multiplier)
    , hardBB_(hardBB.clone())
{
}

PhotoRealPhotonAssumption::PhotoRealPhotonAssumption(const PhotoRealPhotonAssumption& photo)
    : Photonuclear(photo)
    , hardBB_(photo.hardBB_->clone())
{
}

PhotoRealPhotonAssumption::~PhotoRealPhotonAssumption()
{
    delete hardBB_;
}

// ------------------------------------------------------------------------- //
// Bezrukov, Bugaev
// Sov. J. Nucl. Phys. 33 (1981), 635
// eq. 23
// ------------------------------------------------------------------------- //
double PhotoRealPhotonAssumption::DifferentialCrossSection(double energy, double v)
{
    double aux, aum, G, t;

    double nu = v * energy * 1.e-3;

    double sgn = CalculateParametrization(nu);

    const double m1 = 0.54; // in GeV^2
    const double m2 = 1.80; // in GeV^2

    double kappa = 1 - 2 / v + 2 / (v * v);

    // calculate the shadowing factor
    if (components_[component_index_]->GetNucCharge() == 1)
    {
        G = 1;
    }
    else
    {
        // eq. 18
        double tmp = 0.00282 * pow(components_[component_index_]->GetAtomicNum(), 1. / 3) * sgn;
        // eq. 3
        G = (3 / tmp) * (0.5 + ((1 + tmp) * exp(-tmp) - 1) / (tmp * tmp));
    }

    // enhanced formula by Bugaev Shelpin
    // Phys. Rev. D 67 (2003), 034027
    // eq. 4.6
    G *= 3;
    aux = v * particle_def_.mass * 1.e-3;
    t   = aux * aux / (1 - v);
    aum = particle_def_.mass * 1.e-3;
    aum *= aum;
    aux = 2 * aum / t;
    aux = G * ((kappa + 4 * aum / m1) * log(1 + m1 / t) - (kappa * m1) / (m1 + t) - aux)
        + ((kappa + 2 * aum / m2) * log(1 + m2 / t) - aux)
        + aux * (G * (m1 - 4 * t) / (m1 + t) + (m2 / t) * log(1 + t / m2));

    aux *= ALPHA / (8 * PI) * components_[component_index_]->GetAtomicNum() * v * sgn * 1.e-30;

    // hard component by Bugaev, Montaruli, Shelpin, Sokalski
    // Astrop. Phys. 21 (2004), 491
    // in the appendix
    aux += components_[component_index_]->GetAtomicNum() * 1.e-30 * hardBB_->CalculateHardBB(energy, v);

    return multiplier_ * medium_->GetMolDensity() * components_[component_index_]->GetAtomInMolecule()
            * particle_def_.charge * particle_def_.charge * aux;
}

// ------------------------------------------------------------------------- //
// fit of Caldwell at al.
// Phys. Rev Let. 42 (1979), 553
// Table 1
// ------------------------------------------------------------------------- //
double PhotoRealPhotonAssumption::NucleusCrossSectionCaldwell(double nu)
{
    return 49.2 + 11.1 * log(nu) + 151.8 / sqrt(nu);
}

// ------------------------------------------------------------------------- //
// Getter
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
size_t PhotoRealPhotonAssumption::GetHash() const
{
    size_t seed = Parametrization::GetHash();
    boost::hash_combine(seed, hardBB_);

    return seed;
}

// ------------------------------------------------------------------------- //
// Print
// ------------------------------------------------------------------------- //

void PhotoRealPhotonAssumption::print(std::ostream& os) const
{
    os << "hardBB enabled: " << (dynamic_cast<HardBB*>(hardBB_)? 1 : 0) << '\n';
}
/******************************************************************************
*                            Zeus Parametrization                            *
******************************************************************************/

// Signature: (new class, parent class)
PHOTO_PARAM_REAL_IMPL(Zeus, RealPhotonAssumption)

// ------------------------------------------------------------------------- //
double PhotoZeus::CalculateParametrization(double nu)
{
    double aux;

    aux = nu * 2.e-3 * this->components_[this->component_index_]->GetAverageNucleonWeight();
    aux = 63.5 * pow(aux, 0.097) + 145 * pow(aux, -0.5);

    return aux;
}



/******************************************************************************
*                      Bezrukov Bugaev Parametrization                       *
******************************************************************************/

// Signature: (new class, parent class)
PHOTO_PARAM_REAL_IMPL(BezrukovBugaev, RealPhotonAssumption)

// ------------------------------------------------------------------------- //
// Bezrukov, Bugaev
// Sov. J. Nucl. Phys. 33 (1981), 635
// eq. 17
// TODO: look at Sov. J. Nucl. Phys. 32 (1980), 847
// ------------------------------------------------------------------------- //
double PhotoBezrukovBugaev::CalculateParametrization(double nu)
{
    double aux;

    aux = log(0.0213 * nu);
    aux = 114.3 + 1.647 * aux * aux;

    return aux;
}

/******************************************************************************
*                          Kokoulin Parametrization                           *
******************************************************************************/

// Signature: (new class, parent class)
PHOTO_PARAM_REAL_IMPL(Kokoulin, BezrukovBugaev)

// ------------------------------------------------------------------------- //
double PhotoKokoulin::CalculateParametrization(double nu)
{
    if (nu <= 200.)
    {
        if (nu <= 17.)
        {
            return 96.1 + 82. / sqrt(nu);
        } else
        {
            return PhotoBezrukovBugaev::CalculateParametrization(nu);
        }
    } else
    {
        return NucleusCrossSectionCaldwell(nu);
    }
}

/******************************************************************************
*                           Rhode Parametrization                            *
******************************************************************************/

PhotoRhode::PhotoRhode(const ParticleDef& particle_def,
                       const Medium& medium,
                       const EnergyCutSettings& cuts,
                       const RealPhoton& hardBB,
                       double multiplier)
    : PhotoRealPhotonAssumption(particle_def, medium, cuts, hardBB, multiplier)
    , interpolant_(NULL)
{
    double x_aux[] = { 0,           0.1,         0.144544,   0.20893,     0.301995,    0.436516,    0.630957,
                       0.912011,    1.31826,     1.90546,    2.75423,     3.98107,     5.7544,      8.31764,
                       12.0226,     17.378,      25.1189,    36.3078,     52.4807,     75.8577,     109.648,
                       158.489,     229.087,     331.131,    478.63,      691.831,     1000,        1445.44,
                       2089.3,      3019.95,     4365.16,    6309.58,     9120.12,     13182.6,     19054.6,
                       27542.3,     39810.8,     57544,      83176.4,     120226,      173780,      251188,
                       363078,      524807,      758576,     1.09648e+06, 1.58489e+06, 2.29086e+06, 3.3113e+06,
                       4.78628e+06, 6.91828e+06, 9.99996e+06 };

    double y_aux[] = { 0,       0.0666667, 0.0963626, 159.74,  508.103, 215.77,  236.403, 201.919, 151.381,
                       145.407, 132.096,   128.546,   125.046, 121.863, 119.16,  117.022, 115.496, 114.607,
                       114.368, 114.786,   115.864,   117.606, 120.011, 123.08,  126.815, 131.214, 136.278,
                       142.007, 148.401,   155.46,    163.185, 171.574, 180.628, 190.348, 200.732, 211.782,
                       223.497, 235.876,   248.921,   262.631, 277.006, 292.046, 307.751, 324.121, 341.157,
                       358.857, 377.222,   396.253,   415.948, 436.309, 457.334, 479.025 };

    std::vector<double> x(x_aux, x_aux + sizeof(x_aux) / sizeof(double));
    std::vector<double> y(y_aux, y_aux + sizeof(y_aux) / sizeof(double));

    interpolant_ = new Interpolant(x, y, 4, false, false);
}

PhotoRhode::PhotoRhode(const PhotoRhode& photo)
    : PhotoRealPhotonAssumption(photo)
    , interpolant_(new Interpolant(*photo.interpolant_))
{
}

PhotoRhode::~PhotoRhode()
{
    delete interpolant_;
}


Photonuclear* PhotoRhode::create(const ParticleDef& particle_def,
                        const Medium& medium,
                        const EnergyCutSettings& cuts,
                        const RealPhoton& hardBB,
                        double multiplier)
{
    return new PhotoRhode(particle_def, medium, cuts, hardBB, multiplier);
}

// ------------------------------------------------------------------------- //
double PhotoRhode::MeasuredSgN(double e)
{
    return interpolant_->InterpolateArray(e);
}

// ------------------------------------------------------------------------- //
double PhotoRhode::CalculateParametrization(double nu)
{
    if (nu <= 200.)
    {
        return MeasuredSgN(nu);
    } else
    {
        return this->NucleusCrossSectionCaldwell(nu);
    }
}

const std::string PhotoRhode::name_ = "PhotoRhode";

#undef Q2_PHOTO_PARAM_INTEGRAL_IMPL
