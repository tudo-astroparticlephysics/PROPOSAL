
#include <boost/bind.hpp>
#include <boost/functional/hash.hpp>
#include <cmath>

#include "PROPOSAL/crossection/parametrization/PhotoQ2Integration.h"

#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

#define Q2_PHOTO_PARAM_INTEGRAL_IMPL(param)                                                                            \
    Photo##param::Photo##param(const ParticleDef& particle_def,                                                        \
                               const Medium& medium,                                                                   \
                               const EnergyCutSettings& cuts,                                                          \
                               const ShadowEffect& shadow_effect,                                                      \
                               double multiplier)                                                                      \
        : PhotoQ2Integral(particle_def, medium, cuts, shadow_effect, multiplier)                                       \
    {                                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    Photo##param::Photo##param(const Photo##param& photo)                                                              \
        : PhotoQ2Integral(photo)                                                                                       \
    {                                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    Photo##param::~Photo##param() {}\
                                                                                                                       \
    const std::string Photo##param::name_ = "Photo" #param;

/******************************************************************************
*                            Photo Q2 Integration                            *
******************************************************************************/

PhotoQ2Integral::PhotoQ2Integral(const ParticleDef& particle_def,
                                 const Medium& medium,
                                 const EnergyCutSettings& cuts,
                                 const ShadowEffect& shadow_effect,
                                 double multiplier)
    : Photonuclear(particle_def, medium, cuts, multiplier)
    , shadow_effect_(shadow_effect.clone())
    , integral_(IROMB, IMAXS, IPREC)
{
}

PhotoQ2Integral::PhotoQ2Integral(const PhotoQ2Integral& photo)
    : Photonuclear(photo)
    , shadow_effect_(photo.shadow_effect_->clone())
    , integral_(photo.integral_)
{
}

PhotoQ2Integral::~PhotoQ2Integral()
{
    delete shadow_effect_;
}

// ------------------------------------------------------------------------- //
double PhotoQ2Integral::DifferentialCrossSection(double energy, double v)
{
    IntegralLimits limits = GetIntegralLimits(energy);

    double aux, min, max;
    double particle_mass   = particle_def_.mass;
    double particle_charge = particle_def_.charge;

    min = particle_mass * v;
    min *= min / (1 - v);

    if (particle_mass < MPI)
    {
        aux = particle_mass * particle_mass / energy;
        min -= (aux * aux) / (2 * (1 - v));
    }

    max = 2 * components_[component_index_]->GetAverageNucleonWeight() * energy * (v - limits.vMin);

    //  if(form==4) max=Math.min(max, 5.5e6);  // as requested in Butkevich and Mikheyev
    if (min > max)
    {
        return 0;
    }

    aux = integral_.Integrate(min, max, boost::bind(&PhotoQ2Integral::FunctionToQ2Integral, this, energy, v, _1), 4);

    aux *= multiplier_ * medium_->GetMolDensity() * components_[component_index_]->GetAtomInMolecule()
         * particle_charge * particle_charge;

    return aux;
}

// ------------------------------------------------------------------------- //
// Getter
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
size_t PhotoQ2Integral::GetHash() const
{
    size_t seed = Parametrization::GetHash();
    boost::hash_combine(seed, shadow_effect_->GetHash());

    return seed;
}

/******************************************************************************
*                          Specifc Parametrizations                           *
******************************************************************************/

// ------------------------------------------------------------------------- //
// Define the specific parametrizations
// ------------------------------------------------------------------------- //

Q2_PHOTO_PARAM_INTEGRAL_IMPL(AbramowiczLevinLevyMaor91)
Q2_PHOTO_PARAM_INTEGRAL_IMPL(AbramowiczLevinLevyMaor97)
Q2_PHOTO_PARAM_INTEGRAL_IMPL(ButkevichMikhailov)
Q2_PHOTO_PARAM_INTEGRAL_IMPL(RenoSarcevicSu)

// ------------------------------------------------------------------------- //
// Abramowicz Levin Levy Maor 91 Integral
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double PhotoAbramowiczLevinLevyMaor91::FunctionToQ2Integral(double energy, double v, double Q2)
{
    double x, aux, nu, G, F2, R2;

    Components::Component* component = components_[component_index_];

    nu = v * energy;
    x  = Q2 / (2 * component->GetAverageNucleonWeight() * nu);

    G = shadow_effect_->CalculateShadowEffect(*component, x, nu);

    double P, W2;

    aux = x * x;
    P   = 1 - 1.85 * x + 2.45 * aux - 2.35 * aux * x + aux * aux;
    G *= (component->GetNucCharge() + (component->GetAtomicNum() - component->GetNucCharge()) * P);
    W2 = component->GetAverageNucleonWeight() * component->GetAverageNucleonWeight() - Q2 +
         2 * component->GetAverageNucleonWeight() * energy * v;

    double cp1, cp2, cp3, cr1, cr2, cr3, ap1, ap2, ap3, ar1, ar2, ar3;
    double bp1, bp2, bp3, br1, br2, br3, m2o, m2r, L2, m2p, Q2o;

    cp1 = 0.26550;
    cp2 = 0.04856;
    cp3 = 1.04682;
    cr1 = 0.67639;
    cr2 = 0.49027;
    cr3 = 2.66275;
    ap1 = -0.04503;
    ap2 = -0.36407;
    ap3 = 8.17091;
    ar1 = 0.60408;
    ar2 = 0.17353;
    ar3 = 1.61812;
    bp1 = 0.49222;
    bp2 = 0.52116;
    bp3 = 3.55115;
    br1 = 1.26066;
    br2 = 1.83624;
    br3 = 0.81141;
    m2o = 0.30508;
    m2r = 0.20623;
    L2  = 0.06527;
    m2p = 10.67564;
    Q2o = 0.27799;

    // GeV -> MeV conversion
    m2o *= 1e6;
    m2r *= 1e6;
    L2 *= 1e6;
    m2p *= 1e6;
    Q2o *= 1e6;

    // these values are corrected according to the file f2allm.f from Halina Abramowicz
    bp1 *= bp1;
    bp2 *= bp2;
    br1 *= br1;
    br2 *= br2;
    Q2o += L2;

    const double R = 0;

    double cr, ar, cp, ap, br, bp, t;

    t = log(log((Q2 + Q2o) / L2) / log(Q2o / L2));

    if (t < 0)
    {
        t = 0;
    }

    cr = cr1 + cr2 * pow(t, cr3);
    ar = ar1 + ar2 * pow(t, ar3);
    cp = cp1 + (cp1 - cp2) * (1 / (1 + pow(t, cp3)) - 1);
    ap = ap1 + (ap1 - ap2) * (1 / (1 + pow(t, ap3)) - 1);
    br = br1 + br2 * pow(t, br3);
    bp = bp1 + bp2 * pow(t, bp3);

    double xp, xr, F2p, F2r;

    xp  = (Q2 + m2p) / (Q2 + m2p + W2 - component->GetAverageNucleonWeight() * component->GetAverageNucleonWeight());
    xr  = (Q2 + m2r) / (Q2 + m2r + W2 - component->GetAverageNucleonWeight() * component->GetAverageNucleonWeight());
    F2p = cp * pow(xp, ap) * pow(1 - x, bp);
    F2r = cr * pow(xr, ar) * pow(1 - x, br);
    F2  = (Q2 / (Q2 + m2o)) * (F2p + F2r) * G;
    R2  = (2 * (1 + R));

    aux = ME * RE / Q2;
    aux *=
        aux *
        (1 - v - component->GetAverageNucleonWeight() * x * v / (2 * energy) +
         (1 - 2 * particle_def_.mass * particle_def_.mass / Q2) * v * v *
             (1 + 4 * component->GetAverageNucleonWeight() * component->GetAverageNucleonWeight() * x * x / Q2) / R2);

    return (4 * PI * F2 / v) * aux;
}

// ------------------------------------------------------------------------- //
// Abramowicz Levin Levy Maor 97 Integral
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double PhotoAbramowiczLevinLevyMaor97::FunctionToQ2Integral(double energy, double v, double Q2)
{
    double x, aux, nu, G, F2, R2;

    Components::Component* component = components_[component_index_];

    nu = v * energy;
    x  = Q2 / (2 * component->GetAverageNucleonWeight() * nu);

    G = shadow_effect_->CalculateShadowEffect(*component, x, nu);

    double P, W2;

    aux = x * x;
    P   = 1 - 1.85 * x + 2.45 * aux - 2.35 * aux * x + aux * aux;
    G *= (component->GetNucCharge() + (component->GetAtomicNum() - component->GetNucCharge()) * P);
    W2 = component->GetAverageNucleonWeight() * component->GetAverageNucleonWeight() - Q2 +
         2 * component->GetAverageNucleonWeight() * energy * v;

    double cp1, cp2, cp3, cr1, cr2, cr3, ap1, ap2, ap3, ar1, ar2, ar3;
    double bp1, bp2, bp3, br1, br2, br3, m2o, m2r, L2, m2p, Q2o;

    cp1 = 0.28067;
    cp2 = 0.22291;
    cp3 = 2.1979;
    cr1 = 0.80107;
    cr2 = 0.97307;
    cr3 = 3.4942;
    ap1 = -0.0808;
    ap2 = -0.44812;
    ap3 = 1.1709;
    ar1 = 0.58400;
    ar2 = 0.37888;
    ar3 = 2.6063;
    bp1 = 0.60243;
    bp2 = 1.3754;
    bp3 = 1.8439;
    br1 = 0.10711;
    br2 = 1.9386;
    br3 = 0.49338;
    m2o = 0.31985;
    m2r = 0.15052;
    L2  = 0.06527;
    m2p = 49.457;
    Q2o = 0.46017;

    // GeV -> MeV conversion
    m2o *= 1e6;
    m2r *= 1e6;
    L2 *= 1e6;
    m2p *= 1e6;
    Q2o *= 1e6;

    // these values are corrected according to the file f2allm.f from Halina Abramowicz
    bp1 *= bp1;
    bp2 *= bp2;
    br1 *= br1;
    br2 *= br2;
    Q2o += L2;

    const double R = 0;

    double cr, ar, cp, ap, br, bp, t;

    t = log(log((Q2 + Q2o) / L2) / log(Q2o / L2));

    if (t < 0)
    {
        t = 0;
    }

    cr = cr1 + cr2 * pow(t, cr3);
    ar = ar1 + ar2 * pow(t, ar3);
    cp = cp1 + (cp1 - cp2) * (1 / (1 + pow(t, cp3)) - 1);
    ap = ap1 + (ap1 - ap2) * (1 / (1 + pow(t, ap3)) - 1);
    br = br1 + br2 * pow(t, br3);
    bp = bp1 + bp2 * pow(t, bp3);

    double xp, xr, F2p, F2r;

    xp  = (Q2 + m2p) / (Q2 + m2p + W2 - component->GetAverageNucleonWeight() * component->GetAverageNucleonWeight());
    xr  = (Q2 + m2r) / (Q2 + m2r + W2 - component->GetAverageNucleonWeight() * component->GetAverageNucleonWeight());
    F2p = cp * pow(xp, ap) * pow(1 - x, bp);
    F2r = cr * pow(xr, ar) * pow(1 - x, br);
    F2  = (Q2 / (Q2 + m2o)) * (F2p + F2r) * G;
    R2  = (2 * (1 + R));

    aux = ME * RE / Q2;
    aux *=
        aux *
        (1 - v - component->GetAverageNucleonWeight() * x * v / (2 * energy) +
         (1 - 2 * particle_def_.mass * particle_def_.mass / Q2) * v * v *
             (1 + 4 * component->GetAverageNucleonWeight() * component->GetAverageNucleonWeight() * x * x / Q2) / R2);

    return (4 * PI * F2 / v) * aux;
}

// ------------------------------------------------------------------------- //
// Butkevich Mikhailov Integral
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double PhotoButkevichMikhailov::FunctionToQ2Integral(double energy, double v, double Q2)
{
    Components::Component* component = components_[component_index_];

    double x, aux, nu, G, F2, R2;

    nu  =   v*energy;
    x   =   Q2/(2*component->GetAverageNucleonWeight()*nu);

    G = shadow_effect_->CalculateShadowEffect(*component, x, nu);

    const double a  =   0.2513e6;
    const double b  =   0.6186e6;
    const double c  =   3.0292e6;
    const double d  =   1.4817e6;
    const double d0 =   0.0988;
    const double ar =   0.4056;
    const double t  =   1.8152;
    const double As =   0.12;
    const double Bu =   1.2437;
    const double Bd =   0.1853;
    const double R  =   0.25;

    double F2p, F2n, FSp, FNp, FSn, FNn, n, dl, xUv, xDv;

    n   =   1.5*(1 + Q2/(Q2 + c));
    dl  =   d0*(1 + 2*Q2/(Q2 + d));
    aux =   As*pow(x, -dl)*pow(Q2/(Q2 + a), 1 + dl);
    FSp =   aux*pow(1 - x, n + 4);
    FSn =   aux*pow(1 - x, n + t);
    aux =   pow(x, 1 - ar)*pow(1 - x, n)*pow(Q2/(Q2 + b), ar);
    xUv =   Bu*aux;
    xDv =   Bd*aux*(1 - x);
    FNp =   xUv + xDv;
    FNn =   xUv/4 + xDv*4;
    F2p =   FSp + FNp;
    F2n =   FSn + FNn;
    F2  =   G*(component->GetNucCharge()*F2p + (component->GetAtomicNum() - component->GetNucCharge())*F2n);
    R2  =   (2*(1 + R));


    aux =   ME*RE/Q2;
    aux *=  aux*(1 - v - component->GetAverageNucleonWeight()*x*v/(2*energy) +
                 (1 - 2*particle_def_.mass*particle_def_.mass/Q2)
                 *v*v*(1 + 4*component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight()*x*x/Q2)/R2);

    return (4*PI*F2/v)*aux;
}

// ------------------------------------------------------------------------- //
// Reno Sarcevic Su Integral
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double PhotoRenoSarcevicSu::FunctionToQ2Integral(double energy, double v, double Q2)
{
    Components::Component* component = components_[component_index_];

    double x, aux, nu;

    nu  =   v*energy;
    x   =   Q2/(2*component->GetAverageNucleonWeight()*nu);

    // -------------[ Evaluate shadowfactor ]---------------- //

    double a;

    if(component->GetNucCharge()==1)
    {
        a   =   1;
    }
    else
    {
        a = shadow_effect_->CalculateShadowEffect(*component, x, nu);
    }

    double P;

    aux =   x*x;
    P   =   1 - 1.85*x + 2.45*aux - 2.35*aux*x + aux*aux;

    // ---------[ Evaluate ALLM form factor F_2 ]------------ //
    //
    // F_2 = c_i(t) * x_i^{a_i(t)} * (1 - x)^{b_i(t)}; i = P,R
    // ------------------------------------------------------ //

    double cp1, cp2, cp3;
    double cr1, cr2, cr3;

    double ap1, ap2, ap3;
    double ar1, ar2, ar3;

    double bp1, bp2, bp3;
    double br1, br2, br3;

    double m2o, m2r, L2, m2p, Q2o;

    cp1     =   0.28067;
    cp2     =   0.22291;
    cp3     =   2.1979;

    cr1     =   0.80107;
    cr2     =   0.97307;
    cr3     =   3.4942;

    ap1     =   -0.0808;
    ap2     =   -0.44812;
    ap3     =   1.1709;

    ar1     =   0.58400;
    ar2     =   0.37888;
    ar3     =   2.6063;

    bp1     =   0.60243;
    bp2     =   1.3754;
    bp3     =   1.8439;

    br1     =   0.10711;
    br2     =   1.9386;
    br3     =   0.49338;

    m2o     =   0.31985;
    m2r     =   0.15052;
    L2      =   0.06527;
    m2p     =   49.457;
    Q2o     =   0.46017;


    // GeV -> MeV conversion
    m2o     *=  1e6;
    m2r     *=  1e6;
    L2      *=  1e6;
    m2p     *=  1e6;
    Q2o     *=  1e6;

    // these values are corrected according to the file f2allm.f from Halina Abramowicz
    bp1     *=  bp1;
    bp2     *=  bp2;
    br1     *=  br1;
    br2     *=  br2;
    Q2o     +=  L2;

    // R(x, Q^2) is approximated to 0
    const double R = 0;

    double cr, ar, cp, ap, br, bp, t;

    t   =   log(log((Q2 + Q2o)/L2)/log(Q2o/L2));

    if(t<0)
    {
        t=0;
    }

    cr  =   cr1 + cr2*pow(t, cr3);
    ar  =   ar1 + ar2*pow(t, ar3);
    cp  =   cp1 + (cp1 - cp2)*(1/(1 + pow(t, cp3)) - 1);
    ap  =   ap1 + (ap1 - ap2)*(1/(1 + pow(t, ap3)) - 1);
    br  =   br1 + br2*pow(t, br3);
    bp  =   bp1 + bp2*pow(t, bp3);

    double xp, xr, F2p, F2A, F2P, F2R, W2;

    W2  =   component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight()
            - Q2 + 2*component->GetAverageNucleonWeight()*energy*v;
    xp  =   (Q2 + m2p)/(Q2 + m2p + W2 - component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight());
    xr  =   (Q2 + m2r)/(Q2 + m2r + W2 - component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight());
    F2P =   cp*pow(xp, ap)*pow(1 - x, bp);
    F2R =   cr*pow(xr, ar)*pow(1 - x, br);
    F2p =   (Q2/(Q2 + m2o))*(F2P + F2R);
    F2A =   a*(component->GetNucCharge() + (component->GetAtomicNum() - component->GetNucCharge())) *P*F2p;

    // ---------[ Write together cross section ]------------- //

    aux =   ME*RE/Q2;
    aux *=  aux*(1 - v + 0.25*v*v -
                (1 + 4*particle_def_.mass*particle_def_.mass/Q2)*0.25*v*v *
                (1 + 4*component->GetAverageNucleonWeight()*component->GetAverageNucleonWeight()*x*x/Q2) /
                (1 + R));

    return (4*PI*F2A/v)*aux;
}

#undef Q2_PHOTO_PARAM_INTEGRAL_IMPL
