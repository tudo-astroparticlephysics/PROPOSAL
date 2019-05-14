
#include <cmath>

#include "PROPOSAL/crossection/parametrization/PhotoQ2Integration.h"

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/math/Interpolant.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

#define Q2_PHOTO_PARAM_INTEGRAL_IMPL(param)                                                                            \
    Photo##param::Photo##param(const ParticleDef& particle_def,                                                        \
                               const Medium& medium,                                                                   \
                               const EnergyCutSettings& cuts,                                                          \
                               double multiplier,                                                                      \
                               const ShadowEffect& shadow_effect)                                                      \
        : PhotoQ2Integral(particle_def, medium, cuts, multiplier, shadow_effect)                                       \
    {                                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    Photo##param::Photo##param(const Photo##param& photo)                                                              \
        : PhotoQ2Integral(photo)                                                                                       \
    {                                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    Photo##param::~Photo##param() {}                                                                                   \
                                                                                                                       \
    const std::string Photo##param::name_ = "Photo" #param;

/******************************************************************************
 *                            Photo Q2 Integration                            *
 ******************************************************************************/

PhotoQ2Integral::PhotoQ2Integral(const ParticleDef& particle_def,
                                 const Medium& medium,
                                 const EnergyCutSettings& cuts,
                                 double multiplier,
                                 const ShadowEffect& shadow_effect)
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

bool PhotoQ2Integral::compare(const Parametrization& parametrization) const
{
    const PhotoQ2Integral* photo = static_cast<const PhotoQ2Integral*>(&parametrization);

    if (*shadow_effect_ != *photo->shadow_effect_)
        return false;
    if (integral_ != photo->integral_)
        return false;
    else
        return Photonuclear::compare(parametrization);
}

// ------------------------------------------------------------------------- //
double PhotoQ2Integral::DifferentialCrossSection(double energy, double v)
{
    IntegralLimits limits = GetIntegralLimits(energy);

    double aux, q2_min, q2_max;

    q2_min = particle_def_.mass * v;
    q2_min *= q2_min / (1 - v);

    if (particle_def_.mass < MPI)
    {
        aux = particle_def_.mass * particle_def_.mass / energy;
        q2_min -= (aux * aux) / (2 * (1 - v));
    }

    q2_max = 2 * components_[component_index_]->GetAverageNucleonWeight() * energy * (v - limits.vMin);

    //  if(form==4) max=Math.min(max, 5.5e6);  // as requested in Butkevich and Mikheyev
    if (q2_min > q2_max)
    {
        return 0;
    }

    aux = integral_.Integrate(
        q2_min, q2_max, std::bind(&PhotoQ2Integral::FunctionToQ2Integral, this, energy, v, std::placeholders::_1), 4);

    aux *= medium_->GetMolDensity() * components_[component_index_]->GetAtomInMolecule() *
           particle_def_.charge * particle_def_.charge;

    return aux;
}

// ------------------------------------------------------------------------- //
// Getter
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
size_t PhotoQ2Integral::GetHash() const
{
    size_t seed = Parametrization::GetHash();
    hash_combine(seed, shadow_effect_->GetHash());

    return seed;
}

// ------------------------------------------------------------------------- //
// Print
// ------------------------------------------------------------------------- //

void PhotoQ2Integral::print(std::ostream& os) const
{
    os << "shadow effect: "
       << (dynamic_cast<ShadowDuttaRenoSarcevicSeckel*>(shadow_effect_) ? "DuttaRenoSarcevicSeckel"
                                                                        : "ButkevichMikhailov")
       << '\n';
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
// Abramowicz Levin Levy Maor 91
// Phys. Lett. B 269 (1991), 465
// ------------------------------------------------------------------------- //
double PhotoAbramowiczLevinLevyMaor91::FunctionToQ2Integral(double energy, double v, double Q2)
{
    Components::Component* component = components_[component_index_];

    double mass_nucleus = component->GetAverageNucleonWeight();

    // Bjorken x = \frac{Q^2}{2pq}
    double bjorken_x = Q2 / (2 * mass_nucleus * v * energy);

    // Parameter for Pomeron and Reggeon
    const double a1_reggeon = 0.60408;
    const double a2_reggeon = 0.17353;
    const double a3_reggeon = 1.61812;

    double b1_pomeron       = 0.49222;
    double b2_pomeron       = 0.52116;
    const double b3_pomeron = 3.55115;

    const double c1_pomeron = 0.26550;
    const double c2_pomeron = 0.04856;
    const double c3_pomeron = 1.04682;

    const double a1_pomeron = -0.04503;
    const double a2_pomeron = -0.36407;
    const double a3_pomeron = 8.17091;

    double b1_reggeon       = 1.26066;
    double b2_reggeon       = 1.83624;
    const double b3_reggeon = 0.81141;

    const double c1_reggeon = 0.67639;
    const double c2_reggeon = 0.49027;
    const double c3_reggeon = 2.66275;

    // parameter with conversion from GeV^2 to Mev^2
    const double mass_photon_eff = 0.30508 * 1e6;
    const double mass_pomeron    = 10.67564 * 1e6;
    const double mass_reggeon    = 0.20623 * 1e6;
    const double scaleParameter  = 0.06527 * 1e6;
    double Q20_free_param        = 0.27799 * 1e6;

    // these values are corrected according to the file f2allm.f from Halina Abramowicz
    // TODO: do this values really have to be corrected
    //       or is the correection just for ALLM97
    b1_pomeron *= b1_pomeron;
    b2_pomeron *= b2_pomeron;
    b1_reggeon *= b1_reggeon;
    b2_reggeon *= b2_reggeon;
    Q20_free_param += scaleParameter;

    // R(x, Q^2) is approximated to 0
    // relation between structure functions F_1 and F_2
    const double R = 0;

    // t = log( \frac{log \frac{Q^2 + Q_0^2}{\Lambda^2}}{log \frac{Q_0^2}{\Lambda^2}} )
    // eq. 22
    double t = std::log(std::log((Q2 + Q20_free_param) / scaleParameter) / std::log(Q20_free_param / scaleParameter));

    if (t < 0)
        t = 0;

    // parameter, that increase with Q^2
    // eq. 23
    // f(t) = f_1 + f_2 t^{f_3}
    double a_reggeon = a1_reggeon + a2_reggeon * std::pow(t, a3_reggeon);
    double b_reggeon = b1_reggeon + b2_reggeon * std::pow(t, b3_reggeon);
    double c_reggeon = c1_reggeon + c2_reggeon * std::pow(t, c3_reggeon);
    double b_pomeron = b1_pomeron + b2_pomeron * std::pow(t, b3_pomeron);
    // parameter, that decrease with Q^2
    // eq. 24
    // f'(t) = f_1 + (f_1 - f_2) (\frac{1}{1 + t^{f_3}} - 1)
    double a_pomeron = a1_pomeron + (a1_pomeron - a2_pomeron) * (1 / (1 + std::pow(t, a3_pomeron)) - 1);
    double c_pomeron = c1_pomeron + (c1_pomeron - c2_pomeron) * (1 / (1 + std::pow(t, c3_pomeron)) - 1);

    // invariant mass of nucleus and virtual photon
    // W^2 = (p + q)^2 = M^2 + 2MEv - Q^2
    double W2 = mass_nucleus * mass_nucleus + 2 * mass_nucleus * energy * v - Q2;

    // Relation between structure function of proton and structure function of neutron
    // from the BCDMS Collaboration
    // Phys. Lett. B 237 (1990), 599
    // eq. 4
    // P(x) = 1 - 1.85x + 2.45x^2 - 2.35x^3 + x^4
    double relation_proton_neutron = bjorken_x * bjorken_x;
    relation_proton_neutron        = 1 - 1.85 * bjorken_x + 2.45 * relation_proton_neutron -
                              2.35 * relation_proton_neutron * bjorken_x +
                              relation_proton_neutron * relation_proton_neutron;

    // eq. 17 and 18
    // x_i = \frac{Q^2 + m_i}{Q^2 + m_i + W^2 - M^2}
    double bjorken_x_pomeron = (Q2 + mass_pomeron) / (Q2 + mass_pomeron + W2 - mass_nucleus * mass_nucleus);
    double bjorken_x_reggeon = (Q2 + mass_reggeon) / (Q2 + mass_reggeon + W2 - mass_nucleus * mass_nucleus);
    // eq. 20 and 21
    // F_{2, i}(x, Q^2) = c_i(t) x^{a_i(t)} (1 - x)^{b_i(t)}
    double pomeron_contribution = c_pomeron * std::pow(bjorken_x_pomeron, a_pomeron) * std::pow(1 - bjorken_x, b_pomeron);
    double reggeon_controbution = c_reggeon * std::pow(bjorken_x_reggeon, a_reggeon) * std::pow(1 - bjorken_x, b_reggeon);
    // ALLM97 eq. 1
    // F_{2, proton}(x, Q^2) = \frac{Q^2}{Q^2 + m_0^2} (F_{2, Pomeron} + F_{2, Reggeon})
    double structure_function_proton = Q2 / (Q2 + mass_photon_eff) * (pomeron_contribution + reggeon_controbution);
    // structure function from Dutta et al
    // Phys. Rev. D 63 (2001), 094020
    // eq. 3.11
    // F_{2, nucleus} = G(x) (Z + (A - Z)P(x)) F_{2, Proton}
    double structure_function_nucleus =
        structure_function_proton * shadow_effect_->CalculateShadowEffect(*component, bjorken_x, v * energy) *
        (component->GetNucCharge() + (component->GetAtomicNum() - component->GetNucCharge()) * relation_proton_neutron);

    // differential cross section from Dutta et al.
    // Phys. Rev. D 63 (2001), 094020
    // eq. 3.11
    // eq. 3.4
    // \frac{d^2 \sigma}{dQ^2 dx} = \frac{4\pi\alpha^2}{Q^4} \frac{F_{2, Proton}}{v}
    //     ( 1 - v + \frac{M x v}{2E}
    //     - \frac{v^2}{2} (1 - \frac{2m_{particle}}{Q^2})
    //          \frac{1 + \frac{4 M^2 x^2}{Q^2}}{1 + R})
    double result = ME * RE / Q2;
    result *= result * 4 * PI * structure_function_nucleus / v *
              (1 - v - mass_nucleus * bjorken_x * v / (2 * energy) +
               (1 - 2 * particle_def_.mass * particle_def_.mass / Q2) * v * v *
                   (1 + 4 * mass_nucleus * mass_nucleus * bjorken_x * bjorken_x / Q2) / (2 * (1 + R)));

    return result;
}

// ------------------------------------------------------------------------- //
// Abramowicz Levin Levy Maor 97
// arXiv:hep-ph/9712415
// ------------------------------------------------------------------------- //
double PhotoAbramowiczLevinLevyMaor97::FunctionToQ2Integral(double energy, double v, double Q2)
{
    Components::Component* component = components_[component_index_];

    double mass_nucleus = component->GetAverageNucleonWeight();

    // Bjorken x = \frac{Q^2}{2pq}
    double bjorken_x = Q2 / (2 * mass_nucleus * v * energy);

    // --------------------------------------------------------------------- //
    // Evaluate form factor F_2 for the nucleus same like ALLM91
    // but using new parameters
    // the parameters were updated on 19.08.2004 on arXiv
    // because the values of b1_pomeron, b2_pomeron, b1_reggeon, b2_reggeon
    // and the Q20_free_param were misprinted
    // the corrected parameters are used here
    // --------------------------------------------------------------------- //

    // Parameter for Pomeron and Reggeon
    const double a1_pomeron = -0.0808;
    const double a2_pomeron = -0.44812;
    const double a3_pomeron = 1.1709;

    const double b1_pomeron = 0.36292;
    const double b2_pomeron = 1.8917;
    const double b3_pomeron = 1.8439;

    const double c1_pomeron = 0.28067;
    const double c2_pomeron = 0.22291;
    const double c3_pomeron = 2.1979;

    const double a1_reggeon = 0.58400;
    const double a2_reggeon = 0.37888;
    const double a3_reggeon = 2.6063;

    const double b1_reggeon = 0.01147;
    const double b2_reggeon = 3.7582;
    const double b3_reggeon = 0.49338;

    const double c1_reggeon = 0.80107;
    const double c2_reggeon = 0.97307;
    const double c3_reggeon = 3.4942;

    // parameter with conversion from GeV^2 to Mev^2
    const double mass_photon_eff = 0.31985 * 1e6;
    const double mass_reggeon    = 0.15052 * 1e6;
    const double mass_pomeron    = 49.457 * 1e6;
    const double scaleParameter  = 0.06527 * 1e6;
    const double Q20_free_param  = 0.52544 * 1e6;

    // R(x, Q^2) is approximated to 0
    // relation between structure functions F_1 and F_2
    const double R = 0;

    // t = log( \frac{log \frac{Q^2 + Q_0^2}{\Lambda^2}}{log \frac{Q_0^2}{\Lambda^2}} )
    // eq. 22
    double t = std::log(std::log((Q2 + Q20_free_param) / scaleParameter) / std::log(Q20_free_param / scaleParameter));

    if (t < 0)
        t = 0;

    // parameter, that increase with Q^2
    // eq. 23
    // f(t) = f_1 + f_2 t^{f_3}
    double a_reggeon = a1_reggeon + a2_reggeon * std::pow(t, a3_reggeon);
    double b_reggeon = b1_reggeon + b2_reggeon * std::pow(t, b3_reggeon);
    double c_reggeon = c1_reggeon + c2_reggeon * std::pow(t, c3_reggeon);
    double b_pomeron = b1_pomeron + b2_pomeron * std::pow(t, b3_pomeron);
    // parameter, that decrease with Q^2
    // eq. 24
    // f'(t) = f_1 + (f_1 - f_2) (\frac{1}{1 + t^{f_3}} - 1)
    double a_pomeron = a1_pomeron + (a1_pomeron - a2_pomeron) * (1 / (1 + std::pow(t, a3_pomeron)) - 1);
    double c_pomeron = c1_pomeron + (c1_pomeron - c2_pomeron) * (1 / (1 + std::pow(t, c3_pomeron)) - 1);

    // invariant mass of nucleus and virtual photon
    // W^2 = (p + q)^2 = M^2 + 2MEv - Q^2
    double W2 = mass_nucleus * mass_nucleus + 2 * mass_nucleus * energy * v - Q2;

    // Relation between structure function of proton and structure function of neutron
    // from the BCDMS Collaboration
    // Phys. Lett. B 237 (1990), 599
    // eq. 4
    // P(x) = 1 - 1.85x + 2.45x^2 - 2.35x^3 + x^4
    double relation_proton_neutron = bjorken_x * bjorken_x;
    relation_proton_neutron        = 1 - 1.85 * bjorken_x + 2.45 * relation_proton_neutron -
                              2.35 * relation_proton_neutron * bjorken_x +
                              relation_proton_neutron * relation_proton_neutron;

    // eq. 17 and 18
    // x_i = \frac{Q^2 + m_i}{Q^2 + m_i + W^2 - M^2}
    double bjorken_x_pomeron = (Q2 + mass_pomeron) / (Q2 + mass_pomeron + W2 - mass_nucleus * mass_nucleus);
    double bjorken_x_reggeon = (Q2 + mass_reggeon) / (Q2 + mass_reggeon + W2 - mass_nucleus * mass_nucleus);
    // eq. 20 and 21
    // F_{2, i}(x, Q^2) = c_i(t) x^{a_i(t)} (1 - x)^{b_i(t)}
    double pomeron_contribution = c_pomeron * std::pow(bjorken_x_pomeron, a_pomeron) * std::pow(1 - bjorken_x, b_pomeron);
    double reggeon_controbution = c_reggeon * std::pow(bjorken_x_reggeon, a_reggeon) * std::pow(1 - bjorken_x, b_reggeon);
    // ALLM97 eq. 1
    // F_{2, proton}(x, Q^2) = \frac{Q^2}{Q^2 + m_0^2} (F_{2, Pomeron} + F_{2, Reggeon})
    double structure_function_proton = Q2 / (Q2 + mass_photon_eff) * (pomeron_contribution + reggeon_controbution);
    // structure function from Dutta et al
    // Phys. Rev. D 63 (2001), 094020
    // eq. 3.11
    // F_{2, nucleus} = G(x) (Z + (A - Z)P(x)) F_{2, Proton}
    double structure_function_nucleus =
        structure_function_proton * shadow_effect_->CalculateShadowEffect(*component, bjorken_x, v * energy) *
        (component->GetNucCharge() + (component->GetAtomicNum() - component->GetNucCharge()) * relation_proton_neutron);

    // differential cross section from Dutta et al.
    // Phys. Rev. D 63 (2001), 094020
    // eq. 3.11
    // eq. 3.4
    // \frac{d^2 \sigma}{dQ^2 dx} = \frac{4\pi\alpha^2}{Q^4} \frac{F_{2, Proton}}{v}
    //     ( 1 - v + \frac{M x v}{2E}
    //     - \frac{v^2}{2} (1 - \frac{2m_{particle}}{Q^2})
    //          \frac{1 + \frac{4 M^2 x^2}{Q^2}}{1 + R})
    double result = ME * RE / Q2;
    result *= result * 4 * PI * structure_function_nucleus / v *
              (1 - v - mass_nucleus * bjorken_x * v / (2 * energy) +
               (1 - 2 * particle_def_.mass * particle_def_.mass / Q2) * v * v *
                   (1 + 4 * mass_nucleus * mass_nucleus * bjorken_x * bjorken_x / Q2) / (2 * (1 + R)));

    return result;
}

// ------------------------------------------------------------------------- //
// Butkevich Mikheyev Parametrization
// JETP 95 (2002), 11
//
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
double PhotoButkevichMikhailov::FunctionToQ2Integral(double energy, double v, double Q2)
{
    Components::Component* component = components_[component_index_];

    double mass_nucleus = component->GetAverageNucleonWeight();

    // Bjorken x = \frac{Q^2}{2pq}
    double bjorken_x = Q2 / (2 * mass_nucleus * v * energy);

    const double a                 = 0.2513e6;
    const double b                 = 0.6186e6;
    const double c                 = 3.0292e6;
    const double d                 = 1.4817e6;
    const double intercept_pomeron = 0.0988;
    const double intercept_reggeon = 0.4056;
    const double tau               = 1.8152;
    const double A_s               = 0.12;
    const double B_up              = 1.2437;
    const double B_down            = 0.1853;

    // R(x, Q^2) is approximated to 0
    // relation between structure functions F_1 and F_2
    const double R = 0.25;

    double aux;

    // Contribution of Seaquarks and Gluons (called Singlet Term)
    // eq. 26
    // n(Q^2) = 1.5 (1 + \frac{Q^2}{Q^2 + c})
    double n = 1.5 * (1 + Q2 / (Q2 + c));
    // eq. 25
    // \Delta(Q^2) = \Delta_0 (1 + \frac{2 Q^2}{Q^2 + d})
    double dl = intercept_pomeron * (1 + 2 * Q2 / (Q2 + d));
    // eq. 24
    // F_{2, Proton} = A_s x^{-\Delta} (1 - x)^{n + 4} (\frac{Q^2}{Q^2 + a})^{1 + \Delta}
    aux                     = A_s * std::pow(bjorken_x, -dl) * std::pow(Q2 / (Q2 + a), 1 + dl);
    double F_proton_singlet = aux * std::pow(1 - bjorken_x, n + 4);
    // eq. 37
    // F_{2, Neutron} = A_s x^{-\Delta} (1 - x)^{n + \tau} (\frac{Q^2}{Q^2 + a})^{1 + \Delta}
    double F_neutron_singlet = aux * std::pow(1 - bjorken_x, n + tau);
    // Contribution of Valence quarks (called non-Singlet Term)
    // splitted into Up quark and down quark contributions
    // eq. 29
    // xU_v = B_{up} x^{1 - \alpha_{Reggeon}} (1 - x)^{n} (\frac{Q^2}{Q^2 + b})^{\alpha_{Reggeon}}
    aux = std::pow(bjorken_x, 1 - intercept_reggeon) * std::pow(1 - bjorken_x, n) * std::pow(Q2 / (Q2 + b), intercept_reggeon);
    double Up_valence = B_up * aux;
    // eq. 30
    // xD_v = B_{down} x^{1 - \alpha_{Reggeon}} (1 - x)^{n + 1} (\frac{Q^2}{Q^2 + b})^{\alpha_{Reggeon}}
    double Down_valence = B_down * aux * (1 - bjorken_x);
    // eq. 28
    // F_{2, Proton, non-Singlet} = xU_v + xD_v
    double F_proton_non_singlet = Up_valence + Down_valence;
    // eq. 36
    // F_{2, Neutron, non-Singlet} = \frac{1}{4} xU_v + 4 xD_v
    double F_neutron_non_singlet = Up_valence / 4 + Down_valence * 4;
    // eq. 23
    // F_{2, i} = F_{2, i, Singlet} + F_{2, i, non-Singlet}
    double structure_function_proton  = F_proton_singlet + F_proton_non_singlet;
    double structure_function_neutron = F_neutron_singlet + F_neutron_non_singlet;
    // F_{2, nucleus} = G (Z F_{2, Proton} + (A-Z) F_{2, Neutron})
    double structure_function_nucleus =
        shadow_effect_->CalculateShadowEffect(*component, bjorken_x, v * energy) *
        (component->GetNucCharge() * structure_function_proton +
         (component->GetAtomicNum() - component->GetNucCharge()) * structure_function_neutron);

    // differential cross section from Dutta et al.
    // Phys. Rev. D 63 (2001), 094020
    // eq. 3.11
    // eq. 3.4
    // \frac{d^2 \sigma}{dQ^2 dx} = \frac{4\pi\alpha^2}{Q^4} \frac{F_{2, Proton}}{v}
    //     ( 1 - v + \frac{M x v}{2E}
    //     - \frac{v^2}{2} (1 - \frac{2m_{particle}}{Q^2})
    //          \frac{1 + \frac{4 M^2 x^2}{Q^2}}{1 + R})
    double result = ME * RE / Q2;
    result *= result * 4 * PI * structure_function_nucleus / v *
              (1 - v - mass_nucleus * bjorken_x * v / (2 * energy) +
               (1 - 2 * particle_def_.mass * particle_def_.mass / Q2) * v * v *
                   (1 + 4 * mass_nucleus * mass_nucleus * bjorken_x * bjorken_x / Q2) / (2 * (1 + R)));

    return result;
}

// ------------------------------------------------------------------------- //
// Reno Sarcevic Su Integral
// Astrop. Phys. 24 (2005), 107
// this parametrization was calculated for sTaus with spin 0
// the other parametrizations are for charged leptons with spin 1/2
// ------------------------------------------------------------------------- //
double PhotoRenoSarcevicSu::FunctionToQ2Integral(double energy, double v, double Q2)
{
    Components::Component* component = components_[component_index_];

    double mass_nucleus = component->GetAverageNucleonWeight();

    // Bjorken x = \frac{Q^2}{2pq}
    double bjorken_x = Q2 / (2 * mass_nucleus * v * energy);

    // --------------------------------------------------------------------- //
    // Evaluate form factor F_2 for the nucleus same like ALLM91
    // but using new parameters from ALLM97
    // the parameters were updated on 19.08.2004 on arXiv
    // because the values of b1_pomeron, b2_pomeron, b1_reggeon, b2_reggeon
    // and the Q20_free_param were misprinted
    // the corrected parameters are used here
    // --------------------------------------------------------------------- //

    // Parameter for Pomeron and Reggeon
    const double a1_pomeron = -0.0808;
    const double a2_pomeron = -0.44812;
    const double a3_pomeron = 1.1709;

    const double b1_pomeron = 0.36292;
    const double b2_pomeron = 1.8917;
    const double b3_pomeron = 1.8439;

    const double c1_pomeron = 0.28067;
    const double c2_pomeron = 0.22291;
    const double c3_pomeron = 2.1979;

    const double a1_reggeon = 0.58400;
    const double a2_reggeon = 0.37888;
    const double a3_reggeon = 2.6063;

    const double b1_reggeon = 0.01147;
    const double b2_reggeon = 3.7582;
    const double b3_reggeon = 0.49338;

    const double c1_reggeon = 0.80107;
    const double c2_reggeon = 0.97307;
    const double c3_reggeon = 3.4942;

    // parameter with conversion from GeV^2 to Mev^2
    const double mass_photon_eff = 0.31985 * 1e6;
    const double mass_reggeon    = 0.15052 * 1e6;
    const double mass_pomeron    = 49.457 * 1e6;
    const double scaleParameter  = 0.06527 * 1e6;
    const double Q20_free_param  = 0.52544 * 1e6;

    // R(x, Q^2) is approximated to 0
    // relation between structure functions F_1 and F_2
    const double R = 0;

    // t = log( \frac{log \frac{Q^2 + Q_0^2}{\Lambda^2}}{log \frac{Q_0^2}{\Lambda^2}} )
    // eq. 22
    double t = std::log(std::log((Q2 + Q20_free_param) / scaleParameter) / std::log(Q20_free_param / scaleParameter));

    if (t < 0)
        t = 0;

    // parameter, that increase with Q^2
    // eq. 23
    // f(t) = f_1 + f_2 t^{f_3}
    double a_reggeon = a1_reggeon + a2_reggeon * std::pow(t, a3_reggeon);
    double b_reggeon = b1_reggeon + b2_reggeon * std::pow(t, b3_reggeon);
    double c_reggeon = c1_reggeon + c2_reggeon * std::pow(t, c3_reggeon);
    double b_pomeron = b1_pomeron + b2_pomeron * std::pow(t, b3_pomeron);
    // parameter, that decrease with Q^2
    // eq. 24
    // f'(t) = f_1 + (f_1 - f_2) (\frac{1}{1 + t^{f_3}} - 1)
    double a_pomeron = a1_pomeron + (a1_pomeron - a2_pomeron) * (1 / (1 + std::pow(t, a3_pomeron)) - 1);
    double c_pomeron = c1_pomeron + (c1_pomeron - c2_pomeron) * (1 / (1 + std::pow(t, c3_pomeron)) - 1);

    // invariant mass of nucleus and virtual photon
    // W^2 = (p + q)^2 = M^2 + 2MEv - Q^2
    double W2 = mass_nucleus * mass_nucleus + 2 * mass_nucleus * energy * v - Q2;

    // Relation between structure function of proton and structure function of neutron
    // from the BCDMS Collaboration
    // Phys. Lett. B 237 (1990), 599
    // eq. 4
    // P(x) = 1 - 1.85x + 2.45x^2 - 2.35x^3 + x^4
    double relation_proton_neutron = bjorken_x * bjorken_x;
    relation_proton_neutron        = 1 - 1.85 * bjorken_x + 2.45 * relation_proton_neutron -
                              2.35 * relation_proton_neutron * bjorken_x +
                              relation_proton_neutron * relation_proton_neutron;

    // eq. 17 and 18
    // x_i = \frac{Q^2 + m_i}{Q^2 + m_i + W^2 - M^2}
    double bjorken_x_pomeron = (Q2 + mass_pomeron) / (Q2 + mass_pomeron + W2 - mass_nucleus * mass_nucleus);
    double bjorken_x_reggeon = (Q2 + mass_reggeon) / (Q2 + mass_reggeon + W2 - mass_nucleus * mass_nucleus);
    // eq. 20 and 21
    // F_{2, i}(x, Q^2) = c_i(t) x^{a_i(t)} (1 - x)^{b_i(t)}
    double pomeron_contribution = c_pomeron * std::pow(bjorken_x_pomeron, a_pomeron) * std::pow(1 - bjorken_x, b_pomeron);
    double reggeon_controbution = c_reggeon * std::pow(bjorken_x_reggeon, a_reggeon) * std::pow(1 - bjorken_x, b_reggeon);
    // ALLM97 eq. 1
    // F_{2, proton}(x, Q^2) = \frac{Q^2}{Q^2 + m_0^2} (F_{2, Pomeron} + F_{2, Reggeon})
    double structure_function_proton = Q2 / (Q2 + mass_photon_eff) * (pomeron_contribution + reggeon_controbution);
    // structure function from Dutta et al
    // Phys. Rev. D 63 (2001), 094020
    // eq. 3.11
    // F_{2, nucleus} = G(x) (Z + (A - Z)P(x)) F_{2, Proton}
    double structure_function_nucleus =
        structure_function_proton * shadow_effect_->CalculateShadowEffect(*component, bjorken_x, v * energy) *
        (component->GetNucCharge() + (component->GetAtomicNum() - component->GetNucCharge()) * relation_proton_neutron);

    // --------------------------------------------------------------------- //
    // this is the only part, that differs from ALLM parametrization
    // or from the parametrization considering charged leptons with spin 1/2
    // eq. 6 and 8 from Reno et al
    // --------------------------------------------------------------------- //
    // \frac{d^2 \sigma}{dQ^2 dx} = \frac{4\pi\alpha^2}{Q^4} \frac{F_{2, Proton}}{v}
    //     ( 1 - v + \frac{v^2}{4}
    //     - \frac{v^2}{4} (1 - \frac{4m_{particle}}{Q^2})
    //          \frac{1 + \frac{4M^2x^2}{Q^2}}{1 + R})
    double result = ME * RE / Q2;
    result *= result * 4 * PI * structure_function_nucleus / v *
              (1 - v + 0.25 * v * v -
               (1 + 4 * particle_def_.mass * particle_def_.mass / Q2) * 0.25 * v * v *
                   (1 + 4 * mass_nucleus * mass_nucleus * bjorken_x * bjorken_x / Q2) / (1 + R));

    return result;
}

#undef Q2_PHOTO_PARAM_INTEGRAL_IMPL
