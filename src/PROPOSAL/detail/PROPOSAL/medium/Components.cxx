/*! \file   Components.h
 *   \brief  Source file for the components routines.
 *
 *   For more details see the class documentation.
 *
 *   \date   Sat Aug  5 14:47:16 CEST 2017
 *   \author Mario Dunsch
 */

#include <cmath>
#include <functional>
#include <iomanip>
#include <sstream>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;
using namespace PROPOSAL::Components;

std::unique_ptr<std::map<size_t, Component>> Component::component_map = nullptr;

#define COMPONENT_IMPL(cls, SYMBOL, NUCCHARGE, ATOMICNUM)                      \
    cls::cls(double atomInMolecule)                                            \
        : Component(#SYMBOL, NUCCHARGE, ATOMICNUM, atomInMolecule)             \
    {                                                                          \
    }

/******************************************************************************
 *                                  Component                                  *
 ******************************************************************************/

Component::Component(
    std::string name, double nucCharge, double atomicNum, double atomInMolecule)
    : name_(name)
    , nucCharge_(nucCharge)
    , atomicNum_(atomicNum)
    , atomInMolecule_(atomInMolecule)
{
    SetLogConstant();
    SetBPrime();

    averageNucleonWeight_
        = (nucCharge_ * MP + (atomicNum_ - nucCharge_) * MN) / atomicNum_;

    if (nucCharge != 1.0) {
        // see Butkevich, Mikheyev JETP 95 (2002), 11 eq. 45-47
        auto r0 = std::pow(atomicNum, 1.0 / 3.0);
        r0 = 1.12 * r0 - 0.86 / r0;

        wood_saxon_
            = 1.0 - 4.0 * PI * 0.17 * WoodSaxonPotential(r0) / atomicNum_;
    }

    hash = 0;
    hash_combine(hash, nucCharge_, atomicNum_, atomInMolecule_,
                 logConstant_, bPrime_, averageNucleonWeight_, wood_saxon_);

    if (!component_map)
        component_map = std::make_unique<std::map<size_t, Component>>();

    if (component_map->find(hash) == component_map->end())
        component_map->insert({hash, Component(*this)});

}

namespace PROPOSAL {
bool operator==(Component const& lhs, Component const& rhs) noexcept
{
    if (lhs.name_ != rhs.name_)
        return false;
    else if (lhs.nucCharge_ != rhs.nucCharge_)
        return false;
    else if (lhs.atomicNum_ != rhs.atomicNum_)
        return false;
    else if (lhs.atomInMolecule_ != rhs.atomInMolecule_)
        return false;
    else if (lhs.logConstant_ != rhs.logConstant_)
        return false;
    else if (lhs.bPrime_ != rhs.bPrime_)
        return false;
    else if (lhs.averageNucleonWeight_ != rhs.averageNucleonWeight_)
        return false;
    else if (lhs.wood_saxon_ != rhs.wood_saxon_)
        return false;
    else if (lhs.hash != rhs.hash)
        return false;
    return true;
}

std::ostream& operator<<(std::ostream& os, Component const& comp) noexcept
{
    std::stringstream ss;
    ss << " Component (" << &comp << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';
    os << comp.GetName() << std::endl;
    os << "AtomicNuc:"
       << "\t\t" << comp.GetAtomicNum() << std::endl;
    os << "AtomInMolecule:"
       << "\t\t" << comp.GetAtomInMolecule() << std::endl;
    os << "NucCharge:"
       << "\t\t" << comp.GetNucCharge() << std::endl;
    os << "AverageNucleonWeight:"
       << "\t" << comp.GetAverageNucleonWeight() << std::endl;
    os << "Hash:"
       << "\t" << comp.GetHash() << std::endl;
    os << Helper::Centered(60, "");
    return os;
}
}

size_t Component::GetHash() const noexcept
{
    return hash;
}

// ------------------------------------------------------------------------- //
// Methods
// ------------------------------------------------------------------------- //

bool Component::operator!=(const Component& component) const
{
    return !(*this == component);
}

// ------------------------------------------------------------------------- //
void Component::SetLogConstant()
{
    int z = std::round(nucCharge_);
    switch (z) {
    case 1:
        logConstant_ = 202.4;
        break;
    case 2:
        logConstant_ = 151.9;
        break;
    case 3:
        logConstant_ = 159.9;
        break;
    case 4:
        logConstant_ = 172.3;
        break;
    case 5:
        logConstant_ = 177.9;
        break;
    case 6:
        logConstant_ = 178.3;
        break;
    case 7:
        logConstant_ = 176.6;
        break;
    case 8:
        logConstant_ = 173.4;
        break;
    case 9:
        logConstant_ = 170.0;
        break;
    case 10:
        logConstant_ = 165.8;
        break;
    case 11:
        logConstant_ = 165.8;
        break;
    case 12:
        logConstant_ = 167.1;
        break;
    case 13:
        logConstant_ = 169.1;
        break;
    case 14:
        logConstant_ = 170.8;
        break;
    case 15:
        logConstant_ = 172.2;
        break;
    case 16:
        logConstant_ = 173.4;
        break;
    case 17:
        logConstant_ = 174.3;
        break;
    case 18:
        logConstant_ = 174.8;
        break;
    case 19:
        logConstant_ = 175.1;
        break;
    case 20:
        logConstant_ = 175.6;
        break;
    case 21:
        logConstant_ = 176.2;
        break;
    case 22:
        logConstant_ = 176.8;
        break;
    case 26:
        logConstant_ = 175.8;
        break;
    case 29:
        logConstant_ = 173.1;
        break;
    case 32:
        logConstant_ = 173.0;
        break;
    case 35:
        logConstant_ = 173.5;
        break;
    case 42:
        logConstant_ = 175.9;
        break;
    case 50:
        logConstant_ = 177.4;
        break;
    case 53:
        logConstant_ = 178.6;
        break;
    case 74:
        logConstant_ = 177.6;
        break;
    case 82:
        logConstant_ = 178.0;
        break;
    case 92:
        logConstant_ = 179.8;
        break;
    default:
        logConstant_ = 182.7;
    }
}

// ------------------------------------------------------------------------- //
void Component::SetBPrime()
{
    int z = std::round(nucCharge_);
    switch (z) {
    case 1:
        bPrime_ = 446;
        break;
    default:
        bPrime_ = 1429;
    }
}

// ------------------------------------------------------------------------- //
double Component::WoodSaxonPotential(double r0)
{
    // This is the analytical intergral to
    // $\int_{r_0}^{\infty} \frac{r^2}{1+\exp((r-r_0)/a)}dr$
    // with $a=0.54$
    // which can be written to
    // $ar_0^2 \int_0^{\infty} \frac{dx}{1+e^x}$
    // $ + 2a^2r_0 \int_0^{\infty} \frac{x dx}{1+e^x}$
    // $ + a^3 \int_0^{\infty} \frac{x^2 dx}{1+e^x}$
    // and results in
    // $ar_0^2\log(2) + 2a^2r_0\pi^2/12 + 3/2a^3\zeta(3)$
    const double a = 0.54;
    return a
        * (r0 * r0 * std::log(2) + a * r0 * PI * PI / 6 + a * a * 1.5 * ZETA3);
}

namespace PROPOSAL {
double calculate_proton_massnumber_fraction(
    const component_list& comp_list) noexcept
{
    double charge = 0.;
    double nucleons = 0.;
    for (const auto& comp : comp_list) {
        charge += comp.GetAtomInMolecule() * comp.GetNucCharge();
        nucleons += comp.GetAtomInMolecule() * comp.GetAtomicNum();
    }
    return charge / nucleons;
}
} // namespace PROPOSAL

/******************************************************************************
 *                            Different Components                             *
 ******************************************************************************/

COMPONENT_IMPL(Hydrogen, H, 1, 1.00794)
COMPONENT_IMPL(Carbon, C, 6, 12.0011)
COMPONENT_IMPL(Nitrogen, N, 7, 14.0067)
COMPONENT_IMPL(Oxygen, O, 8, 15.9994)
COMPONENT_IMPL(Sodium, Na, 11, 22.989770)
COMPONENT_IMPL(Magnesium, Mg, 12, 24.31)
COMPONENT_IMPL(Sulfur, S, 16, 32.07)
COMPONENT_IMPL(Chlorine, Cl, 17, 35.4527)
COMPONENT_IMPL(Argon, Ar, 18, 39.948)
COMPONENT_IMPL(Potassium, K, 19, 39.10)
COMPONENT_IMPL(Calcium, Ca, 20, 40.08)
COMPONENT_IMPL(Iron, Fe, 26, 55.845)
COMPONENT_IMPL(Copper, Cu, 29, 63.546)
COMPONENT_IMPL(Lead, Pb, 82, 207.2)
COMPONENT_IMPL(Uranium, U, 92, 238.0289)
COMPONENT_IMPL(StandardRock, StandardRock, 11, 22.0)
COMPONENT_IMPL(FrejusRock, FrejusRock, 10.12, 20.34)
