/*! \file   Components.h
 *   \brief  Source file for the components routines.
 *
 *   For more details see the class documentation.
 *
 *   \date   Sat Aug  5 14:47:16 CEST 2017
 *   \author Mario Dunsch
 */

#include <functional>
#include <cmath>
#include <sstream>
#include <iomanip>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/math/Integral.h"
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;
using namespace PROPOSAL::Components;

#define COMPONENT_IMPL(cls, SYMBOL, NUCCHARGE, ATOMICNUM)                                                              \
    cls::cls(double atomInMolecule)                                                                                    \
        : Component(#SYMBOL, NUCCHARGE, ATOMICNUM, atomInMolecule)                                                     \
    {                                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    cls::cls(const cls& component)                                                                                     \
        : Component(component)                                                                                         \
    {                                                                                                                  \
    }

/******************************************************************************
 *                                  OStream                                    *
 ******************************************************************************/

namespace PROPOSAL {

namespace Components {

std::ostream& operator<<(std::ostream& os, Component const& component)
{
    std::stringstream ss;
    ss << " Component (" << &component << ") ";
    os << Helper::Centered(60, ss.str()) << '\n';

    os << std::fixed << std::setprecision(6);

    os << component.GetName() << std::endl;
    os << "AtomicNuc:"
       << "\t\t" << component.GetAtomicNum() << std::endl;
    os << "AtomInMolecule:"
       << "\t\t" << component.GetAtomInMolecule() << std::endl;
    os << "NucCharge:"
       << "\t\t" << component.GetNucCharge() << std::endl;
    os << "AverageNucleonWeight:"
       << "\t" << component.GetAverageNucleonWeight() << std::endl;

    os << Helper::Centered(60, "");
    return os;
}

} // namespace Components

} // namespace PROPOSAL

/******************************************************************************
 *                                  Component                                  *
 ******************************************************************************/

Component::Component(std::string name, double nucCharge, double atomicNum, double atomInMolecule)
    : name_(name)
    , nucCharge_(nucCharge)
    , atomicNum_(atomicNum)
    , atomInMolecule_(atomInMolecule)
    , logConstant_(0.0)
    , bPrime_(0.0)
    , M_(0.0)
    , mN_(0.0)
    , r0_(0.0)
{
    SetLogConstant();
    SetBPrime();

    M_ = (nucCharge_ * MP + (atomicNum_ - nucCharge_) * MN) / atomicNum_;

    if (nucCharge != 1.0)
    {
        Integral integral(IROMB, IMAXS, IPREC);

        r0_ = pow(atomicNum, 1.0 / 3.0);
        r0_ = 1.12 * r0_ - 0.86 / r0_;

        mN_ = 1.0 - 4.0 * PI * 0.17 *
                        integral.Integrate(r0_, -1.0, std::bind(&Component::FunctionToIntegral, this, std::placeholders::_1), 3, 2.0) /
                        atomicNum_;
    }
}

// ------------------------------------------------------------------------- //
// Operators & swap
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
bool Component::operator==(const Component& component) const
{
    if (typeid(*this) != typeid(component))
        return false;
    if (name_ != component.name_)
        return false;
    else if (nucCharge_ != component.nucCharge_)
        return false;
    else if (atomicNum_ != component.atomicNum_)
        return false;
    else if (atomInMolecule_ != component.atomInMolecule_)
        return false;
    else if (logConstant_ != component.logConstant_)
        return false;
    else if (bPrime_ != component.bPrime_)
        return false;
    else if (M_ != component.M_)
        return false;
    else if (mN_ != component.mN_)
        return false;
    else if (r0_ != component.r0_)
        return false;
    else
        return true;
}

bool Component::operator!=(const Component& component) const
{
    return !(*this == component);
}

// ------------------------------------------------------------------------- //
// Methods
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
void Component::SetLogConstant()
{
    int z = std::round(nucCharge_);
    switch (z)
    {
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
    switch (z)
    {
        case 1:
            bPrime_ = 446;
            break;
        default:
            bPrime_ = 1429;
    }
}

// ------------------------------------------------------------------------- //
double Component::FunctionToIntegral(double r)
{
    const double a = 0.54;

    return r * r / (1 + exp((r - r0_) / a));
}

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
