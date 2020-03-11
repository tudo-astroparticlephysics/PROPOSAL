
#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include <cmath>
#include <sstream>
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

Parametrization::Parametrization(const ParticleDef& particle_def,
                                 std::shared_ptr<const Medium> medium,
                                 double multiplier)
    : particle_mass_(particle_def.mass),
      particle_charge_(particle_def.charge),
      particle_low_(particle_def.low),
      medium_(medium),
      components_(medium_->GetComponents()),
      component_index_(0),
      multiplier_(multiplier) {}

Parametrization::Parametrization(const Parametrization& param)
    : particle_mass_(param.particle_mass_),
      particle_charge_(param.particle_charge_),
      particle_low_(param.particle_low_),
      medium_(param.medium_),
      components_(medium_->GetComponents()),
      component_index_(param.component_index_)  // //TODO(mario): Check better
                                                // way Mon 2017/09/04
      ,
      multiplier_(param.multiplier_) {}

Parametrization::~Parametrization() {
}

bool Parametrization::operator==(const Parametrization& parametrization) const {
    if (typeid(*this) != typeid(parametrization))
        return false;
    else
        return this->compare(parametrization);
}

bool Parametrization::operator!=(const Parametrization& parametrization) const {
    return !(*this == parametrization);
}

bool Parametrization::compare(const Parametrization& parametrization) const {
    if (particle_charge_ != parametrization.particle_charge_
        or particle_mass_ != parametrization.particle_mass_
        or particle_low_ != parametrization.particle_low_)
        return false;
    if (*medium_ != *parametrization.medium_)
        return false;
    else if (component_index_ != parametrization.component_index_)
        return false;
    else if (multiplier_ != parametrization.multiplier_)
        return false;
    else
        return true;
}

namespace PROPOSAL{
std::ostream& operator<<(std::ostream& os, Parametrization const& param) {
    std::stringstream ss;
    ss << " Parametrization (" << &param << ") ";
    os << Helper::Centered(80, ss.str()) << '\n';

    os << "name: " << param.GetName() << '\n';

    param.print(os);

    os << "multiplier: " << param.multiplier_ << '\n';
    os << "current component index: " << param.component_index_ << '\n';
    os << "particle_mass: "<< param.particle_mass_ << '\n';
    os << "particle_charge: "<< param.particle_charge_ << '\n';
    os << "particle_low: " << param.particle_low_ << '\n';
    os << *param.medium_ << '\n';
    os << Helper::Centered(80, "");
    return os;
}
} // namespace PROPOSAL

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

//----------------------------------------------------------------------------//
double Parametrization::FunctionToDNdxIntegral(double energy, double variable) {
    return DifferentialCrossSection(energy, variable);
}

// ------------------------------------------------------------------------- //
double Parametrization::FunctionToDEdxIntegral(double energy, double variable) {
    return variable * DifferentialCrossSection(energy, variable);
}

// ------------------------------------------------------------------------- //
double Parametrization::FunctionToDE2dxIntegral(double energy,
                                                double variable) {
    return variable * variable * DifferentialCrossSection(energy, variable);
}

// ------------------------------------------------------------------------- //
// Getter
// ------------------------------------------------------------------------- //

size_t Parametrization::GetHash() const {
    std::size_t seed = 0;
    hash_combine(seed, GetName(), std::abs(particle_charge_),
                 particle_mass_, medium_->GetName());

    return seed;
}
