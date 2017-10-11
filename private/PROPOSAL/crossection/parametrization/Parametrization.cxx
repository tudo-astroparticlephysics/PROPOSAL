
#include <boost/functional/hash.hpp>

#include "PROPOSAL/crossection/parametrization/Parametrization.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/medium/Medium.h"


using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

Parametrization::Parametrization(const ParticleDef& particle_def,
                                 const Medium& medium,
                                 const EnergyCutSettings& cuts,
                                 double multiplier)
    : particle_def_(particle_def)
    , medium_(medium.clone())
    , cut_settings_(cuts)
    , components_(medium_->GetComponents())
    , component_index_(0)
    , multiplier_(multiplier)
{
}

Parametrization::Parametrization(const Parametrization& param)
    : particle_def_(param.particle_def_)
    , medium_(param.medium_->clone())
    , cut_settings_(param.cut_settings_)
    , components_(medium_->GetComponents())
    , component_index_(param.component_index_) // //TODO(mario): Check better way Mon 2017/09/04
    , multiplier_(1.0)
{
}

Parametrization::~Parametrization()
{
    delete medium_;
}

std::ostream& PROPOSAL::operator<<(std::ostream& os, Parametrization const& param)
{
    std::stringstream ss;
    ss << " Parametrization (" << &param << ") ";
    os << Helper::Centered(80, ss.str()) << '\n';

    os << "name: " << param.GetName() << '\n';

    param.print(os);

    os << "multiplier: " << param.multiplier_ << '\n';
    os << "current component index: " << param.component_index_ << '\n';
    os << param.particle_def_ << '\n';
    os << *param.medium_ << '\n';
    os << param.cut_settings_ << '\n';
    os << Helper::Centered(80, "");
    return os;
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

//----------------------------------------------------------------------------//
double Parametrization::FunctionToDNdxIntegral(double energy, double variable)
{
    return DifferentialCrossSection(energy, variable);
}

// ------------------------------------------------------------------------- //
double Parametrization::FunctionToDEdxIntegral(double energy, double variable)
{
    return variable * DifferentialCrossSection(energy, variable);
}

// ------------------------------------------------------------------------- //
double Parametrization::FunctionToDE2dxIntegral(double energy, double variable)
{
    return variable * variable * DifferentialCrossSection(energy, variable);
}

// ------------------------------------------------------------------------- //
// Getter
// ------------------------------------------------------------------------- //

size_t Parametrization::GetHash() const
{
    std::size_t seed = 0;

    boost::hash_combine(seed, GetName());
    boost::hash_combine(seed, particle_def_.name);
    boost::hash_combine(seed, medium_->GetName());
    boost::hash_combine(seed, cut_settings_.GetEcut());
    boost::hash_combine(seed, cut_settings_.GetVcut());
    boost::hash_combine(seed, multiplier_);

    return seed;
}
