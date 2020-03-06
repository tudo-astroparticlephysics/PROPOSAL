
#include <cmath>

#include "PROPOSAL/crossection/parametrization/Annihilation.h"

#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/math/MathMethods.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

// ------------------------------------------------------------------------- //
// Constructor & Destructor
// ------------------------------------------------------------------------- //

Annihilation::Annihilation(const ParticleDef& particle_def,
                                 std::shared_ptr<const Medium> medium,
                                 double multiplier)
        : Parametrization(particle_def, medium, EnergyCutSettings(), multiplier)
{
}

Annihilation::Annihilation(const Annihilation& param)
        : Parametrization(param)
{
}

Annihilation::~Annihilation() {}

bool Annihilation::compare(const Parametrization& parametrization) const
{
    return Parametrization::compare(parametrization);
}

// ------------------------------------------------------------------------- //
// Public methods
// ------------------------------------------------------------------------- //

Parametrization::IntegralLimits Annihilation::GetIntegralLimits(double energy)
{
    // Limits according to simple 2->2 body interactions
    IntegralLimits limits;

    double gamma    = energy / particle_def_.mass;

    limits.vMin     = 0.5 * (1. - std::sqrt( std::max(0., (gamma - 1.)/(gamma + 1.) ) ));
    limits.vMax     = 0.5 * (1. + std::sqrt( std::max(0., (gamma - 1.)/(gamma + 1.) ) ));

    limits.vUp = limits.vMin; //treat interaction as fully stochastic

    return limits;
}

size_t Annihilation::GetHash() const
{
    size_t seed = Parametrization::GetHash();
    hash_combine(seed);

    return seed;
}

// ------------------------------------------------------------------------- //
// Specific implementations
// ------------------------------------------------------------------------- //

AnnihilationHeitler::AnnihilationHeitler(const ParticleDef& particle_def, std::shared_ptr<const Medium> medium, double multiplier)
        : Annihilation(particle_def, medium, multiplier)
{}

AnnihilationHeitler::AnnihilationHeitler(const AnnihilationHeitler& param)
        : Annihilation(param)
{}

AnnihilationHeitler::~AnnihilationHeitler()
{}

double AnnihilationHeitler::DifferentialCrossSection(double energy, double v)
{
    // W. Heitler. The Quantum Theory of Radiation, Clarendon Press, Oxford (1954)
    // Adapted from Geant4 PhysicsReferenceManual

    // v = energy of photon1 / total available energy
    // with the total available energy being the sum of the total positron energy and the electron mass

    double gamma    = energy / particle_def_.mass;
    double aux      = 1. + (2. * gamma) / std::pow(gamma + 1., 2.) - v - 1. / std::pow(gamma + 1., 2.) * 1. / v;
    aux *= medium_->GetMassDensity() * NA * medium_->GetZA() * PI * RE * RE / (gamma - 1.) * 1. / v; // TODO: prefactors

    return aux;
}


const std::string AnnihilationHeitler::name_ = "AnnihilationHeitler";
