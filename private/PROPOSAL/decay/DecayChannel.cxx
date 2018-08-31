
#include <cmath>
#include <iostream>

#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/math/Vector3D.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Output.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

DecayChannel::DecayChannel() {}

DecayChannel::DecayChannel(const DecayChannel&) {}

DecayChannel::~DecayChannel() {}

bool DecayChannel::operator==(const DecayChannel& table) const
{
    return this->compare(table);
}

bool DecayChannel::operator!=(const DecayChannel& def) const
{
    return !(*this == def);
}

// ------------------------------------------------------------------------- //
std::ostream& PROPOSAL::operator<<(std::ostream& os, DecayChannel const& channel)
{
    std::stringstream ss;
    ss << " DecayChannel (" << &channel << ") ";

    os << Helper::Centered(60, ss.str()) << '\n';
    os << "Channel: " << channel.GetName() << '\n';

    channel.print(os);

    os << Helper::Centered(60, "");
    return os;
}

// ------------------------------------------------------------------------- //
void DecayChannel::Boost(Particle& particle, const Vector3D& direction_unnormalized, double beta)
{
    Vector3D direction = direction_unnormalized;
    direction.normalise();

    double gamma = 1.0 / std::sqrt(1.0 - beta * beta);

    Vector3D momentum_vec(particle.GetMomentum() * particle.GetDirection());

    double direction_correction =
        (gamma - 1.0) * scalar_product(momentum_vec, direction) - gamma * beta * particle.GetEnergy();

    momentum_vec = momentum_vec + direction_correction * direction;

    // Energy will be implicit corrected with respect to the mass:
    particle.SetMomentum(momentum_vec.magnitude());

    momentum_vec.normalise();
    particle.SetDirection(momentum_vec);
}

// ------------------------------------------------------------------------- //
void DecayChannel::Boost(DecayProducts& products, const Vector3D& direction, double beta)
{
    for (DecayProducts::const_iterator iter = products.begin(); iter != products.end(); ++iter)
    {
        Boost(**iter, direction, beta);
    }
}

// ------------------------------------------------------------------------- //
double DecayChannel::Momentum(double m1, double m2, double m3)
{
    double kaellen = (m1 - m2 - m3) * (m1 + m2 + m3) * (m1 - m2 + m3) * (m1 + m2 - m3);

    if (kaellen > 0.0)
    {
        return std::sqrt(kaellen) / (2.0 * m1);
    } else
    {
        log_fatal("Kaellen function is negative. Cannot caluclate momentum");
    }
}

// ------------------------------------------------------------------------- //
Vector3D DecayChannel::GenerateRandomDirection()
{
    double phi       = 2.0 * PI * RandomGenerator::Get().RandomDouble();
    double cos_theta = 2.0 * RandomGenerator::Get().RandomDouble() - 1.0;
    double sin_theta = std::sqrt((1.0 - cos_theta) * (1.0 + cos_theta));

    return Vector3D(sin_theta * std::sin(phi), sin_theta * std::cos(phi), cos_theta);
}

// ------------------------------------------------------------------------- //
// Protected Member functions
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
void DecayChannel::CopyParticleProperties(DecayProducts& products, const Particle& particle)
{
    int id = 1;
    for (std::vector<Particle*>::iterator iter = products.begin(); iter != products.end(); ++iter)
    {
        (*iter)->SetPosition(particle.GetPosition());
        (*iter)->SetTime(particle.GetTime());
        (*iter)->SetParentParticleEnergy(particle.GetEnergy());
        (*iter)->SetParticleId(particle.GetParticleId() + id);
        (*iter)->SetParentParticleId(particle.GetParticleId());

        id++;
    }
}
