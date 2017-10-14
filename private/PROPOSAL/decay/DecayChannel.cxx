
#include <iostream>
#include <cmath>

#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/math/Vector3D.h"

#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

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
void DecayChannel::Boost(Particle* particle, const Vector3D& direction_unnormalized, double beta)
{
    Vector3D direction = direction_unnormalized;
    direction.normalise();

    double gamma = 1.0 / std::sqrt(1.0 - beta * beta);

    Vector3D momentum_vec(particle->GetMomentum()*particle->GetDirection());

    double direction_correction =
        (gamma - 1.0) * scalar_product(momentum_vec, direction) - gamma * beta * particle->GetEnergy();

    momentum_vec = momentum_vec + direction_correction * direction;

    // Energy will be implicit corrected with respect to the mass:
    particle->SetMomentum(momentum_vec.magnitude());

    momentum_vec.normalise();
    particle->SetDirection(momentum_vec);
}
