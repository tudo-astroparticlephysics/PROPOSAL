
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>

#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/math/Vector3D.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
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

namespace PROPOSAL {

std::ostream& operator<<(std::ostream& os, DecayChannel const& channel)
{
    std::stringstream ss;
    ss << " DecayChannel (" << &channel << ") ";

    os << Helper::Centered(60, ss.str()) << '\n';
    os << "Channel: " << channel.GetName() << '\n';

    channel.print(os);

    os << Helper::Centered(60, "");
    return os;
}

} // namespace PROPOSAL

// ------------------------------------------------------------------------- //
void DecayChannel::Boost(ParticleState& particle, const Vector3D& direction_unnormalized, double gamma, double betagamma)
{
    Cartesian3D direction = direction_unnormalized;
    direction.normalize();

    Cartesian3D momentum_vec(particle.GetMomentum() * particle.direction);

    double direction_correction =
        (gamma - 1.0) * (momentum_vec * direction) - betagamma * particle.energy;

    momentum_vec = momentum_vec + direction_correction * direction;

    // Energy will be implicit corrected with respect to the mass:
    particle.SetMomentum(momentum_vec.magnitude());

    momentum_vec.normalize();
    particle.direction = momentum_vec;
}

// ------------------------------------------------------------------------- //
void DecayChannel::Boost(std::vector<ParticleState>& secondaries, const Vector3D& direction, double gamma, double betagamma)
{
    for (auto& p : secondaries)
    {
        Boost(p, direction, gamma, betagamma);
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
        // here this term is numerically unstable,
        // to handle also the cases at the edge of the phase space
        // use COMPUTER_PRECISION or std::numeric_limits<double>::epsilon()
        if (std::abs(m1 - m2 - m3) < m1*COMPUTER_PRECISION)
        {
            Logging::Get("proposal.decay")->warn("(m1-m2-m3) in Kaellen function is numerically unstable and slightly negative {}. Now set to zero.",
                (m1 - m2 - m3));
        }
        else
        {
            Logging::Get("proposal.decay")->error("Kaellen function is negative. Cannot caluclate momentum");
        }
        return 0.0;
    }
}

// ------------------------------------------------------------------------- //
Cartesian3D DecayChannel::GenerateRandomDirection()
{
    double phi       = 2.0 * PI * RandomGenerator::Get().RandomDouble();
    double cos_theta = 2.0 * RandomGenerator::Get().RandomDouble() - 1.0;
    double sin_theta = std::sqrt((1.0 - cos_theta) * (1.0 + cos_theta));
    Cartesian3D direction(sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta);
    return direction;
}

// ------------------------------------------------------------------------- //
// Protected Member functions
// ------------------------------------------------------------------------- //
