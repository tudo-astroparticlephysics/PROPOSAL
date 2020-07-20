
#include <cmath>
#include <iostream>

#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/math/Vector3D.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/methods.h"

#include "PROPOSAL/Secondaries.h"

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
void DecayChannel::Boost(DynamicData& particle, const Vector3D& direction_unnormalized, double gamma, double betagamma)
{
    Vector3D direction = direction_unnormalized;
    direction.normalise();

    Vector3D momentum_vec(particle.GetMomentum() * particle.GetDirection());

    double direction_correction =
        (gamma - 1.0) * scalar_product(momentum_vec, direction) - betagamma * particle.GetEnergy();

    momentum_vec = momentum_vec + direction_correction * direction;

    // Energy will be implicit corrected with respect to the mass:
    particle.SetMomentum(momentum_vec.magnitude());

    momentum_vec.normalise();
    momentum_vec.CalculateSphericalCoordinates();
    particle.SetDirection(momentum_vec);
}

// ------------------------------------------------------------------------- //
void DecayChannel::Boost(Secondaries& secondaries, const Vector3D& direction, double gamma, double betagamma)
{
    for (auto& p : secondaries.GetModifyableSecondaries())
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
            log_warn("(m1-m2-m3) in Kaellen function is numerically unstable and slightly negative %f. Now set to zero.",
                (m1 - m2 - m3));
        }
        else
        {
            log_fatal("Kaellen function is negative. Cannot caluclate momentum");
        }
        return 0.0;
    }
}

// ------------------------------------------------------------------------- //
Vector3D DecayChannel::GenerateRandomDirection()
{
    double phi       = 2.0 * PI * RandomGenerator::Get().RandomDouble();
    double cos_theta = 2.0 * RandomGenerator::Get().RandomDouble() - 1.0;
    double sin_theta = std::sqrt((1.0 - cos_theta) * (1.0 + cos_theta));
    Vector3D direction = Vector3D(sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta);
    direction.CalculateSphericalCoordinates();
    // direction.SetSphericalCoordinates(1.0, phi, std::acos(cos_theta));
    return direction;
}

// ------------------------------------------------------------------------- //
// Protected Member functions
// ------------------------------------------------------------------------- //
