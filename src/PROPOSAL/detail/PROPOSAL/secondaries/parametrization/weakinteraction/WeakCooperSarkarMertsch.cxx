#include "PROPOSAL/secondaries/parametrization/weakinteraction/WeakCooperSarkarMertsch.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>
#include <stdexcept>

using namespace PROPOSAL;

size_t secondaries::WeakCooperSarkarMertsch::generate_hash(const ParticleDef& p,
                                                           const Medium& m) {
    size_t hash = 0;
    hash_combine(hash, p.GetHash(), m.GetHash(), "secondaries_coopersarkarmertsch");
    return hash;
}

double secondaries::WeakCooperSarkarMertsch::CalculateRelativeLoss(
        double energy, double rnd, const Component& c) {
    if (!dndx)
        throw std::logic_error("dndx Interpolant for WeakInteraction not defined.");
    for (auto& it : *dndx) {
        if (c.GetHash() == it.first) {
            auto& calc = *std::get<1>(it.second);
            auto rate = rnd * calc.Calculate(energy);
            auto v = calc.GetUpperLimit(energy, rate);
            return v;
        }
    }
    std::ostringstream s;
    s << "Component (" << c.GetName()
    << ") can not be found in the precalculated weak interaction tables.";
    throw std::out_of_range(s.str());
}

std::vector<ParticleState>
secondaries::WeakCooperSarkarMertsch::CalculateSecondaries(StochasticLoss loss,
                                                       const Component& c,
                                                       std::vector<double>& rnd)
{
    // v corresponds to the energy lost to the hadron
    auto v = CalculateRelativeLoss(loss.parent_particle_energy, rnd[0], c);

    auto sec = std::vector<ParticleState>();
    sec.emplace_back(static_cast<ParticleType>(weak_partner_type),
                     loss.position, loss.direction, (1. - v) * loss.energy,
                     loss.time, 0.);
    sec.emplace_back(ParticleType::Hadron,
                     loss.position, loss.direction, v * loss.energy,
                     loss.time, 0.);
    return sec;
}
