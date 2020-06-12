#include "PROPOSAL/secondaries/bremsstrahlung/NaivBremsstrahlung.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>
#include <stdexcept>

using namespace PROPOSAL;

vector<Loss::secondary_t> secondaries::NaivBremsstrahlung::CalculateSecondaries(
    Loss::secondary_t loss,
    array<double, secondaries::NaivBremsstrahlung::n_rnd> rnd)
        {
            std::get<Loss::TYPE>(loss)
                = static_cast<int>(PROPOSAL::ParticleType::Gamma);
            auto sec = vector<Loss::secondary_t>{ move(loss) };
            return sec;
        }
