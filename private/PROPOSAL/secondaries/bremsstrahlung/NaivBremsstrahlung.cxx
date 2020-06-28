#include "PROPOSAL/secondaries/bremsstrahlung/NaivBremsstrahlung.h"
#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"

#include <cmath>
#include <stdexcept>

using namespace PROPOSAL;

vector<Loss::secondary_t> secondaries::NaivBremsstrahlung::CalculateSecondaries(
    double primary_energy, Loss::secondary_t loss, const Component&, vector<double>)
{
    auto primary_lepton = loss;
    std::get<Loss::ENERGY>(primary_lepton) = primary_energy - std::get<Loss::ENERGY>(loss);
    std::get<Loss::TYPE>(primary_lepton) = primary_lepton_type;
    std::get<Loss::TYPE>(loss) = static_cast<int>(PROPOSAL::ParticleType::Gamma);
    auto sec = vector<Loss::secondary_t>{ move(primary_lepton), move(loss) };
    return sec;
}
