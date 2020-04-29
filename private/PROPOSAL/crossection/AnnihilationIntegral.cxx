
#include <functional>
#include <stdexcept>
#include <vector>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crossection/AnnihilationIntegral.h"
#include "PROPOSAL/crossection/parametrization/Annihilation.h"
#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/RandomGenerator.h"

using namespace PROPOSAL;
using std::vector;

double AnnihilationIntegral::CalculatedEdx(double energy)
{
    (void)energy;
    return 0;
}

double AnnihilationIntegral::CalculatedE2dx(double energy)
{
    (void)energy;
    return 0;
}

/* vector<DynamicData> AnnihilationIntegral::CalculateProducedParticles(double energy, Vector3D&& initial_dir, array<double, 2>&& rnd) */
/* { */
/*     auto rates = CalculatedNdx(energy, rnd[0]); */
/*     auto sampled_rate = rnd[1] * accumulate(rates.begin(), rates.end(), 0.); */

/*     vector<DynamicData> particle_list(2); */
/*     particle_list.emplace_back(static_cast<int>(ParticleType::Gamma)); */
/*     particle_list.emplace_back(static_cast<int>(ParticleType::Gamma)); */

/*     auto rate = rates.cbegin(); */
/*     auto integral = dNdx_.begin(); */
/*     for (; rate != rates.end() and integral != dNdx_.end(); */
/*          ++rate, ++integral) { */
/*         sampled_rate -= *rate; */

/*         if (sampled_rate < 0) { */
/*             auto rho = integral->GetUpperLimit(); */

/*             // The available energy is the positron energy plus the mass of the */
/*             // electron */
/*             particle_list[0].SetEnergy((energy + ME) * (1 - rho)); */
/*             particle_list[1].SetEnergy((energy + ME) * rho); */

/*             particle_list[0].SetDirection(initial_dir); */
/*             particle_list[1].SetDirection(initial_dir); */

/*             auto cosphi0 = ((energy + ME) * (1. - rho) - ME) */
/*                 / ((1. - rho) * std::sqrt((energy + ME) * (energy - ME))); */
/*             auto cosphi1 = ((energy + ME) * rho - ME) */
/*                 / (rho * std::sqrt((energy + ME) * (energy - ME))); */

/*             auto rndtheta = RandomGenerator::Get().RandomDouble(); */

/*             particle_list[0].DeflectDirection(cosphi0, rndtheta * 2. * PI); */
/*             particle_list[1].DeflectDirection( */
/*                 cosphi1, std::fmod(rndtheta * 2. * PI + PI, 2. * PI)); */

/*             return particle_list; */
/*         } */
/*     } */
/*     throw std::logic_error( */
/*         "Could not sample ProducedParticles for PhotoPairProduction."); */
/* } */

/* double AnnihilationIntegral::CalculateStochasticLoss(double energy, array<double, 2>) */
/* { */
/*     return energy; // losses are always catastrophic */
/* } */
