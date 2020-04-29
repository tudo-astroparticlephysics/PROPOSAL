
#include <functional>

#include "PROPOSAL/crossection/AnnihilationInterpolant.h"

#include "PROPOSAL/math/InterpolantBuilder.h"

#include "PROPOSAL/particle/Particle.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/methods.h"

using std::pair;
using std::vector;

using namespace PROPOSAL;

double AnnihilationInterpolant::CalculatedEdx(double energy)
{
    (void)energy;
    return 0;
}

double AnnihilationInterpolant::CalculatedE2dx(double energy)
{
    (void)energy;
    return 0;
}

/* double AnnihilationInterpolant::CalculateStochasticLoss( */
/*     double energy, array<double, 2> rnd) */
/* { */
/*     CalculatedNdx(energy, rnd[0]); */
/*     rndc_ = rnd[1]; // save random number for component sampling */
/*     return energy;  // losses are always catastrophic */
/* } */

/* vector<DynamicData> AnnihilationInterpolant::CalculateProducedParticles( */
/*     double energy, Vector3D&& initial_dir, array<double, 2>&& rnd) */
/* { */
/*     auto rates = CalculatedNdx(energy, rnd[0]); */
/*     auto sampled_rate = rnd[1] * accumulate(rates.begin(), rates.end(), 0.); */

/*     vector<DynamicData> particle_list(2); */
/*     particle_list.emplace_back(static_cast<int>(ParticleType::Gamma)); */
/*     particle_list.emplace_back(static_cast<int>(ParticleType::Gamma)); */

/*     auto rate = rates.cbegin(); */
/*     auto dndx = dndx_interpolants_.begin(); */
/*     auto kin_limit = kinematic_limits_.begin(); */
/*     for (; rate != rates.end(); ++rate, ++dndx, ++kin_limit) { */
/*         sampled_rate -= *rate; */

/*         if (sampled_rate < 0) { */
/*             auto v_max = (*kin_limit)(energy).vMax; */
/*             auto v_cut = GetEnergyCut(energy); */

/*             auto v = (*dndx)->FindLimit(energy, rnd[1] * *rate); */

/*             auto rho = logarithm_trafo(v, v_cut, v_max); */

/*             // The available energy is the positron energy plus the mass of the */
/*             // electron */
/*             particle_list[0].SetEnergy((energy + ME) * (1 - rho)); */
/*             particle_list[1].SetEnergy((energy + ME) * rho); */

/*             particle_list[0].SetDirection(initial_dir); */
/*             particle_list[1].SetDirection(initial_dir); */

/*             double cosphi0 = ((energy + ME) * (1. - rho) - ME) */
/*                 / ((1. - rho) * std::sqrt((energy + ME) * (energy - ME))); */
/*             double cosphi1 = ((energy + ME) * rho - ME) */
/*                 / (rho * std::sqrt((energy + ME) * (energy - ME))); */

/*             double rndtheta = RandomGenerator::Get().RandomDouble(); */

/*             particle_list[0].DeflectDirection(cosphi0, rndtheta * 2. * PI); */
/*             particle_list[1].DeflectDirection( */
/*                 cosphi1, std::fmod(rndtheta * 2. * PI + PI, 2. * PI)); */

/*             return particle_list; */
/*         } */
/*     } */

/*     throw std::logic_error( */
/*         "Could not sample ProducedParticles for PhotoPairProduction."); */
/* } */
