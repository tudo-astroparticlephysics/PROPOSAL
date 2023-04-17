#include "gtest/gtest.h"
#include "PROPOSAL/crosssection/ParticleDefaultCrossSectionList.h"
#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/propagation_utility/TimeBuilder.h"
#include "PROPOSAL/propagation_utility/InteractionBuilder.h"
#include "PROPOSAL/propagation_utility/ContRandBuilder.h"
#include "PROPOSAL/propagation_utility/DecayBuilder.h"
#include "PROPOSAL/density_distr/density_homogeneous.h"
#include "PROPOSAL/scattering/ScatteringFactory.h"
#include "PROPOSAL/geometry/Sphere.h"
#include "PROPOSAL/particle/Particle.h"

using namespace PROPOSAL;

TEST(Propagator, min_energy)
{
    auto p_def = MuMinusDef();
    auto medium = Ice();
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 0.05, true);
    auto cross = GetStdCrossSections(p_def, medium, cuts, true);

    auto collection = PropagationUtility::Collection();
    collection.interaction_calc = make_interaction(cross, true);
    collection.displacement_calc = make_displacement(cross, true);
    collection.time_calc = make_time(cross, p_def, true);
    collection.cont_rand = make_contrand(cross, true);
    collection.scattering = make_scattering(MultipleScatteringType::Highland, {}, p_def, medium);

    auto prop_utility = PropagationUtility(collection);

    auto density_distr = std::make_shared<Density_homogeneous>(medium);
    auto world = std::make_shared<Sphere>(Cartesian3D(0, 0, 0), 1e20);

    auto sector = std::make_tuple(world, prop_utility, density_distr);
    std::vector<Sector> sec_vec = {sector};

    auto prop = Propagator(p_def, sec_vec);

    auto init_state = ParticleState();
    init_state.energy = 1e5;
    init_state.position = Cartesian3D(0, 0, 0);
    init_state.direction = Cartesian3D(0, 0, 1);

    double min_energy = 1e4;
    for (size_t i=0; i<1000; i++) {
        auto sec = prop.Propagate(init_state, 1e20, min_energy);
        double final_energy = sec.GetTrack().back().energy;
        // If we propagate stable particles to min_energy, no particles should
        // have an energy bigger than min_energy
        EXPECT_LE(final_energy, min_energy);

        // Particles with a stochastic loss as their last interaction are allowed
        // to have energies smaller than min_energy. Otherwise, they don't.
        if (final_energy < min_energy)
            EXPECT_NE(sec.GetTrackTypes().back(), InteractionType::ContinuousEnergyLoss);
        else
            EXPECT_EQ(sec.GetTrackTypes().back(), InteractionType::ContinuousEnergyLoss);
    }
}

TEST(Propagator, RestMassDecay)
{
    // make sure that all muons decay at the end of their propagation, even if they have reached their rest mass
    // (see issue https://github.com/tudo-astroparticlephysics/PROPOSAL/issues/323)
    auto p_def = MuMinusDef();
    auto medium = Ice();
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 0.05, true);
    auto cross = GetStdCrossSections(p_def, medium, cuts, true);

    auto collection = PropagationUtility::Collection();
    collection.interaction_calc = make_interaction(cross, true);
    collection.displacement_calc = make_displacement(cross, true);
    collection.time_calc = make_time(cross, p_def, true);
    collection.cont_rand = make_contrand(cross, true);
    collection.scattering = make_scattering(MultipleScatteringType::Highland, {}, p_def, medium);
    collection.decay_calc = make_decay(cross, p_def, true);

    auto prop_utility = PropagationUtility(collection);

    auto density_distr = std::make_shared<Density_homogeneous>(medium);
    auto world = std::make_shared<Sphere>(Cartesian3D(0, 0, 0), 1e20);

    auto sector = std::make_tuple(world, prop_utility, density_distr);
    std::vector<Sector> sec_vec = {sector};

    auto prop = Propagator(p_def, sec_vec);

    auto init_state = ParticleState();
    init_state.energy = 150;
    init_state.position = Cartesian3D(0, 0, 0);
    init_state.direction = Cartesian3D(0, 0, 1);

    for (size_t i=0; i<1000; i++) {
        auto sec = prop.Propagate(init_state);
        EXPECT_FALSE(sec.GetDecayProducts().empty()); // check that decay products are produced
    }

}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}