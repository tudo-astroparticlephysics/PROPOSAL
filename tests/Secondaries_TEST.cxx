#include "gtest/gtest.h"

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/geometry/Sphere.h"
#include "PROPOSAL/geometry/Box.h"
#include "PROPOSAL/geometry/Cylinder.h"
#include "PROPOSAL/Secondaries.h"
#include "PROPOSAL/crosssection/ParticleDefaultCrossSectionList.h"
#include "PROPOSAL/propagation_utility/TimeBuilder.h"
#include "PROPOSAL/propagation_utility/InteractionBuilder.h"
#include "PROPOSAL/density_distr/density_homogeneous.h"

using namespace PROPOSAL;

std::shared_ptr<Propagator> GetPropagator() {
    static std::shared_ptr<Propagator> ptr = nullptr;
    if(ptr == nullptr) {
        //create propagator
        auto p_def = MuMinusDef();
        auto medium = Ice();
        auto cuts = std::make_shared<EnergyCutSettings>(INF, 1, false);
        auto cross = GetStdCrossSections(p_def, medium, cuts, true);

        auto collection = PropagationUtility::Collection();
        collection.interaction_calc = make_interaction(cross, true);
        collection.displacement_calc = make_displacement(cross, true);
        collection.time_calc = make_time(cross, p_def, true);
        auto prop_utility = PropagationUtility(collection);

        auto density_distr = std::make_shared<Density_homogeneous>(medium);
        auto world = std::make_shared<Sphere>(Cartesian3D(0, 0, 0), 1e20);
        auto sphere = std::make_shared<Sphere>(Cartesian3D(0, 0, 10000), 1000);

        auto sector1 = std::make_tuple(world, prop_utility, density_distr);
        auto sector2 = std::make_tuple(sphere, prop_utility, density_distr);
        std::vector<Sector> sec_vec = {sector1, sector2};
        auto prop = Propagator(p_def, sec_vec);
        ptr = std::make_shared<Propagator>(prop);
    }
    return ptr;
}

std::shared_ptr<Propagator> GetPropagatorStochastic() {
    static std::shared_ptr<Propagator> ptr = nullptr;
    if(ptr == nullptr) {
        //create propagator
        auto p_def = MuMinusDef();
        auto medium = Ice();
        auto cuts = std::make_shared<EnergyCutSettings>(INF, 0.05, false);
        auto cross = GetStdCrossSections(p_def, medium, cuts, true);

        auto collection = PropagationUtility::Collection();
        collection.interaction_calc = make_interaction(cross, true);
        collection.displacement_calc = make_displacement(cross, true);
        collection.time_calc = make_time(cross, p_def, true);
        auto prop_utility = PropagationUtility(collection);

        auto density_distr = std::make_shared<Density_homogeneous>(medium);
        auto world = std::make_shared<Sphere>(Cartesian3D(0, 0, 0), 1e20);

        auto sector = std::make_tuple(world, prop_utility, density_distr);
        std::vector<Sector> sec_vec = {sector};
        auto prop = Propagator(p_def, sec_vec);
        ptr = std::make_shared<Propagator>(prop);
    }
    return ptr;
}

TEST(SecondaryVector, EntryPointExitPoint)
{
    auto prop = GetPropagator();

    Cartesian3D position(0, 0, 0);
    Cartesian3D direction(0, 0, 1);
    auto energy = 1e8; // MeV
    auto init_state = ParticleState(position, direction, energy, 0., 0.);

    auto secondaries = prop->Propagate(init_state, 50000);

    // Test for geometry in front of track
    auto sphere_infront = Sphere(Cartesian3D(0, 0, -1000), 10);
    EXPECT_TRUE(secondaries.GetEntryPoint(sphere_infront) == nullptr);
    EXPECT_TRUE(secondaries.GetExitPoint(sphere_infront) == nullptr);

    // Test for geometry behin track
    auto sphere_behind = Sphere(Cartesian3D(0, 0, 100000), 10);
    EXPECT_TRUE(secondaries.GetEntryPoint(sphere_behind) == nullptr);
    EXPECT_TRUE(secondaries.GetExitPoint(sphere_behind) == nullptr);

    // Test for track beginning in geometry
    auto sphere_start = Sphere(Cartesian3D(0, 0, 100), 100);
    EXPECT_TRUE(secondaries.GetEntryPoint(sphere_start)->energy == energy);
    EXPECT_TRUE(secondaries.GetEntryPoint(sphere_start)->position == position);

    // Test for track ending in geometry
    auto sphere_end = Sphere(Cartesian3D(0, 0, 49000), 1000);
    EXPECT_TRUE(secondaries.GetExitPoint(sphere_end)->energy == secondaries.back().energy);
    EXPECT_TRUE(secondaries.GetExitPoint(sphere_end)->propagated_distance == secondaries.back().propagated_distance);
    EXPECT_TRUE(secondaries.GetExitPoint(sphere_end)->position == secondaries.back().position);

    // Test points where interaction points should be equal to entry/exit points
    auto sphere = Sphere(Cartesian3D(0, 0, 10000), 1000);
    EXPECT_TRUE(secondaries.GetEntryPoint(sphere)->propagated_distance == 9000);
    EXPECT_TRUE(secondaries.GetExitPoint(sphere)->propagated_distance == 11000);
}

TEST(SecondaryVector, EntryPointExitPointRePropagation)
{
    auto prop = GetPropagator();

    Cartesian3D position(0, 0, 0);
    Cartesian3D direction(0, 0, 1);
    auto energy = 1e8; // MeV
    auto init_state = ParticleState(position, direction, energy, 0., 0.);

    auto secondaries = prop->Propagate(init_state, 1e5);

    auto sphere = Sphere(Cartesian3D(0, 0, 5e4), 500);
    auto entry_point = secondaries.GetEntryPoint(sphere);
    auto exit_point = secondaries.GetExitPoint(sphere);
    int i = 0;
    while (secondaries[i].propagated_distance < entry_point->propagated_distance ) {
        i++;
    }
    auto sec_i = secondaries[i-1];
    auto sec_f = secondaries[i];

    EXPECT_TRUE(sec_i.energy > entry_point->energy);
    EXPECT_TRUE(sec_i.time < entry_point->time);
    EXPECT_TRUE(sec_i.propagated_distance < entry_point->propagated_distance);
    EXPECT_TRUE(sphere.IsInfront(sec_i.position, sec_i.direction));

    EXPECT_TRUE(sec_f.energy < entry_point->energy);
    EXPECT_TRUE(sec_f.time > entry_point->time);
    EXPECT_TRUE(sec_f.propagated_distance > entry_point->propagated_distance);
    EXPECT_FALSE(sphere.IsInfront(sec_f.position, sec_f.direction));
    EXPECT_TRUE(entry_point->propagated_distance == sphere.GetPosition().GetZ() - sphere.GetRadius());

    while (secondaries[i].propagated_distance < exit_point->propagated_distance) {
        i++;
    }

    sec_i = secondaries[i-1]; // point before exit point
    sec_f = secondaries[i]; // points after exit point

    EXPECT_TRUE(sec_i.energy > exit_point->energy);
    EXPECT_TRUE(sec_i.time < exit_point->time);
    EXPECT_TRUE(sec_i.propagated_distance < exit_point->propagated_distance);
    EXPECT_FALSE(sphere.IsBehind(sec_i.position, sec_i.direction));

    EXPECT_TRUE(sec_f.energy < exit_point->energy);
    EXPECT_TRUE(sec_f.time > exit_point->time);
    EXPECT_TRUE(sec_f.propagated_distance > exit_point->propagated_distance);
    EXPECT_TRUE(sphere.IsBehind(sec_f.position, sec_f.direction));
    EXPECT_TRUE(exit_point->propagated_distance == sphere.GetPosition().GetZ() + sphere.GetRadius());
}


TEST(SecondaryVector, HitDetector) {
    // define our dummy particle track
    Secondaries dummy_track(nullptr, std::vector<Sector>{});
    std::vector<Cartesian3D> positions{ {0, 0, 0}, {0, 0, 10}, {0, 0, 20} };

    for (auto p : positions) {
        ParticleState p_state;
        p_state.position = p;
        p_state.direction = Cartesian3D(0, 0, 1);
        dummy_track.push_back(p_state, InteractionType::Undefined);
    }

    // test geometries along propagation axis
    EXPECT_FALSE(dummy_track.HitGeometry(Sphere(Cartesian3D(0, 0, -5), 1)));
    for (double z = -0.5; z <= 20.5; z=z+0.5) {
        EXPECT_TRUE(dummy_track.HitGeometry(Sphere(Cartesian3D(0, 0, z), 1)));
        EXPECT_TRUE(dummy_track.HitGeometry(Box(Cartesian3D(0, 0, z), 2, 2, 2)));
        EXPECT_TRUE(dummy_track.HitGeometry(Cylinder(Cartesian3D(0, 0, z), 5, 1)));
    }
    EXPECT_TRUE(dummy_track.HitGeometry(Sphere(Cartesian3D(0, 0, 20), 1)));
    EXPECT_FALSE(dummy_track.HitGeometry(Sphere(Cartesian3D(0, 0, 30), 1)));

    // test geometries displaced to propagation axis
    EXPECT_FALSE(dummy_track.HitGeometry(Sphere(Cartesian3D(0, 5, 10), 1)));
    EXPECT_FALSE(dummy_track.HitGeometry(Sphere(Cartesian3D(-5, 0, 5), 1)));

    // check border cases
    EXPECT_FALSE(dummy_track.HitGeometry(Sphere(Cartesian3D(0, 0, -1), 1)));
    EXPECT_FALSE(dummy_track.HitGeometry(Sphere(Cartesian3D(0, 1, 10), 1)));
    EXPECT_TRUE(dummy_track.HitGeometry(Sphere(Cartesian3D(0, 0, 21), 1)));
    EXPECT_TRUE(dummy_track.HitGeometry(Box(Cartesian3D(0, 0, 21), 2, 2, 2)));
    EXPECT_TRUE(dummy_track.HitGeometry(Cylinder(Cartesian3D(0, 0, 21), 5, 2)));

}

TEST(SecondaryVector, EnergyConservation) {
    auto prop = GetPropagatorStochastic();

    Cartesian3D position(0, 0, 0);
    Cartesian3D direction(0, 0, 1);
    auto energy = 1e8; // MeV
    auto init_state = ParticleState(position, direction, energy, 0., 0.);

    auto secondaries = prop->Propagate(init_state);

    double sum_continuous_losses = 0.;
    for (auto continuous_loss : secondaries.GetContinuousLosses()) {
        double loss_energy = continuous_loss.energy;
        EXPECT_GT(loss_energy, 0.); // must be positive
        sum_continuous_losses += loss_energy;
    }

    double sum_stochastic_losses = 0.;
    for (auto stochastic_loss : secondaries.GetStochasticLosses()) {
        double loss_energy = stochastic_loss.energy;
        EXPECT_GT(loss_energy, 0); // must be positive
        sum_stochastic_losses += loss_energy;
    }

    // energy needs to be conserved
    EXPECT_DOUBLE_EQ(sum_continuous_losses + sum_stochastic_losses + MuMinusDef().mass, energy);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
