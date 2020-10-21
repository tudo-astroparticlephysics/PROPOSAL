
#include "gtest/gtest.h"

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/Secondaries.h"
#include "PROPOSAL/propagation_utility/InteractionBuilder.h"
#include "PROPOSAL/density_distr/density_homogeneous.h"

using namespace PROPOSAL;

std::shared_ptr<Propagator> GetPropagator() {
    static std::shared_ptr<Propagator> ptr = nullptr;
    if(ptr == nullptr) {
        //create propagator
        InterpolationDef::path_to_tables = "resources/tables";
        auto p_def = MuMinusDef();
        auto medium = Ice();
        auto cuts = std::make_shared<EnergyCutSettings>(500, 1, false);
        auto cross = GetStdCrossSections(p_def, medium, cuts, true);

        auto collection = PropagationUtility::Collection();
        collection.interaction_calc = make_interaction(cross, true);
        collection.displacement_calc = make_displacement(cross, true);
        collection.time_calc = make_time(cross, p_def, true);
        auto prop_utility = PropagationUtility(collection);

        auto density_distr = std::make_shared<Density_homogeneous>(medium);
        auto world = std::make_shared<Sphere>(Vector3D(0, 0, 0), 1e20);
        auto sphere = std::make_shared<Sphere>(Vector3D(0, 0, 10000), 1000);

        auto sector1 = std::make_tuple(world, prop_utility, density_distr);
        auto sector2 = std::make_tuple(sphere, prop_utility, density_distr);
        auto prop = Propagator(p_def, {sector1, sector2});
        ptr = std::make_shared<Propagator>(prop);
    }
    return ptr;
}

TEST(SecondaryVector, EntryPointExitPoint)
{
    auto prop = GetPropagator();

    Vector3D position(0, 0, 0);
    Vector3D direction(0, 0, 1);
    auto energy = 1e8; // MeV
    auto init_state = DynamicData(13, position, direction, energy, energy, 0., 0.);

    auto secondaries = prop->Propagate(init_state, 50000);

    // Test for geometry in front of track
    auto sphere_infront = Sphere(Vector3D(0, 0, -1000), 10);
    EXPECT_TRUE(secondaries.GetEntryPoint(sphere_infront) == nullptr);
    EXPECT_TRUE(secondaries.GetExitPoint(sphere_infront) == nullptr);

    // Test for geometry behin track
    auto sphere_behind = Sphere(Vector3D(0, 0, 100000), 10);
    EXPECT_TRUE(secondaries.GetEntryPoint(sphere_behind) == nullptr);
    EXPECT_TRUE(secondaries.GetExitPoint(sphere_behind) == nullptr);

    // Test for track beginning in geometry
    auto sphere_start = Sphere(Vector3D(0, 0, 100), 100);
    EXPECT_TRUE(secondaries.GetEntryPoint(sphere_start)->GetEnergy() == energy);
    EXPECT_TRUE(secondaries.GetEntryPoint(sphere_start)->GetPropagatedDistance() == 0);
    EXPECT_TRUE(secondaries.GetEntryPoint(sphere_start)->GetPosition() == position);

    // Test for track ending in geometry
    auto sphere_end = Sphere(Vector3D(0, 0, 49000), 1000);
    EXPECT_TRUE(secondaries.GetExitPoint(sphere_end)->GetEnergy() == secondaries.back().GetEnergy());
    EXPECT_TRUE(secondaries.GetExitPoint(sphere_end)->GetPropagatedDistance() == secondaries.back().GetPropagatedDistance());
    EXPECT_TRUE(secondaries.GetExitPoint(sphere_end)->GetPosition() == secondaries.back().GetPosition());

    // Test points where interaction points should be equal to entry/exit points
    auto sphere = Sphere(Vector3D(0, 0, 10000), 1000);
    EXPECT_TRUE(secondaries.GetEntryPoint(sphere)->GetPropagatedDistance() == 9000);
    EXPECT_TRUE(secondaries.GetExitPoint(sphere)->GetPropagatedDistance() == 11000);

}
