
// #include <iostream>
// #include <string>
// #include <cmath>
#include <vector>
#include <cmath>

#include "gtest/gtest.h"

#include "PROPOSAL/PROPOSAL.h"

using namespace std;
using namespace PROPOSAL;

TEST(Propagation , Test_nan) {

    int statistic = 100;
    int EmaxLog10 = 8;

    std::vector<Sector*> collections;
    std::vector<SectorFactory::Definition> sector_defintions;

    // --------------------------------------------------------------------- //
    // Builder Pattern
    // --------------------------------------------------------------------- //

    // Sector* sec = SectorBuilder::Get()
    //                   .Location(1)
    //                   .EnableInterpolation(true)
    //                   .EnableContinousRandomization(true)
    //                   .EnableScattering(true)
    //                   .ScatteringModel(ScatteringFactory::MoliereFirstOrder)
    //                   .DensityCorrection(0.98)
    //                   .ECut(500)
    //                   .VCut(0.01)
    //                   .build();
    //
    // collections.push_back(sec);

    // --------------------------------------------------------------------- //
    // Factory Pattern
    // --------------------------------------------------------------------- //

    PROPOSALParticle particle(MuMinusDef::Get());
    particle.SetDirection(Vector3D(0, 0, -1));

    // ----[ Interpolation ]--------------------------------- //

    InterpolationDef interpolation_def;
    // interpolation_def.path_to_tables =  "../src/resources/tables";
    interpolation_def.raw = true;

    // ----[ Sector ]---------------------------------------- //

    SectorFactory::Definition sec_def;

    sec_def.utility_def.epair_def.lpm_effect = true;
    sec_def.utility_def.brems_def.lpm_effect = true;

    sec_def.location = Sector::ParticleLocation::InsideDetector;

    sec_def.do_continuous_randomization = true;
    sec_def.do_exact_time_calculation = true;

    sec_def.scattering_model = ScatteringFactory::Moliere;

    sec_def.e_cut = 500;
    sec_def.v_cut = 0.05;

    sec_def.geometry_def.shape = GeometryFactory::Sphere;
    sec_def.geometry_def.inner_radius = 0.0;
    sec_def.geometry_def.radius = 1e18;

    sec_def.medium_def.type = MediumFactory::Ice;
    sec_def.medium_def.density_correction = 0.98;

    sector_defintions.push_back(sec_def);

    // Sector* sec2 = SectorFactory::Get().CreateSector(particle, sec_def, interpolation_def);
    //
    // collections.push_back(sec2);

    // --------------------------------------------------------------------- //
    // Create Propagator
    // --------------------------------------------------------------------- //

    // Propagator pr(collections, Sphere(Vector3D(), 1e18, 0));
    Propagator pr(particle, sector_defintions, Sphere(Vector3D(), 1e18, 0), interpolation_def);
    std::cout << "Propagagtor created" << std::endl;

    std::vector<unsigned int> length_sec;

    for(int i =0;i<statistic;i++)
    {
        particle.SetEnergy(pow(10,EmaxLog10));
        particle.SetPropagatedDistance(0);
        particle.SetPosition(Vector3D(0, 0, 0));
        particle.SetDirection(Vector3D(0, 0, -1));

        std::vector<DynamicData*> sec = pr.Propagate();

        length_sec.push_back(sec.size());

        // for (std::vector<DynamicData*>::iterator iter = sec.begin(); iter != sec.end(); ++iter) {
        //     std::cout << (*iter)->GetTypeId() << std::endl;
        // }
        std::cout << particle.GetPropagatedDistance() / 100.0 << std::endl;
    }

    // for(std::vector<unsigned int>::iterator it = length_sec.begin(); it != length_sec.end(); ++it)
    // {
    //     std::cout << "length of secondies: " << *it << std::endl;
    // }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}




