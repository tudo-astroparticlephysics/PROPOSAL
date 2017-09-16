
// #include <iostream>
// #include <string>
// #include <cmath>
#include <vector>
#include <cmath>

#include "gtest/gtest.h"

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/geometry/Sphere.h"
#include "PROPOSAL/geometry/Box.h"
// #include "PROPOSAL/sector/Sector.h"

// #include "PROPOSAL/sector/SectorInterpolant.h"
#include "PROPOSAL/sector/SectorFactory.h"
#include "PROPOSAL/medium/Medium.h"
// #include "PROPOSAL/geometry/Geometry.h"
// #include "PROPOSAL/PROPOSALParticle.h"

#include "PROPOSAL/crossection/parametrization/Bremsstrahlung.h"
#include "PROPOSAL/crossection/BremsstrahlungFactory.h"
#include "PROPOSAL/crossection/BremsInterpolant.h"

#include "PROPOSAL/crossection/parametrization/PhotoRealPhotonAssumption.h"
#include "PROPOSAL/crossection/parametrization/PhotoQ2Integration.h"
#include "PROPOSAL/crossection/PhotonuclearFactory.h"
#include "PROPOSAL/crossection/PhotoInterpolant.h"

using namespace std;
using namespace PROPOSAL;

TEST(Propagation , Test_nan) {

    int statistic = 10;
    int EmaxLog10 = 8;

    std::vector<Sector*> collections;

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

    SectorFactory::Definition sec_def;

    sec_def.location = Sector::ParticleLocation::InsideDetector;

    sec_def.do_interpolation = true;
    sec_def.do_continuous_randomization = true;
    sec_def.do_exact_time_calculation = true;
    sec_def.lpm_effect_enabled = false;

    sec_def.scattering_model = ScatteringFactory::MoliereFirstOrder;
    sec_def.density_correction = 0.98;

    sec_def.e_cut = 500;
    sec_def.v_cut = 0.05;

    sec_def.geometry = GeometryFactory::Sphere;
    sec_def.inner_radius = 0.0;
    sec_def.radius = 1e18;

    sec_def.medium = MediumFactory::Ice;

    sec_def.path_to_tables = "../src/resources/tables";
    sec_def.raw = false;

    Sector* sec2 = SectorFactory::Get().CreateSector(particle, sec_def);

    collections.push_back(sec2);

    // --------------------------------------------------------------------- //
    // Create Propagator
    // --------------------------------------------------------------------- //

    Propagator pr(collections, Sphere(Vector3D(), 1e18, 0));
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

        // for (std::vector<PROPOSALParticle*>::iterator iter = sec.begin(); iter != sec.end(); ++iter) {
        //     std::cout << (*iter)->GetName() << std::endl;
        // }
    }

    for(std::vector<unsigned int>::iterator it = length_sec.begin(); it != length_sec.end(); ++it)
    {
        std::cout << "length of secondies: " << *it << std::endl;
    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}




