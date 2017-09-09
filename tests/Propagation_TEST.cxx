
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
#include "PROPOSAL/crossection/BremsInterpolant.h"

#include "PROPOSAL/crossection/parametrization/PhotoRealPhotonAssumption.h"
#include "PROPOSAL/crossection/parametrization/PhotoQ2Integration.h"

using namespace std;
using namespace PROPOSAL;

TEST(Propagation , Test_nan) {

    int statistic = 10;
    int EmaxLog10 = 8;

    ParticleDef pdef = MuMinusDef::Get();

    PhotoRhode rhode(pdef, Ice(), EnergyCutSettings(500, 0.5), HardBB(pdef));
    PhotoZeus zeus(MuMinusDef::Get(), Ice(), EnergyCutSettings(500, 0.5), HardBB(pdef));
    PhotoBezrukovBugaev bezrukov(MuMinusDef::Get(), Ice(), EnergyCutSettings(500, 0.5), SoftBB());
    PhotoKokoulin kokoulin(MuMinusDef::Get(), Ice(), EnergyCutSettings(500, 0.5), SoftBB());

    // PhotoAbramowiczLevinLevyMaor97Interpolant allm(MuMinusDef::Get(), Ice(), EnergyCutSettings(500, 0.5), ShadowDutta());
    PhotoQ2Interpolant<PhotoAbramowiczLevinLevyMaor97> allm2(MuMinusDef::Get(), Ice(), EnergyCutSettings(500, 0.5), ShadowDutta());

    std::cout << "Result: " << allm2.DifferentialCrossSection(1e4, 0.01) << std::endl;

    // Parametrization::Definition def;
    // def.path_to_tables = "../src/resources/tables";
    // def.raw = false;
    // def.lpm_effect_enabled = false;
    // BremsKelnerKokoulinPetrukhin param(MuMinusDef::Get(), Ice(1.0), EnergyCutSettings(500, 0.05), def);
    // BremsInterpolant brems(param);

    // SectorDef col_def;
    // col_def.location = 1;
    // col_def.do_continuous_randomization_ = true;
    // col_def.do_scattering = true;
    // col_def.scattering_model = ScatteringFactory::ScatteringModel::MoliereFirstOrder;
    //
    //
    // SectorInterpolant col_interpol(Water(0.98), Sphere(Vector3D(), 1e18, 0), EnergyCutSettings(500, 0.5), col_def);
    //
    //
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
    sec_def.do_scattering = true;
    sec_def.lpm_effect_enabled = false;

    sec_def.scattering_model = ScatteringFactory::MoliereFirstOrder;
    sec_def.density_correction = 0.98;

    sec_def.e_cut = 500;
    sec_def.v_cut = 0.05;

    sec_def.geometry = GeometryFactory::Sphere;
    sec_def.inner_radius = 0.0;
    sec_def.radius = 1e18;

    sec_def.medium = MediumFactory::Ice;

    Sector* sec2 = SectorFactory::Get().CreateSector(particle, sec_def);
    std::cout << "Got sector" << std::endl;

    collections.push_back(sec2);

    // --------------------------------------------------------------------- //
    // Create Propagator
    // --------------------------------------------------------------------- //

    // Propagator pr(collections, Sphere(Vector3D(), 1e18, 0));
    // std::cout << "Propagagtor created" << std::endl;

    // Propagator pr("../resources/config_ice.json");
    // pr->EnableInterpolation(*pr->GetParticle());

    // pr->set_seed(seed);

    // std::vector<unsigned int> length_sec;
    //
    // for(int i =0;i<statistic;i++)
    // {
    //     std::cout << "loop: " << i << std::endl;
    //     particle.SetEnergy(pow(10,EmaxLog10));
    //     particle.SetPropagatedDistance(0);
    //     particle.SetPosition(Vector3D(0, 0, 0));
    //     particle.SetDirection(Vector3D(0, 0, -1));
    //
    //     std::vector<PROPOSALParticle*> sec = pr.Propagate();
    //
    //     length_sec.push_back(sec.size());
    //
    //     // for (std::vector<PROPOSALParticle*>::iterator iter = sec.begin(); iter != sec.end(); ++iter) {
    //     //     std::cout << (*iter)->GetName() << std::endl;
    //     // }
    // }
    //
    // for(std::vector<unsigned int>::iterator it = length_sec.begin(); it != length_sec.end(); ++it)
    // {
    //     std::cout << "length of secondies: " << *it << std::endl;
    // }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}




