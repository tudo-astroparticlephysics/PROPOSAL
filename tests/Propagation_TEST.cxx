
// #include <iostream>
// #include <string>
// #include <cmath>
#include <vector>
#include <cmath>

#include "gtest/gtest.h"

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/sector/CollectionInterpolant.h"
#include "PROPOSAL/geometry/Geometry.h"
// #include "PROPOSAL/PROPOSALParticle.h"

using namespace std;
using namespace PROPOSAL;

TEST(Propagation , Test_nan) {

    int statistic = 100;
    int EmaxLog10 = 8;

    // CollectionDef col_def;
    // col_def.location = 1;
    // col_def.do_continuous_randomization_ = true;
    // col_def.do_scattering = true;
    // col_def.scattering_model = ScatteringFactory::ScatteringModel::MoliereFirstOrder;
    //
    // CollectionInterpolant col_interpol(Water(), Sphere(Vector3D(), 1e18, 0), EnergyCutSettings(500, 0.5), col_def);
    //
    // std::vector<Collection*> collections;
    //
    // collections.push_back(&col_interpol);
    //
    // Propagator pr(collections, Sphere(Vector3D(), 1e18, 0));
    Propagator pr("../resources/config_ice.json");

    PROPOSALParticle particle(MuMinusDef::Get());
    particle.SetDirection(Vector3D(0, 0, -1));
    // pr->EnableInterpolation(*pr->GetParticle());

    // pr->set_seed(seed);

    std::vector<unsigned int> length_sec;

    for(int i =0;i<statistic;i++)
    {
        particle.SetEnergy(pow(10,EmaxLog10));
        particle.SetPropagatedDistance(0);

        std::vector<PROPOSALParticle*> sec = pr.Propagate(particle);
        length_sec.push_back(sec.size());
        // std::cout << particle.GetPropagatedDistance() / 100.0 << '\n';

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




