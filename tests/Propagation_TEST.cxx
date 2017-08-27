
// #include <iostream>
// #include <string>
// #include <cmath>
#include <vector>

#include "gtest/gtest.h"

#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/CollectionInterpolant.h"
// #include "PROPOSAL/PROPOSALParticle.h"

using namespace std;
using namespace PROPOSAL;

TEST(Propagation , Test_nan) {

    int statistic = 1;
    int EmaxLog10 = 6;

    CollectionDef col_def;
    col_def.location = 1;
    col_def.do_continuous_randomization_ = true;

    CollectionInterpolant col_interpol(Water(), Sphere(Vector3D(), 1e18, 0), EnergyCutSettings(500, 0.5), col_def);

    std::vector<Collection*> collections;

    collections.push_back(&col_interpol);

    Propagator pr(collections, Sphere(Vector3D(), 1e18, 0));

    PROPOSALParticle particle(MuMinusDef::Get());
    // pr->EnableInterpolation(*pr->GetParticle());

    // pr->set_seed(seed);

    for(int i =0;i<statistic;i++)
    {
        particle.SetEnergy(pow(10,EmaxLog10));
        particle.SetPropagatedDistance(0);

        printf("in propagate: %s\n", particle.GetName().c_str());

        std::vector<PROPOSALParticle*> sec = pr.Propagate(particle);

        // for (std::vector<PROPOSALParticle*>::iterator iter = sec.begin(); iter != sec.end(); ++iter) {
        //     std::cout << (*iter)->GetName() << std::endl;
        // }
    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}




