#include "gtest/gtest.h"
#include "PROPOSAL/Ionization.h"
#include "PROPOSAL/CrossSections.h"
#include <iostream>
#include <string>

using namespace std;


TEST(Ionization , Test_of_dEdx ) {

    ifstream in;
    in.open("bin/Ioniz_dEdx.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dEdx_new;
    double energy;
    double dEdx;
    double ecut;
    double vcut;
    string med;
    string particleName;

    cout.precision(16);


    while(in.good())
    {
        in>>ecut>>vcut>>energy>>med>>particleName>>dEdx;

        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *ioniz = new Ionization(particle, medium, cuts);

        dEdx_new=ioniz->CalculatedEdx();
        ASSERT_NEAR(dEdx_new, dEdx, 1e-14*dEdx);

        delete cuts;
        delete medium;
        delete particle;
        delete ioniz;



    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
