#include "gtest/gtest.h"
#include "PROPOSAL/Epairproduction.h"
#include "PROPOSAL/CrossSections.h"
#include <iostream>
#include <string>

using namespace std;


TEST(Epairproduction , Test_of_dEdx ) {

    ifstream in;
    in.open("bin/Epair_dEdx.txt");

    char firstLine[256];
    in.getline(firstLine,256);
    double dEdx_new;
    double energy;
    double dEdx;
    double ecut;
    double vcut;
    string med;
    string particleName;
    bool lpm;

    cout.precision(16);


    while(in.good())
    {
        in>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;
        Medium *medium = new Medium(med,1.);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        EnergyCutSettings *cuts = new EnergyCutSettings(ecut,vcut);

        CrossSections *epair = new Epairproduction(particle, medium, cuts);

        epair->EnableLpmEffect(lpm);

        dEdx_new=epair->CalculatedEdx();
        ASSERT_NEAR(dEdx_new, dEdx, 1e-5*dEdx);

        delete cuts;
        delete medium;
        delete particle;
        delete epair;



    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
