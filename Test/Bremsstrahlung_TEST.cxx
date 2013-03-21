#include "gtest/gtest.h"
#include "PROPOSAL/Bremsstrahlung.h"
#include "PROPOSAL/CrossSections.h"
#include <iostream>
#include <string>

using namespace std;

TEST(Bremsstrahlung , Test_of_dEdx ) {

    ifstream in;
    in.open("bin/Brems_dEdx.txt");

    char firstLine[256];
    in.getline(firstLine,256);

    double energy;
    double dEdx;
    double ecut;
    double vcut;
    string med;
    string particleName;
    bool lpm;
    int para;

    CrossSections *brems = new Bremsstrahlung();

    while(in.good())
    {
        in>>para>>ecut>>vcut>>lpm>>energy>>med>>particleName>>dEdx;
        Medium *medium = new Medium(med,1.);
        brems->GetEnergyCutSettings()->SetEcut(ecut);
        brems->GetEnergyCutSettings()->SetVcut(vcut);
        brems->SetMedium(medium);
        brems->SetParametrization(para);
        brems->EnableLpmEffect(lpm);
        Particle *particle = new Particle(particleName,1.,1.,1,.20,20,1e5,10);
        particle->SetEnergy(energy);
        brems->SetParticle(particle);
        cout<<lpm<<"\t"<<brems->CalculatedEdx()<<"\t\t"<< dEdx<<endl;
        //ASSERT_FLOAT_EQ(brems->CalculatedEdx(), dEdx);
        delete medium;
        delete particle;



    }
    delete brems;
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
